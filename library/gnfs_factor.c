/*
 * gnfs_factor.c - GNFS factoring driver for 90-digit semiprimes
 *
 * Approach: Use msieve for polynomial selection (or pre-computed poly),
 * GGNFS lattice siever for relation collection, msieve for post-processing.
 *
 * Compile: gcc -O2 -o gnfs_factor library/gnfs_factor.c -lgmp -lm
 * Usage:   timeout 295 ./gnfs_factor <N>
 *
 * Requires: gnfs-lasieve4I12e, msieve in hardcoded paths.
 * Single-threaded: siever runs as single process, msieve with -t 1.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include <gmp.h>

#define MSIEVE_PATH "/tmp/agent-factoring-3/msieve"
#define SIEVER_I12  "/tmp/ggnfs_build/gnfs-lasieve4I12e"
#define SIEVER_I11  "/tmp/ggnfs_build/gnfs-lasieve4I11e"

static struct timeval g_start;
static double elapsed_sec(void) {
    struct timeval now;
    gettimeofday(&now, NULL);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_usec - g_start.tv_usec) / 1e6;
}

/* Create msieve factor base file for polynomial selection */
static int do_poly_select(const char *workdir, const char *N, int deadline) {
    char cmd[4096];
    int remaining = 290 - (int)elapsed_sec();
    if (remaining < deadline) deadline = remaining - 5;
    if (deadline < 5) return -1;

    snprintf(cmd, sizeof(cmd),
        "cd '%s' && timeout %d %s -s msieve.dat -l msieve.log -t 1 -np "
        "\"polydegree=4\" \"poly_deadline=%d\" '%s' 2>/dev/null",
        workdir, deadline + 5, MSIEVE_PATH, deadline, N);

    int ret = system(cmd);
    char fb_path[512];
    snprintf(fb_path, sizeof(fb_path), "%s/msieve.dat.fb", workdir);
    struct stat st;
    return (stat(fb_path, &st) == 0) ? 0 : -1;
}

/* Convert msieve .fb file to GGNFS job format */
static int create_job_file(const char *workdir, const char *N) {
    char fb_path[512], job_path[512];
    snprintf(fb_path, sizeof(fb_path), "%s/msieve.dat.fb", workdir);
    snprintf(job_path, sizeof(job_path), "%s/gnfs.job", workdir);

    FILE *fb = fopen(fb_path, "r");
    if (!fb) {
        fprintf(stderr, "Cannot open %s: %s\n", fb_path, strerror(errno));
        return -1;
    }

    char skew[64] = "", c[6][64], Y[2][64];
    memset(c, 0, sizeof(c));
    memset(Y, 0, sizeof(Y));

    char line[1024];
    while (fgets(line, sizeof(line), fb)) {
        if (strncmp(line, "SKEW ", 5) == 0) sscanf(line+5, "%63s", skew);
        for (int i = 0; i < 6; i++) {
            char tag[8]; snprintf(tag, sizeof(tag), "A%d ", i);
            if (strncmp(line, tag, strlen(tag)) == 0)
                sscanf(line+strlen(tag), "%63s", c[i]);
        }
        if (strncmp(line, "R0 ", 3) == 0) sscanf(line+3, "%63s", Y[0]);
        if (strncmp(line, "R1 ", 3) == 0) sscanf(line+3, "%63s", Y[1]);
    }
    fclose(fb);

    FILE *job = fopen(job_path, "w");
    if (!job) return -1;

    fprintf(job, "n: %s\n", N);
    if (skew[0]) fprintf(job, "skew: %s\n", skew);
    for (int i = 0; i < 6; i++)
        if (c[i][0]) fprintf(job, "c%d: %s\n", i, c[i]);
    for (int i = 0; i < 2; i++)
        if (Y[i][0]) fprintf(job, "Y%d: %s\n", i, Y[i]);

    /* Gimarel-optimized parameters for 90d */
    fprintf(job, "rlim: 350000\n");
    fprintf(job, "alim: 840000\n");
    fprintf(job, "lpbr: 25\n");
    fprintf(job, "lpba: 25\n");
    fprintf(job, "mfbr: 50\n");
    fprintf(job, "mfba: 50\n");
    fprintf(job, "rlambda: 2.400\n");
    fprintf(job, "alambda: 2.400\n");

    fclose(job);
    return 0;
}

/* Write pre-computed polynomial directly as GGNFS job file */
static int write_precomputed_job(const char *workdir, const char *N,
                                  const char *poly_data) {
    char job_path[512];
    snprintf(job_path, sizeof(job_path), "%s/gnfs.job", workdir);

    FILE *job = fopen(job_path, "w");
    if (!job) return -1;
    fprintf(job, "%s", poly_data);
    fclose(job);
    return 0;
}

/* Run GGNFS lattice siever and collect relations */
static int run_sieve(const char *workdir, int sieve_budget_sec,
                      const char *siever_path) {
    int startq = 210000;
    int qrange = 5000;
    int batch = 0;
    int total_rels = 0;
    int target = 1460000;
    char rels_file[512];
    snprintf(rels_file, sizeof(rels_file), "%s/msieve.dat", workdir);

    while (elapsed_sec() < sieve_budget_sec && total_rels < target) {
        int cur_start = startq + batch * qrange;
        int remaining = sieve_budget_sec - (int)elapsed_sec();
        if (remaining < 3) break;

        char batch_file[256];
        snprintf(batch_file, sizeof(batch_file), "%s/rels_%d.dat", workdir, batch);

        char cmd[4096];
        snprintf(cmd, sizeof(cmd),
            "cd '%s' && timeout %d '%s' -f %d -c %d -o rels_%d.dat -a gnfs.job 2>/dev/null",
            workdir, remaining, siever_path, cur_start, qrange, batch);

        int ret = system(cmd);
        int exit_code = WIFEXITED(ret) ? WEXITSTATUS(ret) : -1;

        /* Count new relations */
        struct stat st;
        if (stat(batch_file, &st) == 0 && st.st_size > 0) {
            char count_cmd[512];
            snprintf(count_cmd, sizeof(count_cmd), "wc -l < '%s'", batch_file);
            FILE *p = popen(count_cmd, "r");
            int new_rels = 0;
            if (p) { fscanf(p, "%d", &new_rels); pclose(p); }

            /* Append to msieve data file */
            char cat_cmd[1024];
            snprintf(cat_cmd, sizeof(cat_cmd), "cat '%s' >> '%s'", batch_file, rels_file);
            system(cat_cmd);

            total_rels += new_rels;
            fprintf(stderr, "[%.1fs] q=%d-%d: +%d rels (total: %d/%d)\n",
                    elapsed_sec(), cur_start, cur_start+qrange, new_rels, total_rels, target);
        } else if (exit_code != 0 && exit_code != 124) {
            fprintf(stderr, "Siever error (exit %d) at q=%d\n", exit_code, cur_start);
        }

        batch++;
    }

    fprintf(stderr, "[%.1fs] Sieving done: %d relations\n", elapsed_sec(), total_rels);
    return total_rels;
}

/* Run msieve for filtering + linear algebra + square root */
static int run_postprocess(const char *workdir, const char *N) {
    int remaining = 293 - (int)elapsed_sec();
    if (remaining < 10) {
        fprintf(stderr, "Not enough time for post-processing (%ds left)\n", remaining);
        return -1;
    }

    char cmd[4096];
    snprintf(cmd, sizeof(cmd),
        "cd '%s' && timeout %d %s -s msieve.dat -l msieve.log -t 1 -nc '%s' 2>/dev/null",
        workdir, remaining, MSIEVE_PATH, N);

    fprintf(stderr, "[%.1fs] Post-processing (filter+LA+sqrt)...\n", elapsed_sec());
    int ret = system(cmd);
    fprintf(stderr, "[%.1fs] Post-processing done\n", elapsed_sec());
    return ret;
}

/* Extract factors from msieve log */
static int extract_factors(const char *workdir) {
    char logpath[512];
    snprintf(logpath, sizeof(logpath), "%s/msieve.log", workdir);
    FILE *f = fopen(logpath, "r");
    if (!f) return -1;

    char line[1024];
    int found = 0;
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "prp") || strstr(line, "factor:")) {
            /* Extract just the factor value */
            char *p = strstr(line, "prp");
            if (!p) p = strstr(line, "factor:");
            if (p) {
                /* Find the number */
                while (*p && (*p < '0' || *p > '9')) p++;
                char *end = p;
                while (*end >= '0' && *end <= '9') end++;
                *end = '\0';
                if (strlen(p) > 5) {
                    printf("%s\n", p);
                    found++;
                }
            }
        }
    }
    fclose(f);
    return found > 0 ? 0 : -1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    const char *N = argv[1];
    gettimeofday(&g_start, NULL);

    /* Create work directory */
    char workdir[256];
    snprintf(workdir, sizeof(workdir), "/tmp/gnfs_work_%d", getpid());
    mkdir(workdir, 0755);

    fprintf(stderr, "=== GNFS factoring: %s ===\n", N);

    /* Step 1: Polynomial selection (or use pre-computed) */
    fprintf(stderr, "[%.1fs] Polynomial selection...\n", elapsed_sec());
    int ret = do_poly_select(workdir, N, 25);
    if (ret != 0) {
        fprintf(stderr, "Poly select failed, aborting\n");
        goto cleanup;
    }

    /* Create GGNFS job file */
    ret = create_job_file(workdir, N);
    if (ret != 0) {
        fprintf(stderr, "Failed to create job file\n");
        goto cleanup;
    }
    fprintf(stderr, "[%.1fs] Polynomial ready\n", elapsed_sec());

    /* Step 2: GGNFS sieving */
    /* Budget: leave 20s for post-processing */
    int sieve_budget = 275;
    int total_rels = run_sieve(workdir, sieve_budget, SIEVER_I12);

    if (total_rels < 500000) {
        fprintf(stderr, "Too few relations (%d), unlikely to succeed\n", total_rels);
    }

    /* Step 3: Post-processing */
    ret = run_postprocess(workdir, N);

    /* Step 4: Extract and print factors */
    if (extract_factors(workdir) != 0) {
        fprintf(stderr, "No factors found\n");
        ret = 1;
    }

    fprintf(stderr, "[%.1fs] Total time\n", elapsed_sec());

cleanup:
    {
        char cmd[512];
        snprintf(cmd, sizeof(cmd), "rm -rf '%s'", workdir);
        system(cmd);
    }
    return ret;
}
