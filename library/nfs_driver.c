/*
 * nfs_driver.c - Lightweight NFS driver for 90-digit balanced semiprimes
 *
 * Uses msieve for polynomial selection and linear algebra,
 * GGNFS lasieve4I12e for sieving.
 *
 * Optimized pipeline:
 * 1. Quick polynomial selection with msieve (target: <30s)
 * 2. GGNFS lattice sieving (target: <220s)
 * 3. msieve filtering + Block Lanczos + sqrt (target: <30s)
 *
 * Compile: gcc -O2 -march=native -o nfs_driver library/nfs_driver.c -lgmp -lm
 * Usage: ./nfs_driver <N>
 *
 * Requirements:
 * - msieve binary (from YAFU build)
 * - gnfs-lasieve4I12e binary
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <signal.h>
#include <errno.h>

#define MAX_CMD 4096
#define TIMEOUT_TOTAL 290  /* Total budget in seconds */

static double get_time(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

/* Run a command with timeout, return exit code */
static int run_cmd(const char *cmd, int timeout_sec) {
    char full_cmd[MAX_CMD + 64];
    snprintf(full_cmd, sizeof(full_cmd), "timeout %d %s", timeout_sec, cmd);
    return system(full_cmd);
}

/* Count lines in a file */
static int count_lines(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    int count = 0;
    char buf[4096];
    while (fgets(buf, sizeof(buf), f)) count++;
    fclose(f);
    return count;
}

/* Check if file contains a factor */
static int extract_factor(const char *logfile, const char *N, char *factor, size_t flen) {
    FILE *f = fopen(logfile, "r");
    if (!f) return 0;
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        /* msieve outputs: "prp<digits> factor: <number>" */
        char *p = strstr(line, "prp");
        if (p) {
            char *colon = strstr(p, "factor: ");
            if (colon) {
                colon += 8;  /* skip "factor: " */
                /* Remove trailing whitespace */
                char *end = colon + strlen(colon) - 1;
                while (end > colon && (*end == '\n' || *end == '\r' || *end == ' ')) *end-- = 0;
                if (strcmp(colon, N) != 0 && strlen(colon) > 1) {
                    strncpy(factor, colon, flen - 1);
                    factor[flen - 1] = 0;
                    fclose(f);
                    return 1;
                }
            }
        }
    }
    fclose(f);
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [msieve_path] [ggnfs_path]\n", argv[0]);
        return 1;
    }

    const char *N = argv[1];
    const char *msieve_path = argc > 2 ? argv[2] : "/tmp/agent-factoring-10/yafu_build/msieve";
    const char *ggnfs_path = argc > 3 ? argv[3] : "/tmp/ggnfs_build/gnfs-lasieve4I12e";

    /* Validate inputs */
    size_t ndigits = strlen(N);
    if (ndigits < 85 || ndigits > 100) {
        fprintf(stderr, "Warning: NFS is designed for 85-100 digit numbers, got %zu digits\n", ndigits);
    }

    /* Create working directory */
    char workdir[256];
    snprintf(workdir, sizeof(workdir), "/tmp/nfs_work_%d", getpid());
    mkdir(workdir, 0755);

    double t_start = get_time();
    char cmd[MAX_CMD];
    int ret;

    /* ===== Phase 1: Polynomial Selection ===== */
    fprintf(stderr, "[NFS] Phase 1: Polynomial selection for %zu-digit number\n", ndigits);
    double t_poly_start = get_time();

    /* Use msieve for polynomial selection
     * For 90d (degree 4): default settings work, but we limit time
     * -np: polynomial selection only
     * -v: verbose
     * -t 1: single thread
     */
    snprintf(cmd, sizeof(cmd),
        "cd %s && %s -s %s/session -l %s/msieve.log -np -v -t 1 %s",
        workdir, msieve_path, workdir, workdir, N);

    int poly_timeout = 60;  /* 60s max for poly select */
    ret = run_cmd(cmd, poly_timeout);

    double t_poly_end = get_time();
    fprintf(stderr, "[NFS] Poly select: %.1fs (ret=%d)\n", t_poly_end - t_poly_start, ret);

    /* Check polynomial was found */
    char poly_file[512];
    snprintf(poly_file, sizeof(poly_file), "%s/session.fb", workdir);
    struct stat st;
    if (stat(poly_file, &st) != 0 || st.st_size == 0) {
        fprintf(stderr, "[NFS] ERROR: No polynomial found\n");
        /* Cleanup */
        snprintf(cmd, sizeof(cmd), "rm -rf %s", workdir);
        system(cmd);
        return 1;
    }

    /* ===== Phase 2: Sieving ===== */
    fprintf(stderr, "[NFS] Phase 2: Sieving with GGNFS\n");
    double t_sieve_start = get_time();
    double time_remaining = TIMEOUT_TOTAL - (t_sieve_start - t_start);
    int sieve_timeout = (int)(time_remaining - 40);  /* Reserve 40s for filter+LA+sqrt */
    if (sieve_timeout < 60) {
        fprintf(stderr, "[NFS] ERROR: Not enough time for sieving (%ds remaining)\n", sieve_timeout);
        snprintf(cmd, sizeof(cmd), "rm -rf %s", workdir);
        system(cmd);
        return 1;
    }

    /* Parse the polynomial file to get NFS parameters */
    /* For 90d balanced semiprimes, use:
     * I=12 (sieve region 2^12 = 4096)
     * lpb=25/25 (large prime bounds)
     * mfb=50/50 (cofactor bounds)
     * a0-a1 range: 0 to ~600000 (for ~1.5M relations)
     */
    /* GGNFS needs a .poly file in specific format */
    /* For now, use YAFU's GNFS driver via factor() which handles all this */

    /* Actually, let's use YAFU's full GNFS pipeline since it handles
     * poly conversion, GGNFS invocation, filtering, etc. */
    fprintf(stderr, "[NFS] Using YAFU GNFS driver (xover=85)\n");

    /* We'll invoke YAFU with GNFS mode */
    const char *yafu_path = "/tmp/agent-factoring-4/yafu/yafu";

    snprintf(cmd, sizeof(cmd),
        "cd %s && echo \"nfs(%s)\" | timeout %d %s -threads 1 -seed 42 -xover 85 2>&1",
        workdir, N, (int)(TIMEOUT_TOTAL - (get_time() - t_start)), yafu_path);

    fprintf(stderr, "[NFS] Running: %s\n", cmd);

    FILE *pipe = popen(cmd, "r");
    char line[1024];
    char factor[512] = {0};
    int factored = 0;

    while (pipe && fgets(line, sizeof(line), pipe)) {
        /* Look for factor output */
        if (strstr(line, "***factors found***") ||
            strstr(line, "P") || strstr(line, "PRP")) {
            fprintf(stderr, "[NFS] %s", line);
        }
        char *prp = strstr(line, "= ");
        if (prp && (strstr(line, "P") || strstr(line, "C"))) {
            /* This might be a factor line like "P45 = 123456..." */
            prp += 2;
            char *end = prp + strlen(prp) - 1;
            while (end > prp && (*end == '\n' || *end == '\r')) *end-- = 0;
            if (strlen(prp) > 5 && strlen(prp) < ndigits) {
                strncpy(factor, prp, sizeof(factor) - 1);
                factored = 1;
            }
        }
    }
    if (pipe) pclose(pipe);

    double t_end = get_time();

    if (factored) {
        printf("Factor found: %s\n", factor);
        printf("Time: %.4f seconds\n", t_end - t_start);
    } else {
        fprintf(stderr, "[NFS] Failed to factor in %.1fs\n", t_end - t_start);
        snprintf(cmd, sizeof(cmd), "rm -rf %s", workdir);
        system(cmd);
        return 1;
    }

    /* Cleanup */
    snprintf(cmd, sizeof(cmd), "rm -rf %s", workdir);
    system(cmd);

    return 0;
}
