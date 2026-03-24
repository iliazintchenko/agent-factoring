/*
 * dbrm.c — Aggressive Large-Prime Sieve factorer
 *
 * Uses polynomial evaluation Q(x) = (x+m)^2 - N with aggressive large prime
 * acceptance and graph-based relation combining.
 *
 * Usage: ./dbrm <N>
 * Output: "factor1 factor2" on stdout
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ========== Utility ========== */

static unsigned long power_mod_ul(unsigned long base, unsigned long exp, unsigned long mod) {
    unsigned long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = ((__uint128_t)result * base) % mod;
        base = ((__uint128_t)base * base) % mod;
        exp >>= 1;
    }
    return result;
}

static unsigned long tonelli_shanks(unsigned long n, unsigned long p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;
    if (power_mod_ul(n, (p - 1) / 2, p) != 1) return 0;
    unsigned long Q = p - 1; int S = 0;
    while ((Q & 1) == 0) { Q >>= 1; S++; }
    if (S == 1) return power_mod_ul(n, (p + 1) / 4, p);
    unsigned long z = 2;
    while (power_mod_ul(z, (p - 1) / 2, p) != p - 1) z++;
    unsigned long Mv = S, c = power_mod_ul(z, Q, p);
    unsigned long t = power_mod_ul(n, Q, p), R = power_mod_ul(n, (Q + 1) / 2, p);
    for (;;) {
        if (t == 1) return R;
        unsigned long tt = t; unsigned long i = 0;
        while (tt != 1) { tt = ((__uint128_t)tt * tt) % p; i++; }
        if (i == Mv) return 0;
        unsigned long b = c;
        for (unsigned long j = 0; j < Mv - i - 1; j++) b = ((__uint128_t)b * b) % p;
        Mv = i; c = ((__uint128_t)b * b) % p;
        t = ((__uint128_t)t * c) % p; R = ((__uint128_t)R * b) % p;
    }
}

static int *sieve_primes(int limit, int *count) {
    char *comp = calloc(limit + 1, 1);
    int *primes = malloc(sizeof(int) * (limit / 3 + 100));
    int cnt = 0;
    for (int i = 2; i <= limit; i++) {
        if (!comp[i]) {
            primes[cnt++] = i;
            if ((long long)i * i <= limit)
                for (int j = i * i; j <= limit; j += i) comp[j] = 1;
        }
    }
    free(comp); *count = cnt; return primes;
}

/* ========== Relation types ========== */

/* Exponent vector stored as list of factor-base indices with odd exponent */
typedef struct {
    long x;            /* sieve position */
    int *fb_exp;       /* FB indices with odd exponent (col 0 = sign) */
    int nfb;
    unsigned long lp;  /* large prime (0 = fully smooth) */
} Relation;

/* ========== GF(2) dense matrix ========== */
/* Store matrix as array of uint64_t bitmasks, one row per word-array */

typedef struct {
    int nrows, ncols;
    int words_per_row;
    unsigned long long *data; /* row-major, each row = words_per_row words */
    int *row_origin;          /* which relation(s) contributed to this row */
    int *origin_count;
    int **origins;            /* origins[i] = list of relation indices */
} GF2Matrix;

static void gf2_init(GF2Matrix *m, int nrows, int ncols) {
    m->nrows = nrows;
    m->ncols = ncols;
    m->words_per_row = (ncols + 63) / 64;
    m->data = calloc((size_t)nrows * m->words_per_row, sizeof(unsigned long long));
    m->origins = malloc(sizeof(int*) * nrows);
    m->origin_count = malloc(sizeof(int) * nrows);
    for (int i = 0; i < nrows; i++) {
        m->origins[i] = malloc(sizeof(int));
        m->origins[i][0] = i;
        m->origin_count[i] = 1;
    }
}

static void gf2_set(GF2Matrix *m, int row, int col) {
    m->data[(size_t)row * m->words_per_row + col / 64] ^= 1ULL << (col % 64);
}

static int gf2_get(GF2Matrix *m, int row, int col) {
    return (m->data[(size_t)row * m->words_per_row + col / 64] >> (col % 64)) & 1;
}

static void gf2_xor_rows(GF2Matrix *m, int dst, int src) {
    unsigned long long *d = m->data + (size_t)dst * m->words_per_row;
    unsigned long long *s = m->data + (size_t)src * m->words_per_row;
    for (int w = 0; w < m->words_per_row; w++) d[w] ^= s[w];

    /* Merge origins (symmetric difference for GF(2)) */
    int *oa = m->origins[dst], na = m->origin_count[dst];
    int *ob = m->origins[src], nb = m->origin_count[src];
    int *merged = malloc(sizeof(int) * (na + nb));
    int i = 0, j = 0, k = 0;
    /* Sort both first */
    /* origins should already be sorted */
    while (i < na && j < nb) {
        if (oa[i] < ob[j]) merged[k++] = oa[i++];
        else if (oa[i] > ob[j]) merged[k++] = ob[j++];
        else { i++; j++; }
    }
    while (i < na) merged[k++] = oa[i++];
    while (j < nb) merged[k++] = ob[j++];
    free(m->origins[dst]);
    m->origins[dst] = merged;
    m->origin_count[dst] = k;
}

static int gf2_is_zero(GF2Matrix *m, int row) {
    unsigned long long *r = m->data + (size_t)row * m->words_per_row;
    for (int w = 0; w < m->words_per_row; w++)
        if (r[w]) return 0;
    return 1;
}

/* ========== Parameters ========== */

static void choose_params(int digits, int *B, long *M, unsigned long *lp_bound) {
    /* Balanced parameters: moderate B (fast trial div), large M (many smooth near x=0),
     * lp_bound = B*100 for good LP collision rate */
    /* Continuous L[1/2] scaling, calibrated from empirical testing.
     * B controls factor base size (and smooth probability).
     * M controls sieve range (must be large enough to find fb_size+ relations).
     */
    double ln_N = digits * 2.3026;
    double ln_ln_N = log(ln_N);
    double sqr = sqrt(ln_N * ln_ln_N);

    /* B = L[1/2, 0.45]: conservative to keep FB small and trial div fast */
    *B = (int)(exp(0.45 * sqr));
    if (*B < 200) *B = 200;
    if (*B > 8000000) *B = 8000000;

    /* M = L[1/2, 1.0]: generous sieve range for enough smooth values */
    double M_d = exp(1.0 * sqr);
    *M = (long)M_d;
    if (*M < 10000) *M = 10000;
    /* Cap memory at ~2GB for sieve array (float = 4 bytes) */
    if (*M > 250000000L) *M = 250000000L;

    *lp_bound = (unsigned long)(*B) * 60;
    if (*lp_bound < 10000) *lp_bound = 10000;
}

/* ========== Main ========== */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    #define ELAPSED() ({ struct timespec _t; clock_gettime(CLOCK_MONOTONIC, &_t); \
        (_t.tv_sec - t0.tv_sec) + (_t.tv_nsec - t0.tv_nsec) / 1e9; })

    mpz_t N, sqrtN, tmp, tmp2;
    mpz_inits(N, sqrtN, tmp, tmp2, NULL);
    mpz_set_str(N, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N, 10);

    /* Trial division */
    { int np; int *sp = sieve_primes(1000000, &np);
      for (int i = 0; i < np; i++) {
          if (mpz_divisible_ui_p(N, sp[i])) {
              mpz_divexact_ui(tmp, N, sp[i]);
              gmp_printf("%d %Zd\n", sp[i], tmp); free(sp); return 0;
          }
      } free(sp); }

    int B; long M; unsigned long lp_bound;
    choose_params(digits, &B, &M, &lp_bound);

    mpz_sqrt(sqrtN, N);
    mpz_mul(tmp, sqrtN, sqrtN);
    if (mpz_cmp(tmp, N) < 0) mpz_add_ui(sqrtN, sqrtN, 1);

    /* Build factor base with precomputed x-roots for Q(x) ≡ 0 (mod p) */
    int np;
    int *all_p = sieve_primes(B, &np);
    int *fb = malloc(sizeof(int) * (np + 1));
    unsigned long *fb_root = malloc(sizeof(unsigned long) * (np + 1));
    long *fb_xr1 = malloc(sizeof(long) * (np + 1)); /* x ≡ r - m (mod p) */
    long *fb_xr2 = malloc(sizeof(long) * (np + 1)); /* x ≡ -r - m (mod p) */
    int fb_size = 0;

    for (int i = 0; i < np; i++) {
        int p = all_p[i];
        unsigned long nmod = mpz_fdiv_ui(N, p);
        unsigned long r;
        if (p == 2) { r = nmod; fb[fb_size] = 2; fb_root[fb_size] = r; fb_xr1[fb_size] = 0; fb_xr2[fb_size] = 0; fb_size++; continue; }
        r = tonelli_shanks(nmod, p);
        if (r != 0 || nmod == 0) {
            fb[fb_size] = p;
            fb_root[fb_size] = r;
            unsigned long sqmod = mpz_fdiv_ui(sqrtN, p);
            fb_xr1[fb_size] = ((long)r - (long)sqmod + 2L*(long)p) % p;
            fb_xr2[fb_size] = ((long)(p - r) - (long)sqmod + 2L*(long)p) % p;
            fb_size++;
        }
    }
    free(all_p);

    fprintf(stderr, "[%.1fs] %dd, B=%d (fb=%d), M=%ld, lp=%lu\n",
            ELAPSED(), digits, B, fb_size, M, lp_bound);

    /* ===== SIEVE ===== */
    long sieve_len = 2 * M + 1;
    float *sieve = calloc(sieve_len, sizeof(float));
    if (!sieve) { fprintf(stderr, "OOM sieve %ld\n", sieve_len); return 1; }

    for (int i = 0; i < fb_size; i++) {
        int p = fb[i];
        if (p < 3) continue; /* skip 2 for simplicity, trial div handles it */
        float logp = log2f((float)p);
        unsigned long sqmod = mpz_fdiv_ui(sqrtN, p);
        long r1 = ((long)fb_root[i] - (long)sqmod + 2L*(long)p) % p;
        long r2 = ((long)(p - fb_root[i]) - (long)sqmod + 2L*(long)p) % p;
        long off1 = (r1 + M) % p;
        long off2 = (r2 + M) % p;
        for (long j = off1; j < sieve_len; j += p) sieve[j] += logp;
        if (r1 != r2) for (long j = off2; j < sieve_len; j += p) sieve[j] += logp;

        /* Sieve prime powers */
        long pk = (long)p * p;
        for (int pw = 2; pk <= (long)B && pk > 0; pw++, pk *= p) {
            unsigned long nmod_pk = mpz_fdiv_ui(N, pk);
            /* Find roots of x^2 ≡ N mod p^k via Hensel lift */
            /* Simplified: just use the base root and check */
            unsigned long rr = fb_root[i];
            /* Hensel lift: r_new = r - (r^2 - N) / (2r) mod p^k */
            /* Skip for now, base prime sieve is most important */
            break;
        }
    }

    /* Threshold: sieve accumulates log2(p) for primes dividing Q(x).
     * For a 1LP value: accumulated ≈ log2(Q(x)) - log2(large_prime).
     * Q(x) ranges from ~2*sqrt(N) (at x=0) to ~2*sqrt(N)*M (at x=M).
     * Use conservative threshold based on MINIMUM Q value. */
    float log2_Qmin = (float)(mpz_sizeinbase(N, 2) / 2.0 + 1.0); /* log2(~2*sqrt(N)) */
    float log2_lp = log2f((float)lp_bound);
    /* Generous slack: factor of 2 missing (not sieved), small primes, prime powers */
    float threshold = log2_Qmin - log2_lp - 10.0f;
    if (threshold < 3.0f) threshold = 3.0f;

    fprintf(stderr, "[%.1fs] Sieve done. thresh=%.1f Scanning...\n", ELAPSED(), threshold);

    /* ===== TRIAL DIVISION + RELATION COLLECTION ===== */
    /* Collect smooth (0LP) and single-large-prime (1LP) relations */
    int rel_cap = 4096;
    Relation *smooth_rels = malloc(sizeof(Relation) * rel_cap);
    int n_smooth = 0;
    Relation *partial_rels = malloc(sizeof(Relation) * rel_cap * 4);
    int n_partial = 0;
    int partial_cap = rel_cap * 4;

    mpz_t qx, cofactor, base_val;
    mpz_inits(qx, cofactor, base_val, NULL);

    long checked = 0;
    for (long idx = 0; idx < sieve_len; idx++) {
        if (sieve[idx] < threshold) continue;
        checked++;

        long x = idx - M;
        mpz_set_si(base_val, x);
        mpz_add(base_val, base_val, sqrtN);
        mpz_mul(qx, base_val, base_val);
        mpz_sub(qx, qx, N);
        int is_neg = (mpz_sgn(qx) < 0);
        if (mpz_sgn(qx) == 0) continue;
        mpz_abs(cofactor, qx);

        int fb_buf[512];
        int nfb = 0;
        if (is_neg) fb_buf[nfb++] = 0;

        /* FAST trial division: use precomputed roots to skip non-dividing primes.
         * Q(x) ≡ 0 (mod p) iff x ≡ xr1[i] or xr2[i] (mod p). Only check those. */
        {
            /* Handle p=2 specially */
            if (fb[0] == 2) {
                int exp = 0;
                while (mpz_divisible_ui_p(cofactor, 2)) {
                    mpz_divexact_ui(cofactor, cofactor, 2);
                    exp++;
                }
                if (exp & 1) fb_buf[nfb++] = 1;
            }
            /* Handle odd primes using root pre-screening */
            int start_i = (fb[0] == 2) ? 1 : 0;
            for (int i = start_i; i < fb_size; i++) {
                unsigned long p = fb[i];
                long xmod = x % (long)p;
                if (xmod < 0) xmod += p;
                if (xmod != fb_xr1[i] && xmod != fb_xr2[i]) continue;
                /* p divides Q(x). Trial divide to get full exponent. */
                int exp = 0;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    exp++;
                }
                if (exp & 1) fb_buf[nfb++] = i + 1;
            }
        }

        unsigned long cof = 0;
        if (mpz_cmp_ui(cofactor, 1) == 0) {
            cof = 0; /* smooth */
        } else if (mpz_fits_ulong_p(cofactor)) {
            cof = mpz_get_ui(cofactor);
            if (cof > lp_bound) continue; /* too large */
            /* Quick primality check: skip if obviously composite */
            if (cof > 3 && (cof % 2 == 0 || cof % 3 == 0)) continue;
            mpz_set_ui(tmp, cof);
            if (!mpz_probab_prime_p(tmp, 2)) continue;
        } else {
            continue;
        }

        /* Store relation */
        int *fb_exp = malloc(sizeof(int) * nfb);
        memcpy(fb_exp, fb_buf, sizeof(int) * nfb);

        if (cof == 0) {
            if (n_smooth >= rel_cap) { rel_cap *= 2; smooth_rels = realloc(smooth_rels, sizeof(Relation) * rel_cap); }
            smooth_rels[n_smooth++] = (Relation){x, fb_exp, nfb, 0};
        } else {
            if (n_partial >= partial_cap) { partial_cap *= 2; partial_rels = realloc(partial_rels, sizeof(Relation) * partial_cap); }
            partial_rels[n_partial++] = (Relation){x, fb_exp, nfb, cof};
        }

        if (ELAPSED() > 250) { fprintf(stderr, "Time limit\n"); break; }
    }
    free(sieve);

    fprintf(stderr, "[%.1fs] %d smooth, %d partial (1LP), %ld checked\n",
            ELAPSED(), n_smooth, n_partial, checked);

    /* ===== LARGE PRIME MATCHING ===== */
    /* Sort partials by LP, then combine pairs with same LP */

    /* Sort by LP */
    int cmp_rel_lp(const void *a, const void *b) {
        unsigned long la = ((Relation*)a)->lp, lb = ((Relation*)b)->lp;
        return (la > lb) - (la < lb);
    }
    qsort(partial_rels, n_partial, sizeof(Relation), cmp_rel_lp);

    /* Combine pairs: for each LP appearing 2+ times, take pairs */
    /* A "combined" relation: XOR of two 1LP relations sharing same LP */
    /* The combined relation has the LP cancelled (even exponent) */
    int combined_cap = 4096;
    Relation *combined_rels = malloc(sizeof(Relation) * combined_cap);
    int n_combined = 0;

    int i = 0;
    while (i < n_partial) {
        int j = i;
        while (j < n_partial && partial_rels[j].lp == partial_rels[i].lp) j++;
        /* partial_rels[i..j-1] share the same LP */
        int group_size = j - i;
        /* Combine pairs: take consecutive pairs */
        for (int k = 0; k + 1 < group_size; k += 2) {
            Relation *a = &partial_rels[i + k];
            Relation *b = &partial_rels[i + k + 1];
            /* XOR exponent vectors */
            int *merged = malloc(sizeof(int) * (a->nfb + b->nfb));
            int mi = 0, bi2 = 0, mk = 0;
            /* Sort both (should already be sorted) */
            /* Simple merge-XOR */
            int ai = 0;
            while (ai < a->nfb && bi2 < b->nfb) {
                if (a->fb_exp[ai] < b->fb_exp[bi2]) merged[mk++] = a->fb_exp[ai++];
                else if (a->fb_exp[ai] > b->fb_exp[bi2]) merged[mk++] = b->fb_exp[bi2++];
                else { ai++; bi2++; } /* cancel */
            }
            while (ai < a->nfb) merged[mk++] = a->fb_exp[ai++];
            while (bi2 < b->nfb) merged[mk++] = b->fb_exp[bi2++];

            if (n_combined >= combined_cap) {
                combined_cap *= 2;
                combined_rels = realloc(combined_rels, sizeof(Relation) * combined_cap);
            }
            combined_rels[n_combined].fb_exp = merged;
            combined_rels[n_combined].nfb = mk;
            combined_rels[n_combined].x = a->x; /* store first x; need both for extraction */
            combined_rels[n_combined].lp = b->x; /* abuse: store second x in lp field */
            n_combined++;
        }
        i = j;
    }

    int total_rels = n_smooth + n_combined;
    fprintf(stderr, "[%.1fs] %d combined from LP matching. Total: %d rels for %d cols\n",
            ELAPSED(), n_combined, total_rels, fb_size + 1);

    if (total_rels < fb_size + 2) {
        fprintf(stderr, "Not enough relations (%d < %d). Need more sieving.\n",
                total_rels, fb_size + 2);
        return 1;
    }

    /* ===== BUILD MATRIX ===== */
    int ncols = fb_size + 1; /* col 0 = sign, cols 1..fb_size = factor base */
    GF2Matrix mat;
    gf2_init(&mat, total_rels, ncols);

    /* Add smooth relations */
    for (int r = 0; r < n_smooth; r++) {
        for (int c = 0; c < smooth_rels[r].nfb; c++)
            gf2_set(&mat, r, smooth_rels[r].fb_exp[c]);
    }
    /* Add combined relations */
    for (int r = 0; r < n_combined; r++) {
        int row = n_smooth + r;
        for (int c = 0; c < combined_rels[r].nfb; c++)
            gf2_set(&mat, row, combined_rels[r].fb_exp[c]);
    }

    /* ===== GAUSSIAN ELIMINATION ===== */
    fprintf(stderr, "[%.1fs] Gauss elim on %dx%d...\n", ELAPSED(), total_rels, ncols);

    int *pivot_col = malloc(sizeof(int) * ncols);
    for (int c = 0; c < ncols; c++) pivot_col[c] = -1;

    for (int r = 0; r < total_rels; r++) {
        /* Find leading column */
        int lead = -1;
        unsigned long long *row = mat.data + (size_t)r * mat.words_per_row;
        for (int w = 0; w < mat.words_per_row; w++) {
            if (row[w]) { lead = w * 64 + __builtin_ctzll(row[w]); break; }
        }

        while (lead >= 0 && lead < ncols && pivot_col[lead] >= 0) {
            gf2_xor_rows(&mat, r, pivot_col[lead]);
            lead = -1;
            row = mat.data + (size_t)r * mat.words_per_row;
            for (int w = 0; w < mat.words_per_row; w++) {
                if (row[w]) { lead = w * 64 + __builtin_ctzll(row[w]); break; }
            }
        }

        if (lead < 0 || lead >= ncols) {
            /* Zero row = dependency! Try to extract factor */
            mpz_t x_acc, y_sq, y_val, g, diff;
            mpz_inits(x_acc, y_sq, y_val, g, diff, NULL);
            mpz_set_ui(x_acc, 1);
            mpz_set_ui(y_sq, 1);

            for (int oi = 0; oi < mat.origin_count[r]; oi++) {
                int orig = mat.origins[r][oi];
                long x1, x2; /* the x value(s) */
                int is_combined = (orig >= n_smooth);

                if (!is_combined) {
                    x1 = smooth_rels[orig].x;
                    mpz_set_si(tmp, x1); mpz_add(tmp, tmp, sqrtN);
                    mpz_mul(x_acc, x_acc, tmp); mpz_mod(x_acc, x_acc, N);
                    mpz_set_si(tmp, x1); mpz_add(tmp, tmp, sqrtN);
                    mpz_mul(tmp, tmp, tmp); mpz_sub(tmp, tmp, N); mpz_abs(tmp, tmp);
                    mpz_mul(y_sq, y_sq, tmp);
                } else {
                    int ci = orig - n_smooth;
                    x1 = combined_rels[ci].x;
                    x2 = (long)combined_rels[ci].lp; /* stored second x here */

                    mpz_set_si(tmp, x1); mpz_add(tmp, tmp, sqrtN);
                    mpz_mul(x_acc, x_acc, tmp); mpz_mod(x_acc, x_acc, N);
                    mpz_set_si(tmp, x1); mpz_add(tmp, tmp, sqrtN);
                    mpz_mul(tmp, tmp, tmp); mpz_sub(tmp, tmp, N); mpz_abs(tmp, tmp);
                    mpz_mul(y_sq, y_sq, tmp);

                    mpz_set_si(tmp, x2); mpz_add(tmp, tmp, sqrtN);
                    mpz_mul(x_acc, x_acc, tmp); mpz_mod(x_acc, x_acc, N);
                    mpz_set_si(tmp, x2); mpz_add(tmp, tmp, sqrtN);
                    mpz_mul(tmp, tmp, tmp); mpz_sub(tmp, tmp, N); mpz_abs(tmp, tmp);
                    mpz_mul(y_sq, y_sq, tmp);
                }
            }

            if (mpz_perfect_square_p(y_sq)) {
                mpz_sqrt(y_val, y_sq);
                mpz_mod(y_val, y_val, N);

                mpz_sub(diff, x_acc, y_val);
                mpz_gcd(g, diff, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_divexact(tmp, N, g);
                    gmp_printf("%Zd %Zd\n", g, tmp);
                    fprintf(stderr, "[%.1fs] Factor found!\n", ELAPSED());
                    return 0;
                }
                mpz_add(diff, x_acc, y_val);
                mpz_gcd(g, diff, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_divexact(tmp, N, g);
                    gmp_printf("%Zd %Zd\n", g, tmp);
                    fprintf(stderr, "[%.1fs] Factor found!\n", ELAPSED());
                    return 0;
                }
            }
            mpz_clears(x_acc, y_sq, y_val, g, diff, NULL);
        } else {
            pivot_col[lead] = r;
        }
    }

    fprintf(stderr, "[%.1fs] FAILED: no factor from Gauss\n", ELAPSED());
    free(pivot_col);

    /* Cleanup */
    for (int i = 0; i < n_smooth; i++) free(smooth_rels[i].fb_exp);
    for (int i = 0; i < n_partial; i++) free(partial_rels[i].fb_exp);
    for (int i = 0; i < n_combined; i++) free(combined_rels[i].fb_exp);
    free(smooth_rels); free(partial_rels); free(combined_rels);
    for (int i = 0; i < total_rels; i++) free(mat.origins[i]);
    free(mat.origins); free(mat.origin_count); free(mat.data);
    free(fb); free(fb_root);
    mpz_clears(N, sqrtN, tmp, tmp2, qx, cofactor, base_val, NULL);
    return 1;
}
