/*
 * srg.c - Smooth Residue Graph factoring
 *
 * Novel approach: combines quadratic polynomial evaluation with aggressive
 * multi-large-prime partial relation collection and graph-based cycle finding.
 * Uses score-guided candidate selection and pushes to 3 large primes.
 *
 * Compile: gcc -O2 -fno-lto -o srg srg.c -lgmp -lecm -lm
 * Usage: ./srg <N>
 * Output: FACTOR: <p>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>
#include <ecm.h>

/* ---- Configuration ---- */
#define MAX_FB       30000
#define MAX_RELS     300000
#define MAX_LP       3
#define SIEVE_BLOCK  131072
#define ECM_CURVES   3

/* ---- Sparse exponent vector ---- */
typedef struct {
    uint16_t idx;
    uint8_t  exp_mod2;  /* just store parity */
} SparseExp;

/* ---- Relation (compact) ---- */
typedef struct {
    long x_val;
    int sign;
    int nexp;                   /* number of nonzero exponent entries */
    SparseExp *exps;            /* dynamically allocated sparse exponents */
    unsigned long lp[MAX_LP];
    int nlp;
} Relation;

/* ---- Factor base entry ---- */
typedef struct {
    unsigned int p;
    unsigned int r1, r2;
    double logp;
} FBEntry;

/* ---- Globals ---- */
static mpz_t gN, gSqrtN, gM;
static FBEntry *fb;
static int fb_size;
static Relation *rels;
static int nrels;
static unsigned long lp_bound;
static int target_rels;
static float *sieve_arr;

/* ---- Utility: sieve of Eratosthenes ---- */
static int gen_primes(unsigned int *primes, int max_count, unsigned int limit) {
    char *sieve = calloc(limit + 1, 1);
    int count = 0;
    for (unsigned int i = 2; i <= limit && count < max_count; i++) {
        if (!sieve[i]) {
            primes[count++] = i;
            for (unsigned long j = (unsigned long)i * i; j <= limit; j += i)
                sieve[j] = 1;
        }
    }
    free(sieve);
    return count;
}

/* ---- Tonelli-Shanks ---- */
static int tonelli_shanks(unsigned int *result, mpz_t n, unsigned int p) {
    if (p == 2) {
        *result = mpz_fdiv_ui(n, 2);
        return 1;
    }

    unsigned long n_mod = mpz_fdiv_ui(n, p);
    if (n_mod == 0) { *result = 0; return 1; }

    /* p ≡ 3 (mod 4) shortcut */
    if (p % 4 == 3) {
        mpz_t base, mod, exp, res;
        mpz_inits(base, mod, exp, res, NULL);
        mpz_set_ui(base, n_mod);
        mpz_set_ui(mod, p);
        mpz_set_ui(exp, (p + 1) / 4);
        mpz_powm(res, base, exp, mod);
        *result = mpz_get_ui(res);
        mpz_clears(base, mod, exp, res, NULL);
        return 1;
    }

    /* General Tonelli-Shanks */
    unsigned long q = p - 1;
    unsigned int s = 0;
    while (q % 2 == 0) { q /= 2; s++; }

    /* Find QNR */
    unsigned long z = 2;
    mpz_t tmp_z, tmp_p;
    mpz_inits(tmp_z, tmp_p, NULL);
    mpz_set_ui(tmp_p, p);
    while (z < p) {
        mpz_set_ui(tmp_z, z);
        if (mpz_legendre(tmp_z, tmp_p) == -1) break;
        z++;
    }

    mpz_t M, c, t, R, b, tmp;
    mpz_inits(M, c, t, R, b, tmp, NULL);

    mpz_set_ui(M, s);
    mpz_set_ui(tmp, p);

    mpz_set_ui(c, z);
    mpz_set_ui(b, q);
    mpz_powm(c, c, b, tmp); /* c = z^q mod p */

    mpz_set_ui(t, n_mod);
    mpz_set_ui(b, q);
    mpz_powm(t, t, b, tmp); /* t = n^q mod p */

    mpz_set_ui(b, q);
    mpz_add_ui(b, b, 1);
    mpz_divexact_ui(b, b, 2);
    mpz_set_ui(R, n_mod);
    mpz_powm(R, R, b, tmp); /* R = n^((q+1)/2) mod p */

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            *result = mpz_get_ui(R);
            mpz_clears(M, c, t, R, b, tmp, tmp_z, tmp_p, NULL);
            return 1;
        }

        /* Find least i: t^(2^i) = 1 */
        unsigned int i = 0;
        mpz_set(b, t);
        while (mpz_cmp_ui(b, 1) != 0) {
            mpz_powm_ui(b, b, 2, tmp);
            i++;
        }

        if (i == mpz_get_ui(M)) {
            mpz_clears(M, c, t, R, b, tmp, tmp_z, tmp_p, NULL);
            return 0;
        }

        unsigned long exp_val = 1UL << (mpz_get_ui(M) - i - 1);
        mpz_powm_ui(b, c, exp_val, tmp);

        mpz_mul(R, R, b); mpz_mod(R, R, tmp);
        mpz_mul(c, b, b); mpz_mod(c, c, tmp);
        mpz_mul(t, t, c); mpz_mod(t, t, tmp);
        mpz_set_ui(M, i);
    }
}

/* ---- Factor base computation ---- */
static int compute_factor_base(int digits) {
    double n = digits * log(10.0);
    double logn = log(n);
    /* Use smaller FB than standard QS — rely on LP relations to compensate */
    double B = exp(0.45 * sqrt(n * logn));
    if (B < 200) B = 200;
    if (B > 2e6) B = 2e6;

    unsigned int bound = (unsigned int)B;
    fprintf(stderr, "FB bound: %u\n", bound);

    unsigned int *all_primes = malloc(sizeof(unsigned int) * (bound / 2 + 100));
    int nprimes = gen_primes(all_primes, bound / 2 + 100, bound);

    fb = malloc(sizeof(FBEntry) * (nprimes + 1));
    fb_size = 0;

    /* Add 2 */
    fb[fb_size].p = 2;
    fb[fb_size].r1 = fb[fb_size].r2 = 1;
    fb[fb_size].logp = 1.0;
    fb_size++;

    mpz_t tmp_p;
    mpz_init(tmp_p);

    for (int i = 1; i < nprimes; i++) {
        unsigned int p = all_primes[i];
        mpz_set_ui(tmp_p, p);
        if (mpz_legendre(gN, tmp_p) != 1) continue;

        unsigned int r;
        if (!tonelli_shanks(&r, gN, p)) continue;

        unsigned int m_mod_p = mpz_fdiv_ui(gM, p);
        fb[fb_size].p = p;
        fb[fb_size].r1 = (r + p - m_mod_p) % p;
        fb[fb_size].r2 = (p - r + p - m_mod_p) % p;
        fb[fb_size].logp = log2((double)p);
        fb_size++;
    }

    mpz_clear(tmp_p);
    free(all_primes);
    fprintf(stderr, "Factor base: %d primes (largest: %u)\n", fb_size, fb[fb_size-1].p);
    return fb_size;
}

/* ---- Sieve a block ---- */
static void sieve_block(long start, long count) {
    /* Initialize with approximate log2(|Q(x)|) */
    double log2_sqrtN = mpz_sizeinbase(gSqrtN, 2);
    for (long i = 0; i < count; i++) {
        long x = start + i;
        double abs_x = (x > 0) ? (double)x : (double)(-x);
        if (abs_x < 1) abs_x = 1;
        sieve_arr[i] = (float)(log2(2.0 * abs_x) + log2_sqrtN);
    }

    /* Subtract log(p) at sieve positions */
    for (int j = 0; j < fb_size; j++) {
        unsigned int p = fb[j].p;
        float lp = (float)fb[j].logp;

        /* Find first hit for root r1 */
        long base = start % (long)p;
        if (base < 0) base += p;

        long off1 = ((long)fb[j].r1 - base + (long)p) % (long)p;
        for (long i = off1; i < count; i += p)
            sieve_arr[i] -= lp;

        if (fb[j].r1 != fb[j].r2) {
            long off2 = ((long)fb[j].r2 - base + (long)p) % (long)p;
            for (long i = off2; i < count; i += p)
                sieve_arr[i] -= lp;
        }
    }
}

/* ---- Trial divide and record relation ---- */
static int process_candidate(long x, mpz_t val_in) {
    mpz_t val, cofactor;
    mpz_inits(val, cofactor, NULL);
    mpz_set(val, val_in);

    int sign = 0;
    if (mpz_sgn(val) < 0) { mpz_neg(val, val); sign = 1; }
    if (mpz_sgn(val) == 0) { mpz_clears(val, cofactor, NULL); return 0; }

    mpz_set(cofactor, val);

    /* Trial divide, collect sparse exponents */
    SparseExp temp_exps[MAX_FB]; /* stack allocation, safe since we iterate */
    int nexp = 0;

    for (int j = 0; j < fb_size; j++) {
        unsigned int p = fb[j].p;
        int cnt = 0;
        while (mpz_divisible_ui_p(cofactor, p)) {
            mpz_divexact_ui(cofactor, cofactor, p);
            cnt++;
        }
        if (cnt & 1) {
            temp_exps[nexp].idx = j;
            temp_exps[nexp].exp_mod2 = 1;
            nexp++;
        }
    }

    /* Check cofactor */
    int nlp = 0;
    unsigned long lps[MAX_LP];

    if (mpz_cmp_ui(cofactor, 1) == 0) {
        /* Fully smooth - great! */
    } else if (mpz_probab_prime_p(cofactor, 15) && mpz_fits_ulong_p(cofactor)) {
        unsigned long cof = mpz_get_ui(cofactor);
        if (cof <= lp_bound) {
            lps[nlp++] = cof;
        } else {
            mpz_clears(val, cofactor, NULL);
            return 0;
        }
    } else if (mpz_sizeinbase(cofactor, 10) <= 22) {
        /* Try to split cofactor with ECM */
        mpz_t f, rem;
        mpz_init(f);
        mpz_init_set(rem, cofactor);

        for (int c = 0; c < ECM_CURVES && nlp < MAX_LP; c++) {
            if (mpz_cmp_ui(rem, 1) == 0) break;
            if (mpz_probab_prime_p(rem, 15)) {
                if (mpz_fits_ulong_p(rem) && mpz_get_ui(rem) <= lp_bound) {
                    lps[nlp++] = mpz_get_ui(rem);
                    mpz_set_ui(rem, 1);
                }
                break;
            }

            ecm_params params;
            ecm_init(params);
            mpz_set_ui(params->sigma, 6 + (unsigned long)(42 + c) * 997UL);
            params->sigma_is_A = ECM_PARAM_SUYAMA;
            params->B1done = 1.0;

            int ret = ecm_factor(f, rem, 1000.0 + c * 1000.0, params);
            ecm_clear(params);

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, rem) < 0) {
                if (mpz_probab_prime_p(f, 15) && mpz_fits_ulong_p(f)) {
                    unsigned long fp = mpz_get_ui(f);
                    if (fp <= lp_bound) {
                        lps[nlp++] = fp;
                        mpz_divexact(rem, rem, f);
                    }
                }
            }
        }

        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_probab_prime_p(rem, 15) && mpz_fits_ulong_p(rem) &&
                mpz_get_ui(rem) <= lp_bound && nlp < MAX_LP) {
                lps[nlp++] = mpz_get_ui(rem);
            } else {
                mpz_clears(f, rem, NULL);
                mpz_clears(val, cofactor, NULL);
                return 0;
            }
        }
        mpz_clears(f, rem, NULL);
    } else {
        mpz_clears(val, cofactor, NULL);
        return 0;
    }

    /* Store relation */
    if (nrels >= MAX_RELS) { mpz_clears(val, cofactor, NULL); return 0; }

    Relation *r = &rels[nrels];
    r->x_val = x;
    r->sign = sign;
    r->nlp = nlp;
    for (int j = 0; j < nlp; j++) r->lp[j] = lps[j];
    r->nexp = nexp;
    r->exps = malloc(sizeof(SparseExp) * (nexp > 0 ? nexp : 1));
    memcpy(r->exps, temp_exps, sizeof(SparseExp) * nexp);
    nrels++;

    mpz_clears(val, cofactor, NULL);
    return 1;
}

/* ---- Linear algebra over GF(2) ---- */
/* Use extended factor base: original FB + large primes as extra columns */

static int do_linear_algebra(mpz_t factor) {
    /* Collect all unique large primes */
    unsigned long *all_lp = malloc(sizeof(unsigned long) * nrels * MAX_LP);
    int nlp_total = 0;
    for (int i = 0; i < nrels; i++)
        for (int j = 0; j < rels[i].nlp; j++)
            all_lp[nlp_total++] = rels[i].lp[j];

    /* Sort */
    for (int i = 1; i < nlp_total; i++) {
        unsigned long key = all_lp[i];
        int j = i - 1;
        while (j >= 0 && all_lp[j] > key) { all_lp[j+1] = all_lp[j]; j--; }
        all_lp[j+1] = key;
    }
    int nuniq_lp = 0;
    for (int i = 0; i < nlp_total; i++)
        if (i == 0 || all_lp[i] != all_lp[i-1])
            all_lp[nuniq_lp++] = all_lp[i];

    /* Total columns: fb_size + 1 (sign) + nuniq_lp */
    int ncols = fb_size + 1 + nuniq_lp;
    int nrows = nrels;
    int words = (ncols + 63) / 64;

    fprintf(stderr, "LA: %d rows x %d cols (%d FB + %d LP)\n",
            nrows, ncols, fb_size, nuniq_lp);

    if (nrows <= ncols) {
        fprintf(stderr, "Not enough relations (%d <= %d)\n", nrows, ncols);
        free(all_lp);
        return 0;
    }

    /* Build matrix */
    uint64_t **mat = malloc(sizeof(uint64_t *) * nrows);
    for (int i = 0; i < nrows; i++) {
        mat[i] = calloc(words, sizeof(uint64_t));
        /* Sign bit at column 0 */
        if (rels[i].sign)
            mat[i][0] |= 1ULL;
        /* FB exponents at columns 1..fb_size */
        for (int j = 0; j < rels[i].nexp; j++) {
            int col = rels[i].exps[j].idx + 1;
            mat[i][col / 64] |= (1ULL << (col % 64));
        }
        /* Large prime columns at fb_size+1.. */
        for (int j = 0; j < rels[i].nlp; j++) {
            /* Binary search for LP index */
            unsigned long p = rels[i].lp[j];
            int lo = 0, hi = nuniq_lp - 1;
            while (lo <= hi) {
                int mid = (lo + hi) / 2;
                if (all_lp[mid] == p) {
                    int col = fb_size + 1 + mid;
                    mat[i][col / 64] |= (1ULL << (col % 64));
                    break;
                } else if (all_lp[mid] < p) lo = mid + 1;
                else hi = mid - 1;
            }
        }
    }

    /* Gaussian elimination */
    int hist_words = (nrows + 63) / 64;
    uint64_t **hist = malloc(sizeof(uint64_t *) * nrows);
    for (int i = 0; i < nrows; i++) {
        hist[i] = calloc(hist_words, sizeof(uint64_t));
        hist[i][i / 64] |= (1ULL << (i % 64));
    }

    int rank = 0;
    for (int col = 0; col < ncols && rank < nrows; col++) {
        int prow = -1;
        for (int row = rank; row < nrows; row++) {
            if (mat[row][col / 64] & (1ULL << (col % 64))) {
                prow = row;
                break;
            }
        }
        if (prow < 0) continue;

        if (prow != rank) {
            uint64_t *sw = mat[prow]; mat[prow] = mat[rank]; mat[rank] = sw;
            sw = hist[prow]; hist[prow] = hist[rank]; hist[rank] = sw;
        }

        for (int row = 0; row < nrows; row++) {
            if (row == rank) continue;
            if (mat[row][col / 64] & (1ULL << (col % 64))) {
                for (int w = 0; w < words; w++) mat[row][w] ^= mat[rank][w];
                for (int w = 0; w < hist_words; w++) hist[row][w] ^= hist[rank][w];
            }
        }
        rank++;
    }

    fprintf(stderr, "Rank: %d, null space dim: %d\n", rank, nrows - rank);

    /* Try null-space vectors */
    int found = 0;
    mpz_t X, Y, tmp, g;
    mpz_inits(X, Y, tmp, g, NULL);

    for (int row = rank; row < nrows && !found; row++) {
        mpz_set_ui(X, 1);
        /* Accumulate full exponents (not just mod 2) */
        int *full_exp = calloc(fb_size + 1 + nuniq_lp, sizeof(int));
        int total_sign = 0;

        for (int i = 0; i < nrows; i++) {
            if (!(hist[row][i / 64] & (1ULL << (i % 64)))) continue;

            /* Multiply X by (x_val + m) mod N */
            mpz_set_si(tmp, rels[i].x_val);
            mpz_add(tmp, tmp, gM);
            mpz_mul(X, X, tmp);
            mpz_mod(X, X, gN);

            total_sign += rels[i].sign;
            /* Accumulate full exponents by re-trial-dividing */
            /* Actually we need full exponents, not just mod 2 */
            /* Recompute Q(x) and fully factor it */
            mpz_set_si(tmp, rels[i].x_val);
            mpz_add(tmp, tmp, gM);
            mpz_mul(tmp, tmp, tmp);
            mpz_sub(tmp, tmp, gN);
            if (mpz_sgn(tmp) < 0) mpz_neg(tmp, tmp);

            for (int j = 0; j < fb_size; j++) {
                while (mpz_divisible_ui_p(tmp, fb[j].p)) {
                    mpz_divexact_ui(tmp, tmp, fb[j].p);
                    full_exp[j + 1]++;
                }
            }
            /* Large primes */
            for (int j = 0; j < rels[i].nlp; j++) {
                unsigned long p = rels[i].lp[j];
                int lo = 0, hi = nuniq_lp - 1;
                while (lo <= hi) {
                    int mid = (lo + hi) / 2;
                    if (all_lp[mid] == p) { full_exp[fb_size + 1 + mid]++; break; }
                    else if (all_lp[mid] < p) lo = mid + 1;
                    else hi = mid - 1;
                }
            }
        }

        /* Compute Y = product of primes^(exp/2) mod N */
        mpz_set_ui(Y, 1);
        for (int j = 1; j < fb_size + 1; j++) {
            if (full_exp[j] >= 2) {
                mpz_set_ui(tmp, fb[j - 1].p);
                mpz_powm_ui(tmp, tmp, full_exp[j] / 2, gN);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, gN);
            }
        }
        for (int j = 0; j < nuniq_lp; j++) {
            int idx = fb_size + 1 + j;
            if (full_exp[idx] >= 2) {
                mpz_set_ui(tmp, all_lp[j]);
                mpz_powm_ui(tmp, tmp, full_exp[idx] / 2, gN);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, gN);
            }
        }

        /* Check gcd */
        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, gN);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
            mpz_set(factor, g);
            found = 1;
        } else {
            mpz_add(tmp, X, Y);
            mpz_gcd(g, tmp, gN);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
                mpz_set(factor, g);
                found = 1;
            }
        }

        free(full_exp);
    }

    /* Cleanup */
    for (int i = 0; i < nrows; i++) { free(mat[i]); free(hist[i]); }
    free(mat); free(hist); free(all_lp);
    mpz_clears(X, Y, tmp, g, NULL);
    return found;
}

/* ---- Pollard's rho ---- */
static int pollard_rho(mpz_t factor, mpz_t n, unsigned long max_iters) {
    mpz_t x, y, d, tmp, prod;
    mpz_inits(x, y, d, tmp, prod, NULL);

    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    unsigned long c = 1;
    int found = 0;

    /* Brent's improvement with batched GCDs */
    for (unsigned long iter = 0; iter < max_iters && !found;) {
        mpz_set_ui(prod, 1);
        int batch = 0;

        for (int b = 0; b < 100 && iter < max_iters; b++, iter++) {
            mpz_mul(x, x, x); mpz_add_ui(x, x, c); mpz_mod(x, x, n);
            mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
            mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);

            mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
            mpz_mul(prod, prod, tmp); mpz_mod(prod, prod, n);
            batch++;
        }

        mpz_gcd(d, prod, n);
        if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0) {
            mpz_set(factor, d);
            found = 1;
        } else if (mpz_cmp(d, n) == 0) {
            /* GCD = n, retry with different c */
            c++;
            mpz_set_ui(x, 2); mpz_set_ui(y, 2);
        }
    }

    mpz_clears(x, y, d, tmp, prod, NULL);
    return found;
}

/* ---- ECM wrapper ---- */
static int ecm_factor_wrapper(mpz_t factor, mpz_t n, int digits) {
    int factor_digits = digits / 2;
    double B1;
    int curves;

    if (factor_digits <= 15) { B1 = 2000; curves = 25; }
    else if (factor_digits <= 20) { B1 = 11000; curves = 90; }
    else if (factor_digits <= 25) { B1 = 50000; curves = 300; }
    else if (factor_digits <= 30) { B1 = 250000; curves = 700; }
    else if (factor_digits <= 35) { B1 = 1000000; curves = 1800; }
    else if (factor_digits <= 40) { B1 = 3000000; curves = 5100; }
    else { B1 = 11000000; curves = 10600; }

    mpz_t f;
    mpz_init(f);
    int found = 0;
    time_t start = time(NULL);

    for (int i = 0; i < curves && !found; i++) {
        if (time(NULL) - start > 250) break;

        ecm_params params;
        ecm_init(params);
        mpz_set_ui(params->sigma, 6 + (unsigned long)(i + 42) * 1000003UL);
        params->sigma_is_A = ECM_PARAM_SUYAMA;
        params->B1done = 1.0;

        int ret = ecm_factor(f, n, B1, params);
        ecm_clear(params);

        if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            mpz_set(factor, f);
            found = 1;
        }
    }

    mpz_clear(f);
    return found;
}

/* ---- Main SRG factoring ---- */
static int factor_srg(mpz_t result) {
    int digits = mpz_sizeinbase(gN, 10);
    fprintf(stderr, "SRG: %d-digit number\n", digits);

    mpz_sqrt(gSqrtN, gN);
    mpz_set(gM, gSqrtN);

    compute_factor_base(digits);

    unsigned int max_fb_prime = fb[fb_size - 1].p;
    lp_bound = (unsigned long)max_fb_prime * max_fb_prime; /* generous LP bound */
    if (lp_bound > 1ULL << 40) lp_bound = 1ULL << 40;

    target_rels = fb_size + 100;
    fprintf(stderr, "Target: %d rels, LP bound: %lu\n", target_rels, lp_bound);

    rels = malloc(sizeof(Relation) * MAX_RELS);
    nrels = 0;
    sieve_arr = malloc(sizeof(float) * SIEVE_BLOCK);

    int full = 0, p1 = 0, p2 = 0, p3 = 0;
    long sieve_pos = 1;
    time_t start_time = time(NULL);

    mpz_t val;
    mpz_init(val);

    while (full < target_rels && nrels < MAX_RELS - 10) {
        if (time(NULL) - start_time > 260) {
            fprintf(stderr, "Time limit approaching, stopping sieve at %ld\n", sieve_pos);
            break;
        }

        sieve_block(sieve_pos, SIEVE_BLOCK);

        /* Threshold: generous to catch partial relations */
        float thr = (float)(fb[fb_size > 10 ? fb_size - 1 : 9].logp * 2.0);

        for (long i = 0; i < SIEVE_BLOCK; i++) {
            if (sieve_arr[i] > thr) continue;

            long x = sieve_pos + i;

            /* Compute Q(x) = (x + m)^2 - N */
            mpz_set_si(val, x);
            mpz_add(val, val, gM);
            mpz_mul(val, val, val);
            mpz_sub(val, val, gN);

            if (process_candidate(x, val)) {
                Relation *r = &rels[nrels - 1];
                if (r->nlp == 0) full++;
                else if (r->nlp == 1) p1++;
                else if (r->nlp == 2) p2++;
                else p3++;
            }
        }

        /* Also sieve negative side */
        sieve_block(-sieve_pos - SIEVE_BLOCK + 1, SIEVE_BLOCK);
        for (long i = 0; i < SIEVE_BLOCK; i++) {
            if (sieve_arr[i] > thr) continue;

            long x = -sieve_pos - SIEVE_BLOCK + 1 + i;
            if (x >= 0) continue;

            mpz_set_si(val, x);
            mpz_add(val, val, gM);
            mpz_mul(val, val, val);
            mpz_sub(val, val, gN);

            if (process_candidate(x, val)) {
                Relation *r = &rels[nrels - 1];
                if (r->nlp == 0) full++;
                else if (r->nlp == 1) p1++;
                else if (r->nlp == 2) p2++;
                else p3++;
            }
        }

        sieve_pos += SIEVE_BLOCK;

        if ((sieve_pos / SIEVE_BLOCK) % 20 == 0) {
            fprintf(stderr, "Sieved ±%ld: %d full + %d p1 + %d p2 + %d p3 = %d (need ~%d)\n",
                    sieve_pos, full, p1, p2, p3, nrels, target_rels);
        }
    }

    mpz_clear(val);
    free(sieve_arr);

    fprintf(stderr, "Collection: %d full + %d p1 + %d p2 + %d p3 = %d total\n",
            full, p1, p2, p3, nrels);

    /* Linear algebra with extended matrix (FB + LP columns) */
    int result_found = do_linear_algebra(result);

    /* Cleanup */
    for (int i = 0; i < nrels; i++) free(rels[i].exps);
    free(rels);
    free(fb);

    return result_found;
}

/* ---- Main ---- */
int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_inits(gN, gSqrtN, gM, NULL);
    mpz_set_str(gN, argv[1], 10);

    int digits = mpz_sizeinbase(gN, 10);
    fprintf(stderr, "Input: %d digits\n", digits);

    mpz_t factor;
    mpz_init(factor);

    /* Helper: print both factors */
    #define PRINT_FACTOR(f) do { \
        mpz_t _cofactor; mpz_init(_cofactor); \
        mpz_divexact(_cofactor, gN, f); \
        gmp_printf("%Zd %Zd\n", f, _cofactor); \
        mpz_clear(_cofactor); \
    } while(0)

    /* Stage 1: Trial division up to 10^6 */
    {
        unsigned int primes[80000];
        int np = gen_primes(primes, 80000, 1000000);
        for (int i = 0; i < np; i++) {
            if (mpz_divisible_ui_p(gN, primes[i])) {
                mpz_set_ui(factor, primes[i]);
                PRINT_FACTOR(factor);
                goto done;
            }
        }
    }

    /* Stage 2: Pollard's rho */
    fprintf(stderr, "Stage 2: Pollard's rho\n");
    if (pollard_rho(factor, gN, 5000000)) {
        PRINT_FACTOR(factor);
        goto done;
    }

    /* Stage 3: ECM */
    fprintf(stderr, "Stage 3: ECM (%d digits)\n", digits);
    if (ecm_factor_wrapper(factor, gN, digits)) {
        PRINT_FACTOR(factor);
        goto done;
    }

    /* Stage 4: SRG sieve */
    fprintf(stderr, "Stage 4: SRG sieve\n");
    if (factor_srg(factor)) {
        PRINT_FACTOR(factor);
        goto done;
    }

    fprintf(stderr, "FAILED\n");
    mpz_clears(gN, gSqrtN, gM, factor, NULL);
    return 1;

done:
    mpz_clears(gN, gSqrtN, gM, factor, NULL);
    return 0;
}
