/*
 * srg.c - Smooth Residue Graph factoring
 *
 * Novel approach: combines quadratic polynomial evaluation with aggressive
 * multi-large-prime partial relation collection and graph-based cycle finding.
 * Uses SIQS-style multiple polynomials for smaller Q(x) values.
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
#define SIEVE_HALF   65536
#define SIEVE_BLOCK  (2 * SIEVE_HALF)
#define ECM_CURVES   3
#define MAX_APRIME   3      /* max number of factor base primes in 'a' */

/* ---- Sparse exponent vector ---- */
typedef struct {
    uint16_t idx;
    uint8_t  exp_mod2;  /* just store parity */
} SparseExp;

/* ---- Relation (compact) ---- */
typedef struct {
    mpz_t ax_plus_b;            /* a*x + b for this polynomial */
    int sign;
    int nexp;                   /* number of nonzero exponent entries */
    SparseExp *exps;            /* dynamically allocated sparse exponents */
    unsigned long lp[MAX_LP];
    int nlp;
    /* indices of factor base primes composing 'a' and their parity */
    int a_nfactors;
    uint16_t a_fb_idx[MAX_APRIME];
    uint8_t  a_exp_mod2[MAX_APRIME]; /* parity of exponent of a-prime in Q/a */
} Relation;

/* ---- Factor base entry ---- */
typedef struct {
    unsigned int p;
    unsigned int sqrt_n_mod_p;  /* sqrt(N) mod p, the raw Tonelli-Shanks root */
    unsigned int r1, r2;        /* sieve roots for current polynomial */
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

    /* p == 3 (mod 4) shortcut */
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

    mpz_t Mv, c, t, R, b, tmp;
    mpz_inits(Mv, c, t, R, b, tmp, NULL);

    mpz_set_ui(Mv, s);
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
            mpz_clears(Mv, c, t, R, b, tmp, tmp_z, tmp_p, NULL);
            return 1;
        }

        /* Find least i: t^(2^i) = 1 */
        unsigned int i = 0;
        mpz_set(b, t);
        while (mpz_cmp_ui(b, 1) != 0) {
            mpz_powm_ui(b, b, 2, tmp);
            i++;
        }

        if (i == mpz_get_ui(Mv)) {
            mpz_clears(Mv, c, t, R, b, tmp, tmp_z, tmp_p, NULL);
            return 0;
        }

        unsigned long exp_val = 1UL << (mpz_get_ui(Mv) - i - 1);
        mpz_powm_ui(b, c, exp_val, tmp);

        mpz_mul(R, R, b); mpz_mod(R, R, tmp);
        mpz_mul(c, b, b); mpz_mod(c, c, tmp);
        mpz_mul(t, t, c); mpz_mod(t, t, tmp);
        mpz_set_ui(Mv, i);
    }
}

/* ---- Factor base computation ---- */
static int compute_factor_base(int digits) {
    double n = digits * log(10.0);
    double logn = log(n);
    /* Use smaller FB than standard QS -- rely on LP relations to compensate */
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
    fb[fb_size].sqrt_n_mod_p = 1;
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

        fb[fb_size].p = p;
        fb[fb_size].sqrt_n_mod_p = r;  /* store raw root for MPQS use */
        /* r1/r2 will be computed per-polynomial */
        fb[fb_size].r1 = 0;
        fb[fb_size].r2 = 0;
        fb[fb_size].logp = log2((double)p);
        fb_size++;
    }

    mpz_clear(tmp_p);
    free(all_primes);
    fprintf(stderr, "Factor base: %d primes (largest: %u)\n", fb_size, fb[fb_size-1].p);
    return fb_size;
}

/* ---- Sieve a block for multi-polynomial ---- */
/* Sieve interval [-M, M] given precomputed roots in fb[].r1, fb[].r2 */
static void sieve_block_mpoly(float *sarr, long count, double log2_approx_size) {
    /* Initialize sieve with approximate log2(|Q_a(x)/a|) */
    for (long i = 0; i < count; i++) {
        long x = i - SIEVE_HALF;
        double abs_x = (x > 0) ? (double)x : (double)(-x);
        if (abs_x < 1) abs_x = 1;
        /* |Q_a(x)/a| ~ (2*|x|*sqrt(N)/a + ...) but we approximate */
        sarr[i] = (float)(log2(abs_x + 1.0) + log2_approx_size);
    }

    /* Subtract log(p) at sieve positions */
    for (int j = 0; j < fb_size; j++) {
        unsigned int p = fb[j].p;
        float lp = (float)fb[j].logp;
        unsigned int r1 = fb[j].r1;
        unsigned int r2 = fb[j].r2;

        /* Sieve interval is [0, count) representing x in [-M, M) */
        /* We need positions where (a*x + b)^2 - N == 0 mod p */
        /* These are stored as offsets from position 0 in r1, r2 */

        if (r1 < (unsigned int)count) {
            for (long i = r1; i < count; i += p)
                sarr[i] -= lp;
        }

        if (r1 != r2 && r2 < (unsigned int)count) {
            for (long i = r2; i < count; i += p)
                sarr[i] -= lp;
        }
    }
}

/* ---- Trial divide and record relation for MPQS ---- */
/* val_in = Q_a(x)/a, ax_b = a*x + b, a_fb_idx/a_nfactors describe the 'a' primes */
static int process_candidate_mpqs(mpz_t val_in, mpz_t ax_b,
                                   int a_nfactors, uint16_t *a_fb_idx) {
    mpz_t val, cofactor;
    mpz_inits(val, cofactor, NULL);
    mpz_set(val, val_in);

    int sign = 0;
    if (mpz_sgn(val) < 0) { mpz_neg(val, val); sign = 1; }
    if (mpz_sgn(val) == 0) { mpz_clears(val, cofactor, NULL); return 0; }

    mpz_set(cofactor, val);

    /* Trial divide, collect sparse exponents */
    SparseExp temp_exps[MAX_FB];
    int nexp = 0;

    /* Track exponents for the 'a' primes separately (they appear in Q/a with
       one fewer power than in Q, since we divided by a) */
    int a_prime_parity[MAX_APRIME];
    memset(a_prime_parity, 0, sizeof(a_prime_parity));

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
        /* Check if this is one of the 'a' primes */
        for (int k = 0; k < a_nfactors; k++) {
            if (a_fb_idx[k] == j) {
                a_prime_parity[k] = cnt & 1;
            }
        }
    }

    /* Check cofactor */
    int nlp = 0;
    unsigned long lps[MAX_LP];

    if (mpz_cmp_ui(cofactor, 1) == 0) {
        /* Fully smooth */
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
    mpz_init_set(r->ax_plus_b, ax_b);
    r->sign = sign;
    r->nlp = nlp;
    for (int j = 0; j < nlp; j++) r->lp[j] = lps[j];
    r->nexp = nexp;
    r->exps = malloc(sizeof(SparseExp) * (nexp > 0 ? nexp : 1));
    memcpy(r->exps, temp_exps, sizeof(SparseExp) * nexp);
    r->a_nfactors = a_nfactors;
    for (int k = 0; k < a_nfactors; k++) {
        r->a_fb_idx[k] = a_fb_idx[k];
        r->a_exp_mod2[k] = a_prime_parity[k];
    }
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
    mpz_t X, Y, tmp, g, qa_div_a;
    mpz_inits(X, Y, tmp, g, qa_div_a, NULL);

    for (int row = rank; row < nrows && !found; row++) {
        mpz_set_ui(X, 1);
        /* Accumulate full exponents (not just mod 2) */
        int *full_exp = calloc(fb_size + 1 + nuniq_lp, sizeof(int));

        for (int i = 0; i < nrows; i++) {
            if (!(hist[row][i / 64] & (1ULL << (i % 64)))) continue;

            /* Multiply X by (a*x + b) mod N */
            mpz_mul(X, X, rels[i].ax_plus_b);
            mpz_mod(X, X, gN);

            /* Recompute Q_a(x)/a = (a*x+b)^2/a - N/a ... actually:
               (a*x+b)^2 - N = a * (Q_a(x)/a)
               So Q_a(x)/a = ((a*x+b)^2 - N) / a
               We need to fully factor this to get exponents. */
            mpz_mul(tmp, rels[i].ax_plus_b, rels[i].ax_plus_b);
            mpz_sub(tmp, tmp, gN);
            /* tmp = (a*x+b)^2 - N, which should be divisible by a */
            /* Divide by a (product of the a-primes) */
            for (int k = 0; k < rels[i].a_nfactors; k++) {
                unsigned int ap = fb[rels[i].a_fb_idx[k]].p;
                mpz_divexact_ui(tmp, tmp, ap);
            }
            /* tmp is now Q_a(x)/a, possibly negative */
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

        /* Compute Y = product of primes^(exp/2) mod N
         * The relation is: product((a_i*x_i + b_i))^2 = product(a_i) * product(Q_a(x_i)/a_i) mod N
         * The LHS is X^2, the RHS must be a perfect square.
         * Actually: (a*x+b)^2 = a * (Q/a) + N, so (a*x+b)^2 == a*(Q/a) mod N
         * Product of (a*x+b)^2 = product(a) * product(Q/a) mod N
         * So X^2 = product(a) * product(Q/a) mod N
         * We need product(a) * product(Q/a) to be a perfect square.
         * Since a = q1*q2*... (each prime appearing once), each a-prime appears
         * once per relation in which it's used.
         *
         * The exponent vector we stored (for GF(2) purposes) accounts for the
         * Q/a factorization plus the a-primes (handled separately).
         * But for the full exponents we need to also count the a-prime contributions.
         */
        /* Add a-prime contributions: each relation contributes 1 to each a-prime */
        for (int i = 0; i < nrows; i++) {
            if (!(hist[row][i / 64] & (1ULL << (i % 64)))) continue;
            for (int k = 0; k < rels[i].a_nfactors; k++) {
                full_exp[rels[i].a_fb_idx[k] + 1]++;
            }
        }

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
    mpz_clears(X, Y, tmp, g, qa_div_a, NULL);
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
            c++;
            mpz_set_ui(x, 2); mpz_set_ui(y, 2);
        }
    }

    mpz_clears(x, y, d, tmp, prod, NULL);
    return found;
}

/* ---- ECM wrapper with progressive bounds ---- */
static int ecm_factor_wrapper(mpz_t factor, mpz_t n, int digits) {
    static const struct { double B1; int curves; } levels[] = {
        {     2000,   25 },
        {    11000,   90 },
        {    50000,  300 },
        {   250000,  700 },
        {  1000000, 1800 },
        {  3000000, 5100 },
        { 11000000, 10600 },
        { 43000000, 19300 },
    };
    int nlevels = sizeof(levels) / sizeof(levels[0]);

    int factor_digits = digits / 2;
    int start_level = 0;
    if (factor_digits > 15) start_level = 1;
    if (factor_digits > 20) start_level = 2;
    if (factor_digits > 25) start_level = 3;
    if (factor_digits > 30) start_level = 4;

    mpz_t f;
    mpz_init(f);
    int found = 0;
    time_t start = time(NULL);
    int total_curves = 0;

    for (int level = start_level; level < nlevels && !found; level++) {
        double B1 = levels[level].B1;
        int ncurves = levels[level].curves;

        for (int i = 0; i < ncurves && !found; i++) {
            /* Budget: for 65+ digit numbers, cap ECM at 120s for sieve time */
            int ecm_budget = (digits > 64) ? 120 : 270;
            if (time(NULL) - start > ecm_budget) goto done;

            ecm_params params;
            ecm_init(params);
            mpz_set_ui(params->sigma, 6 + (unsigned long)(total_curves + 42) * 1000003UL);
            params->sigma_is_A = ECM_PARAM_SUYAMA;
            params->B1done = 1.0;

            int ret = ecm_factor(f, n, B1, params);
            ecm_clear(params);
            total_curves++;

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
                mpz_set(factor, f);
                found = 1;
            }
        }

        if (!found && level < nlevels - 1) {
            fprintf(stderr, "ECM level %d (B1=%.0f, %d curves) done, trying next...\n",
                    level, B1, ncurves);
        }
    }

done:
    mpz_clear(f);
    return found;
}

/* ---- Modular inverse mod p (unsigned int) ---- */
static unsigned int modinv_ui(unsigned int a, unsigned int p) {
    long g0 = (long)p, g1 = (long)a;
    long u0 = 0, u1 = 1;
    while (g1 != 0) {
        long q = g0 / g1;
        long t = g0 - q * g1; g0 = g1; g1 = t;
        t = u0 - q * u1; u0 = u1; u1 = t;
    }
    if (u0 < 0) u0 += (long)p;
    return (unsigned int)u0;
}

/* ---- Main SRG factoring with MPQS ---- */
static int factor_srg(mpz_t result) {
    int digits = mpz_sizeinbase(gN, 10);
    fprintf(stderr, "SRG-MPQS: %d-digit number\n", digits);

    mpz_sqrt(gSqrtN, gN);
    mpz_set(gM, gSqrtN);

    compute_factor_base(digits);

    unsigned int max_fb_prime = fb[fb_size - 1].p;
    unsigned long lp_mult = 50;
    if (digits > 60) lp_mult = 20;
    if (digits > 70) lp_mult = 10;
    lp_bound = (unsigned long)max_fb_prime * lp_mult;
    if (lp_bound > 1ULL << 40) lp_bound = 1ULL << 40;

    target_rels = fb_size + 100;
    fprintf(stderr, "Target: %d rels, LP bound: %lu\n", target_rels, lp_bound);

    rels = malloc(sizeof(Relation) * MAX_RELS);
    nrels = 0;
    sieve_arr = malloc(sizeof(float) * SIEVE_BLOCK);

    int full = 0, p1 = 0, p2 = 0, p3 = 0;
    time_t start_time = time(NULL);

    /* Compute optimal 'a' size:
     * We want |Q_a(x)/a| ~ x * 2*sqrt(N)/a for x in [-M,M]
     * At x=M: |Q_a(x)/a| ~ M * 2*sqrt(N)/a
     * For optimal smoothness: a ~ sqrt(2*sqrt(N)/M)
     * This makes max |Q_a(x)/a| ~ 2*M*sqrt(2*sqrt(N)/M) = 2*sqrt(2*M*sqrt(N)) */
    mpz_t target_a, a_val, b_val, tmp, tmp2, ainv;
    mpz_inits(target_a, a_val, b_val, tmp, tmp2, ainv, NULL);

    /* target_a = floor(sqrt(2 * sqrt(N) / M)) */
    mpz_mul_ui(tmp, gSqrtN, 2);
    mpz_tdiv_q_ui(tmp, tmp, SIEVE_HALF);
    mpz_sqrt(target_a, tmp);

    char ta_str[256];
    mpz_get_str(ta_str, 10, target_a);
    fprintf(stderr, "Target a ~ %s (%zu digits)\n", ta_str,
            mpz_sizeinbase(target_a, 10));

    /* Determine how many factor base primes to multiply to reach target_a.
     * We want primes near (target_a)^(1/k) for k primes.
     * Try k=2 or k=3. */
    int a_nprimes = 2;
    double log_target_a = mpz_sizeinbase(target_a, 2) * log(2.0);
    double ideal_log_per_prime = log_target_a / 2.0;
    if (ideal_log_per_prime > log((double)max_fb_prime)) {
        a_nprimes = 3;
        ideal_log_per_prime = log_target_a / 3.0;
    }
    if (a_nprimes == 2 && ideal_log_per_prime < log(3.0)) {
        a_nprimes = 1;
        ideal_log_per_prime = log_target_a;
    }

    /* Find the FB index range near the ideal prime size */
    double ideal_prime_size = exp(ideal_log_per_prime);
    fprintf(stderr, "Using %d primes per 'a', ideal prime ~ %.0f\n",
            a_nprimes, ideal_prime_size);

    /* Find center index in FB near ideal_prime_size */
    int center_idx = 5; /* skip tiny primes */
    for (int j = 5; j < fb_size; j++) {
        if (fb[j].p >= (unsigned int)ideal_prime_size) {
            center_idx = j;
            break;
        }
        center_idx = j;
    }

    /* We'll pick primes from a window around center_idx */
    int window_lo = center_idx - 30;
    int window_hi = center_idx + 30;
    if (window_lo < 3) window_lo = 3;  /* skip 2 and very small primes */
    if (window_hi >= fb_size) window_hi = fb_size - 1;

    fprintf(stderr, "A-prime window: FB[%d..%d] = [%u..%u]\n",
            window_lo, window_hi, fb[window_lo].p, fb[window_hi].p);

    /* Simple deterministic polynomial generation using seed 42 */
    unsigned long rng_state = 42;
    #define SIMPLE_RNG() (rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL)

    mpz_t qa_val, ax_b;
    mpz_inits(qa_val, ax_b, NULL);

    int poly_count = 0;
    int max_polys = 1000000;  /* safety limit */

    /* Pre-allocate b solutions for CRT */
    mpz_t b_solutions[MAX_APRIME * 2];  /* up to 2 solutions per prime */
    for (int i = 0; i < MAX_APRIME * 2; i++)
        mpz_init(b_solutions[i]);

    while (nrels < target_rels && nrels < MAX_RELS - 10 && poly_count < max_polys) {
        if (time(NULL) - start_time > 260) {
            fprintf(stderr, "Time limit approaching, stopping sieve\n");
            break;
        }

        /* Select a_nprimes factor base indices for this polynomial */
        uint16_t a_idx[MAX_APRIME];
        int ok = 1;

        if (a_nprimes == 1) {
            SIMPLE_RNG();
            a_idx[0] = window_lo + (int)(rng_state % (unsigned long)(window_hi - window_lo + 1));
        } else if (a_nprimes == 2) {
            for (int attempt = 0; attempt < 20; attempt++) {
                SIMPLE_RNG();
                a_idx[0] = window_lo + (int)(rng_state % (unsigned long)(window_hi - window_lo + 1));
                SIMPLE_RNG();
                a_idx[1] = window_lo + (int)(rng_state % (unsigned long)(window_hi - window_lo + 1));
                if (a_idx[0] != a_idx[1]) break;
                if (attempt == 19) ok = 0;
            }
            /* Sort indices */
            if (a_idx[0] > a_idx[1]) { uint16_t t = a_idx[0]; a_idx[0] = a_idx[1]; a_idx[1] = t; }
        } else { /* a_nprimes == 3 */
            for (int attempt = 0; attempt < 50; attempt++) {
                SIMPLE_RNG();
                a_idx[0] = window_lo + (int)(rng_state % (unsigned long)(window_hi - window_lo + 1));
                SIMPLE_RNG();
                a_idx[1] = window_lo + (int)(rng_state % (unsigned long)(window_hi - window_lo + 1));
                SIMPLE_RNG();
                a_idx[2] = window_lo + (int)(rng_state % (unsigned long)(window_hi - window_lo + 1));
                if (a_idx[0] != a_idx[1] && a_idx[0] != a_idx[2] && a_idx[1] != a_idx[2]) break;
                if (attempt == 49) ok = 0;
            }
            /* Sort */
            for (int i = 0; i < 3; i++)
                for (int j = i+1; j < 3; j++)
                    if (a_idx[i] > a_idx[j]) { uint16_t t = a_idx[i]; a_idx[i] = a_idx[j]; a_idx[j] = t; }
        }

        if (!ok) { poly_count++; continue; }

        /* Compute a = product of selected primes */
        mpz_set_ui(a_val, 1);
        for (int k = 0; k < a_nprimes; k++) {
            mpz_mul_ui(a_val, a_val, fb[a_idx[k]].p);
        }

        /* Compute b such that b^2 == N (mod a) using CRT.
         * For each prime q_k, we know sqrt(N) mod q_k = fb[a_idx[k]].sqrt_n_mod_p
         * We combine using CRT to get b mod a. */

        if (a_nprimes == 1) {
            unsigned int q = fb[a_idx[0]].p;
            unsigned int sq = fb[a_idx[0]].sqrt_n_mod_p;
            /* Choose root closer to a/2 for better centering */
            mpz_set_ui(b_val, sq);
        } else if (a_nprimes == 2) {
            unsigned int q0 = fb[a_idx[0]].p;
            unsigned int q1 = fb[a_idx[1]].p;
            unsigned int s0 = fb[a_idx[0]].sqrt_n_mod_p;
            unsigned int s1 = fb[a_idx[1]].sqrt_n_mod_p;

            /* CRT: b = s0 * q1 * (q1^-1 mod q0) + s1 * q0 * (q0^-1 mod q1) mod (q0*q1) */
            unsigned int q1_inv_q0 = modinv_ui(q1 % q0, q0);
            unsigned int q0_inv_q1 = modinv_ui(q0 % q1, q1);

            mpz_set_ui(b_val, 0);
            mpz_set_ui(tmp, s0);
            mpz_mul_ui(tmp, tmp, q1);
            mpz_mul_ui(tmp, tmp, q1_inv_q0);
            mpz_add(b_val, b_val, tmp);

            mpz_set_ui(tmp, s1);
            mpz_mul_ui(tmp, tmp, q0);
            mpz_mul_ui(tmp, tmp, q0_inv_q1);
            mpz_add(b_val, b_val, tmp);

            mpz_mod(b_val, b_val, a_val);
        } else { /* a_nprimes == 3 */
            unsigned int q[3], s[3];
            for (int k = 0; k < 3; k++) {
                q[k] = fb[a_idx[k]].p;
                s[k] = fb[a_idx[k]].sqrt_n_mod_p;
            }

            /* CRT for 3 moduli */
            mpz_set_ui(b_val, 0);
            for (int k = 0; k < 3; k++) {
                /* M_k = a / q[k] */
                mpz_divexact_ui(tmp, a_val, q[k]);
                /* M_k^-1 mod q[k] */
                unsigned int mk_mod_qk = mpz_fdiv_ui(tmp, q[k]);
                unsigned int mk_inv = modinv_ui(mk_mod_qk, q[k]);
                /* term = s[k] * M_k * mk_inv */
                mpz_mul_ui(tmp, tmp, (unsigned long)s[k] * mk_inv % q[k]);
                /* Actually do it more carefully to avoid overflow */
                mpz_divexact_ui(tmp, a_val, q[k]);
                mpz_mul_ui(tmp2, tmp, s[k]);
                mpz_mul_ui(tmp2, tmp2, mk_inv);
                mpz_add(b_val, b_val, tmp2);
            }
            mpz_mod(b_val, b_val, a_val);
        }

        /* Verify b^2 == N mod a */
        mpz_mul(tmp, b_val, b_val);
        mpz_mod(tmp, tmp, a_val);
        mpz_mod(tmp2, gN, a_val);
        if (mpz_cmp(tmp, tmp2) != 0) {
            /* Try the other root: a - b */
            mpz_sub(b_val, a_val, b_val);
            mpz_mul(tmp, b_val, b_val);
            mpz_mod(tmp, tmp, a_val);
            if (mpz_cmp(tmp, tmp2) != 0) {
                /* This shouldn't happen, but skip if CRT failed */
                poly_count++;
                continue;
            }
        }

        /* Adjust b so that b == sqrt(N) mod a AND b is close to sqrt(N)
         * We want b such that b^2 == N (mod a) and |b| is close to sqrt(N).
         * Set b = b_val + k*a where k = round((sqrt(N) - b_val) / a) */
        mpz_sub(tmp, gSqrtN, b_val);
        mpz_tdiv_q(tmp, tmp, a_val);
        mpz_mul(tmp, tmp, a_val);
        mpz_add(b_val, b_val, tmp);

        /* Now b_val ~ sqrt(N) and b_val^2 == N (mod a) */
        /* Verify */
        mpz_mul(tmp, b_val, b_val);
        mpz_sub(tmp, tmp, gN);
        if (!mpz_divisible_p(tmp, a_val)) {
            poly_count++;
            continue;
        }

        /* Compute sieve roots for this polynomial.
         * Q_a(x) = (a*x + b)^2 - N
         * For prime p in factor base (p not dividing a):
         *   (a*x + b)^2 == N (mod p)
         *   a*x + b == +-sqrt(N) (mod p)
         *   x == (+-sqrt(N) - b) * a^{-1} (mod p)
         *
         * The sieve array represents x in [-M, M], stored at index x + M.
         * So sieve position = x + M, and root offset = ((+-r - b) * a^{-1} mod p) + M mod p.
         */
        for (int j = 0; j < fb_size; j++) {
            unsigned int p = fb[j].p;
            if (p == 2) {
                /* Handle 2 specially */
                unsigned long b_mod2 = mpz_fdiv_ui(b_val, 2);
                unsigned long a_mod2 = mpz_fdiv_ui(a_val, 2);
                if (a_mod2 == 0) {
                    /* a is even -- skip (shouldn't happen since a is product of odd FB primes) */
                    fb[j].r1 = fb[j].r2 = SIEVE_BLOCK + 1;
                    continue;
                }
                /* (a*x+b)^2 == N mod 2, need N mod 2 == 0 or always true mod 2 */
                /* For p=2, just sieve every position (crude but fine) */
                fb[j].r1 = 0;
                fb[j].r2 = SIEVE_BLOCK + 1; /* disable second root */
                continue;
            }

            /* Check if p divides a */
            int p_divides_a = 0;
            for (int k = 0; k < a_nprimes; k++) {
                if (fb[a_idx[k]].p == p) { p_divides_a = 1; break; }
            }

            if (p_divides_a) {
                /* When p | a: (a*x+b)^2 - N == 0 mod p means b^2 == N mod p
                 * which is true by construction. So Q_a(x) == 0 mod p for all x,
                 * but Q_a(x)/a may or may not be divisible by p.
                 * For Q_a(x)/a: since p|a, Q_a(x) = (ax+b)^2 - N, and p | (ax+b)^2 - N means
                 * p | b^2 - N (since p|a), which is true. But (ax+b)^2 - N = a*(Q/a), and p|a,
                 * so Q/a = ((ax+b)^2 - N)/a. We need to find when p | Q/a.
                 * Q/a = ((ax+b)^2 - N)/a. Since p|a, write a = p*a'. Then:
                 * ((p*a'*x + b)^2 - N) / (p*a') = (p^2*a'^2*x^2 + 2*p*a'*x*b + b^2 - N)/(p*a')
                 * = p*a'*x^2 + 2*x*b/1 + (b^2-N)/(p*a')  -- not exactly...
                 * Actually just solve: Q/a == 0 mod p.
                 * ((ax+b)^2 - N)/a == 0 mod p
                 * (ax+b)^2 - N == 0 mod (a*p)
                 * Since a = p*a', need (p*a'*x+b)^2 == N mod p^2*a'
                 * For just mod p: we need (b^2 - N)/a + 2*b*x == 0 mod p (linear approx)
                 * x == -(b^2-N)/(2*a*b) mod p
                 */
                unsigned long b_mod_p = mpz_fdiv_ui(b_val, p);
                unsigned long a_mod_pp = mpz_fdiv_ui(a_val, (unsigned long)p * p);
                mpz_mul(tmp, b_val, b_val);
                mpz_sub(tmp, tmp, gN);
                mpz_divexact(tmp, tmp, a_val);
                long qba_mod_p = (long)mpz_fdiv_ui(tmp, p);
                /* x == -qba / (2*b) mod p (but only if b != 0 mod p) */
                if (b_mod_p == 0) {
                    /* Single root at x=0 offset; just sieve every p-th */
                    fb[j].r1 = SIEVE_HALF % p;
                    fb[j].r2 = SIEVE_BLOCK + 1;
                } else {
                    unsigned int two_b = (2 * b_mod_p) % p;
                    unsigned int inv_2b = modinv_ui(two_b, p);
                    long x_root = (-(long)qba_mod_p * (long)inv_2b) % (long)p;
                    if (x_root < 0) x_root += p;
                    /* Convert to sieve position: x + M */
                    unsigned int pos = ((long)x_root + SIEVE_HALF) % p;
                    fb[j].r1 = pos;
                    fb[j].r2 = SIEVE_BLOCK + 1; /* only one root when p|a */
                }
                continue;
            }

            /* Normal case: p does not divide a */
            unsigned int r = fb[j].sqrt_n_mod_p;  /* sqrt(N) mod p */
            unsigned long a_mod_p = mpz_fdiv_ui(a_val, p);
            unsigned long b_mod_p = mpz_fdiv_ui(b_val, p);
            unsigned int a_inv = modinv_ui((unsigned int)a_mod_p, p);

            /* root1: x = (r - b) * a^{-1} mod p */
            long x1 = ((long)r - (long)b_mod_p % (long)p + (long)p) % (long)p;
            x1 = (x1 * (long)a_inv) % (long)p;
            /* root2: x = (-r - b) * a^{-1} mod p */
            long x2 = ((long)(p - r) - (long)b_mod_p % (long)p + (long)p) % (long)p;
            x2 = (x2 * (long)a_inv) % (long)p;

            /* Convert to sieve array positions: sieve[i] represents x = i - M */
            /* So position = (x_root + M) mod p, then sieve at pos, pos+p, pos+2p, ... */
            unsigned int pos1 = ((x1 + SIEVE_HALF) % p + p) % p;
            unsigned int pos2 = ((x2 + SIEVE_HALF) % p + p) % p;

            fb[j].r1 = pos1;
            fb[j].r2 = pos2;
        }

        /* Compute approximate log2 of |Q_a(x)/a| at x=M:
         * |Q_a(x)/a| ~ 2*M*sqrt(N)/a */
        double log2_approx = log2(2.0) + log2((double)SIEVE_HALF)
                            + mpz_sizeinbase(gSqrtN, 2)
                            - mpz_sizeinbase(a_val, 2);
        /* Subtract the x-dependent part (added per-position in sieve_block_mpoly) */
        double log2_base = log2_approx - log2((double)SIEVE_HALF);

        /* Sieve the block */
        sieve_block_mpoly(sieve_arr, SIEVE_BLOCK, log2_base);

        /* Threshold */
        float thr = (float)(fb[fb_size > 10 ? fb_size - 1 : 9].logp * 2.0);

        /* Scan for smooth candidates */
        for (long i = 0; i < SIEVE_BLOCK; i++) {
            if (sieve_arr[i] > thr) continue;

            long x = i - SIEVE_HALF;

            /* Compute (a*x + b) */
            mpz_mul_si(ax_b, a_val, x);
            mpz_add(ax_b, ax_b, b_val);

            /* Compute Q_a(x) = (a*x+b)^2 - N */
            mpz_mul(qa_val, ax_b, ax_b);
            mpz_sub(qa_val, qa_val, gN);

            /* Divide by a to get Q_a(x)/a */
            if (!mpz_divisible_p(qa_val, a_val)) continue;
            mpz_divexact(qa_val, qa_val, a_val);

            if (process_candidate_mpqs(qa_val, ax_b, a_nprimes, a_idx)) {
                Relation *r = &rels[nrels - 1];
                if (r->nlp == 0) full++;
                else if (r->nlp == 1) p1++;
                else if (r->nlp == 2) p2++;
                else p3++;
            }
        }

        poly_count++;

        if (poly_count % 50 == 0) {
            fprintf(stderr, "Poly %d: %d full + %d p1 + %d p2 + %d p3 = %d (need ~%d)\n",
                    poly_count, full, p1, p2, p3, nrels, target_rels);
        }
    }

    mpz_clears(qa_val, ax_b, NULL);
    for (int i = 0; i < MAX_APRIME * 2; i++)
        mpz_clear(b_solutions[i]);
    mpz_clears(target_a, a_val, b_val, tmp, tmp2, ainv, NULL);

    free(sieve_arr);

    fprintf(stderr, "Collection done (%d polys): %d full + %d p1 + %d p2 + %d p3 = %d total\n",
            poly_count, full, p1, p2, p3, nrels);

    /* Linear algebra with extended matrix (FB + LP columns) */
    int result_found = do_linear_algebra(result);

    /* Cleanup */
    for (int i = 0; i < nrels; i++) {
        mpz_clear(rels[i].ax_plus_b);
        free(rels[i].exps);
    }
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

    /* Stage 4: SRG sieve (MPQS) */
    fprintf(stderr, "Stage 4: SRG-MPQS sieve\n");
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
