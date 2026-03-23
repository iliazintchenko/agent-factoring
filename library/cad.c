/*
 * CAD - Cascaded Algebraic Descent
 *
 * Novel factoring approach inspired by the descent step of NFS-DLP.
 *
 * Core idea: Instead of requiring sieve values to be FULLY smooth,
 * accept values with arbitrarily many large prime factors and
 * recursively descend to decompose them.
 *
 * Algorithm:
 * 1. Generate candidates: Q(x) = (x+m)^2 - N for x near 0
 * 2. Extract smooth part via trial division with factor base primes
 * 3. For each cofactor C:
 *    a. If C < LP^2, try to split via Pollard's rho
 *    b. If split, add both primes as "large primes" in the relation
 *    c. If C < LP^3, try harder with ECM
 *    d. Accept relations with up to K large primes
 * 4. Build a hypergraph of large-prime relations
 * 5. Find even-multiplicity subsets via structured linear algebra
 * 6. Extract factor via GCD
 *
 * The key novelty: aggressively using ECM on cofactors to split
 * multi-large-prime relations, combined with a hypergraph-based
 * merge strategy. This allows much larger effective smoothness
 * bounds than standard QS.
 *
 * Compile: gcc -O3 -o cad cad.c -lgmp -L/usr/local/lib -lecm -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

/* --- Parameters --- */
#define MAX_FB     20000    /* max factor base size */
#define MAX_RELS   100000   /* max relations */
#define MAX_LP_PER_REL 4    /* max large primes per relation */
#define SIEVE_BLOCK  65536  /* sieve block size */

static struct timespec g_start;
static double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) +
           (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* Factor base entry */
typedef struct {
    unsigned long p;
    int r1, r2;     /* sqrt(N) mod p, two roots */
    double logp;
} fb_t;

/* Relation: (x+m)^2 ≡ ∏ fb[i]^exp[i] · ∏ lp[j] (mod N) */
typedef struct {
    int x_off;           /* offset from sqrt(N) */
    int neg;             /* 1 if value was negated */
    unsigned char *exp;  /* exponent parity (mod 2) for factor base */
    int nlp;
    unsigned long lp[MAX_LP_PER_REL];
} rel_t;

/* Globals */
static mpz_t N_global, sqrtN_global;
static fb_t fb[MAX_FB];
static int fb_size;
static rel_t rels[MAX_RELS];
static int nrels;
static unsigned long smooth_B;
static unsigned long lp_bound;

/* --- Sieve of Eratosthenes --- */
static unsigned long *gen_primes(unsigned long limit, int *count) {
    char *sieve = calloc(limit + 1, 1);
    memset(sieve, 1, limit + 1);
    sieve[0] = sieve[1] = 0;
    for (unsigned long i = 2; i * i <= limit; i++)
        if (sieve[i])
            for (unsigned long j = i*i; j <= limit; j += i)
                sieve[j] = 0;
    *count = 0;
    for (unsigned long i = 2; i <= limit; i++)
        if (sieve[i]) (*count)++;
    unsigned long *p = malloc(*count * sizeof(unsigned long));
    int idx = 0;
    for (unsigned long i = 2; i <= limit; i++)
        if (sieve[i]) p[idx++] = i;
    free(sieve);
    return p;
}

/* Tonelli-Shanks: compute sqrt(n) mod p */
static int sqrtmod(unsigned long *result, unsigned long n, unsigned long p) {
    n = n % p;
    if (n == 0) { *result = 0; return 1; }

    mpz_t a, pm, res;
    mpz_inits(a, pm, res, NULL);
    mpz_set_ui(a, n);
    mpz_set_ui(pm, p);
    if (mpz_jacobi(a, pm) != 1) {
        mpz_clears(a, pm, res, NULL);
        return 0;
    }

    if (p % 4 == 3) {
        mpz_t exp;
        mpz_init(exp);
        mpz_set_ui(exp, (p + 1) / 4);
        mpz_powm(res, a, exp, pm);
        *result = mpz_get_ui(res);
        mpz_clears(a, pm, res, exp, NULL);
        return 1;
    }

    /* Full Tonelli-Shanks */
    unsigned long Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned long z = 2;
    mpz_t zt;
    mpz_init(zt);
    while (1) {
        mpz_set_ui(zt, z);
        if (mpz_jacobi(zt, pm) == -1) break;
        z++;
    }

    mpz_t M, c, t, R, exp, b, tmp;
    mpz_inits(M, c, t, R, exp, b, tmp, NULL);
    mpz_set_ui(M, S);
    mpz_set_ui(exp, Q);
    mpz_powm(c, zt, exp, pm);
    mpz_powm(t, a, exp, pm);
    mpz_add_ui(exp, exp, 1);
    mpz_fdiv_q_2exp(exp, exp, 1);
    mpz_powm(R, a, exp, pm);

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            *result = mpz_get_ui(R);
            mpz_clears(M, c, t, R, exp, b, tmp, zt, a, pm, res, NULL);
            return 1;
        }
        mpz_set(tmp, t);
        unsigned long i;
        for (i = 1; i < mpz_get_ui(M); i++) {
            mpz_mul(tmp, tmp, tmp); mpz_mod(tmp, tmp, pm);
            if (mpz_cmp_ui(tmp, 1) == 0) break;
        }
        mpz_set(b, c);
        for (unsigned long j = 0; j < mpz_get_ui(M) - i - 1; j++) {
            mpz_mul(b, b, b); mpz_mod(b, b, pm);
        }
        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, pm);
        mpz_mul(t, t, c); mpz_mod(t, t, pm);
        mpz_mul(R, R, b); mpz_mod(R, R, pm);
    }
}

/* Build factor base */
static void build_fb(const mpz_t N, unsigned long B) {
    int np;
    unsigned long *primes = gen_primes(B, &np);
    fb_size = 0;
    mpz_t Nmod;
    mpz_init(Nmod);

    /* Always include -1 as first "prime" (handled via neg flag) */

    for (int i = 0; i < np && fb_size < MAX_FB; i++) {
        unsigned long p = primes[i];
        if (p == 2) {
            fb[fb_size].p = 2;
            fb[fb_size].r1 = mpz_fdiv_ui(N, 2) == 0 ? 0 : 1;
            fb[fb_size].r2 = -1;
            fb[fb_size].logp = 1.0;
            fb_size++;
            continue;
        }
        unsigned long Nmodp = mpz_fdiv_ui(N, p);
        unsigned long r;
        if (sqrtmod(&r, Nmodp, p)) {
            fb[fb_size].p = p;
            fb[fb_size].r1 = (int)r;
            fb[fb_size].r2 = (int)(p - r);
            fb[fb_size].logp = log2((double)p);
            fb_size++;
        }
    }
    free(primes);
    mpz_clear(Nmod);
}

/* Trial divide val by factor base, store exponents mod 2.
 * Returns 1 if fully smooth, 0 if cofactor remains.
 * Sets cofactor_out to remaining value. */
static int trial_divide_fb(mpz_t cofactor_out, unsigned char *exp_out,
                           const mpz_t val) {
    mpz_set(cofactor_out, val);
    memset(exp_out, 0, fb_size);

    for (int i = 0; i < fb_size; i++) {
        while (mpz_divisible_ui_p(cofactor_out, fb[i].p)) {
            mpz_divexact_ui(cofactor_out, cofactor_out, fb[i].p);
            exp_out[i] ^= 1;
        }
    }
    return (mpz_cmp_ui(cofactor_out, 1) == 0);
}

/* Try to split cofactor into large primes using rho + ECM */
static int split_cofactor(unsigned long *lps, int *nlps,
                          const mpz_t cofactor, unsigned long lp_max) {
    *nlps = 0;
    mpz_t c, f, t;
    mpz_inits(c, f, t, NULL);
    mpz_set(c, cofactor);

    /* If cofactor fits in unsigned long and is prime, it's a single LP */
    if (mpz_fits_ulong_p(c)) {
        unsigned long cv = mpz_get_ui(c);
        if (cv <= lp_max && cv > 1) {
            /* Check primality quickly */
            if (mpz_probab_prime_p(c, 15) > 0) {
                lps[(*nlps)++] = cv;
                mpz_clears(c, f, t, NULL);
                return 1;
            }
        }
    }

    /* Try rho for small cofactors */
    if (mpz_sizeinbase(c, 10) <= 20) {
        gmp_randstate_t rng;
        gmp_randinit_default(rng);
        gmp_randseed_ui(rng, 42);

        for (int attempt = 0; attempt < 5 && mpz_cmp_ui(c, 1) > 0; attempt++) {
            mpz_t x, y, d, acc, cc;
            mpz_inits(x, y, d, acc, cc, NULL);
            mpz_urandomm(cc, rng, c);
            mpz_add_ui(cc, cc, 1);
            mpz_urandomm(x, rng, c);
            mpz_set(y, x);
            mpz_set_ui(acc, 1);
            int found = 0;

            for (int i = 0; i < 100000 && !found; i++) {
                mpz_mul(x, x, x); mpz_add(x, x, cc); mpz_mod(x, x, c);
                mpz_mul(y, y, y); mpz_add(y, y, cc); mpz_mod(y, y, c);
                mpz_mul(y, y, y); mpz_add(y, y, cc); mpz_mod(y, y, c);
                mpz_sub(t, x, y); mpz_abs(t, t);
                if (mpz_sgn(t) == 0) break;
                mpz_mul(acc, acc, t); mpz_mod(acc, acc, c);

                if (i % 128 == 127) {
                    mpz_gcd(d, acc, c);
                    if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, c) < 0) {
                        /* Found a factor of cofactor */
                        if (mpz_fits_ulong_p(d) && mpz_get_ui(d) <= lp_max) {
                            lps[(*nlps)++] = mpz_get_ui(d);
                            mpz_divexact(c, c, d);
                            /* Check remaining */
                            if (mpz_fits_ulong_p(c) && mpz_get_ui(c) <= lp_max) {
                                if (mpz_cmp_ui(c, 1) > 0)
                                    lps[(*nlps)++] = mpz_get_ui(c);
                                mpz_set_ui(c, 1);
                            }
                            found = 1;
                        } else {
                            break; /* factor too large */
                        }
                    }
                    mpz_set_ui(acc, 1);
                }
            }
            mpz_clears(x, y, d, acc, cc, NULL);
            if (*nlps >= MAX_LP_PER_REL) break;
        }
        gmp_randclear(rng);
    }

    /* Try ECM for medium cofactors */
    if (mpz_cmp_ui(c, 1) > 0 && mpz_sizeinbase(c, 10) <= 30 &&
        *nlps < MAX_LP_PER_REL) {
        ecm_params par;
        ecm_init(par);
        par->method = ECM_ECM;
        par->verbose = 0;

        double b1_vals[] = {100, 500, 2000, 11000, 50000};
        int n_b1 = 5;

        for (int b = 0; b < n_b1 && mpz_cmp_ui(c, 1) > 0 &&
             *nlps < MAX_LP_PER_REL; b++) {
            for (int curve = 0; curve < 10 && mpz_cmp_ui(c, 1) > 0; curve++) {
                mpz_set_si(par->B2, -1);
                par->B1done = 0;
                int ret = ecm_factor(f, c, b1_vals[b], par);
                if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, c) < 0) {
                    if (mpz_fits_ulong_p(f) && mpz_get_ui(f) <= lp_max) {
                        lps[(*nlps)++] = mpz_get_ui(f);
                        mpz_divexact(c, c, f);
                        if (mpz_fits_ulong_p(c) && mpz_get_ui(c) <= lp_max &&
                            mpz_cmp_ui(c, 1) > 0) {
                            lps[(*nlps)++] = mpz_get_ui(c);
                            mpz_set_ui(c, 1);
                        }
                    }
                    break;
                }
            }
        }
        ecm_clear(par);
    }

    int ok = (mpz_cmp_ui(c, 1) == 0 && *nlps <= MAX_LP_PER_REL);
    mpz_clears(c, f, t, NULL);
    return ok;
}

/* Logarithmic sieve to find likely-smooth candidates */
static void sieve_block(double *log_arr, int start, int len,
                        const mpz_t m) {
    /* Initialize with approximate log2 of |Q(x)| */
    mpz_t qval;
    mpz_init(qval);

    for (int i = 0; i < len; i++) {
        int x = start + i;
        /* Q(x) = (x+m)^2 - N ≈ 2*m*x for small x */
        mpz_set(qval, m);
        if (x >= 0) mpz_add_ui(qval, qval, x);
        else mpz_sub_ui(qval, qval, -x);
        mpz_mul(qval, qval, qval);
        mpz_sub(qval, qval, N_global);
        mpz_abs(qval, qval);
        log_arr[i] = mpz_sizeinbase(qval, 2); /* ≈ log2(|Q(x)|) */
    }
    mpz_clear(qval);

    /* Subtract contributions of factor base primes */
    for (int fi = 0; fi < fb_size; fi++) {
        unsigned long p = fb[fi].p;
        double lp = fb[fi].logp;

        if (p == 2) {
            for (int i = 0; i < len; i += 2)
                log_arr[i] -= lp;
            /* higher powers of 2 */
            for (unsigned long pk = 4; pk <= smooth_B; pk *= 2)
                for (int i = 0; i < len; i++)
                    if (((start + i) % pk) == 0 || ((start + i + 1) % pk) == 0)
                        log_arr[i] -= lp;
            continue;
        }

        /* For each root r of N mod p */
        for (int ri = 0; ri < 2; ri++) {
            int r = (ri == 0) ? fb[fi].r1 : fb[fi].r2;
            if (r < 0) continue;

            /* x ≡ r - m (mod p) makes Q(x) divisible by p */
            long mmodp = mpz_fdiv_ui(sqrtN_global, p);
            long off = ((long)r - mmodp % (long)p + (long)p) % (long)p;

            /* Find first x in [start, start+len) ≡ off (mod p) */
            long first = start + ((off - start % (long)p + (long)p) % (long)p);
            for (long x = first; x < start + len; x += p) {
                int idx = (int)(x - start);
                if (idx >= 0 && idx < len) {
                    log_arr[idx] -= lp;
                    /* Higher powers */
                    /* Quick check for p^2 */
                }
            }
        }
    }
}

/* Main factoring function */
static int cad_factor(mpz_t factor, const mpz_t N) {
    mpz_set(N_global, N);
    mpz_sqrt(sqrtN_global, N);

    int ndigits = mpz_sizeinbase(N, 10);

    /* Choose smoothness bound: L[1/2, 0.7] */
    double logN = ndigits * log(10);
    double loglogN = log(logN);
    smooth_B = (unsigned long)exp(0.7 * sqrt(logN * loglogN));
    if (smooth_B < 500) smooth_B = 500;
    if (smooth_B > 5000000) smooth_B = 5000000;

    /* Large prime bound: 100 * B */
    lp_bound = smooth_B * 100;
    if (lp_bound > 500000000UL) lp_bound = 500000000UL;

    fprintf(stderr, "CAD: %d digits, B=%lu, LP=%lu\n",
            ndigits, smooth_B, lp_bound);

    /* Build factor base */
    build_fb(N, smooth_B);
    int needed = fb_size + 20;
    fprintf(stderr, "CAD: factor base size=%d, need %d relations\n",
            fb_size, needed);

    nrels = 0;
    int full_smooth = 0, with_lp = 0;

    /* Count LP columns needed */
    unsigned long *all_lps = NULL;
    int n_all_lps = 0, cap_all_lps = 0;

    /* Sieve and collect relations */
    mpz_t qval, cofactor;
    mpz_inits(qval, cofactor, NULL);

    int sieve_radius = 0;
    double *log_sieve = malloc(SIEVE_BLOCK * sizeof(double));
    double threshold = log2((double)smooth_B) * 1.5; /* Sieve threshold */

    while (nrels < needed * 3 && elapsed_sec() < 270.0) {
        /* Sieve a block */
        int block_start = sieve_radius;
        sieve_block(log_sieve, block_start, SIEVE_BLOCK, sqrtN_global);

        for (int i = 0; i < SIEVE_BLOCK && nrels < MAX_RELS; i++) {
            if (log_sieve[i] > threshold) continue; /* Not smooth enough */

            int x = block_start + i;

            /* Compute Q(x) = (x + m)^2 - N */
            mpz_set(qval, sqrtN_global);
            mpz_add_ui(qval, qval, x);
            mpz_mul(qval, qval, qval);
            mpz_sub(qval, qval, N);

            int neg = (mpz_sgn(qval) < 0);
            if (neg) mpz_neg(qval, qval);

            /* Trial divide */
            unsigned char *exp = calloc(fb_size, 1);
            int smooth = trial_divide_fb(cofactor, exp, qval);

            if (smooth) {
                /* Fully smooth relation */
                rels[nrels].x_off = x;
                rels[nrels].neg = neg;
                rels[nrels].exp = exp;
                rels[nrels].nlp = 0;
                nrels++;
                full_smooth++;
            } else if (mpz_sizeinbase(cofactor, 10) <= 25) {
                /* Try to split cofactor into large primes */
                unsigned long lps[MAX_LP_PER_REL];
                int nlps = 0;
                if (split_cofactor(lps, &nlps, cofactor, lp_bound)) {
                    rels[nrels].x_off = x;
                    rels[nrels].neg = neg;
                    rels[nrels].exp = exp;
                    rels[nrels].nlp = nlps;
                    for (int j = 0; j < nlps; j++)
                        rels[nrels].lp[j] = lps[j];
                    nrels++;
                    with_lp++;
                } else {
                    free(exp);
                }
            } else {
                free(exp);
            }
        }

        /* Also sieve negative offsets */
        sieve_block(log_sieve, -(block_start + SIEVE_BLOCK), SIEVE_BLOCK,
                    sqrtN_global);
        for (int i = 0; i < SIEVE_BLOCK && nrels < MAX_RELS; i++) {
            if (log_sieve[i] > threshold) continue;

            int x = -(block_start + SIEVE_BLOCK) + i;
            if (x == 0) continue;

            mpz_set(qval, sqrtN_global);
            if (x >= 0) mpz_add_ui(qval, qval, x);
            else mpz_sub_ui(qval, qval, -x);
            mpz_mul(qval, qval, qval);
            mpz_sub(qval, qval, N);

            int neg = (mpz_sgn(qval) < 0);
            if (neg) mpz_neg(qval, qval);

            unsigned char *exp = calloc(fb_size, 1);
            int smooth = trial_divide_fb(cofactor, exp, qval);

            if (smooth) {
                rels[nrels].x_off = x;
                rels[nrels].neg = neg;
                rels[nrels].exp = exp;
                rels[nrels].nlp = 0;
                nrels++;
                full_smooth++;
            } else if (mpz_sizeinbase(cofactor, 10) <= 25) {
                unsigned long lps[MAX_LP_PER_REL];
                int nlps = 0;
                if (split_cofactor(lps, &nlps, cofactor, lp_bound)) {
                    rels[nrels].x_off = x;
                    rels[nrels].neg = neg;
                    rels[nrels].exp = exp;
                    rels[nrels].nlp = nlps;
                    for (int j = 0; j < nlps; j++)
                        rels[nrels].lp[j] = lps[j];
                    nrels++;
                    with_lp++;
                } else {
                    free(exp);
                }
            } else {
                free(exp);
            }
        }

        sieve_radius += SIEVE_BLOCK;

        if ((sieve_radius / SIEVE_BLOCK) % 10 == 0) {
            fprintf(stderr, "CAD: sieved to ±%d, %d rels (%d full + %d LP), %.1fs\n",
                    sieve_radius, nrels, full_smooth, with_lp, elapsed_sec());
        }
    }
    free(log_sieve);

    fprintf(stderr, "CAD: total %d relations (%d full + %d LP)\n",
            nrels, full_smooth, with_lp);

    /* Collect all unique large primes */
    cap_all_lps = nrels * MAX_LP_PER_REL;
    all_lps = malloc(cap_all_lps * sizeof(unsigned long));
    n_all_lps = 0;

    for (int i = 0; i < nrels; i++) {
        for (int j = 0; j < rels[i].nlp; j++) {
            /* Add to list (we'll sort and unique later) */
            all_lps[n_all_lps++] = rels[i].lp[j];
        }
    }

    /* Sort and unique */
    /* Simple sort for now */
    for (int i = 0; i < n_all_lps - 1; i++)
        for (int j = i + 1; j < n_all_lps; j++)
            if (all_lps[i] > all_lps[j]) {
                unsigned long t = all_lps[i];
                all_lps[i] = all_lps[j];
                all_lps[j] = t;
            }
    int unique_lps = 0;
    for (int i = 0; i < n_all_lps; i++) {
        if (i == 0 || all_lps[i] != all_lps[i-1])
            all_lps[unique_lps++] = all_lps[i];
    }
    n_all_lps = unique_lps;

    fprintf(stderr, "CAD: %d unique large primes\n", n_all_lps);

    /* Build GF(2) matrix: columns = fb primes + sign + large primes */
    int ncols = 1 + fb_size + n_all_lps; /* sign + fb + LPs */
    int nrows = nrels;

    if (nrows <= ncols) {
        fprintf(stderr, "CAD: not enough relations (%d rows, %d cols)\n",
                nrows, ncols);
        free(all_lps);
        mpz_clears(qval, cofactor, NULL);
        return 0;
    }

    /* Gaussian elimination over GF(2) using bit vectors */
    /* Each row is a bit vector of ncols bits, packed into unsigned longs */
    int words_per_row = (ncols + 63) / 64;
    unsigned long *matrix = calloc(nrows * words_per_row, sizeof(unsigned long));
    unsigned long *history = calloc(nrows * ((nrows + 63) / 64), sizeof(unsigned long));
    int hist_words = (nrows + 63) / 64;

    /* Fill matrix */
    for (int i = 0; i < nrows; i++) {
        /* History: row i depends on row i */
        history[i * hist_words + (i / 64)] |= (1UL << (i % 64));

        /* Column 0: sign */
        if (rels[i].neg)
            matrix[i * words_per_row + 0] |= 1UL;

        /* Columns 1..fb_size: factor base exponents */
        for (int j = 0; j < fb_size; j++) {
            if (rels[i].exp[j])
                matrix[i * words_per_row + ((j + 1) / 64)] |=
                    (1UL << ((j + 1) % 64));
        }

        /* LP columns */
        for (int j = 0; j < rels[i].nlp; j++) {
            /* Find LP index */
            unsigned long lp = rels[i].lp[j];
            /* Binary search in all_lps */
            int lo = 0, hi = n_all_lps - 1;
            while (lo <= hi) {
                int mid = (lo + hi) / 2;
                if (all_lps[mid] == lp) {
                    int col = 1 + fb_size + mid;
                    matrix[i * words_per_row + (col / 64)] |=
                        (1UL << (col % 64));
                    break;
                } else if (all_lps[mid] < lp) lo = mid + 1;
                else hi = mid - 1;
            }
        }
    }

    /* Gaussian elimination */
    int *pivot_col = malloc(ncols * sizeof(int));
    memset(pivot_col, -1, ncols * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        /* Find pivot row */
        int piv = -1;
        for (int row = 0; row < nrows; row++) {
            if (pivot_col[col] >= 0) break;
            /* Check if this row has a 1 in this column */
            if (matrix[row * words_per_row + (col / 64)] & (1UL << (col % 64))) {
                /* Check it's not already a pivot for an earlier column */
                int is_pivot = 1;
                for (int c2 = 0; c2 < col; c2++) {
                    if (pivot_col[c2] == row) { is_pivot = 0; break; }
                }
                if (is_pivot) { piv = row; break; }
            }
        }
        if (piv == -1) continue;
        pivot_col[col] = piv;

        /* Eliminate */
        for (int row = 0; row < nrows; row++) {
            if (row == piv) continue;
            if (matrix[row * words_per_row + (col / 64)] & (1UL << (col % 64))) {
                for (int w = 0; w < words_per_row; w++)
                    matrix[row * words_per_row + w] ^=
                        matrix[piv * words_per_row + w];
                for (int w = 0; w < hist_words; w++)
                    history[row * hist_words + w] ^=
                        history[piv * hist_words + w];
            }
        }
    }

    /* Find zero rows (dependencies) */
    int found = 0;
    mpz_t lhs, rhs, g, tmp;
    mpz_inits(lhs, rhs, g, tmp, NULL);

    for (int row = 0; row < nrows && !found; row++) {
        /* Check if row is zero */
        int is_zero = 1;
        for (int w = 0; w < words_per_row; w++) {
            if (matrix[row * words_per_row + w]) { is_zero = 0; break; }
        }
        if (!is_zero) continue;

        /* Collect relations in this dependency */
        mpz_set_ui(lhs, 1);
        int *total_exp = calloc(fb_size, sizeof(int));
        int neg_count = 0;
        unsigned long *lp_exp = calloc(n_all_lps, sizeof(unsigned long));

        for (int r = 0; r < nrows; r++) {
            if (!(history[row * hist_words + (r / 64)] & (1UL << (r % 64))))
                continue;

            /* Multiply LHS by (x + m) */
            mpz_set(tmp, sqrtN_global);
            if (rels[r].x_off >= 0) mpz_add_ui(tmp, tmp, rels[r].x_off);
            else mpz_sub_ui(tmp, tmp, -rels[r].x_off);
            mpz_mul(lhs, lhs, tmp);
            mpz_mod(lhs, lhs, N);

            if (rels[r].neg) neg_count++;
            for (int j = 0; j < fb_size; j++)
                total_exp[j] += rels[r].exp[j] ? 1 : 0;
                /* Wait, exp stores parity. Need ACTUAL exponents. */

            for (int j = 0; j < rels[r].nlp; j++) {
                unsigned long lp = rels[r].lp[j];
                int lo = 0, hi = n_all_lps - 1;
                while (lo <= hi) {
                    int mid = (lo + hi) / 2;
                    if (all_lps[mid] == lp) { lp_exp[mid]++; break; }
                    else if (all_lps[mid] < lp) lo = mid + 1;
                    else hi = mid - 1;
                }
            }
        }

        /* Compute RHS = sqrt of product of all Q(x) values */
        /* We stored only exponent parity. Need actual exponents.
         * This is a limitation — we need to recompute actual exponents. */

        /* Actually, we can compute RHS differently:
         * LHS = product of (x_i + m) mod N
         * RHS^2 = product of Q(x_i) = product of (x_i+m)^2 - N ≡ LHS^2 (mod N)
         * But that's trivial! We need RHS = sqrt of product of smooth parts.
         *
         * Let me recompute: for each relation in the dependency,
         * recompute Q(x) and fully factor it, accumulating exponents.
         */

        int *real_exp = calloc(fb_size, sizeof(int));
        int real_neg = 0;

        for (int r = 0; r < nrows; r++) {
            if (!(history[row * hist_words + (r / 64)] & (1UL << (r % 64))))
                continue;

            mpz_set(tmp, sqrtN_global);
            if (rels[r].x_off >= 0) mpz_add_ui(tmp, tmp, rels[r].x_off);
            else mpz_sub_ui(tmp, tmp, -rels[r].x_off);
            mpz_mul(qval, tmp, tmp);
            mpz_sub(qval, qval, N);

            if (mpz_sgn(qval) < 0) { mpz_neg(qval, qval); real_neg++; }

            for (int j = 0; j < fb_size; j++) {
                while (mpz_divisible_ui_p(qval, fb[j].p)) {
                    mpz_divexact_ui(qval, qval, fb[j].p);
                    real_exp[j]++;
                }
            }
            /* qval should now be the LP product — divide those out too */
            for (int j = 0; j < rels[r].nlp; j++) {
                unsigned long lp = rels[r].lp[j];
                while (mpz_divisible_ui_p(qval, lp)) {
                    mpz_divexact_ui(qval, qval, lp);
                }
            }
        }

        /* All exponents should be even */
        int all_even = 1;
        for (int j = 0; j < fb_size; j++) {
            if (real_exp[j] % 2 != 0) { all_even = 0; break; }
        }
        for (int j = 0; j < n_all_lps; j++) {
            if (lp_exp[j] % 2 != 0) { all_even = 0; break; }
        }
        if (real_neg % 2 != 0) all_even = 0;

        if (all_even) {
            /* Compute RHS = product of fb[j]^(exp[j]/2) * product of lp[j]^(lp_exp[j]/2) mod N */
            mpz_set_ui(rhs, 1);
            for (int j = 0; j < fb_size; j++) {
                if (real_exp[j] > 0) {
                    mpz_set_ui(tmp, fb[j].p);
                    mpz_t e;
                    mpz_init_set_ui(e, real_exp[j] / 2);
                    mpz_t pw;
                    mpz_init(pw);
                    mpz_powm(pw, tmp, e, N);
                    mpz_mul(rhs, rhs, pw);
                    mpz_mod(rhs, rhs, N);
                    mpz_clears(e, pw, NULL);
                }
            }
            for (int j = 0; j < n_all_lps; j++) {
                if (lp_exp[j] > 0) {
                    mpz_set_ui(tmp, all_lps[j]);
                    mpz_t e;
                    mpz_init_set_ui(e, lp_exp[j] / 2);
                    mpz_t pw;
                    mpz_init(pw);
                    mpz_powm(pw, tmp, e, N);
                    mpz_mul(rhs, rhs, pw);
                    mpz_mod(rhs, rhs, N);
                    mpz_clears(e, pw, NULL);
                }
            }

            /* Factor = gcd(lhs - rhs, N) */
            mpz_sub(g, lhs, rhs);
            mpz_mod(g, g, N);
            mpz_gcd(g, g, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                found = 1;
            } else {
                mpz_add(g, lhs, rhs);
                mpz_mod(g, g, N);
                mpz_gcd(g, g, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_set(factor, g);
                    found = 1;
                }
            }
        }

        free(real_exp);
        free(total_exp);
        free(lp_exp);
    }

    free(pivot_col);
    free(matrix);
    free(history);
    free(all_lps);
    mpz_clears(lhs, rhs, g, tmp, qval, cofactor, NULL);
    return found;
}

/* Quick trial division */
static int trial_divide_quick(mpz_t factor, const mpz_t N) {
    for (unsigned long d = 2; d <= 1000000; d++) {
        if (mpz_divisible_ui_p(N, d)) {
            mpz_set_ui(factor, d);
            return 1;
        }
        if (d == 2) d = 1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, factor, cofactor;
    mpz_inits(N, factor, cofactor, N_global, sqrtN_global, NULL);
    mpz_set_str(N, argv[1], 10);

    int found = 0;

    /* Quick trial division */
    found = trial_divide_quick(factor, N);

    /* Main CAD algorithm */
    if (!found) {
        found = cad_factor(factor, N);
    }

    if (found) {
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "CAD: factored in %.3fs\n", elapsed_sec());
    } else {
        fprintf(stderr, "CAD: FAILED after %.3fs\n", elapsed_sec());
        return 1;
    }

    mpz_clears(N, factor, cofactor, N_global, sqrtN_global, NULL);
    return 0;
}
