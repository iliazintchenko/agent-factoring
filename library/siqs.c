/*
 * Self-Initializing Quadratic Sieve (SIQS)
 * Usage: ./siqs <N>
 * Targets: 50-100 digit balanced semiprimes
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>

/* ==================== Parameters ==================== */

typedef struct {
    int fb_size;
    int sieve_size;    /* total sieve array = 2*sieve_size */
    int lp_mult;       /* large prime multiplier (max LP = lp_mult * largest FB prime) */
    int num_a_factors;
    double thresh_adj; /* threshold adjustment (lower = more candidates) */
} params_t;

/* Parameters tuned for digit count.
 * Key insight: for SIQS, smaller FB + lower threshold = faster relation finding.
 * Formula: L = exp(sqrt(ln(N)*ln(ln(N)))), FB ≈ L^(1/sqrt(2)), M ≈ L^(1/sqrt(2))
 */
static params_t get_params(int digits) {
    if (digits <= 40) return (params_t){100, 32768, 30, 3, 0.75};
    if (digits <= 45) return (params_t){200, 32768, 40, 4, 0.78};
    if (digits <= 50) return (params_t){350, 65536, 40, 5, 0.80};
    if (digits <= 55) return (params_t){650, 65536, 50, 6, 0.82};
    if (digits <= 60) return (params_t){1200, 65536, 50, 7, 0.84};
    if (digits <= 65) return (params_t){2500, 131072, 60, 8, 0.86};
    if (digits <= 70) return (params_t){4500, 131072, 60, 9, 0.88};
    if (digits <= 75) return (params_t){8000, 196608, 70, 10, 0.90};
    if (digits <= 80) return (params_t){15000, 262144, 80, 11, 0.91};
    if (digits <= 85) return (params_t){25000, 262144, 90, 11, 0.92};
    if (digits <= 90) return (params_t){42000, 393216, 100, 12, 0.93};
    if (digits <= 95) return (params_t){65000, 524288, 110, 12, 0.94};
    return (params_t){100000, 524288, 120, 13, 0.95};
}

/* ==================== Modular Arithmetic ==================== */

/* Tonelli-Shanks: compute sqrt(n_mod_p) mod p.  Returns 0 if no sqrt exists. */
static unsigned int tonelli_shanks(unsigned long n_mod_p, unsigned int p) {
    if (n_mod_p == 0) return 0;
    if (p == 2) return n_mod_p & 1;

    /* Check Euler criterion */
    mpz_t base, exp, mod, result;
    mpz_inits(base, exp, mod, result, NULL);
    mpz_set_ui(base, n_mod_p);
    mpz_set_ui(mod, p);
    mpz_set_ui(exp, (p - 1) / 2);
    mpz_powm(result, base, exp, mod);
    if (mpz_cmp_ui(result, 1) != 0) {
        mpz_clears(base, exp, mod, result, NULL);
        return 0; /* not a QR */
    }

    if (p % 4 == 3) {
        mpz_set_ui(exp, (p + 1) / 4);
        mpz_powm(result, base, exp, mod);
        unsigned int r = mpz_get_ui(result);
        mpz_clears(base, exp, mod, result, NULL);
        return r;
    }

    /* General case */
    unsigned long Q = p - 1;
    int S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned long z = 2;
    while (1) {
        mpz_set_ui(base, z);
        mpz_set_ui(exp, (p - 1) / 2);
        mpz_powm(result, base, exp, mod);
        if (mpz_cmp_ui(result, p - 1) == 0) break;
        z++;
    }

    mpz_t M, c, t, R, b, tmp;
    mpz_inits(M, c, t, R, b, tmp, NULL);
    mpz_set_ui(M, S);
    mpz_set_ui(c, z); mpz_set_ui(exp, Q); mpz_powm(c, c, exp, mod);
    mpz_set_ui(t, n_mod_p); mpz_powm(t, t, exp, mod);
    mpz_set_ui(R, n_mod_p); mpz_set_ui(exp, (Q + 1) / 2); mpz_powm(R, R, exp, mod);

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            unsigned int r = mpz_get_ui(R);
            mpz_clears(base, exp, mod, result, M, c, t, R, b, tmp, NULL);
            return r;
        }
        int i = 0;
        mpz_set(tmp, t);
        while (mpz_cmp_ui(tmp, 1) != 0) {
            mpz_mul(tmp, tmp, tmp); mpz_mod(tmp, tmp, mod);
            i++;
        }
        mpz_set(b, c);
        for (int j = 0; j < (int)mpz_get_ui(M) - i - 1; j++) {
            mpz_mul(b, b, b); mpz_mod(b, b, mod);
        }
        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, mod);
        mpz_mul(t, t, c); mpz_mod(t, t, mod);
        mpz_mul(R, R, b); mpz_mod(R, R, mod);
    }
}

/* ==================== Factor Base ==================== */

typedef struct {
    unsigned int *prime;
    unsigned int *root1; /* sqrt(N) mod p */
    unsigned int *root2; /* p - root1 */
    unsigned char *logp;
    int size;
    int alloc;
} fb_t;

static fb_t *fb_create(mpz_t N, int target) {
    fb_t *fb = calloc(1, sizeof(fb_t));
    fb->alloc = target + 10;
    fb->prime = malloc(fb->alloc * sizeof(unsigned int));
    fb->root1 = malloc(fb->alloc * sizeof(unsigned int));
    fb->root2 = malloc(fb->alloc * sizeof(unsigned int));
    fb->logp  = malloc(fb->alloc * sizeof(unsigned char));

    /* Entry 0: -1 (sign) */
    fb->prime[0] = 1; fb->root1[0] = fb->root2[0] = 0; fb->logp[0] = 0;
    fb->size = 1;

    /* Entry 1: 2 */
    fb->prime[1] = 2; fb->root1[1] = fb->root2[1] = 1; fb->logp[1] = 1;
    fb->size = 2;

    /* Sieve for primes up to some bound */
    int bound = target * 15 + 10000;
    char *sieve = calloc(bound, 1);
    for (int i = 2; (long)i * i < bound; i++)
        if (!sieve[i]) for (int j = i*i; j < bound; j += i) sieve[j] = 1;

    for (int p = 3; p < bound && fb->size < target; p += 2) {
        if (sieve[p]) continue;
        unsigned long n_mod_p = mpz_fdiv_ui(N, p);
        unsigned int r = tonelli_shanks(n_mod_p, p);
        if (r == 0) continue; /* N is not QR mod p */

        fb->prime[fb->size] = p;
        fb->root1[fb->size] = r;
        fb->root2[fb->size] = p - r;
        fb->logp[fb->size]  = (unsigned char)(log2(p) + 0.5);
        fb->size++;
    }
    free(sieve);
    return fb;
}

/* ==================== Relations ==================== */

typedef struct {
    mpz_t ax_b;          /* (ax+b) value for sqrt step */
    unsigned char *exps;  /* exponent vector mod 2, packed */
    unsigned long lp;     /* 0 = full, else large prime */
} rel_t;

/* ==================== Large prime hash table ==================== */

#define LP_HASH_SIZE 1048576
#define LP_HASH_MASK (LP_HASH_SIZE - 1)

typedef struct lp_entry {
    unsigned long lp;
    int rel_idx;  /* index into partial relations array */
    struct lp_entry *next;
} lp_entry_t;

static lp_entry_t **lp_hash;
static lp_entry_t *lp_pool;
static int lp_pool_used;

static void lp_init(int max_entries) {
    lp_hash = calloc(LP_HASH_SIZE, sizeof(lp_entry_t*));
    lp_pool = malloc(max_entries * sizeof(lp_entry_t));
    lp_pool_used = 0;
}

static int lp_find(unsigned long lp) {
    unsigned int h = (unsigned int)(lp * 2654435761UL) & LP_HASH_MASK;
    for (lp_entry_t *e = lp_hash[h]; e; e = e->next)
        if (e->lp == lp) return e->rel_idx;
    return -1;
}

static void lp_insert(unsigned long lp, int idx) {
    unsigned int h = (unsigned int)(lp * 2654435761UL) & LP_HASH_MASK;
    lp_entry_t *e = &lp_pool[lp_pool_used++];
    e->lp = lp;
    e->rel_idx = idx;
    e->next = lp_hash[h];
    lp_hash[h] = e;
}

/* ==================== Gaussian Elimination ==================== */

typedef unsigned long long u64;

/*
 * Row-based GF(2) Gaussian elimination.
 * Matrix: nrels rows x fb_size columns (each row = relation's exponent vector mod 2).
 * Augmented with identity to track which relations combine.
 * After reduction, rows with zero left-part give dependencies.
 */
static int find_dependencies(unsigned char **rel_exps, int nrels, int fb_size,
                              int ***deps, int **deps_len, int max_deps) {
    int nrows = nrels;
    int ncols = fb_size;
    int fb_words = (ncols + 63) / 64;
    int id_words = (nrows + 63) / 64;
    int total_words = fb_words + id_words;

    /* Each row: [fb_size bits | nrels identity bits] */
    u64 **mat = malloc(nrows * sizeof(u64*));
    for (int r = 0; r < nrows; r++) {
        mat[r] = calloc(total_words, sizeof(u64));
        /* Fill FB exponents */
        for (int c = 0; c < ncols; c++)
            if (rel_exps[r][c] & 1)
                mat[r][c / 64] |= (1ULL << (c % 64));
        /* Identity bit */
        int id_bit = ncols + r; /* bit position in the augmented part */
        /* But we store identity in words starting at fb_words */
        mat[r][fb_words + r / 64] |= (1ULL << (r % 64));
    }

    /* Row-reduce the left part (FB columns) */
    int *pivot_row = malloc(ncols * sizeof(int));
    for (int c = 0; c < ncols; c++) pivot_row[c] = -1;

    for (int col = 0; col < ncols; col++) {
        /* Find pivot: first row with bit set in this column (among non-pivoted rows) */
        int pr = -1;
        for (int r = 0; r < nrows; r++) {
            if ((mat[r][col / 64] >> (col % 64)) & 1) {
                /* Check this row hasn't been used as a pivot for an earlier column */
                int used = 0;
                for (int cc = 0; cc < col; cc++)
                    if (pivot_row[cc] == r) { used = 1; break; }
                if (!used) { pr = r; break; }
            }
        }
        if (pr < 0) continue;
        pivot_row[col] = pr;

        /* Eliminate: XOR pivot row into all other rows that have this column set */
        for (int r = 0; r < nrows; r++) {
            if (r == pr) continue;
            if ((mat[r][col / 64] >> (col % 64)) & 1) {
                for (int w = 0; w < total_words; w++)
                    mat[r][w] ^= mat[pr][w];
            }
        }
    }

    /* Find rows where the left part (FB bits) is all zero */
    int ndeps = 0;
    *deps = malloc(max_deps * sizeof(int*));
    *deps_len = malloc(max_deps * sizeof(int));

    for (int r = 0; r < nrows && ndeps < max_deps; r++) {
        /* Check if FB part is all zero */
        int zero = 1;
        for (int w = 0; w < fb_words && zero; w++) {
            u64 mask = ~0ULL;
            if (w == fb_words - 1 && ncols % 64 != 0)
                mask = (1ULL << (ncols % 64)) - 1;
            if (mat[r][w] & mask) zero = 0;
        }
        if (!zero) continue;

        /* Extract dependency from identity part */
        int *dep = malloc(nrows * sizeof(int));
        int dep_len = 0;
        for (int w = 0; w < id_words; w++) {
            u64 bits = mat[r][fb_words + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w * 64 + bit;
                if (idx < nrows) dep[dep_len++] = idx;
                bits &= bits - 1;
            }
        }
        if (dep_len > 0) {
            (*deps)[ndeps] = dep;
            (*deps_len)[ndeps] = dep_len;
            ndeps++;
        } else {
            free(dep);
        }
    }

    for (int r = 0; r < nrows; r++) free(mat[r]);
    free(mat);
    free(pivot_row);
    return ndeps;
}

/* ==================== Main ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, argv[1], 10);

    /* Check for perfect square */
    if (mpz_perfect_square_p(N)) {
        mpz_t sq; mpz_init(sq);
        mpz_sqrt(sq, N);
        gmp_printf("%Zd\n", sq);
        mpz_clears(N, sq, NULL);
        return 0;
    }

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);
    params_t P = get_params(digits);

    fprintf(stderr, "SIQS: %d digits (%d bits), FB=%d, M=%d, a_factors=%d\n",
            digits, bits, P.fb_size, P.sieve_size, P.num_a_factors);

    /* Multiplier selection - try small multipliers to increase smoothness */
    /* TODO: implement Knuth-Schroeppel multiplier selection */

    /* Generate factor base */
    fb_t *fb = fb_create(N, P.fb_size);
    fprintf(stderr, "FB: %d primes, largest=%u\n", fb->size, fb->prime[fb->size-1]);

    int M = P.sieve_size;
    int target_rels = fb->size + 64;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;

    /* Sieve threshold:
     * We sieve g(x) = ((ax+b)^2 - N) / a = a*x^2 + 2*b*x + c
     * Max |g(x)| ≈ a*M^2 ≈ sqrt(2N)*M
     * So log2(|g(x)|) ≈ bits/2 + 1 + log2(M)
     */
    int expected_logG = bits / 2 + 1 + (int)(log2(M));
    int threshold = (int)(expected_logG * P.thresh_adj);

    fprintf(stderr, "Target rels=%d, threshold=%d, LP bound=%lu\n",
            target_rels, threshold, lp_bound);

    /* Allocate sieve */
    unsigned char *sieve = malloc(2 * M);

    /* Relations: full relations stored with exponent vectors */
    int max_rels = target_rels + 1000;
    int max_partials = 200000;
    int total_alloc = max_rels + max_partials;

    mpz_t *rel_ax_b = malloc(total_alloc * sizeof(mpz_t));
    unsigned char **rel_exps = malloc(total_alloc * sizeof(unsigned char*));
    /* For sqrt step we need full (not mod 2) exponents. */
    int **rel_full = malloc(total_alloc * sizeof(int*));
    /* Track large primes for combined relations (for sqrt step) */
    unsigned long *rel_lp = calloc(total_alloc, sizeof(unsigned long));
    for (int i = 0; i < total_alloc; i++) {
        mpz_init(rel_ax_b[i]);
        rel_exps[i] = calloc(fb->size, 1);
        rel_full[i] = calloc(fb->size, sizeof(int));
    }

    int num_full = 0;
    int num_partial_stored = 0;
    int partial_base = max_rels; /* partials stored starting here */

    lp_init(max_partials);

    /* Polynomial variables */
    mpz_t poly_a, poly_b;
    mpz_inits(poly_a, poly_b, NULL);

    /* target_a = sqrt(2N) / M */
    mpz_t target_a, sqrt2N;
    mpz_inits(target_a, sqrt2N, NULL);
    mpz_mul_ui(sqrt2N, N, 2);
    mpz_sqrt(target_a, sqrt2N);
    mpz_tdiv_q_ui(target_a, target_a, M);

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, time(NULL) ^ getpid());
    srand(time(NULL) ^ getpid());

    /* Precompute sieve solutions workspace */
    unsigned int *soln1 = malloc(fb->size * sizeof(unsigned int));
    unsigned int *soln2 = malloc(fb->size * sizeof(unsigned int));

    mpz_t ax_b, Qx, residue, tmp, tmp2;
    mpz_inits(ax_b, Qx, residue, tmp, tmp2, NULL);

    int poly_count = 0;

    while (num_full < target_rels) {
        /* Generate polynomial a = product of num_a_factors primes from middle of FB */
        int s = P.num_a_factors;
        int lo = fb->size / 4;
        int hi = 3 * fb->size / 4;
        if (lo < 2) lo = 2;
        if (hi <= lo) hi = fb->size - 1;

        int a_idx[20]; /* indices of primes making up a */
        mpz_set_ui(poly_a, 1);
        for (int i = 0; i < s; i++) {
            int idx, ok;
            do {
                idx = lo + (rand() % (hi - lo));
                ok = 1;
                for (int j = 0; j < i; j++) if (a_idx[j] == idx) { ok = 0; break; }
            } while (!ok);
            a_idx[i] = idx;
            mpz_mul_ui(poly_a, poly_a, fb->prime[idx]);
        }

        /* Compute b via CRT: b^2 ≡ N (mod a) */
        mpz_set_ui(poly_b, 0);
        int crt_ok = 1;
        for (int i = 0; i < s; i++) {
            unsigned int p = fb->prime[a_idx[i]];
            unsigned int r = fb->root1[a_idx[i]];
            mpz_t ai, ai_inv, mod_p;
            mpz_inits(ai, ai_inv, mod_p, NULL);
            mpz_divexact_ui(ai, poly_a, p);
            mpz_set_ui(mod_p, p);
            if (!mpz_invert(ai_inv, ai, mod_p)) { crt_ok = 0; mpz_clears(ai, ai_inv, mod_p, NULL); break; }
            mpz_mul(ai_inv, ai_inv, ai);
            mpz_mul_ui(ai_inv, ai_inv, r);
            mpz_add(poly_b, poly_b, ai_inv);
            mpz_clears(ai, ai_inv, mod_p, NULL);
        }
        if (!crt_ok) { poly_count++; continue; }
        mpz_mod(poly_b, poly_b, poly_a);

        /* Verify b^2 ≡ N mod a */
        mpz_mul(tmp, poly_b, poly_b);
        mpz_sub(tmp, tmp, N);
        mpz_mod(tmp, tmp, poly_a);
        if (mpz_sgn(tmp) != 0) {
            mpz_sub(poly_b, poly_a, poly_b);
            mpz_mul(tmp, poly_b, poly_b);
            mpz_sub(tmp, tmp, N);
            mpz_mod(tmp, tmp, poly_a);
            if (mpz_sgn(tmp) != 0) { poly_count++; continue; }
        }

        /* Compute sieve solutions for each FB prime */
        for (int i = 2; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            unsigned long a_mod = mpz_fdiv_ui(poly_a, p);
            if (a_mod == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }

            /* a_inv = modular inverse of a mod p */
            mpz_set_ui(tmp, a_mod);
            mpz_set_ui(tmp2, p);
            if (!mpz_invert(tmp, tmp, tmp2)) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned long a_inv = mpz_get_ui(tmp);

            unsigned long b_mod = mpz_fdiv_ui(poly_b, p);

            /* x = a_inv * (root - b) mod p, shifted by M */
            unsigned long r1 = fb->root1[i];
            unsigned long r2 = fb->root2[i];

            unsigned long x1 = (a_inv * ((r1 + p - b_mod) % p)) % p;
            unsigned long x2 = (a_inv * ((r2 + p - b_mod) % p)) % p;

            /* Offset into sieve array [0, 2M) */
            soln1[i] = (unsigned int)((x1 + M) % p);
            soln2[i] = (unsigned int)((x2 + M) % p);
        }

        /* Sieve */
        memset(sieve, 0, 2 * M);

        /* Skip small primes for sieving (they contribute little per hit) */
        for (int i = 2; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            unsigned char lp = fb->logp[i];
            if (soln1[i] == 0xFFFFFFFF) continue;

            /* Sieve with first root */
            for (unsigned int j = soln1[i]; j < (unsigned int)(2 * M); j += p)
                sieve[j] += lp;

            /* Sieve with second root (if different) */
            if (soln1[i] != soln2[i]) {
                for (unsigned int j = soln2[i]; j < (unsigned int)(2 * M); j += p)
                    sieve[j] += lp;
            }
        }

        /* Scan for smooth candidates */
        for (int j = 0; j < 2 * M && num_full < target_rels; j++) {
            if (sieve[j] < threshold) continue;

            long x = (long)j - M;
            if (x == 0) continue;

            /* Compute ax+b (for sqrt step) and g(x) = ((ax+b)^2 - N)/a (for smoothness) */
            mpz_mul_si(ax_b, poly_a, x);
            mpz_add(ax_b, ax_b, poly_b);
            mpz_mul(Qx, ax_b, ax_b);
            mpz_sub(Qx, Qx, N);
            /* g(x) = Qx / a - this should be exact since b^2 ≡ N (mod a) */
            mpz_divexact(residue, Qx, poly_a);

            if (mpz_sgn(residue) == 0) continue;

            /* Trial divide g(x) by FB.
             * Relation: (ax+b)^2 = a * g(x) + N => (ax+b)^2 ≡ a * g(x) (mod N)
             * So exponent vector includes: a-factor primes (from 'a') + g(x) factors
             */
            int neg = (mpz_sgn(residue) < 0);
            mpz_abs(residue, residue);

            int store_idx = num_full;
            unsigned char *exps = rel_exps[store_idx];
            int *fexp = rel_full[store_idx];
            memset(exps, 0, fb->size);
            memset(fexp, 0, fb->size * sizeof(int));
            if (neg) { exps[0] = 1; fexp[0] = 1; }

            /* Include the 'a' factor primes in the exponent vector.
             * a = product of a_idx primes, each with exponent 1 */
            for (int k = 0; k < s; k++) {
                int fi = a_idx[k];
                fexp[fi] += 1;
                exps[fi] = fexp[fi] & 1;
            }

            /* Divide g(x) by FB primes */
            for (int i = 1; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                if (mpz_divisible_ui_p(residue, p)) {
                    do {
                        mpz_divexact_ui(residue, residue, p);
                        fexp[i]++;
                    } while (mpz_divisible_ui_p(residue, p));
                    exps[i] = fexp[i] & 1;
                }
            }

            if (mpz_cmp_ui(residue, 1) == 0) {
                /* Full relation */
                mpz_set(rel_ax_b[num_full], ax_b);
                num_full++;
            } else if (mpz_cmp_ui(residue, lp_bound) <= 0 && mpz_fits_ulong_p(residue)) {
                unsigned long lp = mpz_get_ui(residue);
                /* Check if we have a matching partial */
                int match = lp_find(lp);
                if (match >= 0) {
                    /* Combine partials: add exponent vectors, multiply ax_b values */
                    int pi = match;
                    for (int k = 0; k < fb->size; k++) {
                        fexp[k] += rel_full[pi][k];
                        exps[k] = fexp[k] & 1;
                    }
                    /* lp^2 is a perfect square → even exponent in mod 2 vector.
                     * But for sqrt step, we need lp as factor in Y. Track it. */
                    mpz_mul(rel_ax_b[num_full], ax_b, rel_ax_b[pi]);
                    rel_lp[num_full] = lp; /* store LP for sqrt step */
                    num_full++;
                } else {
                    /* Store as partial */
                    int pi = partial_base + num_partial_stored;
                    if (pi < total_alloc) {
                        mpz_set(rel_ax_b[pi], ax_b);
                        memcpy(rel_exps[pi], exps, fb->size);
                        memcpy(rel_full[pi], fexp, fb->size * sizeof(int));
                        lp_insert(lp, pi);
                        num_partial_stored++;
                    }
                }
            }
        }

        poly_count++;
        if (poly_count % 200 == 0) {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - t0.tv_sec) + (now.tv_nsec - t0.tv_nsec) / 1e9;
            double rels_per_sec = num_full / (elapsed > 0 ? elapsed : 1);
            int remaining = target_rels - num_full;
            fprintf(stderr, "  poly=%d rels=%d/%d (partials=%d) %.1f rels/s ETA=%.0fs elapsed=%.1fs\n",
                    poly_count, num_full, target_rels, num_partial_stored,
                    rels_per_sec, remaining / (rels_per_sec > 0 ? rels_per_sec : 1), elapsed);
            if (elapsed > 280) {
                fprintf(stderr, "TIMEOUT approaching\n");
                break;
            }
        }
    }

    struct timespec t_sieve;
    clock_gettime(CLOCK_MONOTONIC, &t_sieve);
    double sieve_time = (t_sieve.tv_sec - t0.tv_sec) + (t_sieve.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Sieving done: %d rels from %d polys in %.1fs\n", num_full, poly_count, sieve_time);

    if (num_full < fb->size + 1) {
        fprintf(stderr, "Not enough relations (%d < %d)\n", num_full, fb->size + 1);
        return 1;
    }

    /* Gaussian elimination */
    fprintf(stderr, "Linear algebra: %d x %d matrix\n", fb->size, num_full);

    int **deps;
    int *deps_len;
    int ndeps = find_dependencies(rel_exps, num_full, fb->size, &deps, &deps_len, 64);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* Try each dependency */
    for (int d = 0; d < ndeps; d++) {
        /* Compute X = product of (ax+b) mod N */
        /* Compute Y = sqrt(product of Q(x)) using exponents */
        mpz_t X, Y, g;
        mpz_inits(X, Y, g, NULL);
        mpz_set_ui(X, 1);
        mpz_set_ui(Y, 1);

        int *total_exp = calloc(fb->size, sizeof(int));
        for (int k = 0; k < deps_len[d]; k++) {
            int ri = deps[d][k];
            mpz_mul(X, X, rel_ax_b[ri]);
            mpz_mod(X, X, N);
            for (int f = 0; f < fb->size; f++)
                total_exp[f] += rel_full[ri][f];
        }

        /* Y = product of p^(total_exp[p]/2) mod N */
        int valid = 1;
        for (int f = 0; f < fb->size; f++) {
            if (total_exp[f] & 1) { valid = 0; break; }
        }
        if (!valid) { free(total_exp); mpz_clears(X, Y, g, NULL); continue; }

        for (int f = 1; f < fb->size; f++) {
            int e = total_exp[f] / 2;
            if (e > 0) {
                mpz_set_ui(tmp, fb->prime[f]);
                mpz_powm_ui(tmp, tmp, e, N);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }
        }
        /* Include large primes from combined partials */
        for (int k = 0; k < deps_len[d]; k++) {
            int ri = deps[d][k];
            if (rel_lp[ri] > 0) {
                mpz_set_ui(tmp, rel_lp[ri]);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }
        }

        /* Verify: X^2 ≡ Y^2 (mod N) */
        mpz_t X2, Y2;
        mpz_inits(X2, Y2, NULL);
        mpz_mul(X2, X, X); mpz_mod(X2, X2, N);
        mpz_mul(Y2, Y, Y); mpz_mod(Y2, Y2, N);
        if (mpz_cmp(X2, Y2) != 0) {
            fprintf(stderr, "  dep %d: X^2 != Y^2 mod N (VERIFICATION FAILED)\n", d);
            mpz_clears(X2, Y2, NULL);
            free(total_exp); mpz_clears(X, Y, g, NULL); continue;
        }
        mpz_clears(X2, Y2, NULL);

        /* gcd(X - Y, N) */
        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other);
            mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stderr, "SIQS done: %.3fs (dep %d/%d)\n",
                    (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9, d+1, ndeps);
            return 0;
        }

        mpz_add(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other);
            mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stderr, "SIQS done: %.3fs (dep %d/%d)\n",
                    (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9, d+1, ndeps);
            return 0;
        }

        free(total_exp);
        mpz_clears(X, Y, g, NULL);
    }

    fprintf(stderr, "SIQS FAILED: no dependency produced a factor\n");
    return 1;
}
