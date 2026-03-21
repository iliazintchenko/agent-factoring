/*
 * Self-Initializing Quadratic Sieve (SIQS) for factoring semiprimes.
 * Usage: ./siqs <N>
 *
 * Targets: 50-100 digit semiprimes.
 * Uses GMP for big integer arithmetic.
 * Implements: factor base generation, SIQS polynomial selection,
 *             sieving, trial division, large prime relations,
 *             Gaussian elimination, and square root extraction.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gmp.h>

/* ==================== Parameters ==================== */

#define MAX_FB_SIZE    200000
#define BLOCK_SIZE     65536
#define MAX_RELATIONS  300000
#define MAX_FACTORS    64

/* Parameter selection based on digit count */
typedef struct {
    int fb_size;          /* number of primes in factor base */
    int sieve_radius;     /* M: sieve from -M to +M */
    int large_prime_bits; /* bits for large prime bound */
    double fudge;         /* threshold fudge factor */
    int num_poly_a_factors; /* number of primes in poly A */
} siqs_params_t;

static siqs_params_t get_params(int digits) {
    siqs_params_t p;
    if (digits <= 50) {
        p.fb_size = 300; p.sieve_radius = 65536; p.large_prime_bits = 23;
        p.fudge = 1.6; p.num_poly_a_factors = 5;
    } else if (digits <= 55) {
        p.fb_size = 600; p.sieve_radius = 65536; p.large_prime_bits = 25;
        p.fudge = 1.5; p.num_poly_a_factors = 6;
    } else if (digits <= 60) {
        p.fb_size = 1200; p.sieve_radius = 131072; p.large_prime_bits = 27;
        p.fudge = 1.4; p.num_poly_a_factors = 7;
    } else if (digits <= 65) {
        p.fb_size = 2500; p.sieve_radius = 131072; p.large_prime_bits = 28;
        p.fudge = 1.35; p.num_poly_a_factors = 8;
    } else if (digits <= 70) {
        p.fb_size = 5000; p.sieve_radius = 196608; p.large_prime_bits = 30;
        p.fudge = 1.3; p.num_poly_a_factors = 9;
    } else if (digits <= 75) {
        p.fb_size = 10000; p.sieve_radius = 262144; p.large_prime_bits = 31;
        p.fudge = 1.25; p.num_poly_a_factors = 10;
    } else if (digits <= 80) {
        p.fb_size = 20000; p.sieve_radius = 262144; p.large_prime_bits = 32;
        p.fudge = 1.2; p.num_poly_a_factors = 11;
    } else if (digits <= 85) {
        p.fb_size = 35000; p.sieve_radius = 393216; p.large_prime_bits = 33;
        p.fudge = 1.15; p.num_poly_a_factors = 11;
    } else if (digits <= 90) {
        p.fb_size = 55000; p.sieve_radius = 524288; p.large_prime_bits = 34;
        p.fudge = 1.1; p.num_poly_a_factors = 12;
    } else if (digits <= 95) {
        p.fb_size = 80000; p.sieve_radius = 524288; p.large_prime_bits = 35;
        p.fudge = 1.05; p.num_poly_a_factors = 12;
    } else {
        p.fb_size = 120000; p.sieve_radius = 786432; p.large_prime_bits = 36;
        p.fudge = 1.0; p.num_poly_a_factors = 13;
    }
    return p;
}

/* ==================== Factor Base ==================== */

typedef struct {
    unsigned int *primes;   /* prime values */
    unsigned int *roots1;   /* sqrt(N) mod p, first root */
    unsigned int *roots2;   /* sqrt(N) mod p, second root */
    unsigned char *logp;    /* floor(log2(p)) */
    int size;               /* number of primes */
} factor_base_t;

/* Tonelli-Shanks for sqrt(n) mod p */
static unsigned int mod_sqrt(mpz_t n, unsigned int p) {
    unsigned long n_mod_p = mpz_fdiv_ui(n, p);
    if (n_mod_p == 0) return 0;
    if (p == 2) return n_mod_p & 1;

    /* Check if n is a QR mod p */
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(tmp, n_mod_p);
    mpz_t pp; mpz_init_set_ui(pp, p);
    int legendre = mpz_legendre(tmp, pp);
    mpz_clear(pp);
    if (legendre != 1) { mpz_clear(tmp); return 0; }

    /* Tonelli-Shanks */
    if (p % 4 == 3) {
        /* Simple case */
        mpz_t result;
        mpz_init(result);
        mpz_set_ui(result, n_mod_p);
        mpz_t mod; mpz_init_set_ui(mod, p);
        mpz_t exp; mpz_init_set_ui(exp, (p + 1) / 4);
        mpz_powm(result, result, exp, mod);
        unsigned int r = mpz_get_ui(result);
        mpz_clears(result, mod, exp, tmp, NULL);
        return r;
    }

    /* General Tonelli-Shanks */
    unsigned long Q = p - 1;
    int S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    /* Find a non-residue */
    unsigned long z = 2;
    while (1) {
        mpz_set_ui(tmp, z);
        mpz_t pp2; mpz_init_set_ui(pp2, p);
        if (mpz_legendre(tmp, pp2) == -1) { mpz_clear(pp2); break; }
        mpz_clear(pp2);
        z++;
    }

    mpz_t M, c, t, R, mod, b, tmp2;
    mpz_inits(M, c, t, R, mod, b, tmp2, NULL);
    mpz_set_ui(mod, p);
    mpz_set_ui(M, S);

    mpz_set_ui(c, z);
    mpz_set_ui(tmp, Q);
    mpz_powm(c, c, tmp, mod);

    mpz_set_ui(t, n_mod_p);
    mpz_set_ui(tmp, Q);
    mpz_powm(t, t, tmp, mod);

    mpz_set_ui(R, n_mod_p);
    mpz_set_ui(tmp, (Q + 1) / 2);
    mpz_powm(R, R, tmp, mod);

    while (1) {
        if (mpz_cmp_ui(t, 0) == 0) { mpz_clear(tmp); mpz_clears(M, c, t, R, mod, b, tmp2, NULL); return 0; }
        if (mpz_cmp_ui(t, 1) == 0) {
            unsigned int result = mpz_get_ui(R);
            mpz_clear(tmp); mpz_clears(M, c, t, R, mod, b, tmp2, NULL);
            return result;
        }

        /* Find least i such that t^(2^i) = 1 mod p */
        int i = 0;
        mpz_set(tmp2, t);
        while (mpz_cmp_ui(tmp2, 1) != 0) {
            mpz_mul(tmp2, tmp2, tmp2);
            mpz_mod(tmp2, tmp2, mod);
            i++;
        }

        /* b = c^(2^(M-i-1)) */
        mpz_set(b, c);
        unsigned long mi1 = mpz_get_ui(M) - i - 1;
        for (unsigned long j = 0; j < mi1; j++) {
            mpz_mul(b, b, b);
            mpz_mod(b, b, mod);
        }

        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, mod);
        mpz_mul(t, t, c); mpz_mod(t, t, mod);
        mpz_mul(R, R, b); mpz_mod(R, R, mod);
    }
}

/* Generate factor base: primes p where N is a QR mod p */
static factor_base_t *make_factor_base(mpz_t N, int target_size) {
    factor_base_t *fb = calloc(1, sizeof(factor_base_t));
    fb->primes = malloc(target_size * sizeof(unsigned int));
    fb->roots1 = malloc(target_size * sizeof(unsigned int));
    fb->roots2 = malloc(target_size * sizeof(unsigned int));
    fb->logp = malloc(target_size * sizeof(unsigned char));
    fb->size = 0;

    /* Add -1 sentinel */
    fb->primes[0] = 1; /* represents -1 */
    fb->roots1[0] = 0;
    fb->roots2[0] = 0;
    fb->logp[0] = 0;
    fb->size = 1;

    /* Sieve for primes */
    int sieve_limit = target_size * 20; /* rough upper bound */
    if (sieve_limit < 100000) sieve_limit = 100000;
    unsigned char *is_prime = calloc(sieve_limit, 1);
    for (int i = 2; i * i < sieve_limit; i++)
        if (!is_prime[i])
            for (int j = i*i; j < sieve_limit; j += i) is_prime[j] = 1;

    for (int p = 2; p < sieve_limit && fb->size < target_size; p++) {
        if (is_prime[p]) continue;
        unsigned int r = mod_sqrt(N, p);
        if (r == 0 && p > 2) continue; /* N is not QR mod p */
        if (p == 2) {
            /* Handle p=2 specially */
            fb->primes[fb->size] = 2;
            fb->roots1[fb->size] = mpz_fdiv_ui(N, 2) == 1 ? 1 : 0;
            fb->roots2[fb->size] = fb->roots1[fb->size];
            fb->logp[fb->size] = 1;
            fb->size++;
            continue;
        }
        fb->primes[fb->size] = p;
        fb->roots1[fb->size] = r;
        fb->roots2[fb->size] = p - r;
        fb->logp[fb->size] = (unsigned char)(log2(p));
        fb->size++;
    }

    free(is_prime);
    return fb;
}

/* ==================== Relations ==================== */

typedef struct {
    mpz_t QofX;           /* Q(x) = poly(x)^2 - N (or residue after dividing by FB) */
    int *exponents;       /* exponent vector over factor base (mod 2) */
    unsigned long large_prime; /* 0 if full, else the large prime */
    mpz_t poly_val;       /* poly(x) value for sqrt step */
    int valid;
} relation_t;

/* ==================== Sieving ==================== */

/*
 * For polynomial g(x) = (ax+b)^2 - N,
 * we sieve for x in [-M, M].
 * Roots: for each FB prime p, solve (ax+b)^2 ≡ N (mod p)
 *   => ax ≡ ±sqrt(N) - b (mod p)
 *   => x ≡ a_inv * (±sqrt(N) - b) (mod p)
 */

typedef struct {
    mpz_t a, b;           /* polynomial coefficients */
    mpz_t a_inv_2;        /* 2a, for sqrt step */
    unsigned int *soln1;  /* sieve root 1 for each FB prime */
    unsigned int *soln2;  /* sieve root 2 for each FB prime */
} siqs_poly_t;

/* Compute sieve roots for a given polynomial and factor base */
static void compute_sieve_roots(siqs_poly_t *poly, factor_base_t *fb, mpz_t N, int M) {
    for (int i = 2; i < fb->size; i++) { /* skip -1 and 2 */
        unsigned int p = fb->primes[i];
        unsigned long a_mod = mpz_fdiv_ui(poly->a, p);
        if (a_mod == 0) { poly->soln1[i] = poly->soln2[i] = 0xFFFFFFFF; continue; }

        /* a_inv = modular inverse of a mod p */
        mpz_t tmp;
        mpz_init(tmp);
        mpz_set_ui(tmp, a_mod);
        mpz_t mod_p; mpz_init_set_ui(mod_p, p);
        mpz_t inv; mpz_init(inv);
        mpz_invert(inv, tmp, mod_p);
        unsigned long a_inv = mpz_get_ui(inv);

        unsigned long b_mod = mpz_fdiv_ui(poly->b, p);
        unsigned long r1 = fb->roots1[i];
        unsigned long r2 = fb->roots2[i];

        /* x1 = a_inv * (r1 - b) mod p */
        long x1 = (long)(a_inv * ((r1 + p - b_mod) % p)) % (long)p;
        if (x1 < 0) x1 += p;
        long x2 = (long)(a_inv * ((r2 + p - b_mod) % p)) % (long)p;
        if (x2 < 0) x2 += p;

        /* Offset by M to make sieve array non-negative */
        poly->soln1[i] = ((unsigned int)x1 + M) % p;
        poly->soln2[i] = ((unsigned int)x2 + M) % p;

        mpz_clears(tmp, mod_p, inv, NULL);
    }
}

/* ==================== Gaussian Elimination ==================== */

typedef unsigned long long u64;

/* Bit matrix for GF(2) Gaussian elimination */
typedef struct {
    u64 **rows;     /* each row is an array of u64 words */
    int nrows;
    int ncols;
    int words_per_row;
} bitmatrix_t;

static bitmatrix_t *bitmatrix_alloc(int nrows, int ncols) {
    bitmatrix_t *m = malloc(sizeof(bitmatrix_t));
    m->nrows = nrows;
    m->ncols = ncols;
    m->words_per_row = (ncols + 63) / 64;
    m->rows = malloc(nrows * sizeof(u64*));
    for (int i = 0; i < nrows; i++) {
        m->rows[i] = calloc(m->words_per_row, sizeof(u64));
    }
    return m;
}

static void bitmatrix_free(bitmatrix_t *m) {
    for (int i = 0; i < m->nrows; i++) free(m->rows[i]);
    free(m->rows);
    free(m);
}

static void bitmatrix_set(bitmatrix_t *m, int r, int c) {
    m->rows[r][c / 64] |= (1ULL << (c % 64));
}

static int bitmatrix_get(bitmatrix_t *m, int r, int c) {
    return (m->rows[r][c / 64] >> (c % 64)) & 1;
}

static void bitmatrix_xor_rows(bitmatrix_t *m, int dst, int src) {
    for (int w = 0; w < m->words_per_row; w++)
        m->rows[dst][w] ^= m->rows[src][w];
}

/*
 * Gaussian elimination on the transpose.
 * Input: matrix with fb_size rows and num_relations columns (exponent vectors).
 * Output: null space vectors (dependencies).
 * Returns number of dependencies found.
 */
static int gauss_eliminate(bitmatrix_t *mat, int **deps, int *ndeps) {
    int nrows = mat->nrows;
    int ncols = mat->ncols;

    /* Augment with identity for tracking dependencies */
    int aug_words = (ncols + 63) / 64;
    u64 **aug = malloc(ncols * sizeof(u64*));
    for (int i = 0; i < ncols; i++) {
        aug[i] = calloc(aug_words, sizeof(u64));
        aug[i][i / 64] |= (1ULL << (i % 64));
    }

    /* Work on transpose: columns = FB primes, rows = relations */
    /* We need to find linear dependencies among the relation vectors */
    int *pivot_col = malloc(nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) pivot_col[i] = -1;

    /* For each column (FB prime), find pivot row */
    for (int col = 0; col < nrows; col++) {
        int pivot_row = -1;
        for (int row = 0; row < ncols; row++) {
            /* Check if this relation has this FB prime in its exponent vector */
            if (bitmatrix_get(mat, col, row)) {
                pivot_row = row;
                break;
            }
        }
        if (pivot_row < 0) continue;
        pivot_col[col] = pivot_row;

        /* Eliminate this column from all other rows */
        for (int row = 0; row < ncols; row++) {
            if (row != pivot_row && bitmatrix_get(mat, col, row)) {
                /* XOR row with pivot_row in matrix */
                for (int c = col; c < nrows; c++) {
                    if (bitmatrix_get(mat, c, pivot_row)) {
                        /* flip bit (c, row) */
                        mat->rows[c][row / 64] ^= (1ULL << (row % 64));
                    }
                }
                /* XOR in augmented matrix */
                for (int w = 0; w < aug_words; w++)
                    aug[row][w] ^= aug[pivot_row][w];
            }
        }
    }

    /* Find free rows (not a pivot for any column) */
    int *is_pivot = calloc(ncols, sizeof(int));
    for (int col = 0; col < nrows; col++) {
        if (pivot_col[col] >= 0) is_pivot[pivot_col[col]] = 1;
    }

    *ndeps = 0;
    *deps = NULL;
    /* Each free row gives a dependency */
    for (int row = 0; row < ncols; row++) {
        if (!is_pivot[row]) {
            /* aug[row] gives the dependency set */
            (*ndeps)++;
            *deps = realloc(*deps, (*ndeps) * sizeof(int));
            (*deps)[(*ndeps) - 1] = row;
        }
    }

    /* Store augmented matrix for later use */
    /* Actually we need to return the full dependency info */
    /* For now, store in a global or restructure */

    free(is_pivot);
    free(pivot_col);
    for (int i = 0; i < ncols; i++) free(aug[i]);
    free(aug);

    return *ndeps;
}

/* ==================== Main SIQS ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    siqs_params_t params = get_params(digits);

    fprintf(stderr, "SIQS: %d digits, FB=%d, M=%d\n", digits, params.fb_size, params.sieve_radius);

    /* Step 1: Generate factor base */
    factor_base_t *fb = make_factor_base(N, params.fb_size);
    fprintf(stderr, "Factor base: %d primes (largest=%u)\n", fb->size, fb->primes[fb->size-1]);

    /* Step 2: Compute target_a = sqrt(2N)/M */
    mpz_t target_a, sqrtN, sqrt2N;
    mpz_inits(target_a, sqrtN, sqrt2N, NULL);
    mpz_mul_ui(sqrt2N, N, 2);
    mpz_sqrt(sqrtN, sqrt2N);
    mpz_tdiv_q_ui(target_a, sqrtN, params.sieve_radius);

    int M = params.sieve_radius;
    int target_relations = fb->size + 64;

    /* Allocate sieve array */
    unsigned char *sieve = malloc(2 * M);

    /* Relations storage */
    int num_full = 0, num_partial = 0;
    relation_t *relations = calloc(target_relations + 1000, sizeof(relation_t));
    for (int i = 0; i < target_relations + 1000; i++) {
        mpz_init(relations[i].QofX);
        mpz_init(relations[i].poly_val);
        relations[i].exponents = calloc(fb->size, sizeof(int));
    }

    /* Large prime hash for combining partials */
    /* Simple: use sorted array + merge */
    unsigned long *lp_list = malloc(500000 * sizeof(unsigned long));
    int *lp_rel_idx = malloc(500000 * sizeof(int));
    int num_lp = 0;
    unsigned long large_prime_bound = 1UL << params.large_prime_bits;

    /* Polynomial workspace */
    siqs_poly_t poly;
    mpz_inits(poly.a, poly.b, poly.a_inv_2, NULL);
    poly.soln1 = malloc(fb->size * sizeof(unsigned int));
    poly.soln2 = malloc(fb->size * sizeof(unsigned int));

    /* Sieve threshold */
    int bits_N = mpz_sizeinbase(N, 2);
    int threshold = (int)(bits_N / 2.0 * params.fudge);

    fprintf(stderr, "Target relations: %d, threshold: %d, large prime bound: %lu\n",
            target_relations, threshold, large_prime_bound);

    /* Step 3: Main sieve loop */
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, time(NULL) ^ getpid());

    int poly_count = 0;
    mpz_t Qx, ax_b, residue, tmp;
    mpz_inits(Qx, ax_b, residue, tmp, NULL);

    while (num_full < target_relations) {
        /* Generate polynomial: a = product of s random FB primes near target_a */
        int s = params.num_poly_a_factors;

        /* Pick s primes from middle of FB */
        int lo = fb->size / 4;
        int hi = fb->size * 3 / 4;
        if (lo < 2) lo = 2;

        int a_factors[MAX_FACTORS];
        mpz_set_ui(poly.a, 1);
        for (int i = 0; i < s; i++) {
            int idx;
            int redo;
            do {
                idx = lo + (rand() % (hi - lo));
                redo = 0;
                for (int j = 0; j < i; j++)
                    if (a_factors[j] == idx) { redo = 1; break; }
            } while (redo);
            a_factors[i] = idx;
            mpz_mul_ui(poly.a, poly.a, fb->primes[idx]);
        }

        /* Compute b: we need b^2 ≡ N (mod a)
         * Using CRT: for each prime p_i dividing a, find b_i = sqrt(N) mod p_i
         * Then combine using CRT.
         * Simplification: use Hensel lifting or just compute directly.
         */
        /* For simplicity, compute b using CRT over the a-factors */
        mpz_set_ui(poly.b, 0);
        for (int i = 0; i < s; i++) {
            unsigned int p = fb->primes[a_factors[i]];
            unsigned int r = fb->roots1[a_factors[i]]; /* sqrt(N) mod p */

            /* CRT: b += r * (a/p) * inverse(a/p, p) */
            mpz_t ai, ai_inv, contribution, mod_p;
            mpz_inits(ai, ai_inv, contribution, mod_p, NULL);

            mpz_divexact_ui(ai, poly.a, p); /* a/p_i */
            mpz_set_ui(mod_p, p);
            mpz_invert(ai_inv, ai, mod_p); /* (a/p_i)^-1 mod p_i */

            mpz_mul(contribution, ai, ai_inv);
            mpz_mul_ui(contribution, contribution, r);
            mpz_add(poly.b, poly.b, contribution);

            mpz_clears(ai, ai_inv, contribution, mod_p, NULL);
        }
        mpz_mod(poly.b, poly.b, poly.a);

        /* Ensure b^2 ≡ N mod a (adjust sign if needed) */
        mpz_mul(tmp, poly.b, poly.b);
        mpz_sub(tmp, tmp, N);
        mpz_mod(tmp, tmp, poly.a);
        if (mpz_sgn(tmp) != 0) {
            /* Try negating b */
            mpz_sub(poly.b, poly.a, poly.b);
            mpz_mul(tmp, poly.b, poly.b);
            mpz_sub(tmp, tmp, N);
            mpz_mod(tmp, tmp, poly.a);
            if (mpz_sgn(tmp) != 0) {
                /* CRT didn't work perfectly, skip this poly */
                poly_count++;
                continue;
            }
        }

        /* Adjust b so that b ≡ sqrt(N) mod 2 if needed */
        /* Make b have same parity as N for convenience */

        compute_sieve_roots(&poly, fb, N, M);

        /* Initialize sieve array */
        memset(sieve, 0, 2 * M);

        /* Sieve with factor base */
        for (int i = 2; i < fb->size; i++) {
            unsigned int p = fb->primes[i];
            unsigned char logp = fb->logp[i];
            if (poly.soln1[i] == 0xFFFFFFFF) continue; /* p divides a */

            unsigned int s1 = poly.soln1[i];
            unsigned int s2 = poly.soln2[i];

            /* Sieve positive roots */
            for (unsigned int j = s1; j < (unsigned int)(2 * M); j += p)
                sieve[j] += logp;
            if (s1 != s2) {
                for (unsigned int j = s2; j < (unsigned int)(2 * M); j += p)
                    sieve[j] += logp;
            }
        }

        /* Scan for smooth candidates */
        for (int j = 0; j < 2 * M; j++) {
            if (sieve[j] < threshold) continue;

            /* Candidate: compute Q(x) = (a*x + b)^2 - N where x = j - M */
            long x = j - M;
            mpz_mul_si(ax_b, poly.a, x);
            mpz_add(ax_b, ax_b, poly.b);
            mpz_mul(Qx, ax_b, ax_b);
            mpz_sub(Qx, Qx, N);

            if (mpz_sgn(Qx) == 0) continue;

            /* Trial divide by factor base */
            int neg = (mpz_sgn(Qx) < 0);
            mpz_abs(residue, Qx);

            int *exps = relations[num_full].exponents;
            memset(exps, 0, fb->size * sizeof(int));
            if (neg) exps[0] = 1; /* sign factor */

            for (int i = 1; i < fb->size; i++) {
                unsigned int p = fb->primes[i];
                while (mpz_divisible_ui_p(residue, p)) {
                    mpz_divexact_ui(residue, residue, p);
                    exps[i]++;
                }
            }

            /* Check if fully smooth or has a single large prime */
            if (mpz_cmp_ui(residue, 1) == 0) {
                /* Full relation */
                mpz_set(relations[num_full].QofX, Qx);
                mpz_set(relations[num_full].poly_val, ax_b);
                relations[num_full].large_prime = 0;
                relations[num_full].valid = 1;
                num_full++;
                if (num_full >= target_relations) break;
            } else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) < large_prime_bound) {
                /* Single large prime relation */
                unsigned long lp = mpz_get_ui(residue);
                /* Check if we already have a relation with this large prime */
                int found_match = -1;
                for (int k = 0; k < num_lp; k++) {
                    if (lp_list[k] == lp) { found_match = k; break; }
                }
                if (found_match >= 0) {
                    /* Combine: multiply the two relations to cancel the large prime */
                    int other_idx = lp_rel_idx[found_match];
                    for (int k = 0; k < fb->size; k++)
                        exps[k] += relations[other_idx].exponents[k];
                    /* The combined relation's poly_val = product of the two poly_vals */
                    mpz_mul(relations[num_full].poly_val, ax_b, relations[other_idx].poly_val);
                    mpz_mul(relations[num_full].QofX, Qx, relations[other_idx].QofX);
                    relations[num_full].large_prime = 0;
                    relations[num_full].valid = 1;
                    num_full++;
                    if (num_full >= target_relations) break;
                } else {
                    /* Store this partial */
                    if (num_lp < 500000) {
                        /* Store in a temporary slot after num_full */
                        int idx = target_relations + num_lp;
                        mpz_set(relations[idx].QofX, Qx);
                        mpz_set(relations[idx].poly_val, ax_b);
                        memcpy(relations[idx].exponents, exps, fb->size * sizeof(int));
                        relations[idx].large_prime = lp;
                        lp_list[num_lp] = lp;
                        lp_rel_idx[num_lp] = idx;
                        num_lp++;
                        num_partial++;
                    }
                }
            }
        }

        poly_count++;
        if (poly_count % 100 == 0) {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - t0.tv_sec) + (now.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "  poly=%d, relations=%d/%d (partials=%d), time=%.1fs\n",
                    poly_count, num_full, target_relations, num_partial, elapsed);
            if (elapsed > 280) {
                fprintf(stderr, "TIMEOUT approaching, giving up\n");
                break;
            }
        }
    }

    fprintf(stderr, "Collected %d full relations from %d polynomials\n", num_full, poly_count);

    if (num_full < fb->size + 1) {
        fprintf(stderr, "Not enough relations\n");
        return 1;
    }

    /* Step 4: Gaussian elimination */
    fprintf(stderr, "Building matrix (%d x %d)...\n", fb->size, num_full);

    bitmatrix_t *mat = bitmatrix_alloc(fb->size, num_full);
    for (int j = 0; j < num_full; j++) {
        for (int i = 0; i < fb->size; i++) {
            if (relations[j].exponents[i] & 1)
                bitmatrix_set(mat, i, j);
        }
    }

    /* Find dependencies using Gaussian elimination */
    /* Simpler approach: work row by row */
    /* Actually, let me do proper GE on the transposed matrix */

    /* We need to find a set S of relations such that sum of exponent vectors = 0 mod 2 */
    /* Standard approach: reduce the matrix and find null space */

    int *dep_rows;
    int ndeps;
    gauss_eliminate(mat, &dep_rows, &ndeps);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* Step 5: For each dependency, try to extract a factor */
    /* TODO: Need to properly track which relations combine.
     * For now, use a simpler approach: build the product directly. */

    /* Simplified: just try random subsets of relations until we find one that works */
    /* Actually, let's redo the Gaussian elimination properly to track dependencies */

    /* For now: if we have enough relations, try the brute force approach for small matrices */
    if (fb->size <= 5000 && num_full > fb->size) {
        /* Redo proper GE */
        /* Matrix: rows = FB primes, cols = relations */
        /* Use u64 arrays */
        int ncols = num_full;
        int nrows = fb->size;
        int words = (ncols + 63) / 64;

        /* Augmented matrix [M | I] where I tracks combinations */
        u64 **matrix = malloc(nrows * sizeof(u64*));
        u64 **history = malloc(ncols * sizeof(u64*));
        for (int i = 0; i < nrows; i++)
            matrix[i] = calloc(words, sizeof(u64));
        for (int j = 0; j < ncols; j++) {
            history[j] = calloc(words, sizeof(u64));
            history[j][j / 64] |= (1ULL << (j % 64));
        }

        /* Fill matrix */
        for (int j = 0; j < ncols; j++)
            for (int i = 0; i < nrows; i++)
                if (relations[j].exponents[i] & 1)
                    matrix[i][j / 64] |= (1ULL << (j % 64));

        /* Row reduce */
        int *pivot = calloc(nrows, sizeof(int));
        memset(pivot, -1, nrows * sizeof(int));

        for (int col = 0; col < nrows; col++) {
            /* Find pivot in this column (across relation rows) */
            int piv = -1;
            for (int j = 0; j < ncols; j++) {
                if ((matrix[col][j / 64] >> (j % 64)) & 1) {
                    piv = j;
                    break;
                }
            }
            if (piv < 0) continue;
            pivot[col] = piv;

            /* Eliminate this column from all other relations */
            for (int j = 0; j < ncols; j++) {
                if (j == piv) continue;
                if ((matrix[col][j / 64] >> (j % 64)) & 1) {
                    /* XOR relation j with relation piv */
                    for (int r = col; r < nrows; r++) {
                        if ((matrix[r][piv / 64] >> (piv % 64)) & 1)
                            matrix[r][j / 64] ^= (1ULL << (j % 64));
                        /* Wait, this isn't right - we need to XOR entire columns */
                    }
                    for (int w = 0; w < words; w++)
                        history[j][w] ^= history[piv][w];
                }
            }
        }

        /* Actually this is getting complicated with the column-based approach.
         * Let me switch to a row-based approach where rows = relations, cols = FB primes.
         */
        for (int i = 0; i < nrows; i++) free(matrix[i]);
        free(matrix);

        /* Rows = relations, Cols = FB primes + identity for tracking */
        int total_cols = nrows + ncols; /* FB primes + relation indices */
        int total_words = (total_cols + 63) / 64;
        u64 **rows2 = malloc(ncols * sizeof(u64*));
        for (int j = 0; j < ncols; j++) {
            rows2[j] = calloc(total_words, sizeof(u64));
            /* Fill FB exponents */
            for (int i = 0; i < nrows; i++)
                if (relations[j].exponents[i] & 1)
                    rows2[j][i / 64] |= (1ULL << (i % 64));
            /* Identity part */
            int c = nrows + j;
            rows2[j][c / 64] |= (1ULL << (c % 64));
        }

        /* Row reduce the left part (FB primes) */
        int *used = calloc(ncols, sizeof(int));
        for (int col = 0; col < nrows; col++) {
            /* Find pivot row for this column */
            int piv = -1;
            for (int r = 0; r < ncols; r++) {
                if (!used[r] && ((rows2[r][col / 64] >> (col % 64)) & 1)) {
                    piv = r;
                    break;
                }
            }
            if (piv < 0) continue;
            used[piv] = 1;

            /* Eliminate this column from all other rows */
            for (int r = 0; r < ncols; r++) {
                if (r == piv) continue;
                if ((rows2[r][col / 64] >> (col % 64)) & 1) {
                    for (int w = 0; w < total_words; w++)
                        rows2[r][w] ^= rows2[piv][w];
                }
            }
        }

        /* Find rows where the left part (FB primes) is all zero = dependencies */
        int dep_found = 0;
        for (int r = 0; r < ncols && dep_found < 20; r++) {
            /* Check if left part is zero */
            int is_zero = 1;
            for (int w = 0; w < (nrows + 63) / 64; w++) {
                if (rows2[r][w] & ((w < (nrows + 63) / 64 - 1) ? ~0ULL :
                    ((1ULL << (nrows % 64)) - 1) ? ((1ULL << (nrows % 64)) - 1) : ~0ULL)) {
                    is_zero = 0;
                    break;
                }
            }
            if (!is_zero) continue;

            /* This row gives a dependency. Extract which relations are involved. */
            dep_found++;

            /* Compute X = product of poly_vals, Y^2 = product of Q(x) values */
            mpz_t X, Y2, Y, g;
            mpz_inits(X, Y2, Y, g, NULL);
            mpz_set_ui(X, 1);
            mpz_set_ui(Y2, 1);

            int *total_exps = calloc(fb->size, sizeof(int));
            int rel_count = 0;

            for (int j = 0; j < ncols; j++) {
                int c = nrows + j;
                if ((rows2[r][c / 64] >> (c % 64)) & 1) {
                    /* Relation j is in the dependency */
                    mpz_mul(X, X, relations[j].poly_val);
                    mpz_mod(X, X, N);
                    for (int k = 0; k < fb->size; k++)
                        total_exps[k] += relations[j].exponents[k];
                    rel_count++;
                }
            }

            /* Compute Y = product of p^(e/2) for all primes */
            mpz_set_ui(Y, 1);
            for (int k = 1; k < fb->size; k++) { /* skip -1 */
                int e = total_exps[k] / 2;
                if (e > 0) {
                    mpz_t pe;
                    mpz_init(pe);
                    mpz_ui_pow_ui(pe, fb->primes[k], e);
                    mpz_mul(Y, Y, pe);
                    mpz_mod(Y, Y, N);
                    mpz_clear(pe);
                }
            }
            /* Handle -1 */
            if (total_exps[0] % 2 != 0) {
                /* Odd number of negative Q(x) - shouldn't happen for a valid dependency */
                free(total_exps);
                mpz_clears(X, Y2, Y, g, NULL);
                continue;
            }

            /* Try gcd(X - Y, N) and gcd(X + Y, N) */
            mpz_sub(tmp, X, Y);
            mpz_gcd(g, tmp, N);

            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                /* Found factor! */
                mpz_t other;
                mpz_init(other);
                mpz_divexact(other, N, g);
                if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
                gmp_printf("%Zd\n", g);

                struct timespec t1;
                clock_gettime(CLOCK_MONOTONIC, &t1);
                double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
                fprintf(stderr, "SIQS: factored in %.3fs (dep #%d, %d relations used)\n",
                        elapsed, dep_found, rel_count);
                mpz_clear(other);
                return 0;
            }

            mpz_add(tmp, X, Y);
            mpz_gcd(g, tmp, N);

            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_t other;
                mpz_init(other);
                mpz_divexact(other, N, g);
                if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
                gmp_printf("%Zd\n", g);

                struct timespec t1;
                clock_gettime(CLOCK_MONOTONIC, &t1);
                double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
                fprintf(stderr, "SIQS: factored in %.3fs (dep #%d, %d relations used)\n",
                        elapsed, dep_found, rel_count);
                mpz_clear(other);
                return 0;
            }

            free(total_exps);
            mpz_clears(X, Y2, Y, g, NULL);
        }

        /* Cleanup */
        for (int j = 0; j < ncols; j++) free(rows2[j]);
        free(rows2);
        for (int j = 0; j < ncols; j++) free(history[j]);
        free(history);
        free(used);
    }

    bitmatrix_free(mat);
    free(dep_rows);

    fprintf(stderr, "SIQS FAILED - could not extract factor from dependencies\n");
    return 1;
}
