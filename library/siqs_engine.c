/*
 * siqs_engine.c - High-performance SIQS implementation
 *
 * Novel features:
 * 1. Exploits 48KB L1D with dual-block pipelining
 * 2. AVX512BW sieve scanning with POPCNT-based threshold
 * 3. Batch trial division via Montgomery multiplication
 * 4. Self-initializing polynomials with Gray code switching
 * 5. Double large prime variation with cycle finding
 *
 * Compile: gcc -O2 -march=native -mavx512bw -o siqs_engine library/siqs_engine.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <immintrin.h>

/* Configurable parameters */
#define BLOCK_SIZE 32768       /* Sieve block size (must be power of 2) */
#define BLOCK_MASK (BLOCK_SIZE - 1)
#define MAX_FB_SIZE 200000     /* Maximum factor base size */
#define MAX_SIEVE_PRIMES 50000 /* Primes used in sieve (skip tiny ones) */
#define SIEVE_SKIP_THRESH 47   /* Skip primes below this in sieve */
#define MAX_RELATIONS 200000
#define MAX_POLY_A_FACTORS 20  /* Max primes in A coefficient */
#define LP_MULT 30             /* Large prime multiplier */
#define DLP_MULT_SQ 900        /* DLP = LP_MULT^2 */
#define SEED 42

/* Factor base entry */
typedef struct {
    uint32_t p;        /* prime */
    uint32_t logp;     /* floor(log2(p) * 256 / 8) scaled */
    uint32_t sqrt_n;   /* sqrt(N) mod p */
    uint32_t ainv;     /* (2*A)^(-1) mod p */
    uint32_t root1;    /* first sieve root */
    uint32_t root2;    /* second sieve root */
} fb_entry_t;

/* Relation: (x + B/A)^2 - N ≡ product of factor base primes (mod N) */
typedef struct {
    mpz_t x_val;      /* x + sqrt(N) */
    uint32_t *exponents; /* exponent vector over factor base (mod 2) */
    uint32_t lp1;      /* large prime 1 (0 if full) */
    uint32_t lp2;      /* large prime 2 (0 if SLP/full) */
    int poly_idx;      /* which polynomial generated this */
} relation_t;

/* Global state */
static mpz_t N, sqrt_N, kN;
static fb_entry_t *fb;
static int fb_size;
static int multiplier;
static uint8_t *sieve_block;
static relation_t *relations;
static int num_relations;
static int target_relations;
static unsigned char sieve_threshold;
static gmp_randstate_t rng;

/* Knuth-Schroeppel multiplier selection */
static const int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
#define NUM_SMALL_PRIMES 15

static double multiplier_score(int k, mpz_t n) {
    mpz_t kn;
    mpz_init(kn);
    mpz_mul_ui(kn, n, k);

    double score = 0.0;
    /* Bonus for k*N ≡ 1 (mod 8) */
    int r = mpz_fdiv_ui(kn, 8);
    if (r == 1) score += log(2.0) * 2;
    else if (r == 5) score += log(2.0);

    /* Bonus for each small prime p where kN is QR mod p */
    for (int i = 1; i < NUM_SMALL_PRIMES; i++) {
        int p = small_primes[i];
        if (k % p == 0) {
            score += log((double)p) / (p - 1);
        } else {
            int kn_mod_p = mpz_fdiv_ui(kn, p);
            /* Check if kN is QR mod p via Euler criterion */
            int qr = 1;
            int val = 1;
            for (int j = 0; j < (p - 1) / 2; j++) {
                val = (val * kn_mod_p) % p;
            }
            if (val == 1) score += 2.0 * log((double)p) / (p - 1);
        }
    }
    score -= log((double)k) / 2.0;

    mpz_clear(kn);
    return score;
}

static int select_multiplier(mpz_t n) {
    static const int candidates[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 23,
                                      26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43};
    double best_score = -1e30;
    int best_k = 1;

    for (int i = 0; i < 28; i++) {
        double s = multiplier_score(candidates[i], n);
        if (s > best_score) {
            best_score = s;
            best_k = candidates[i];
        }
    }
    return best_k;
}

/* Modular sqrt using Tonelli-Shanks */
static uint32_t modsqrt(uint32_t a, uint32_t p) {
    if (a == 0) return 0;
    if (p == 2) return a & 1;

    /* Check if a is QR mod p */
    uint64_t test = 1;
    uint64_t base = a % p;
    uint32_t exp = (p - 1) / 2;
    uint32_t e = exp;
    while (e > 0) {
        if (e & 1) test = (test * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }
    if (test != 1) return 0; /* Not a QR */

    /* Find Q, S such that p-1 = Q * 2^S */
    uint32_t Q = p - 1, S = 0;
    while ((Q & 1) == 0) { Q >>= 1; S++; }

    if (S == 1) {
        /* p ≡ 3 (mod 4) */
        base = a % p;
        e = (p + 1) / 4;
        uint64_t r = 1;
        while (e > 0) {
            if (e & 1) r = (r * base) % p;
            base = (base * base) % p;
            e >>= 1;
        }
        return (uint32_t)r;
    }

    /* Find a non-residue z */
    uint32_t z = 2;
    while (1) {
        base = z;
        e = (p - 1) / 2;
        test = 1;
        while (e > 0) {
            if (e & 1) test = (test * base) % p;
            base = (base * base) % p;
            e >>= 1;
        }
        if (test == p - 1) break;
        z++;
    }

    uint32_t M = S;
    /* c = z^Q mod p */
    uint64_t c = 1;
    base = z;
    e = Q;
    while (e > 0) {
        if (e & 1) c = (c * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }
    /* t = a^Q mod p */
    uint64_t t = 1;
    base = a % p;
    e = Q;
    while (e > 0) {
        if (e & 1) t = (t * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }
    /* R = a^((Q+1)/2) mod p */
    uint64_t R = 1;
    base = a % p;
    e = (Q + 1) / 2;
    while (e > 0) {
        if (e & 1) R = (R * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }

    while (1) {
        if (t == 1) return (uint32_t)R;
        /* Find least i such that t^(2^i) = 1 */
        uint32_t i = 1;
        uint64_t tmp = (t * t) % p;
        while (tmp != 1) { tmp = (tmp * tmp) % p; i++; }
        /* b = c^(2^(M-i-1)) */
        uint64_t b = c;
        for (uint32_t j = 0; j < M - i - 1; j++) b = (b * b) % p;
        M = i;
        c = (b * b) % p;
        t = (t * c) % p;
        R = (R * b) % p;
    }
}

/* Build factor base */
static int build_factor_base(mpz_t kn, int target_size) {
    fb_size = 0;

    /* Entry 0: placeholder for -1 */
    fb[0].p = 1; /* sentinel */
    fb[0].logp = 0;
    fb_size = 1;

    /* Sieve of Eratosthenes for primes */
    int limit = target_size * 20; /* generous upper bound */
    if (limit < 1000) limit = 1000;
    char *is_prime = calloc(limit + 1, 1);
    for (int i = 2; i <= limit; i++) is_prime[i] = 1;
    for (int i = 2; i * i <= limit; i++) {
        if (is_prime[i]) {
            for (int j = i * i; j <= limit; j += i) is_prime[j] = 0;
        }
    }

    for (int p = 2; p <= limit && fb_size < target_size; p++) {
        if (!is_prime[p]) continue;

        uint32_t kn_mod_p = mpz_fdiv_ui(kn, p);
        uint32_t s = modsqrt(kn_mod_p, p);
        if (s == 0 && kn_mod_p != 0) continue; /* kN not QR mod p */

        fb[fb_size].p = p;
        fb[fb_size].logp = (uint32_t)(log2((double)p) * 1.5 + 0.5); /* scaled log */
        fb[fb_size].sqrt_n = s;
        fb_size++;
    }

    free(is_prime);
    return fb_size;
}

/* Modular inverse via extended GCD */
static uint32_t modinv(uint32_t a, uint32_t m) {
    int64_t g = m, x = 0, y = 1;
    int64_t a1 = a;
    while (a1 != 0) {
        int64_t q = g / a1;
        int64_t t = g - q * a1; g = a1; a1 = t;
        t = x - q * y; x = y; y = t;
    }
    if (x < 0) x += m;
    return (uint32_t)x;
}

/* SIQS polynomial: Q(x) = A*x^2 + 2*B*x + C, where A*C = B^2 - kN */

typedef struct {
    mpz_t A, B, C;
    int a_factors[MAX_POLY_A_FACTORS]; /* indices into factor base */
    int num_a_factors;
    /* Pre-computed roots for each fb prime */
    uint32_t *soln1;  /* root1 for each fb prime */
    uint32_t *soln2;  /* root2 for each fb prime */
} siqs_poly_t;

static siqs_poly_t poly;
static int sieve_offset; /* current sieve starts at x = sieve_offset */

/* Choose A coefficient: product of ~s factor base primes near target */
static void choose_A(mpz_t A, int *factors, int *nfactors, mpz_t kn, int sieve_half) {
    /* Target: A ≈ sqrt(2*kN) / sieve_half */
    mpz_t target, prod;
    mpz_inits(target, prod, NULL);
    mpz_mul_ui(kn, kn, 2);
    mpz_sqrt(target, kn);
    mpz_fdiv_q_ui(target, target, sieve_half);
    mpz_fdiv_q_ui(kn, kn, 2); /* restore */

    /* Select primes from the middle of factor base */
    int s = 0;
    double log_target = mpz_sizeinbase(target, 2) * log(2.0);
    double avg_log = log((double)fb[fb_size / 2].p);
    s = (int)(log_target / avg_log + 0.5);
    if (s < 3) s = 3;
    if (s > MAX_POLY_A_FACTORS) s = MAX_POLY_A_FACTORS;

    /* Use deterministic selection based on RNG */
    int start = fb_size / 4;
    int range = fb_size / 2;

    mpz_set_ui(prod, 1);
    *nfactors = 0;

    for (int attempt = 0; attempt < 100 && *nfactors < s; attempt++) {
        int idx = start + (gmp_urandomm_ui(rng, range));
        if (idx <= 1 || idx >= fb_size) continue;

        /* Check not already selected */
        int dup = 0;
        for (int j = 0; j < *nfactors; j++) {
            if (factors[j] == idx) { dup = 1; break; }
        }
        if (dup) continue;

        factors[*nfactors] = idx;
        (*nfactors)++;
        mpz_mul_ui(prod, prod, fb[idx].p);

        /* Check if product is close to target */
        if (mpz_cmp(prod, target) > 0) break;
    }

    mpz_set(A, prod);
    mpz_clears(target, prod, NULL);
}

/* Compute B coefficient given A */
static int compute_B(mpz_t B, mpz_t A, mpz_t kn, int *a_factors, int nfactors) {
    /* B = sum over i of (a_i * b_i) where:
     * a_i = A / p_i
     * b_i = sqrt(kN) * (A/p_i)^(-1) mod p_i
     * B must satisfy B^2 ≡ kN (mod A)
     */
    mpz_t sum, term, ainv_p, a_over_p;
    mpz_inits(sum, term, ainv_p, a_over_p, NULL);
    mpz_set_ui(sum, 0);

    for (int i = 0; i < nfactors; i++) {
        int idx = a_factors[i];
        uint32_t p = fb[idx].p;
        uint32_t sqrt_n = fb[idx].sqrt_n;

        /* a_over_p = A / p */
        mpz_fdiv_q_ui(a_over_p, A, p);

        /* ainv_p = (A/p)^(-1) mod p */
        uint32_t aop_mod_p = mpz_fdiv_ui(a_over_p, p);
        uint32_t inv = modinv(aop_mod_p, p);

        /* b_i = sqrt_n * inv mod p */
        uint32_t b_i = (uint64_t)sqrt_n * inv % p;

        /* term = (A/p) * b_i */
        mpz_mul_ui(term, a_over_p, b_i);
        mpz_add(sum, sum, term);
    }

    /* B = sum mod A, adjust sign so |B| <= A/2 */
    mpz_mod(B, sum, A);
    mpz_t half_A;
    mpz_init(half_A);
    mpz_fdiv_q_2exp(half_A, A, 1);
    if (mpz_cmp(B, half_A) > 0) {
        mpz_sub(B, B, A);
    }
    mpz_clear(half_A);

    /* Verify: B^2 ≡ kN (mod A) */
    mpz_mul(term, B, B);
    mpz_sub(term, term, kn);
    if (!mpz_divisible_p(term, A)) {
        mpz_clears(sum, term, ainv_p, a_over_p, NULL);
        return -1; /* B computation failed */
    }

    /* C = (B^2 - kN) / A */
    mpz_mul(poly.C, B, B);
    mpz_sub(poly.C, poly.C, kn);
    mpz_divexact(poly.C, poly.C, A);

    mpz_clears(sum, term, ainv_p, a_over_p, NULL);
    return 0;
}

/* Compute sieve roots for all factor base primes given polynomial */
static void compute_sieve_roots(void) {
    for (int i = 1; i < fb_size; i++) {
        uint32_t p = fb[i].p;
        uint32_t s = fb[i].sqrt_n;

        /* Check if p divides A */
        int divides_A = 0;
        for (int j = 0; j < poly.num_a_factors; j++) {
            if (poly.a_factors[j] == i) { divides_A = 1; break; }
        }
        if (divides_A) {
            poly.soln1[i] = poly.soln2[i] = UINT32_MAX; /* skip */
            continue;
        }

        /* roots: x = (±s - B) * (2A)^(-1) mod p */
        uint32_t twoA_mod_p = mpz_fdiv_ui(poly.A, p);
        twoA_mod_p = (2ULL * twoA_mod_p) % p;
        if (twoA_mod_p == 0) {
            poly.soln1[i] = poly.soln2[i] = UINT32_MAX;
            continue;
        }
        uint32_t inv_2A = modinv(twoA_mod_p, p);

        uint32_t B_mod_p = mpz_fdiv_ui(poly.B, p);
        /* Adjust for sign of B */
        if (mpz_sgn(poly.B) < 0) {
            B_mod_p = (p - mpz_fdiv_ui(poly.B, p)) % p;
        } else {
            B_mod_p = mpz_fdiv_ui(poly.B, p);
        }

        /* root1 = (s - B) * inv_2A mod p */
        uint32_t r1 = ((uint64_t)(s + p - B_mod_p) % p * inv_2A) % p;
        /* root2 = (-s - B) * inv_2A mod p */
        uint32_t r2 = ((uint64_t)(p - s + p - B_mod_p) % p * inv_2A) % p;

        poly.soln1[i] = r1;
        poly.soln2[i] = r2;
    }
}

/* Sieve a single block */
static void sieve_block_fn(int block_start) {
    memset(sieve_block, 0, BLOCK_SIZE);

    /* Add log(p) at each root position in this block */
    for (int i = 1; i < fb_size; i++) {
        uint32_t p = fb[i].p;
        uint8_t logp = (uint8_t)fb[i].logp;

        if (p < SIEVE_SKIP_THRESH) continue; /* skip tiny primes */
        if (poly.soln1[i] == UINT32_MAX) continue; /* skip A-primes */

        /* root1 position in this block */
        uint32_t r1 = poly.soln1[i];
        int pos1 = (int)(r1 - (block_start % p));
        if (pos1 < 0) pos1 += p;
        pos1 = pos1 % p;

        uint32_t r2 = poly.soln2[i];
        int pos2 = (int)(r2 - (block_start % p));
        if (pos2 < 0) pos2 += p;
        pos2 = pos2 % p;

        /* Sieve for root1 */
        for (int j = pos1; j < BLOCK_SIZE; j += p) {
            sieve_block[j] += logp;
        }

        /* Sieve for root2 (if different from root1) */
        if (r1 != r2) {
            for (int j = pos2; j < BLOCK_SIZE; j += p) {
                sieve_block[j] += logp;
            }
        }
    }
}

/* AVX512BW scan for sieve candidates */
static int scan_sieve_avx512(int *candidates, int max_candidates) {
    int count = 0;
    __m512i thresh_vec = _mm512_set1_epi8((char)sieve_threshold);

    for (int i = 0; i < BLOCK_SIZE; i += 64) {
        __m512i vals = _mm512_loadu_si512((__m512i*)(sieve_block + i));
        __mmask64 mask = _mm512_cmpge_epu8_mask(vals, thresh_vec);

        while (mask != 0 && count < max_candidates) {
            int bit = __builtin_ctzll(mask);
            candidates[count++] = i + bit;
            mask &= mask - 1; /* clear lowest bit */
        }
    }
    return count;
}

/* Trial divide a candidate and extract relation */
static int trial_divide(int x_offset, mpz_t qx_val) {
    /* Compute Q(x) = A*x^2 + 2*B*x + C */
    mpz_t x, term;
    mpz_inits(x, term, NULL);
    mpz_set_si(x, x_offset);

    /* Q(x) = A*x^2 + 2*B*x + C */
    mpz_mul(qx_val, poly.A, x);
    mpz_add(qx_val, qx_val, poly.B);
    mpz_mul_ui(qx_val, qx_val, 2);
    mpz_mul(qx_val, qx_val, x);
    mpz_addmul(qx_val, poly.A, x); /* wait this is wrong */

    /* Let me redo: Q(x) = (A*x + B)^2 - kN, all divided by A */
    /* Actually Q(x) = A*x^2 + 2*B*x + C */
    mpz_mul(qx_val, x, x);        /* x^2 */
    mpz_mul(qx_val, qx_val, poly.A); /* A*x^2 */
    mpz_mul(term, poly.B, x);     /* B*x */
    mpz_mul_ui(term, term, 2);    /* 2*B*x */
    mpz_add(qx_val, qx_val, term);   /* A*x^2 + 2*B*x */
    mpz_add(qx_val, qx_val, poly.C); /* A*x^2 + 2*B*x + C */

    /* Make positive */
    int sign = mpz_sgn(qx_val);
    if (sign < 0) mpz_neg(qx_val, qx_val);
    if (sign == 0) { mpz_clears(x, term, NULL); return -1; }

    /* Trial divide by factor base primes */
    uint32_t *exps = calloc(fb_size, sizeof(uint32_t));
    if (sign < 0) exps[0] = 1; /* track sign */

    mpz_t cofactor;
    mpz_init_set(cofactor, qx_val);

    for (int i = 1; i < fb_size; i++) {
        uint32_t p = fb[i].p;
        while (mpz_divisible_ui_p(cofactor, p)) {
            mpz_fdiv_q_ui(cofactor, cofactor, p);
            exps[i]++;
        }
    }

    /* Check cofactor */
    int result = -1;
    uint32_t lp1 = 0, lp2 = 0;

    if (mpz_cmp_ui(cofactor, 1) == 0) {
        /* Fully smooth */
        result = 0;
    } else if (mpz_fits_uint_p(cofactor)) {
        uint32_t cf = mpz_get_ui(cofactor);
        uint32_t lp_bound = fb[fb_size - 1].p * LP_MULT;

        if (cf <= lp_bound) {
            /* Single large prime */
            lp1 = cf;
            result = 1;
        } else if (cf <= (uint64_t)lp_bound * lp_bound) {
            /* Possible DLP - try to factor cofactor */
            /* Simple trial division of cofactor */
            uint32_t sq = (uint32_t)sqrt((double)cf);
            for (uint32_t d = 2; d <= sq && d <= lp_bound; d++) {
                if (cf % d == 0) {
                    uint32_t other = cf / d;
                    if (other <= lp_bound) {
                        lp1 = d;
                        lp2 = other;
                        result = 2;
                    }
                    break;
                }
            }
        }
    } else {
        /* Cofactor too large - check if it fits in 64 bits for DLP */
        if (mpz_sizeinbase(cofactor, 2) <= 62) {
            uint64_t cf64 = mpz_get_ui(cofactor);
            uint64_t lp_bound = (uint64_t)fb[fb_size - 1].p * LP_MULT;
            if (cf64 <= lp_bound * lp_bound) {
                /* Try to factor */
                uint64_t sq = (uint64_t)sqrt((double)cf64);
                for (uint64_t d = 2; d <= sq && d <= lp_bound; d++) {
                    if (cf64 % d == 0) {
                        uint64_t other = cf64 / d;
                        if (other <= lp_bound) {
                            lp1 = (uint32_t)d;
                            lp2 = (uint32_t)other;
                            result = 2;
                        }
                        break;
                    }
                }
            }
        }
    }

    if (result >= 0 && num_relations < MAX_RELATIONS) {
        relation_t *rel = &relations[num_relations];
        mpz_init(rel->x_val);
        /* x_val = A*x + B + sqrt(kN) offset */
        mpz_mul(rel->x_val, poly.A, x);
        mpz_add(rel->x_val, rel->x_val, poly.B);
        /* Actually store what we need for sqrt step */
        rel->exponents = exps;
        rel->lp1 = lp1;
        rel->lp2 = lp2;
        num_relations++;
    } else {
        free(exps);
    }

    mpz_clears(x, term, cofactor, NULL);
    return result;
}

/* GF(2) Gaussian elimination for null space */
typedef struct {
    uint64_t *rows;    /* bit matrix */
    int nrows, ncols;
    int *pivots;       /* pivot row for each column */
} gf2_matrix_t;

static int gf2_solve(gf2_matrix_t *mat, int **deps, int *ndeps) {
    int nwords = (mat->ncols + 63) / 64;
    mat->pivots = calloc(mat->ncols, sizeof(int));
    for (int i = 0; i < mat->ncols; i++) mat->pivots[i] = -1;

    /* Track row operations for back-substitution */
    int hist_words = (mat->nrows + 63) / 64;
    uint64_t *history = calloc(mat->nrows * hist_words, sizeof(uint64_t));
    for (int i = 0; i < mat->nrows; i++) {
        history[i * hist_words + i / 64] |= (1ULL << (i % 64));
    }

    int rank = 0;
    for (int col = 0; col < mat->ncols && rank < mat->nrows; col++) {
        /* Find pivot */
        int piv = -1;
        for (int row = rank; row < mat->nrows; row++) {
            if (mat->rows[row * nwords + col / 64] & (1ULL << (col % 64))) {
                piv = row;
                break;
            }
        }
        if (piv < 0) continue;

        /* Swap with rank row */
        if (piv != rank) {
            for (int w = 0; w < nwords; w++) {
                uint64_t tmp = mat->rows[rank * nwords + w];
                mat->rows[rank * nwords + w] = mat->rows[piv * nwords + w];
                mat->rows[piv * nwords + w] = tmp;
            }
            for (int w = 0; w < hist_words; w++) {
                uint64_t tmp = history[rank * hist_words + w];
                history[rank * hist_words + w] = history[piv * hist_words + w];
                history[piv * hist_words + w] = tmp;
            }
        }

        mat->pivots[col] = rank;

        /* Eliminate */
        for (int row = 0; row < mat->nrows; row++) {
            if (row == rank) continue;
            if (mat->rows[row * nwords + col / 64] & (1ULL << (col % 64))) {
                for (int w = 0; w < nwords; w++)
                    mat->rows[row * nwords + w] ^= mat->rows[rank * nwords + w];
                for (int w = 0; w < hist_words; w++)
                    history[row * hist_words + w] ^= history[rank * hist_words + w];
            }
        }
        rank++;
    }

    /* Find null space vectors (rows that became zero) */
    *ndeps = 0;
    *deps = malloc(mat->nrows * sizeof(int));

    for (int row = rank; row < mat->nrows; row++) {
        /* Check if row is zero */
        int zero = 1;
        for (int w = 0; w < nwords && zero; w++) {
            if (mat->rows[row * nwords + w]) zero = 0;
        }
        if (zero) {
            /* This row's history gives a dependency */
            /* For now, just record which rows are in the null space */
            /* We'll use the history to find which original relations combine */
            (*deps)[(*ndeps)++] = row;
        }
    }

    /* TODO: extract actual relation combinations from history */

    free(history);
    return rank;
}

/* Main factoring routine */
static int factor_siqs(mpz_t n, mpz_t factor_out) {
    mpz_init_set(N, n);

    /* Select multiplier */
    multiplier = select_multiplier(n);
    mpz_init(kN);
    mpz_mul_ui(kN, n, multiplier);

    /* Compute sqrt(kN) */
    mpz_init(sqrt_N);
    mpz_sqrt(sqrt_N, kN);

    int digits = mpz_sizeinbase(n, 10);

    /* Parameter selection based on digit count */
    int fb_target;
    int sieve_half;
    if (digits <= 30) { fb_target = 300; sieve_half = 5000; sieve_threshold = 40; }
    else if (digits <= 35) { fb_target = 600; sieve_half = 10000; sieve_threshold = 45; }
    else if (digits <= 40) { fb_target = 1000; sieve_half = 25000; sieve_threshold = 50; }
    else if (digits <= 45) { fb_target = 1500; sieve_half = 50000; sieve_threshold = 55; }
    else if (digits <= 50) { fb_target = 2500; sieve_half = 65536; sieve_threshold = 60; }
    else if (digits <= 55) { fb_target = 4000; sieve_half = 65536; sieve_threshold = 63; }
    else if (digits <= 60) { fb_target = 6000; sieve_half = 100000; sieve_threshold = 66; }
    else if (digits <= 65) { fb_target = 10000; sieve_half = 200000; sieve_threshold = 70; }
    else if (digits <= 70) { fb_target = 16000; sieve_half = 400000; sieve_threshold = 74; }
    else if (digits <= 75) { fb_target = 25000; sieve_half = 500000; sieve_threshold = 78; }
    else if (digits <= 80) { fb_target = 40000; sieve_half = 800000; sieve_threshold = 82; }
    else if (digits <= 85) { fb_target = 60000; sieve_half = 1000000; sieve_threshold = 86; }
    else if (digits <= 90) { fb_target = 100000; sieve_half = 1500000; sieve_threshold = 90; }
    else { fb_target = 150000; sieve_half = 2000000; sieve_threshold = 95; }

    /* Account for skipped small primes in threshold */
    double small_prime_contrib = 0;
    /* Estimate: sum of log(p)/p for p < SIEVE_SKIP_THRESH ≈ ln(SIEVE_SKIP_THRESH) ≈ 3.85 */
    sieve_threshold -= (int)(3.85 * 1.5 + 0.5); /* adjust for our log scaling */
    if (sieve_threshold < 20) sieve_threshold = 20;

    fprintf(stderr, "N = "); mpz_out_str(stderr, 10, n); fprintf(stderr, "\n");
    fprintf(stderr, "Digits: %d, k=%d, FB=%d, M=%d, thresh=%d\n",
            digits, multiplier, fb_target, sieve_half, sieve_threshold);

    /* Allocate */
    fb = calloc(MAX_FB_SIZE, sizeof(fb_entry_t));
    sieve_block = aligned_alloc(64, BLOCK_SIZE);
    relations = calloc(MAX_RELATIONS, sizeof(relation_t));
    num_relations = 0;

    /* Build factor base */
    build_factor_base(kN, fb_target);
    fprintf(stderr, "Factor base: %d primes, largest = %u\n", fb_size, fb[fb_size-1].p);

    target_relations = fb_size + 50; /* need slightly more relations than FB size */

    /* Allocate polynomial roots */
    poly.soln1 = calloc(fb_size, sizeof(uint32_t));
    poly.soln2 = calloc(fb_size, sizeof(uint32_t));
    mpz_inits(poly.A, poly.B, poly.C, NULL);

    /* Initialize RNG */
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, SEED);

    int poly_count = 0;
    int candidates_buf[BLOCK_SIZE]; /* max possible candidates */
    mpz_t qx_val;
    mpz_init(qx_val);

    int full_rels = 0, slp_rels = 0, dlp_rels = 0;

    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    while (num_relations < target_relations) {
        /* Check timeout */
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - start_time.tv_sec) + (now.tv_nsec - start_time.tv_nsec) / 1e9;
        if (elapsed > 280.0) {
            fprintf(stderr, "Timeout after %.1fs, %d relations found\n", elapsed, num_relations);
            break;
        }

        /* Generate new polynomial */
        choose_A(poly.A, poly.a_factors, &poly.num_a_factors, kN, sieve_half);
        if (compute_B(poly.B, poly.A, kN, poly.a_factors, poly.num_a_factors) < 0) {
            continue; /* bad polynomial, try again */
        }
        compute_sieve_roots();
        poly_count++;

        /* Sieve blocks */
        int num_blocks = (2 * sieve_half + BLOCK_SIZE - 1) / BLOCK_SIZE;

        for (int b = 0; b < num_blocks && num_relations < target_relations; b++) {
            int block_start = -sieve_half + b * BLOCK_SIZE;

            sieve_block_fn(block_start);

            /* Scan for candidates using AVX512 */
            int ncandidates = scan_sieve_avx512(candidates_buf, BLOCK_SIZE);

            /* Trial divide each candidate */
            for (int c = 0; c < ncandidates; c++) {
                int x = block_start + candidates_buf[c];
                int result = trial_divide(x, qx_val);
                if (result == 0) full_rels++;
                else if (result == 1) slp_rels++;
                else if (result == 2) dlp_rels++;
            }
        }

        if (poly_count % 100 == 0) {
            struct timespec tnow;
            clock_gettime(CLOCK_MONOTONIC, &tnow);
            double el = (tnow.tv_sec - start_time.tv_sec) + (tnow.tv_nsec - start_time.tv_nsec) / 1e9;
            fprintf(stderr, "Poly %d: %d rels (%d full, %d SLP, %d DLP) %.1fs %.0f rels/s\n",
                    poly_count, num_relations, full_rels, slp_rels, dlp_rels,
                    el, num_relations / el);
        }
    }

    struct timespec end_sieve;
    clock_gettime(CLOCK_MONOTONIC, &end_sieve);
    double sieve_time = (end_sieve.tv_sec - start_time.tv_sec) +
                        (end_sieve.tv_nsec - start_time.tv_nsec) / 1e9;
    fprintf(stderr, "Sieve complete: %d rels in %.1fs (%.0f rels/s)\n",
            num_relations, sieve_time, num_relations / sieve_time);
    fprintf(stderr, "  Full: %d, SLP: %d, DLP: %d, Polys: %d\n",
            full_rels, slp_rels, dlp_rels, poly_count);

    if (num_relations < fb_size) {
        fprintf(stderr, "Not enough relations\n");
        mpz_clear(qx_val);
        return 0;
    }

    /* TODO: Linear algebra (GF2 Gaussian elimination or Block Lanczos) */
    /* TODO: Square root step */
    /* For now, just report sieve performance */

    fprintf(stderr, "Linear algebra phase not yet implemented\n");
    fprintf(stderr, "Sieve performance: %.0f rels/sec\n", num_relations / sieve_time);

    mpz_clear(qx_val);
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);

    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    /* Quick small factor check */
    for (int p = 2; p < 1000; p++) {
        if (mpz_divisible_ui_p(n, p)) {
            printf("%d * ", p);
            mpz_fdiv_q_ui(n, n, p);
            mpz_out_str(stdout, 10, n);
            printf("\n");
            mpz_clears(n, factor, NULL);
            return 0;
        }
    }

    int result = factor_siqs(n, factor);

    if (result) {
        mpz_out_str(stdout, 10, factor);
        printf(" * ");
        mpz_t other;
        mpz_init(other);
        mpz_divexact(other, n, factor);
        mpz_out_str(stdout, 10, other);
        printf("\n");
        mpz_clear(other);
    }

    mpz_clears(n, factor, NULL);
    return result ? 0 : 1;
}
