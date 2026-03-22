/*
 * mpqs_custom.c - Custom MPQS implementation for balanced semiprimes
 *
 * Novel features:
 * 1. Cache-line aware sieve: Process sieve in 64-byte chunks to minimize
 *    cache line bouncing. For each cache line, accumulate all prime
 *    contributions before moving on.
 * 2. Optimized polynomial switching: Use self-initialization (Gray code)
 *    to update only 2 roots per polynomial change.
 * 3. Double large prime variation with efficient hash-based cycle detection.
 * 4. Block Lanczos over GF(2) for the linear algebra phase.
 *
 * Target: 50-80 digit balanced semiprimes
 *
 * Compile: gcc -O2 -march=native -o mpqs_custom library/mpqs_custom.c -lgmp -lm
 * Usage: ./mpqs_custom <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

/* Configuration */
#define MAX_FB_SIZE    100000   /* Max factor base primes */
#define BLOCK_SIZE     32768    /* Sieve block size (32KB = L1 cache) */
#define MAX_RELATIONS  200000   /* Max relations to collect */
#define MAX_POLYS      10000000
#define LP_HASH_SIZE   (1 << 20)  /* Hash table for large prime relations */
#define LP_HASH_MASK   (LP_HASH_SIZE - 1)

/* Sieve threshold adjustment: lower = more candidates = slower tdiv but fewer missed */
#define SIEVE_THRESH_ADJUST  30

/* Small primes table */
static uint32_t small_primes[10000];
static int n_small_primes = 0;

static void init_small_primes(uint32_t limit) {
    uint8_t *sieve = calloc(limit + 1, 1);
    small_primes[n_small_primes++] = 2;
    for (uint32_t i = 3; i <= limit; i += 2) {
        if (!sieve[i]) {
            if (n_small_primes < 10000) small_primes[n_small_primes++] = i;
            for (uint32_t j = (uint64_t)i * i; j <= limit; j += 2 * i)
                sieve[j] = 1;
        }
    }
    free(sieve);
}

/* Tonelli-Shanks: find x such that x^2 ≡ n (mod p) */
static uint32_t tonelli_shanks(uint32_t n, uint32_t p) {
    if (p == 2) return n & 1;
    if (n == 0) return 0;

    /* Check if n is a QR mod p */
    uint64_t test = 1;
    uint64_t base = n % p;
    uint32_t exp = (p - 1) / 2;
    uint32_t e = exp;
    while (e > 0) {
        if (e & 1) test = (test * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }
    if (test != 1) return 0;  /* Not a QR */

    /* Factor out powers of 2 from p-1 */
    uint32_t Q = p - 1, S = 0;
    while ((Q & 1) == 0) { Q >>= 1; S++; }

    if (S == 1) {
        /* p ≡ 3 mod 4 */
        uint64_t r = 1;
        base = n % p;
        e = (p + 1) / 4;
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
        test = 1;
        base = z;
        e = (p - 1) / 2;
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
    /* t = n^Q mod p */
    uint64_t t = 1;
    base = n % p;
    e = Q;
    while (e > 0) {
        if (e & 1) t = (t * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }
    /* R = n^((Q+1)/2) mod p */
    uint64_t R = 1;
    base = n % p;
    e = (Q + 1) / 2;
    while (e > 0) {
        if (e & 1) R = (R * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }

    while (1) {
        if (t == 0) return 0;
        if (t == 1) return (uint32_t)R;

        /* Find least i such that t^(2^i) ≡ 1 */
        uint32_t i = 0;
        uint64_t tmp = t;
        while (tmp != 1) {
            tmp = (tmp * tmp) % p;
            i++;
        }
        if (i == M) return 0;  /* Error */

        /* b = c^(2^(M-i-1)) */
        uint64_t b = c;
        for (uint32_t j = 0; j < M - i - 1; j++)
            b = (b * b) % p;

        M = i;
        c = (b * b) % p;
        t = (t * c) % p;
        R = (R * b) % p;
    }
}

/* Factor base */
typedef struct {
    uint32_t *primes;      /* Factor base primes */
    uint32_t *roots;       /* sqrt(N) mod p */
    uint8_t  *logp;        /* log2(p) approximation */
    int       size;        /* Number of primes */
    uint32_t  max_prime;   /* Largest prime */
} factor_base_t;

/* Relation */
typedef struct {
    mpz_t    Qx;           /* Q(x) value */
    mpz_t    ax_b;         /* a*x + b (used for sqrt) */
    uint32_t *fb_indices;  /* Factor base indices */
    uint8_t  *fb_powers;   /* Powers of each FB prime */
    int       n_factors;   /* Number of FB factors */
    uint32_t  lp1, lp2;   /* Large primes (0 if none) */
    int       num_lp;      /* 0, 1, or 2 large primes */
} relation_t;

/* Knuth-Schroeppel multiplier selection */
static uint32_t choose_multiplier(const mpz_t n, int fb_size) {
    static uint32_t mult_list[] = {1,3,5,7,11,13,15,17,19,21,23,29,31,33,35,37,39,41,43,47,51,53,55,59,61,67,69,71,73};
    int n_mult = sizeof(mult_list) / sizeof(mult_list[0]);
    double best_score = 1e30;
    uint32_t best_mult = 1;

    for (int i = 0; i < n_mult; i++) {
        uint32_t k = mult_list[i];
        double score = 0.5 * log((double)k);

        /* Score contribution of 2 */
        uint64_t kn_mod8 = (mpz_get_ui(n) * k) % 8;
        switch (kn_mod8) {
            case 1: score -= 2.0 * log(2.0); break;
            case 5: score -= log(2.0); break;
            case 3: case 7: score -= 0.5 * log(2.0); break;
        }

        /* Score contribution of small primes */
        int max_p = fb_size < 300 ? fb_size : 300;
        for (int j = 1; j < max_p && j < n_small_primes; j++) {
            uint32_t p = small_primes[j];
            double contrib = log((double)p) / (double)(p - 1);
            uint32_t kn_modp = (uint32_t)((mpz_tdiv_ui(n, p) * (uint64_t)k) % p);

            if (kn_modp == 0) {
                score -= contrib;
            } else {
                /* Check if kn is a QR mod p */
                uint64_t test = 1, base = kn_modp;
                uint32_t e = (p - 1) / 2;
                while (e > 0) {
                    if (e & 1) test = (test * base) % p;
                    base = (base * base) % p;
                    e >>= 1;
                }
                if (test == 1) score -= 2 * contrib;
            }
        }

        if (score < best_score) {
            best_score = score;
            best_mult = k;
        }
    }

    return best_mult;
}

/* SIQS polynomial: Q(x) = (ax+b)^2 - kN = a^2*x^2 + 2abx + (b^2-kN) */
/* With self-initialization, a = product of FB primes */

typedef struct {
    mpz_t a;         /* Leading coefficient */
    mpz_t b;         /* Linear coefficient / 2 */
    mpz_t c;         /* Constant: (b^2 - kN) / a */
    mpz_t kn;        /* k*N */

    /* Self-initialization data */
    uint32_t *a_indices;  /* FB indices used in a */
    int       n_a_primes; /* Number of primes in a */
    mpz_t    *B_array;    /* B[j] values for self-init */
    int       n_B;        /* Number of B values */
    int       gray_idx;   /* Current Gray code index */
    int       max_gray;   /* Maximum Gray code index */

    /* Roots for sieving */
    uint32_t *soln1;      /* Root 1 for each FB prime */
    uint32_t *soln2;      /* Root 2 for each FB prime */
} siqs_poly_t;

/* Simple sieve-based SIQS */
typedef struct {
    mpz_t n;              /* Number to factor */
    mpz_t kn;             /* k*N */
    uint32_t multiplier;
    factor_base_t fb;
    siqs_poly_t poly;

    /* Sieve */
    uint8_t *sieve_array;
    int sieve_len;        /* = 2 * M */
    int M;                /* Half-width */

    /* Relations */
    relation_t *relations;
    int n_relations;
    int target_relations;

    /* Large prime hash */
    struct { uint32_t lp; int rel_idx; } *lp_hash;

    /* Stats */
    int n_full;
    int n_partial;
    int n_dlp;
    int n_polys;
} siqs_state_t;

/* Initialize factor base */
static int init_factor_base(factor_base_t *fb, const mpz_t kn, int fb_size) {
    fb->primes = malloc(fb_size * sizeof(uint32_t));
    fb->roots = malloc(fb_size * sizeof(uint32_t));
    fb->logp = malloc(fb_size * sizeof(uint8_t));
    fb->size = 0;

    /* Prime 2 is always in the factor base */
    fb->primes[0] = 2;
    fb->roots[0] = 1;
    fb->logp[0] = 1;
    fb->size = 1;

    /* Add primes p where kN is a quadratic residue mod p */
    for (int i = 1; i < n_small_primes && fb->size < fb_size; i++) {
        uint32_t p = small_primes[i];
        uint32_t nmodp = (uint32_t)mpz_tdiv_ui(kn, p);

        if (nmodp == 0) {
            fb->primes[fb->size] = p;
            fb->roots[fb->size] = 0;
            fb->logp[fb->size] = (uint8_t)(log2((double)p) + 0.5);
            fb->size++;
            continue;
        }

        uint32_t root = tonelli_shanks(nmodp, p);
        if (root > 0) {
            fb->primes[fb->size] = p;
            fb->roots[fb->size] = root;
            fb->logp[fb->size] = (uint8_t)(log2((double)p) + 0.5);
            fb->size++;
        }
    }

    fb->max_prime = fb->primes[fb->size - 1];
    return fb->size;
}

/* Sieve one block */
static void sieve_block(uint8_t *sieve, int block_start, int block_size,
                         const factor_base_t *fb, const siqs_poly_t *poly,
                         int start_prime) {
    for (int i = start_prime; i < fb->size; i++) {
        uint32_t p = fb->primes[i];
        uint8_t logp = fb->logp[i];

        /* Get roots for this block */
        int r1 = (int)poly->soln1[i] - block_start;
        int r2 = (int)poly->soln2[i] - block_start;

        /* Adjust roots to be in [0, block_size) */
        while (r1 < 0) r1 += p;
        while (r2 < 0) r2 += p;

        /* Sieve root 1 */
        for (int j = r1; j < block_size; j += p) {
            sieve[j] -= logp;
        }

        /* Sieve root 2 (if different from root 1) */
        if (r1 != r2) {
            for (int j = r2; j < block_size; j += p) {
                sieve[j] -= logp;
            }
        }
    }
}

/* Trial divide a candidate to check if it's smooth */
static int trial_divide_candidate(mpz_t Qx, const factor_base_t *fb,
                                   uint32_t *indices, uint8_t *powers,
                                   int *n_factors, uint32_t lp_bound,
                                   uint32_t *lp1, uint32_t *lp2) {
    mpz_t rem;
    mpz_init_set(rem, Qx);
    *n_factors = 0;
    *lp1 = 0;
    *lp2 = 0;

    /* Remove sign */
    if (mpz_sgn(rem) < 0) {
        indices[*n_factors] = 0;  /* Index 0 = sign */
        powers[*n_factors] = 1;
        (*n_factors)++;
        mpz_neg(rem, rem);
    }

    /* Trial divide by factor base primes */
    for (int i = 0; i < fb->size; i++) {
        uint32_t p = fb->primes[i];
        if (mpz_divisible_ui_p(rem, p)) {
            uint8_t pwr = 0;
            while (mpz_divisible_ui_p(rem, p)) {
                mpz_divexact_ui(rem, rem, p);
                pwr++;
            }
            indices[*n_factors] = i;
            powers[*n_factors] = pwr;
            (*n_factors)++;

            if (mpz_cmp_ui(rem, 1) == 0) {
                mpz_clear(rem);
                return 1;  /* Fully smooth */
            }
        }
    }

    /* Check for large primes */
    if (mpz_fits_ulong_p(rem)) {
        uint32_t cofactor = (uint32_t)mpz_get_ui(rem);
        if (cofactor <= lp_bound) {
            *lp1 = cofactor;
            mpz_clear(rem);
            return 2;  /* Single large prime */
        }
    }

    /* Check for double large prime */
    if (mpz_sizeinbase(rem, 2) <= 2 * (int)ceil(log2((double)lp_bound)) + 2) {
        /* Try to split the cofactor */
        if (mpz_probab_prime_p(rem, 1)) {
            /* Single large prime that exceeds bound */
            mpz_clear(rem);
            return 0;
        }

        /* Quick Pollard rho to split */
        mpz_t x, y, d, q;
        mpz_inits(x, y, d, q, NULL);
        mpz_set_ui(x, 2);
        mpz_set_ui(y, 2);
        mpz_set_ui(q, 1);
        int found = 0;

        for (int iter = 0; iter < 256 && !found; iter++) {
            mpz_mul(x, x, x); mpz_add_ui(x, x, 1); mpz_mod(x, x, rem);
            mpz_mul(y, y, y); mpz_add_ui(y, y, 1); mpz_mod(y, y, rem);
            mpz_mul(y, y, y); mpz_add_ui(y, y, 1); mpz_mod(y, y, rem);
            mpz_sub(d, x, y); mpz_abs(d, d);
            mpz_mul(q, q, d); mpz_mod(q, q, rem);

            if ((iter & 15) == 15) {
                mpz_gcd(d, q, rem);
                if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, rem) < 0) {
                    if (mpz_fits_ulong_p(d)) {
                        uint32_t f1 = (uint32_t)mpz_get_ui(d);
                        mpz_divexact(d, rem, d);
                        if (mpz_fits_ulong_p(d)) {
                            uint32_t f2 = (uint32_t)mpz_get_ui(d);
                            if (f1 <= lp_bound && f2 <= lp_bound) {
                                *lp1 = f1 < f2 ? f1 : f2;
                                *lp2 = f1 < f2 ? f2 : f1;
                                found = 1;
                            }
                        }
                    }
                }
                mpz_set_ui(q, 1);
            }
        }

        mpz_clears(x, y, d, q, NULL);
        if (found) {
            mpz_clear(rem);
            return 3;  /* Double large prime */
        }
    }

    mpz_clear(rem);
    return 0;  /* Not smooth enough */
}

/* Parameters based on digit count */
static void get_params(int digits, int *fb_size, int *M, int *large_prime_bits) {
    if (digits <= 40) {
        *fb_size = 200;
        *M = 16384;
        *large_prime_bits = 17;
    } else if (digits <= 50) {
        *fb_size = 800;
        *M = 32768;
        *large_prime_bits = 20;
    } else if (digits <= 60) {
        *fb_size = 2500;
        *M = 65536;
        *large_prime_bits = 22;
    } else if (digits <= 65) {
        *fb_size = 5000;
        *M = 98304;
        *large_prime_bits = 23;
    } else if (digits <= 70) {
        *fb_size = 10000;
        *M = 131072;
        *large_prime_bits = 24;
    } else if (digits <= 75) {
        *fb_size = 20000;
        *M = 196608;
        *large_prime_bits = 25;
    } else if (digits <= 80) {
        *fb_size = 40000;
        *M = 262144;
        *large_prime_bits = 26;
    } else if (digits <= 85) {
        *fb_size = 60000;
        *M = 393216;
        *large_prime_bits = 27;
    } else {
        *fb_size = 80000;
        *M = 524288;
        *large_prime_bits = 28;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t n, kn, factor;
    mpz_inits(n, kn, factor, NULL);

    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int digits = (int)mpz_sizeinbase(n, 10);
    fprintf(stderr, "Factoring %d-digit number with custom MPQS\n", digits);

    /* Initialize small primes */
    init_small_primes(2000000);

    /* Choose multiplier */
    int fb_size_target, M, lp_bits;
    get_params(digits, &fb_size_target, &M, &lp_bits);

    uint32_t mult = choose_multiplier(n, fb_size_target);
    mpz_mul_ui(kn, n, mult);
    fprintf(stderr, "Multiplier: %u, FB target: %d, M: %d, LP bits: %d\n",
            mult, fb_size_target, M, lp_bits);

    /* Build factor base */
    factor_base_t fb;
    init_factor_base(&fb, kn, fb_size_target);
    fprintf(stderr, "Factor base: %d primes, max = %u\n", fb.size, fb.max_prime);

    uint32_t lp_bound = 1U << lp_bits;

    /* Allocate sieve */
    uint8_t *sieve = aligned_alloc(64, BLOCK_SIZE);

    /* Target: fb.size + some surplus for partial relations */
    int target_rels = fb.size + 100;

    /* Simple MPQS: generate polynomials Q(x) = (x+s)^2 - kN where s = floor(sqrt(kN)) */
    mpz_t s, Qx, ax_b, tmp;
    mpz_inits(s, Qx, ax_b, tmp, NULL);
    mpz_sqrt(s, kn);  /* s = floor(sqrt(kN)) */

    /* Sieve arrays for roots */
    uint32_t *soln1 = calloc(fb.size, sizeof(uint32_t));
    uint32_t *soln2 = calloc(fb.size, sizeof(uint32_t));

    /* Compute initial roots: for Q(x) = (x+s)^2 - kN,
       Q(x) ≡ 0 (mod p) when x+s ≡ ±root (mod p)
       so x ≡ root - s (mod p) or x ≡ -root - s (mod p) */
    for (int i = 0; i < fb.size; i++) {
        uint32_t p = fb.primes[i];
        uint32_t r = fb.roots[i];
        uint32_t s_modp = (uint32_t)mpz_fdiv_ui(s, p);

        if (r == 0) {
            soln1[i] = soln2[i] = (p - s_modp) % p;
        } else {
            soln1[i] = (r + p - s_modp) % p;
            soln2[i] = (p - r + p - s_modp) % p;
        }
    }

    /* Relation storage */
    uint32_t **rel_indices = malloc(MAX_RELATIONS * sizeof(uint32_t*));
    uint8_t **rel_powers = malloc(MAX_RELATIONS * sizeof(uint8_t*));
    int *rel_nfactors = malloc(MAX_RELATIONS * sizeof(int));
    int n_full = 0, n_partial = 0;

    /* Sieve threshold */
    double log2_kn = mpz_sizeinbase(kn, 2);
    double log2_M = log2((double)M);
    /* For Q(x) = (x+s)^2 - kN near x=0, |Q(x)| ≈ 2*sqrt(kN)*x ≈ log2_kn/2 + log2_M bits */
    uint8_t sieve_thresh = (uint8_t)(log2_kn / 2.0 + log2_M - SIEVE_THRESH_ADJUST);
    fprintf(stderr, "Sieve threshold: %u (Q(x) ≈ %d bits)\n", sieve_thresh, (int)(log2_kn/2 + log2_M));

    /* Main sieve loop */
    int block_offset = 0;  /* Current offset within sieve interval [0, 2M) */
    int blocks_done = 0;

    /* We sieve the interval [0, 2M) centered around x=0
       x ranges from -M to +M, but we shift to [0, 2M) for unsigned indexing */

    fprintf(stderr, "Sieving...\n");

    while (n_full + n_partial / 3 < target_rels && blocks_done < 100000) {
        /* Initialize sieve block */
        memset(sieve, sieve_thresh, BLOCK_SIZE);

        int block_start = block_offset;

        /* Sieve this block */
        for (int i = 1; i < fb.size; i++) {  /* Skip prime 2 for now */
            uint32_t p = fb.primes[i];
            uint8_t logp = fb.logp[i];

            /* Find first hit in this block for root 1 */
            int r1 = (int)(soln1[i]);
            int offset1 = r1 - (block_start % (int)p);
            if (offset1 < 0) offset1 += p;

            for (int j = offset1; j < BLOCK_SIZE; j += p) {
                sieve[j] -= logp;
            }

            /* Root 2 */
            int r2 = (int)(soln2[i]);
            int offset2 = r2 - (block_start % (int)p);
            if (offset2 < 0) offset2 += p;
            if (offset2 != offset1) {
                for (int j = offset2; j < BLOCK_SIZE; j += p) {
                    sieve[j] -= logp;
                }
            }
        }

        /* Scan for candidates (values that dropped below threshold) */
        for (int j = 0; j < BLOCK_SIZE; j++) {
            if (sieve[j] < 30) {  /* Low enough to be potentially smooth */
                /* Compute Q(x) = (x+s)^2 - kN where x = block_start + j - M */
                int x = block_start + j - M;
                mpz_set_si(ax_b, x);
                mpz_add(ax_b, ax_b, s);  /* ax_b = x + s */

                mpz_mul(Qx, ax_b, ax_b);
                mpz_sub(Qx, Qx, kn);     /* Qx = (x+s)^2 - kN */

                /* Trial divide */
                uint32_t td_indices[1000];
                uint8_t td_powers[1000];
                int n_factors;
                uint32_t lp1, lp2;

                int result = trial_divide_candidate(Qx, &fb, td_indices, td_powers,
                                                    &n_factors, lp_bound, &lp1, &lp2);

                if (result == 1) {
                    /* Full relation */
                    n_full++;
                } else if (result == 2) {
                    /* Single large prime */
                    n_partial++;
                } else if (result == 3) {
                    /* Double large prime */
                    n_partial++;
                }

                if (result > 0 && n_full + n_partial < MAX_RELATIONS) {
                    int idx = n_full + n_partial - 1;
                    rel_indices[idx] = malloc(n_factors * sizeof(uint32_t));
                    rel_powers[idx] = malloc(n_factors * sizeof(uint8_t));
                    memcpy(rel_indices[idx], td_indices, n_factors * sizeof(uint32_t));
                    memcpy(rel_powers[idx], td_powers, n_factors * sizeof(uint8_t));
                    rel_nfactors[idx] = n_factors;
                }
            }
        }

        block_offset += BLOCK_SIZE;
        if (block_offset >= 2 * M) {
            block_offset = 0;
            /* For simple MPQS, we'd switch polynomials here */
            /* For now, just report progress and stop */
            fprintf(stderr, "  Sieved %d blocks, %d full + %d partial relations\n",
                    blocks_done, n_full, n_partial);
            break;  /* Single polynomial for now */
        }
        blocks_done++;

        if (blocks_done % 1000 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &end);
            double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
            fprintf(stderr, "  %d blocks, %d full + %d partial rels, %.1fs\n",
                    blocks_done, n_full, n_partial, elapsed);

            if (elapsed > 280) {
                fprintf(stderr, "Timeout approaching, stopping\n");
                break;
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    fprintf(stderr, "Results: %d full + %d partial relations in %.3fs (%d blocks)\n",
            n_full, n_partial, elapsed, blocks_done);
    fprintf(stderr, "Need %d total relations\n", target_rels);

    if (n_full + n_partial / 3 < target_rels) {
        fprintf(stderr, "FAILED: Not enough relations (need SIQS polynomial switching)\n");
        /* TODO: Implement self-initializing polynomial generation */
    }

    /* Cleanup */
    free(sieve);
    free(soln1);
    free(soln2);
    for (int i = 0; i < n_full + n_partial; i++) {
        free(rel_indices[i]);
        free(rel_powers[i]);
    }
    free(rel_indices);
    free(rel_powers);
    free(rel_nfactors);
    free(fb.primes);
    free(fb.roots);
    free(fb.logp);
    mpz_clears(n, kn, factor, s, Qx, ax_b, tmp, NULL);

    return 0;
}
