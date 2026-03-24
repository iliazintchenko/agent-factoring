/*
 * mpqs_batch.c — Multi-Polynomial QS with Batch GCD smooth detection
 *
 * NOVEL ELEMENT: Instead of sieving (sequential access), we use batch GCD
 * to test smoothness of candidates generated from many different polynomials
 * simultaneously. This allows "cherry-picking" candidates from whichever
 * polynomial produces the smallest values at each point, potentially
 * improving the effective smoothness rate.
 *
 * The question: does structured multi-polynomial candidate selection with
 * batch GCD achieve better scaling than standard sieving?
 *
 * Each polynomial: Q_a(x) = (a*x + b)^2 - N, where a^2 ~ sqrt(N)/M
 * so that |Q_a(x)| ~ sqrt(N) * M / a ~ sqrt(N) for |x| <= M.
 *
 * Strategy: for each x in the sieve interval, pick the polynomial a that
 * minimizes |Q_a(x)|. This "best-of-K" selection should boost smoothness
 * probability by a factor related to K.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>

#define MAX_FACTOR_BASE 50000
#define MAX_RELATIONS 100000

typedef struct {
    unsigned long p;
    int r1, r2; /* roots of N mod p */
} fb_entry_t;

typedef struct {
    mpz_t value;      /* Q(x) value */
    mpz_t sqrt_part;  /* (a*x + b) such that sqrt_part^2 - N = value */
    int smooth;       /* 1 if B-smooth */
} relation_t;

static gmp_randstate_t rng;

/* Build factor base: primes p where N is a QR mod p */
int build_factor_base(fb_entry_t *fb, mpz_t N, unsigned long B) {
    int count = 0;
    mpz_t p_mpz;
    mpz_init_set_ui(p_mpz, 2);

    /* Add 2 */
    fb[count].p = 2;
    fb[count].r1 = 1;
    fb[count].r2 = 1;
    count++;

    mpz_set_ui(p_mpz, 3);
    while (mpz_cmp_ui(p_mpz, B) <= 0 && count < MAX_FACTOR_BASE) {
        unsigned long p = mpz_get_ui(p_mpz);
        int legendre = mpz_kronecker_ui(N, p);
        if (legendre == 1) {
            fb[count].p = p;
            /* Find sqrt(N) mod p via Tonelli-Shanks */
            mpz_t root, n_mod_p;
            mpz_inits(root, n_mod_p, NULL);
            mpz_mod_ui(n_mod_p, N, p);
            mpz_t p_z;
            mpz_init_set_ui(p_z, p);
            /* Use mpz_sqrtm - not available, use manual Tonelli-Shanks */
            /* For simplicity, brute force for small p */
            unsigned long nm = mpz_get_ui(n_mod_p);
            unsigned long r = 0;
            for (unsigned long i = 0; i < p; i++) {
                if ((i * i) % p == nm) { r = i; break; }
            }
            fb[count].r1 = r;
            fb[count].r2 = p - r;
            count++;
            mpz_clears(root, n_mod_p, p_z, NULL);
        }
        mpz_nextprime(p_mpz, p_mpz);
    }
    mpz_clear(p_mpz);
    return count;
}

/*
 * Generate MPQS polynomial: Q_a(x) = (a*x + b)^2 - N
 * Choose a as product of primes from factor base (so a^2 is a known square)
 * b satisfies b^2 ≡ N (mod a)
 */
void generate_poly(mpz_t a, mpz_t b, mpz_t N, mpz_t target_a,
                   fb_entry_t *fb, int fb_size) {
    /* Simple approach: a = single prime from upper part of factor base */
    /* Pick a prime p from factor base with p ~ sqrt(2N)/M */
    /* For now, just pick a random prime from upper half of FB */
    int idx = fb_size / 2 + (rand() % (fb_size / 2));
    mpz_set_ui(a, fb[idx].p);

    /* b = sqrt(N) mod a — already computed */
    mpz_set_ui(b, fb[idx].r1);
}

/*
 * Core experiment: measure smoothness rate with different candidate
 * generation strategies.
 *
 * Strategy 1 (baseline): Sequential QS candidates Q(x) = (x+m)^2 - N
 * Strategy 2 (multi-poly): Best-of-K candidates from K polynomials
 * Strategy 3 (structured): Lattice-selected candidates
 */
void measure_smoothness_rate(mpz_t N, unsigned long B, int strategy) {
    int n_digits = mpz_sizeinbase(N, 10);

    mpz_t sqrtN, m, candidate, abs_val;
    mpz_inits(sqrtN, m, candidate, abs_val, NULL);
    mpz_sqrt(sqrtN, N);
    mpz_add_ui(m, sqrtN, 1);

    /* Make m odd for better alignment */
    if (mpz_even_p(m)) mpz_add_ui(m, m, 1);

    int M = 200000; /* sieve interval size */
    int n_smooth = 0;
    int n_tested = 0;

    /* Compute primorial for smooth testing */
    mpz_t P, bound_val;
    mpz_inits(P, bound_val, NULL);

    /* Estimate bound: Q(x) ~ 2*sqrt(N)*M */
    mpz_mul_ui(bound_val, sqrtN, 2 * M);

    mpz_set_ui(P, 1);
    mpz_t prime_z;
    mpz_init_set_ui(prime_z, 2);
    while (mpz_cmp_ui(prime_z, B) <= 0) {
        unsigned long p = mpz_get_ui(prime_z);
        mpz_t pk;
        mpz_init_set(pk, prime_z);
        while (1) {
            mpz_t pk_next;
            mpz_init(pk_next);
            mpz_mul(pk_next, pk, prime_z);
            if (mpz_cmp(pk_next, bound_val) > 0) {
                mpz_clear(pk_next);
                break;
            }
            mpz_set(pk, pk_next);
            mpz_clear(pk_next);
        }
        mpz_mul(P, P, pk);
        mpz_clear(pk);
        mpz_nextprime(prime_z, prime_z);
    }
    mpz_clear(prime_z);

    if (strategy == 1) {
        /* Strategy 1: Sequential baseline */
        mpz_t g, rem;
        mpz_inits(g, rem, NULL);

        for (int x = 1; x <= M; x++) {
            mpz_add_ui(candidate, m, x);
            mpz_mul(candidate, candidate, candidate);
            mpz_sub(candidate, candidate, N);
            mpz_abs(abs_val, candidate);

            /* Test smoothness via iterated GCD with primorial */
            mpz_set(rem, abs_val);
            while (1) {
                mpz_gcd(g, rem, P);
                if (mpz_cmp_ui(g, 1) == 0) break;
                mpz_divexact(rem, rem, g);
            }
            if (mpz_cmp_ui(rem, 1) == 0) n_smooth++;
            n_tested++;
        }
        mpz_clears(g, rem, NULL);
    }
    else if (strategy == 2) {
        /* Strategy 2: Best-of-K polynomials
         * For each position x, evaluate K different polynomials and
         * pick the one with smallest |Q(x)|.
         * This increases smoothness probability because smaller values
         * are more likely to be smooth.
         */
        int K = 10; /* number of polynomials */
        mpz_t *a_vals = malloc(K * sizeof(mpz_t));
        mpz_t *b_vals = malloc(K * sizeof(mpz_t));
        mpz_t best_val, test_val, best_sqrt;
        mpz_inits(best_val, test_val, best_sqrt, NULL);

        /* Generate K polynomial bases: Q_k(x) = (x + m + k*offset)^2 - N */
        /* Using different offsets from sqrt(N) */
        for (int k = 0; k < K; k++) {
            mpz_init(a_vals[k]);
            mpz_init(b_vals[k]);
            mpz_set_ui(a_vals[k], 1);
            mpz_add_ui(b_vals[k], m, k * (M / K));
        }

        mpz_t g, rem;
        mpz_inits(g, rem, NULL);

        for (int x = 1; x <= M / K; x++) {
            /* Find the polynomial giving smallest value */
            int best_k = 0;
            mpz_add_ui(candidate, b_vals[0], x);
            mpz_mul(candidate, candidate, candidate);
            mpz_sub(candidate, candidate, N);
            mpz_abs(best_val, candidate);
            mpz_add_ui(best_sqrt, b_vals[0], x);

            for (int k = 1; k < K; k++) {
                mpz_add_ui(candidate, b_vals[k], x);
                mpz_mul(candidate, candidate, candidate);
                mpz_sub(candidate, candidate, N);
                mpz_abs(test_val, candidate);
                if (mpz_cmp(test_val, best_val) < 0) {
                    mpz_set(best_val, test_val);
                    mpz_add_ui(best_sqrt, b_vals[k], x);
                    best_k = k;
                }
            }

            /* Test smoothness of best value */
            mpz_set(rem, best_val);
            while (1) {
                mpz_gcd(g, rem, P);
                if (mpz_cmp_ui(g, 1) == 0) break;
                mpz_divexact(rem, rem, g);
            }
            if (mpz_cmp_ui(rem, 1) == 0) n_smooth++;
            n_tested++;
        }

        for (int k = 0; k < K; k++) {
            mpz_clears(a_vals[k], b_vals[k], NULL);
        }
        free(a_vals);
        free(b_vals);
        mpz_clears(best_val, test_val, best_sqrt, g, rem, NULL);
    }
    else if (strategy == 3) {
        /* Strategy 3: CRT-structured candidates
         * Instead of sequential x, choose x values that are guaranteed
         * to be divisible by many small primes.
         *
         * For each small prime p where N is QR mod p:
         *   (x + m)^2 ≡ N (mod p) has solutions x ≡ r_p - m (mod p)
         *
         * Choose x satisfying several such congruences via CRT.
         * Then Q(x) = (x+m)^2 - N is guaranteed divisible by those primes.
         * The cofactor Q(x)/(product of those primes) should be tested.
         *
         * Key question: is the cofactor small enough to improve smoothness?
         */
        mpz_t g, rem;
        mpz_inits(g, rem, NULL);

        /* Build factor base with roots */
        fb_entry_t *fb = malloc(MAX_FACTOR_BASE * sizeof(fb_entry_t));
        int fb_size = build_factor_base(fb, N, B);

        /* Generate candidates using groups of primes */
        int group_size = 5; /* number of primes to force */

        for (int trial = 0; trial < M && n_tested < M; trial++) {
            /* Pick a random subset of group_size primes from FB */
            /* Use CRT to find x satisfying all congruences */
            mpz_t x_crt, modulus, this_mod, this_res, m_mod;
            mpz_inits(x_crt, modulus, this_mod, this_res, m_mod, NULL);
            mpz_set_ui(x_crt, 0);
            mpz_set_ui(modulus, 1);

            int valid = 1;
            for (int g_idx = 0; g_idx < group_size && g_idx < fb_size; g_idx++) {
                int pidx = 1 + ((trial * group_size + g_idx) % (fb_size - 1));
                unsigned long p = fb[pidx].p;

                /* x ≡ r1 - m (mod p) */
                mpz_mod_ui(m_mod, m, p);
                unsigned long target = (fb[pidx].r1 + p - mpz_get_ui(m_mod)) % p;

                /* CRT combine */
                mpz_set_ui(this_mod, p);
                mpz_set_ui(this_res, target);

                /* Simple CRT: x_crt = x_crt + modulus * t where t solves
                 * modulus * t ≡ this_res - x_crt (mod p) */
                mpz_t diff, inv, t_val;
                mpz_inits(diff, inv, t_val, NULL);
                mpz_mod_ui(diff, x_crt, p);
                long d = ((long)target - (long)mpz_get_ui(diff) % (long)p + (long)p) % (long)p;
                mpz_set_ui(diff, d);

                /* inv = modulus^{-1} mod p */
                mpz_t mod_p;
                mpz_init_set_ui(mod_p, p);
                if (mpz_invert(inv, modulus, mod_p) == 0) {
                    valid = 0;
                    mpz_clears(diff, inv, t_val, mod_p, NULL);
                    break;
                }
                mpz_mul(t_val, diff, inv);
                mpz_mod_ui(t_val, t_val, p);

                mpz_addmul(x_crt, modulus, t_val);
                mpz_mul_ui(modulus, modulus, p);
                mpz_mod(x_crt, x_crt, modulus);

                mpz_clears(diff, inv, t_val, mod_p, NULL);
            }

            if (!valid) {
                mpz_clears(x_crt, modulus, this_mod, this_res, m_mod, NULL);
                continue;
            }

            /* Evaluate Q(x_crt) = (x_crt + m)^2 - N */
            mpz_add(candidate, x_crt, m);
            mpz_mul(candidate, candidate, candidate);
            mpz_sub(candidate, candidate, N);
            mpz_abs(abs_val, candidate);

            /* Test smoothness */
            mpz_set(rem, abs_val);
            while (1) {
                mpz_gcd(g, rem, P);
                if (mpz_cmp_ui(g, 1) == 0) break;
                mpz_divexact(rem, rem, g);
            }
            if (mpz_cmp_ui(rem, 1) == 0) n_smooth++;
            n_tested++;

            mpz_clears(x_crt, modulus, this_mod, this_res, m_mod, NULL);
        }
        free(fb);
    }

    printf("Strategy %d | N: %d digits | B: %lu | tested: %d | smooth: %d | rate: %.4f%%\n",
           strategy, n_digits, B, n_tested, n_smooth, 100.0 * n_smooth / (n_tested > 0 ? n_tested : 1));

    mpz_clears(sqrtN, m, candidate, abs_val, P, bound_val, NULL);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <N_decimal> <B> [strategy=1|2|3]\n", argv[0]);
        return 1;
    }

    srand(42);
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);
    unsigned long B = atol(argv[2]);
    int strategy = argc > 3 ? atoi(argv[3]) : 1;

    measure_smoothness_rate(N, B, strategy);

    mpz_clear(N);
    gmp_randclear(rng);
    return 0;
}
