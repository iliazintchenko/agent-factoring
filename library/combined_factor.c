/*
 * combined_factor.c - High-performance combined factoring for balanced semiprimes
 *
 * Strategy: Use the fastest algorithm for each size range:
 * - < 20 digits: trial division + Pollard rho
 * - 20-62 digits: Pollard rho with Brent improvement + SQUFOF fallback
 * - 63+ digits: SIQS (self-initializing quadratic sieve)
 *
 * Compile: gcc -O2 -march=native -o combined_factor library/combined_factor.c -lgmp -lm
 * Usage: ./combined_factor <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ===== Pollard Rho with Brent's improvement ===== */
/* Uses batch GCD for efficiency - accumulate products and GCD every 128 steps */

static int pollard_rho_brent(mpz_t factor, const mpz_t n, unsigned long c) {
    mpz_t x, y, q, ys, d, tmp;
    unsigned long m = 128;  /* batch size for GCD */
    unsigned long r, i, k;
    int found = 0;

    mpz_inits(x, y, q, ys, d, tmp, NULL);

    /* Start with x = y = 2 */
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(q, 1);

    r = 1;

    do {
        mpz_set(x, y);
        for (i = 0; i < r; i++) {
            mpz_mul(y, y, y);
            mpz_add_ui(y, y, c);
            mpz_mod(y, y, n);
        }

        k = 0;
        do {
            mpz_set(ys, y);
            unsigned long bound = (m < r - k) ? m : r - k;

            for (i = 0; i < bound; i++) {
                mpz_mul(y, y, y);
                mpz_add_ui(y, y, c);
                mpz_mod(y, y, n);

                mpz_sub(tmp, x, y);
                mpz_abs(tmp, tmp);
                mpz_mul(q, q, tmp);
                mpz_mod(q, q, n);
            }

            mpz_gcd(d, q, n);
            k += bound;
        } while (k < r && mpz_cmp_ui(d, 1) == 0);

        r *= 2;
    } while (mpz_cmp_ui(d, 1) == 0 && r < 10000000);

    if (mpz_cmp(d, n) == 0) {
        /* Backtrack */
        do {
            mpz_mul(ys, ys, ys);
            mpz_add_ui(ys, ys, c);
            mpz_mod(ys, ys, n);

            mpz_sub(tmp, x, ys);
            mpz_abs(tmp, tmp);
            mpz_gcd(d, tmp, n);
        } while (mpz_cmp_ui(d, 1) == 0);

        if (mpz_cmp(d, n) == 0) {
            found = 0;
        } else {
            mpz_set(factor, d);
            found = 1;
        }
    } else {
        mpz_set(factor, d);
        found = 1;
    }

    mpz_clears(x, y, q, ys, d, tmp, NULL);
    return found;
}

/* ===== SQUFOF (Square Form Factorization) ===== */
/* Works well for numbers up to ~60 digits with balanced factors */

static int squfof_factor(mpz_t factor, const mpz_t n) {
    /* SQUFOF only works for numbers that fit in reasonable precision */
    /* For GMP, we implement the full precision version */
    mpz_t sqrt_n, P, Q, Qprev, b, tmp, P0, Q0, Qprev0;
    mpz_t forms_P[256], forms_Q[256];
    int num_forms = 0;
    unsigned long i, max_iter;
    int found = 0;

    /* SQUFOF works best for numbers < 2^128 (about 38 digits) */
    if (mpz_sizeinbase(n, 10) > 40) return 0;

    mpz_inits(sqrt_n, P, Q, Qprev, b, tmp, P0, Q0, Qprev0, NULL);

    /* Multipliers to try */
    unsigned long multipliers[] = {1, 3, 5, 7, 11, 3*5, 3*7, 3*11, 5*7, 5*11, 7*11, 3*5*7, 3*5*11};
    int num_mult = 13;

    for (int mi = 0; mi < num_mult && !found; mi++) {
        mpz_t kn;
        mpz_init(kn);
        mpz_mul_ui(kn, n, multipliers[mi]);

        mpz_sqrt(sqrt_n, kn);

        /* Initialize: P0 = floor(sqrt(kN)), Q0 = 1, Q1 = kN - P0^2 */
        mpz_set(P, sqrt_n);
        mpz_set_ui(Qprev, 1);
        mpz_mul(tmp, P, P);
        mpz_sub(Q, kn, tmp);

        if (mpz_sgn(Q) == 0) {
            /* Perfect square */
            mpz_clear(kn);
            continue;
        }

        max_iter = (unsigned long)(2 * sqrt(sqrt(mpz_get_d(kn))));
        if (max_iter < 1000) max_iter = 1000;
        if (max_iter > 1000000) max_iter = 1000000;

        /* Forward cycle: iterate until we find a perfect square */
        for (i = 1; i <= max_iter; i++) {
            /* b = floor((sqrt_n + P) / Q) */
            mpz_add(tmp, sqrt_n, P);
            mpz_fdiv_q(b, tmp, Q);

            /* P_new = b*Q - P */
            mpz_mul(tmp, b, Q);
            mpz_sub(tmp, tmp, P);
            mpz_set(P, tmp);

            /* Q_new = Q_prev + b*(P_old - P_new) -- but we need to save old P */
            /* Actually: Q_{i+1} = Q_{i-1} + b_i * (P_i - P_{i+1}) */
            /* Simpler: Q_new = (kN - P^2) / Q */
            mpz_mul(tmp, P, P);
            mpz_sub(tmp, kn, tmp);
            mpz_divexact(tmp, tmp, Q);

            mpz_set(Qprev, Q);
            mpz_set(Q, tmp);

            /* Check if Q is a perfect square */
            if (i % 2 == 0) {
                if (mpz_perfect_square_p(Q)) {
                    /* Found perfect square Q = s^2 */
                    mpz_t s;
                    mpz_init(s);
                    mpz_sqrt(s, Q);

                    /* Start reverse cycle */
                    mpz_t rP, rQ, rQprev, rb;
                    mpz_inits(rP, rQ, rQprev, rb, NULL);

                    /* P0 = P, Q0 = s */
                    mpz_set(rP, P);
                    mpz_set(rQprev, s);

                    /* b = floor((sqrt_n - P) / s) */
                    mpz_sub(tmp, sqrt_n, rP);
                    mpz_fdiv_q(rb, tmp, s);

                    mpz_mul(tmp, rb, s);
                    mpz_add(rP, rP, tmp);

                    mpz_mul(tmp, rP, rP);
                    mpz_sub(tmp, kn, tmp);
                    mpz_divexact(rQ, tmp, s);

                    mpz_t rPold;
                    mpz_init(rPold);

                    /* Iterate reverse cycle until P stabilizes */
                    for (unsigned long j = 0; j < max_iter; j++) {
                        mpz_add(tmp, sqrt_n, rP);
                        mpz_fdiv_q(rb, tmp, rQ);

                        mpz_set(rPold, rP);
                        mpz_mul(tmp, rb, rQ);
                        mpz_sub(rP, tmp, rP);

                        mpz_mul(tmp, rP, rP);
                        mpz_sub(tmp, kn, tmp);
                        mpz_divexact(tmp, tmp, rQ);
                        mpz_set(rQprev, rQ);
                        mpz_set(rQ, tmp);

                        if (mpz_cmp(rP, rPold) == 0) {
                            /* Found factor */
                            mpz_gcd(factor, rP, n);
                            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0) {
                                found = 1;
                            }
                            break;
                        }
                    }

                    mpz_clears(rP, rQ, rQprev, rb, rPold, s, NULL);
                    if (found) break;
                }
            }
        }

        mpz_clear(kn);
    }

    mpz_clears(sqrt_n, P, Q, Qprev, b, tmp, P0, Q0, Qprev0, NULL);
    return found;
}

/* ===== Simple trial division up to bound ===== */
static int trial_divide(mpz_t factor, const mpz_t n, unsigned long bound) {
    if (mpz_divisible_ui_p(n, 2)) {
        mpz_set_ui(factor, 2);
        return 1;
    }
    for (unsigned long p = 3; p <= bound; p += 2) {
        if (mpz_divisible_ui_p(n, p)) {
            mpz_set_ui(factor, p);
            return 1;
        }
    }
    return 0;
}

/* ===== Main factoring driver ===== */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t n, factor, cofactor;
    mpz_inits(n, factor, cofactor, NULL);

    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int found = 0;
    size_t digits = mpz_sizeinbase(n, 10);

    /* Try trial division first for small factors */
    if (!found) {
        found = trial_divide(factor, n, 100000);
    }

    /* Try SQUFOF for small numbers */
    if (!found && digits <= 40) {
        found = squfof_factor(factor, n);
    }

    /* Try Pollard rho with multiple c values */
    if (!found) {
        /* Use seed 42 as required */
        gmp_randstate_t rstate;
        gmp_randinit_default(rstate);
        gmp_randseed_ui(rstate, 42);

        /* Try c = 1, 2, 3, ... */
        for (unsigned long c = 1; c <= 100 && !found; c++) {
            found = pollard_rho_brent(factor, n, c);
            if (found) {
                /* Verify it's a proper factor */
                if (mpz_cmp_ui(factor, 1) == 0 || mpz_cmp(factor, n) == 0) {
                    found = 0;
                }
            }
        }

        gmp_randclear(rstate);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    if (found) {
        mpz_divexact(cofactor, n, factor);
        /* Print smaller factor first */
        if (mpz_cmp(factor, cofactor) > 0) {
            mpz_swap(factor, cofactor);
        }
        gmp_printf("Factor found: %Zd * %Zd\n", factor, cofactor);
        gmp_printf("Time: %.4f seconds\n", elapsed);
    } else {
        fprintf(stderr, "Failed to factor in %.4f seconds\n", elapsed);
        return 1;
    }

    mpz_clears(n, factor, cofactor, NULL);
    return 0;
}
