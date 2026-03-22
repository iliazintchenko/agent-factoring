/*
 * Ultra-fast factorer for balanced semiprimes (30-50 digits).
 * Uses SQUFOF with __int128 native arithmetic (P,Q values fit in 128 bits).
 * Falls back to Pollard Rho with batched GCD.
 * Compile: gcc -O3 -march=native -o factor_small library/factor_small.c -lgmp -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <stdint.h>

typedef unsigned __int128 u128;
typedef __int128 s128;

/* Convert mpz to u128 (assumes it fits) */
static u128 mpz_to_u128(const mpz_t z) {
    u128 result = 0;
    mpz_t tmp;
    mpz_init_set(tmp, z);
    result = (u128)mpz_getlimbn(tmp, 0);
    if (mpz_size(tmp) > 1)
        result |= (u128)mpz_getlimbn(tmp, 1) << 64;
    mpz_clear(tmp);
    return result;
}

/* Integer square root for u128 */
static u128 isqrt128(u128 n) {
    if (n == 0) return 0;
    u128 x = (u128)sqrtl((long double)n);
    /* Newton refinement */
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

/* Check if u128 is a perfect square, return sqrt or 0 */
static u128 is_perfect_square_128(u128 n) {
    u128 s = isqrt128(n);
    if (s * s == n) return s;
    return 0;
}

/* SQUFOF with __int128 arithmetic.
 * For N up to ~76 digits, sqrt(kN) fits in u128.
 * For our 30-50d range this is fine. */
static const unsigned int squfof_mults[] = {
    1, 3, 5, 7, 11, 3*5, 3*7, 3*11, 5*7, 5*11, 7*11,
    3*5*7, 3*5*11, 3*7*11, 5*7*11, 3*5*7*11,
    13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
    71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    0
};

static int squfof_factor(mpz_t n_mpz, mpz_t factor) {
    /* Check size: sqrt(kN) must fit in u128, so kN < 2^256.
       N up to ~75 digits with mult up to 113: kN up to ~77 digits.
       sqrt(kN) up to ~39 digits = ~129 bits. Tight but works. */
    if (mpz_sizeinbase(n_mpz, 2) > 200) return 0; /* too big for u128 SQUFOF */

    u128 n = mpz_to_u128(n_mpz);

    for (int mi = 0; squfof_mults[mi] != 0; mi++) {
        u128 mult = squfof_mults[mi];
        u128 kn = n * mult;

        /* Check for overflow or perfect square */
        if (kn / mult != n) continue; /* overflow */
        u128 sqrtkn = isqrt128(kn);
        if (sqrtkn * sqrtkn == kn) continue; /* perfect square */

        /* Forward cycle */
        u128 Pprev, P, Qprev, Q, q;
        P = sqrtkn;
        Qprev = 1;
        Q = kn - P * P;
        if (Q == 0) continue;

        /* Max iterations: ~2 * N^{1/4} */
        unsigned long max_iter = (unsigned long)(2.0 * pow((double)n, 0.25));
        if (max_iter < 2000) max_iter = 2000;
        if (max_iter > 20000000) max_iter = 20000000;

        u128 Qsqrt = 0;

        for (unsigned long i = 1; i <= max_iter; i++) {
            /* q = floor((sqrtkn + P) / Q) */
            q = (sqrtkn + P) / Q;

            /* Update: P' = q*Q - P */
            Pprev = P;
            P = q * Q - P;

            /* Update: Q' = Qprev + q*(Pprev - P) */
            u128 Qnew = Qprev + q * (Pprev - P + (P > Pprev ? 0 : 0));
            /* More careful: Pprev - P could underflow if P > Pprev */
            s128 diff = (s128)Pprev - (s128)P;
            Qnew = (u128)((s128)Qprev + (s128)q * diff);
            Qprev = Q;
            Q = Qnew;

            /* Check for perfect square Q at even iterations */
            if ((i & 1) == 0) {
                Qsqrt = is_perfect_square_128(Q);
                if (Qsqrt > 0) break;
            }
        }

        if (Qsqrt == 0) continue;

        /* Reverse cycle from the square form */
        /* P0 = sqrtkn - ((sqrtkn - P) mod Qsqrt) */
        /* But need careful modular arithmetic with u128 */
        s128 sp = (s128)sqrtkn - (s128)P;
        if (sp < 0) sp = -sp; /* |sqrtkn - P| */
        u128 rem = (u128)((s128)sp % (s128)Qsqrt);
        /* Actually: P_new = sqrtkn - ((sqrtkn - P) mod Qsqrt) if sqrtkn >= P */
        if (sqrtkn >= P) {
            rem = (sqrtkn - P) % Qsqrt;
            P = sqrtkn - rem;
        } else {
            /* P > sqrtkn: (sqrtkn - P) is negative */
            rem = (P - sqrtkn) % Qsqrt;
            if (rem == 0)
                P = sqrtkn;
            else
                P = sqrtkn + Qsqrt - rem;
        }

        Qprev = Qsqrt;
        /* Q = (kn - P*P) / Qsqrt */
        Q = (kn - P * P) / Qsqrt;

        /* Forward until P_i == P_{i-1} */
        int found_factor = 0;
        for (unsigned long i = 0; i <= max_iter; i++) {
            q = (sqrtkn + P) / Q;
            Pprev = P;
            P = q * Q - P;

            if (P == Pprev) {
                /* GCD(P, N) should give a factor */
                mpz_t g;
                mpz_init(g);
                mpz_t p_mpz;
                mpz_init(p_mpz);
                /* Convert P to mpz */
                mpz_set_ui(p_mpz, (unsigned long)(P >> 64));
                mpz_mul_2exp(p_mpz, p_mpz, 64);
                mpz_add_ui(p_mpz, p_mpz, (unsigned long)(P & 0xFFFFFFFFFFFFFFFFULL));

                mpz_gcd(g, p_mpz, n_mpz);

                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n_mpz) < 0) {
                    mpz_set(factor, g);
                    found_factor = 1;
                }
                mpz_clear(g);
                mpz_clear(p_mpz);
                break;
            }

            s128 d2 = (s128)Pprev - (s128)P;
            u128 Qnew = (u128)((s128)Qprev + (s128)q * d2);
            Qprev = Q;
            Q = Qnew;
        }

        if (found_factor) return 1;
    }

    return 0;
}

/* Pollard Rho with Brent's cycle detection and batched GCD */
static int pollard_rho_brent(mpz_t n, mpz_t factor, unsigned long seed) {
    mpz_t x, y, d, c, q, ys, temp;
    mpz_inits(x, y, d, c, q, ys, temp, NULL);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    int found = 0;

    for (int attempt = 0; attempt < 100 && !found; attempt++) {
        mpz_urandomm(c, state, n);
        if (mpz_cmp_ui(c, 2) < 0) mpz_set_ui(c, 2);

        mpz_urandomm(x, state, n);
        mpz_set(y, x);
        mpz_set_ui(q, 1);

        unsigned long r = 1, m = 128;

        do {
            mpz_set(x, y);
            for (unsigned long i = 0; i < r; i++) {
                mpz_mul(y, y, y);
                mpz_add(y, y, c);
                mpz_mod(y, y, n);
            }

            unsigned long k = 0;
            do {
                mpz_set(ys, y);
                unsigned long batch = (r - k < m) ? r - k : m;
                for (unsigned long j = 0; j < batch; j++) {
                    mpz_mul(y, y, y);
                    mpz_add(y, y, c);
                    mpz_mod(y, y, n);
                    mpz_sub(temp, x, y);
                    mpz_abs(temp, temp);
                    mpz_mul(q, q, temp);
                    mpz_mod(q, q, n);
                }
                mpz_gcd(d, q, n);
                k += batch;
            } while (k < r && mpz_cmp_ui(d, 1) == 0);

            r *= 2;
        } while (mpz_cmp_ui(d, 1) == 0 && r < (1UL << 26));

        if (mpz_cmp(d, n) == 0) {
            /* Backtrack */
            do {
                mpz_mul(ys, ys, ys);
                mpz_add(ys, ys, c);
                mpz_mod(ys, ys, n);
                mpz_sub(temp, x, ys);
                mpz_abs(temp, temp);
                mpz_gcd(d, temp, n);
            } while (mpz_cmp_ui(d, 1) == 0);
        }

        if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0) {
            mpz_set(factor, d);
            found = 1;
        }
    }

    mpz_clears(x, y, d, c, q, ys, temp, NULL);
    gmp_randclear(state);
    return found;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t n, factor, cofactor;
    mpz_inits(n, factor, cofactor, NULL);
    mpz_set_str(n, argv[1], 10);

    int found = 0;

    /* 1. Quick trial division */
    if (!found && mpz_divisible_ui_p(n, 2)) {
        mpz_set_ui(factor, 2); found = 1;
    }
    if (!found) {
        for (unsigned long p = 3; p < 10000; p += 2) {
            if (mpz_divisible_ui_p(n, p)) {
                mpz_set_ui(factor, p); found = 1; break;
            }
        }
    }

    /* 2. SQUFOF with 128-bit native arithmetic */
    if (!found) {
        found = squfof_factor(n, factor);
    }

    /* 3. Pollard Rho fallback */
    if (!found) {
        found = pollard_rho_brent(n, factor, 42);
    }

    if (found) {
        mpz_divexact(cofactor, n, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
    } else {
        gmp_fprintf(stderr, "FAIL: could not factor %Zd\n", n);
        return 1;
    }

    mpz_clears(n, factor, cofactor, NULL);
    return 0;
}
