/* Debug: generate relations and check one dependency in detail */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* (... reuse the simple QS code but add MPQS and verify deps ...) */
/* Actually, let me just use the LGSH binary with more verbose output */
/* Instead, let me test: take two fully smooth QS relations and combine them */

int main() {
    mpz_t N, m;
    mpz_inits(N, m, NULL);
    mpz_set_str(N, "202930257038936745993812575957", 10);
    mpz_sqrt(m, N); mpz_add_ui(m, m, 1);

    /* Manually find two smooth values and verify extraction */
    /* f(x) = (x+m)^2 - N. Check x = 1, 2, 3, ... for small smooth values */
    for (long x = -100000; x < 100000; x++) {
        mpz_t val, res;
        mpz_inits(val, res, NULL);
        mpz_set_si(val, x); mpz_add(val, val, m);
        mpz_mul(res, val, val); mpz_sub(res, res, N);
        if (mpz_sgn(res) < 0) mpz_neg(res, res);

        /* Check if res is 1000-smooth */
        mpz_t tmp; mpz_init_set(tmp, res);
        for (int p = 2; p <= 1000; p++) {
            while (mpz_divisible_ui_p(tmp, p)) mpz_divexact_ui(tmp, tmp, p);
        }
        if (mpz_cmp_ui(tmp, 1) == 0) {
            gmp_printf("x=%ld: sv=%Zd, f(x)=%Zd (smooth)\n", x, val, res);
        }
        mpz_clears(val, res, tmp, NULL);
        
        /* Just find first few */
        static int count = 0;
        if (mpz_cmp_ui(val, 1) == 0) count = 0; /* reset count check - wrong logic */
    }

    mpz_clears(N, m, NULL);
    return 0;
}
