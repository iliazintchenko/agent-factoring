/* Quick test: verify relation correctness on the failing semiprime */
#include <gmp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main() {
    mpz_t N, m, val, res, sq, prod;
    mpz_inits(N, m, val, res, sq, prod, NULL);
    mpz_set_str(N, "2821457049646560553430084268822185147537933587", 10);
    mpz_sqrt(m, N); mpz_add_ui(m, m, 1);

    /* Test: f(0) = m^2 - N */
    mpz_mul(sq, m, m);
    mpz_sub(res, sq, N);
    gmp_printf("m = %Zd\n", m);
    gmp_printf("m^2 - N = %Zd\n", res);
    gmp_printf("|m^2 - N| bits = %d\n", (int)mpz_sizeinbase(res, 2));

    /* Verify: (x+m)^2 mod N for x=1 */
    mpz_set_ui(val, 1);
    mpz_add(val, val, m);
    mpz_mul(sq, val, val);
    mpz_mod(sq, sq, N);
    gmp_printf("(1+m)^2 mod N = %Zd\n", sq);

    mpz_sub(res, val, m); /* actually val = 1+m, so val - m = 1... */
    /* f(1) = (1+m)^2 - N */
    mpz_set_ui(val, 1);
    mpz_add(val, val, m);
    mpz_mul(res, val, val);
    mpz_sub(res, res, N);
    gmp_printf("f(1) = %Zd\n", res);

    /* Verify f(1) mod N == (1+m)^2 mod N */
    mpz_mod(prod, res, N);
    gmp_printf("f(1) mod N = %Zd\n", prod);
    gmp_printf("(1+m)^2 mod N = %Zd\n", sq);
    if (mpz_cmp(prod, sq) == 0) printf("MATCH!\n"); else printf("MISMATCH!\n");

    /* Factor N with GMP's built-in */
    /* Check if N is actually semiprime */
    printf("N has %d digits, %d bits\n", (int)strlen("2821457049646560553430084268822185147537933587"),
           (int)mpz_sizeinbase(N, 2));

    mpz_clears(N, m, val, res, sq, prod, NULL);
    return 0;
}
