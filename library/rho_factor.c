#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

// Pollard's rho with Brent's improvement
// Fast for small numbers where one factor is << sqrt(N)
// For balanced semiprimes, this isn't ideal but for <60d it's fast enough

static void rho_factor(mpz_t n, mpz_t factor) {
    mpz_t x, y, d, c, temp;
    mpz_inits(x, y, d, c, temp, NULL);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 42);
    
    while (1) {
        mpz_urandomm(c, state, n);
        if (mpz_cmp_ui(c, 0) == 0) mpz_set_ui(c, 1);
        
        mpz_set_ui(x, 2);
        mpz_set_ui(y, 2);
        mpz_set_ui(d, 1);
        
        unsigned long m = 1;
        while (mpz_cmp_ui(d, 1) == 0) {
            mpz_set(x, y);
            for (unsigned long i = 0; i < m; i++) {
                mpz_mul(y, y, y);
                mpz_add(y, y, c);
                mpz_mod(y, y, n);
            }
            
            mpz_set(temp, y);
            unsigned long k = 0;
            while (k < m && mpz_cmp_ui(d, 1) == 0) {
                unsigned long batch = (m - k < 128) ? m - k : 128;
                mpz_set_ui(d, 1);
                for (unsigned long j = 0; j < batch; j++) {
                    mpz_mul(temp, temp, temp);
                    mpz_add(temp, temp, c);
                    mpz_mod(temp, temp, n);
                    
                    mpz_sub(factor, temp, x);
                    mpz_abs(factor, factor);
                    mpz_mul(d, d, factor);
                    mpz_mod(d, d, n);
                }
                mpz_gcd(d, d, n);
                k += batch;
                
                if (mpz_cmp(d, n) == 0) break;
            }
            
            if (mpz_cmp(d, n) == 0) {
                mpz_set_ui(d, 1);
                break; // retry with different c
            }
            
            m *= 2;
        }
        
        if (mpz_cmp_ui(d, 1) != 0 && mpz_cmp(d, n) != 0) {
            mpz_set(factor, d);
            mpz_clears(x, y, d, c, temp, NULL);
            gmp_randclear(state);
            return;
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    
    mpz_t n, factor, cofactor;
    mpz_inits(n, factor, cofactor, NULL);
    mpz_set_str(n, argv[1], 10);
    
    rho_factor(n, factor);
    mpz_divexact(cofactor, n, factor);
    
    if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
    
    gmp_printf("%Zd %Zd\n", factor, cofactor);
    
    mpz_clears(n, factor, cofactor, NULL);
    return 0;
}
