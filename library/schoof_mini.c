#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    mpz_t N; mpz_init_set_ui(N, 143);
    
    /* Compute x^143 mod (x^3 + 115*x + 130) mod 143 */
    /* Work in Z/143Z[x]/(x^3+115x+130) */
    /* Elements: a0 + a1*x + a2*x^2 with a0,a1,a2 in Z/143Z */
    
    /* x^2 mod poly = x^2 (degree < 3) */
    /* x^3 mod poly = -115x - 130 (since x^3 = -115x - 130 mod poly) */
    /* x^4 = x * x^3 = x*(-115x - 130) = -115x^2 - 130x mod 143 */
    /* = 28x^2 + 13x mod 143 */
    
    /* Let's just do repeated squaring manually */
    /* State: (c0, c1, c2) representing c0 + c1*x + c2*x^2 */
    long a = 115, b = 130;
    long c0, c1, c2; /* current polynomial */
    long r0, r1, r2; /* result polynomial */
    
    /* result = 1 */
    r0 = 1; r1 = 0; r2 = 0;
    /* base = x */
    c0 = 0; c1 = 1; c2 = 0;
    
    int n = 143;
    for (int bit = 0; bit < 8; bit++) {
        if ((n >> bit) & 1) {
            /* result *= base mod poly mod 143 */
            long t0 = (r0*c0 + r1*c2*(-a) + r2*c1*(-a) + r2*c2*(-b)) % 143;
            long t1 = (r0*c1 + r1*c0 + r1*c2*(-b) + r2*c2*(-a)) % 143; /* Hmm this is getting complex */
            /* Let me just do full multiplication and reduce */
            /* r * c = polynomial of degree up to 4, reduce mod x^3+ax+b */
            long p0 = r0*c0;
            long p1 = r0*c1 + r1*c0;
            long p2 = r0*c2 + r1*c1 + r2*c0;
            long p3 = r1*c2 + r2*c1;
            long p4 = r2*c2;
            /* Reduce x^4: x^4 = x * x^3 = x*(-ax-b) = -ax^2 - bx */
            p2 += p4 * (-a); p1 += p4 * (-b); p4 = 0;
            /* Reduce x^3: x^3 = -ax - b */
            p1 += p3 * (-a); p0 += p3 * (-b); p3 = 0;
            r0 = ((p0 % 143) + 143) % 143;
            r1 = ((p1 % 143) + 143) % 143;
            r2 = ((p2 % 143) + 143) % 143;
        }
        /* base *= base */
        long p0 = c0*c0;
        long p1 = 2*c0*c1;
        long p2 = c1*c1 + 2*c0*c2;
        long p3 = 2*c1*c2;
        long p4 = c2*c2;
        p2 += p4 * (-a); p1 += p4 * (-b); p4 = 0;
        p1 += p3 * (-a); p0 += p3 * (-b); p3 = 0;
        c0 = ((p0 % 143) + 143) % 143;
        c1 = ((p1 % 143) + 143) % 143;
        c2 = ((p2 % 143) + 143) % 143;
        
        printf("bit %d: result = %ld + %ld*x + %ld*x^2, base = %ld + %ld*x + %ld*x^2\n",
               bit, r0, r1, r2, c0, c1, c2);
    }
    
    printf("\nx^143 mod (x^3+115x+130, 143) = %ld + %ld*x + %ld*x^2\n", r0, r1, r2);
    printf("x^143 - x = %ld + %ld*x + %ld*x^2\n", r0, (r1 - 1 + 143) % 143, r2);
    
    /* Check gcd of each coefficient with 143 */
    for (int i = 0; i < 3; i++) {
        long coef[] = {r0, (r1-1+143)%143, r2};
        long g = 143;
        long c_val = coef[i];
        while (c_val > 0) { long t = g % c_val; g = c_val; c_val = t; }
        printf("gcd(%ld, 143) = %ld\n", coef[i], g);
    }
    
    mpz_clear(N);
    return 0;
}
