/*
 * nfs_hensel_sqrt.c - NFS algebraic square root via Hensel lifting
 *
 * Algorithm:
 * 1. Find prime p where f(x) has d distinct roots mod p
 * 2. For each root r_j mod p:
 *    a. Lift root: r_j mod p -> r_j mod p^2 -> ... -> r_j mod p^e (Hensel)
 *    b. Compute S_j = prod(a_i - b_i*r_j) mod p^e
 *    c. Compute T_j = sqrt(S_j) mod p^e (Hensel lift of sqrt from mod p)
 * 3. Lagrange interpolation with lifted roots: T(m) mod p^e
 * 4. When p^e > N, reduce T(m) mod N to get Y
 * 5. gcd(X - Y, N) or gcd(X + Y, N) gives factor
 *
 * The 2^d sign choices are handled by trying all at the base level.
 *
 * This is a standalone test; integrate with gnfs_work.c for full NFS.
 *
 * Compile: gcc -O3 -o nfs_hensel library/nfs_hensel_sqrt.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#define MAX_DEG 6

/* Hensel lift a root: given r with f(r) ≡ 0 (mod p^k),
 * compute r' with f(r') ≡ 0 (mod p^(2k)).
 * r' = r - f(r) * f'(r)^(-1) mod p^(2k) */
static void hensel_lift_root(mpz_t r_out, const mpz_t r_in,
                             mpz_t *coeff, int degree,
                             const mpz_t modulus) {
    mpz_t f_val, fp_val, tmp, inv;
    mpz_inits(f_val, fp_val, tmp, inv, NULL);

    /* Evaluate f(r) mod modulus */
    mpz_set_ui(f_val, 0);
    for (int i = degree; i >= 0; i--) {
        mpz_mul(f_val, f_val, r_in);
        mpz_add(f_val, f_val, coeff[i]);
        mpz_mod(f_val, f_val, modulus);
    }

    /* Evaluate f'(r) mod modulus */
    mpz_set_ui(fp_val, 0);
    for (int i = degree; i >= 1; i--) {
        mpz_mul(fp_val, fp_val, r_in);
        mpz_mul_ui(tmp, coeff[i], i);
        mpz_add(fp_val, fp_val, tmp);
        mpz_mod(fp_val, fp_val, modulus);
    }

    /* r' = r - f(r) * f'(r)^(-1) mod modulus */
    if (mpz_invert(inv, fp_val, modulus) == 0) {
        /* f'(r) not invertible - this shouldn't happen for simple roots */
        mpz_set(r_out, r_in);
    } else {
        mpz_mul(tmp, f_val, inv);
        mpz_sub(r_out, r_in, tmp);
        mpz_mod(r_out, r_out, modulus);
    }

    mpz_clears(f_val, fp_val, tmp, inv, NULL);
}

/* Hensel lift a square root: given t with t^2 ≡ s (mod p^k),
 * compute t' with t'^2 ≡ s (mod p^(2k)).
 * t' = t + (s - t^2) / (2t) mod p^(2k)
 * = (t + s/t) / 2 mod p^(2k)
 * Or: t' = t - (t^2 - s) * (2t)^(-1) mod p^(2k) */
static void hensel_lift_sqrt(mpz_t t_out, const mpz_t t_in,
                             const mpz_t s, const mpz_t modulus) {
    mpz_t t2, diff, inv, two_t;
    mpz_inits(t2, diff, inv, two_t, NULL);

    mpz_mul(t2, t_in, t_in);       /* t^2 */
    mpz_sub(diff, t2, s);          /* t^2 - s */
    mpz_mod(diff, diff, modulus);

    mpz_mul_ui(two_t, t_in, 2);    /* 2t */
    mpz_mod(two_t, two_t, modulus);

    if (mpz_invert(inv, two_t, modulus) == 0) {
        mpz_set(t_out, t_in);
    } else {
        mpz_mul(diff, diff, inv);
        mpz_sub(t_out, t_in, diff);
        mpz_mod(t_out, t_out, modulus);
    }

    mpz_clears(t2, diff, inv, two_t, NULL);
}

/* Tonelli-Shanks for sqrt(n) mod p */
static unsigned long sqrt_mod(unsigned long n, unsigned long p) {
    if (n == 0) return 0;
    unsigned long long nn = n%p, m = p;
    unsigned long long r=1, b=nn, e=(p-1)/2;
    while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; }
    if (r != 1) return 0;
    if (p%4==3) { r=1; b=nn; e=(p+1)/4; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } return (unsigned long)r; }
    unsigned long Q2=p-1, S=0; while (Q2%2==0) { Q2/=2; S++; }
    unsigned long z=2;
    for (;;) { r=1; b=z; e=(p-1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } if (r==m-1) break; z++; }
    unsigned long long M2=S,c,t,R2;
    r=1; b=z; e=Q2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } c=r;
    r=1; b=nn; e=Q2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } t=r;
    r=1; b=nn; e=(Q2+1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } R2=r;
    for (;;) { if (t==1) return (unsigned long)R2; int i2=0; unsigned long long tt=t;
        while (tt!=1) { tt=tt*tt%p; i2++; }
        unsigned long long bb=c; for (int j=0; j<(int)M2-i2-1; j++) bb=bb*bb%p;
        M2=i2; c=bb*bb%p; t=t*c%p; R2=R2*bb%p; }
}

/* Compute the algebraic square root Y = T(m) mod N
 * Given: polynomial f of degree d, m with f(m) = N,
 * relations (a_i, b_i), rational square root X */
int algebraic_sqrt(mpz_t Y_out, mpz_t *coeff, int degree, mpz_t m,
                   mpz_t N, mpz_t X_rat,
                   long *a_vals, unsigned long *b_vals, int nrels,
                   int *sign_combo /* which of 2^d sign combos to use */) {
    /* Find a prime p where f has d distinct roots */
    unsigned long p = 0;
    unsigned long roots[MAX_DEG];
    int nroots = 0;

    for (unsigned long q = 100003; q < 200000; q += 2) {
        mpz_t qz; mpz_init_set_ui(qz, q);
        if (!mpz_probab_prime_p(qz, 3)) { mpz_clear(qz); continue; }
        mpz_clear(qz);

        int nr = 0;
        for (unsigned long x = 0; x < q && nr <= degree; x++) {
            unsigned long long val = 0;
            for (int j = degree; j >= 0; j--) {
                unsigned long cj = mpz_fdiv_ui(coeff[j], q);
                val = (val * x + cj) % q;
            }
            if (val == 0) roots[nr++] = x;
        }
        if (nr == degree) { p = q; nroots = nr; break; }
    }

    if (p == 0) {
        fprintf(stderr, "Could not find suitable prime for Hensel lifting\n");
        return -1;
    }
    fprintf(stderr, "  Hensel prime p=%lu with %d roots\n", p, nroots);

    /* For each root r_j:
     * 1. Compute S_j = prod(a_i - b_i*r_j) mod p
     * 2. Compute T_j = sqrt(S_j) mod p (with sign choice)
     * 3. Hensel lift both r_j and T_j to mod p^e where p^e > N */

    /* Determine lifting target */
    int bits_needed = mpz_sizeinbase(N, 2) + 10;
    int bits_per_lift = 0;
    { unsigned long pp = p; while (pp > 0) { bits_per_lift++; pp >>= 1; } }
    int num_lifts = 0;
    int current_bits = bits_per_lift;
    while (current_bits < bits_needed) { current_bits *= 2; num_lifts++; }

    fprintf(stderr, "  Need %d Hensel lifts (p=%lu, %d bits -> %d bits)\n",
            num_lifts, p, bits_per_lift, bits_needed);

    mpz_t lifted_roots[MAX_DEG];
    mpz_t lifted_sqrts[MAX_DEG];
    mpz_t S_vals[MAX_DEG];
    mpz_t modulus, tmp;
    mpz_inits(modulus, tmp, NULL);

    for (int j = 0; j < degree; j++) {
        mpz_init_set_ui(lifted_roots[j], roots[j]);
        mpz_init(lifted_sqrts[j]);
        mpz_init(S_vals[j]);
    }

    /* Step 1: Compute S_j = prod(a_i - b_i*r_j) mod p and initial sqrt */
    for (int j = 0; j < degree; j++) {
        unsigned long long prod_val = 1;
        for (int i = 0; i < nrels; i++) {
            long long term = ((a_vals[i] % (long long)p) + p) % p;
            term = (term + p - ((unsigned long long)b_vals[i] % p * roots[j]) % p) % p;
            prod_val = prod_val * (unsigned long long)term % p;
        }
        mpz_set_ui(S_vals[j], (unsigned long)prod_val);

        unsigned long sq = sqrt_mod((unsigned long)prod_val, p);
        if (sq == 0 && prod_val != 0) {
            fprintf(stderr, "  sqrt(%lu) mod %lu does not exist for root %d!\n",
                    (unsigned long)prod_val, p, j);
            return -1;
        }

        /* Apply sign choice */
        if (*sign_combo & (1 << j))
            mpz_set_ui(lifted_sqrts[j], p - sq);
        else
            mpz_set_ui(lifted_sqrts[j], sq);
    }

    /* Step 2: Hensel lift roots and square roots */
    mpz_set_ui(modulus, p);

    for (int lift = 0; lift < num_lifts; lift++) {
        mpz_t new_modulus;
        mpz_init(new_modulus);
        mpz_mul(new_modulus, modulus, modulus); /* p^(2^(lift+1)) */

        /* Lift each root r_j from mod modulus to mod new_modulus */
        for (int j = 0; j < degree; j++) {
            hensel_lift_root(lifted_roots[j], lifted_roots[j],
                            coeff, degree, new_modulus);
        }

        /* Recompute S_j mod new_modulus using lifted roots */
        for (int j = 0; j < degree; j++) {
            mpz_set_ui(S_vals[j], 1);
            for (int i = 0; i < nrels; i++) {
                mpz_set_si(tmp, a_vals[i]);
                mpz_submul_ui(tmp, lifted_roots[j], b_vals[i]);
                mpz_mod(tmp, tmp, new_modulus);
                mpz_mul(S_vals[j], S_vals[j], tmp);
                mpz_mod(S_vals[j], S_vals[j], new_modulus);
            }
        }

        /* Lift each sqrt T_j from mod modulus to mod new_modulus */
        for (int j = 0; j < degree; j++) {
            hensel_lift_sqrt(lifted_sqrts[j], lifted_sqrts[j],
                            S_vals[j], new_modulus);
        }

        mpz_set(modulus, new_modulus);
        mpz_clear(new_modulus);
    }

    /* Step 3: Lagrange interpolation to compute T(m) mod modulus */
    /* T(m) = sum_j T_j * prod_{k!=j} (m - r_k) / (r_j - r_k) */
    mpz_t m_lifted;
    mpz_init(m_lifted);
    mpz_mod(m_lifted, m, modulus);

    mpz_set_ui(Y_out, 0);
    for (int j = 0; j < degree; j++) {
        mpz_t num, den, inv;
        mpz_inits(num, den, inv, NULL);

        mpz_set(num, lifted_sqrts[j]);
        mpz_set_ui(den, 1);

        for (int k = 0; k < degree; k++) {
            if (k == j) continue;
            /* num *= (m - r_k) */
            mpz_sub(tmp, m_lifted, lifted_roots[k]);
            mpz_mod(tmp, tmp, modulus);
            mpz_mul(num, num, tmp);
            mpz_mod(num, num, modulus);

            /* den *= (r_j - r_k) */
            mpz_sub(tmp, lifted_roots[j], lifted_roots[k]);
            mpz_mod(tmp, tmp, modulus);
            mpz_mul(den, den, tmp);
            mpz_mod(den, den, modulus);
        }

        mpz_invert(inv, den, modulus);
        mpz_mul(num, num, inv);
        mpz_mod(num, num, modulus);

        mpz_add(Y_out, Y_out, num);
        mpz_mod(Y_out, Y_out, modulus);

        mpz_clears(num, den, inv, NULL);
    }

    /* Reduce mod N */
    mpz_mod(Y_out, Y_out, N);

    /* Cleanup */
    for (int j = 0; j < degree; j++) {
        mpz_clear(lifted_roots[j]);
        mpz_clear(lifted_sqrts[j]);
        mpz_clear(S_vals[j]);
    }
    mpz_clears(modulus, tmp, m_lifted, NULL);

    return 0;
}

/* Test with a simple example */
int main(int argc, char *argv[]) {
    printf("NFS Hensel lifting algebraic square root module\n");
    printf("Use algebraic_sqrt() function with gnfs_work.c\n");
    return 0;
}
