/*
 * schoof_factor2.c — Factor via Schoof-like torsion analysis.
 *
 * For E: y^2 = x^3 + ax + b mod N, compute x^N mod (x^3+ax+b, N).
 * The coefficients of the result, minus x, reveal factors via GCD.
 *
 * Cost: O(log N) polynomial multiplications in degree-2 ring = O(log^3 N) bits.
 * Success probability per curve: ~1/2 for ℓ=2.
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Element of Z/NZ[x]/(x^3+ax+b): represented as c0 + c1*x + c2*x^2 */
typedef struct { mpz_t c[3]; } elem_t;

void elem_init(elem_t *e) { for (int i=0;i<3;i++) mpz_init(e->c[i]); }
void elem_clear(elem_t *e) { for (int i=0;i<3;i++) mpz_clear(e->c[i]); }
void elem_set(elem_t *d, const elem_t *s) { for(int i=0;i<3;i++) mpz_set(d->c[i],s->c[i]); }

/* Multiply two elements in Z/NZ[x]/(x^3+ax+b), where x^3 = -ax - b */
void elem_mul(elem_t *res, const elem_t *u, const elem_t *v,
              const mpz_t a, const mpz_t b, const mpz_t N) {
    mpz_t p[5], tmp;
    for (int i=0;i<5;i++) mpz_init(p[i]);
    mpz_init(tmp);

    /* Polynomial multiply: degree up to 4 */
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++) {
            mpz_mul(tmp, u->c[i], v->c[j]);
            mpz_add(p[i+j], p[i+j], tmp);
        }

    /* Reduce x^4 = x * x^3 = x*(-a*x - b) = -a*x^2 - b*x */
    /* p[4]*x^4 → -a*p[4]*x^2 - b*p[4]*x */
    mpz_mul(tmp, p[4], a); mpz_sub(p[2], p[2], tmp);
    mpz_mul(tmp, p[4], b); mpz_sub(p[1], p[1], tmp);

    /* Reduce x^3 = -a*x - b */
    /* p[3]*x^3 → -a*p[3]*x - b*p[3] */
    mpz_mul(tmp, p[3], a); mpz_sub(p[1], p[1], tmp);
    mpz_mul(tmp, p[3], b); mpz_sub(p[0], p[0], tmp);

    for (int i=0;i<3;i++) { mpz_mod(res->c[i], p[i], N); }

    for (int i=0;i<5;i++) mpz_clear(p[i]);
    mpz_clear(tmp);
}

/* Compute x^exp in Z/NZ[x]/(x^3+ax+b) by repeated squaring */
void elem_powx(elem_t *res, const mpz_t exp, const mpz_t a, const mpz_t b, const mpz_t N) {
    elem_t base, tmp;
    elem_init(&base); elem_init(&tmp);

    /* res = 1 */
    mpz_set_ui(res->c[0], 1); mpz_set_ui(res->c[1], 0); mpz_set_ui(res->c[2], 0);
    /* base = x */
    mpz_set_ui(base.c[0], 0); mpz_set_ui(base.c[1], 1); mpz_set_ui(base.c[2], 0);

    size_t bits = mpz_sizeinbase(exp, 2);
    for (size_t i = 0; i < bits; i++) {
        if (mpz_tstbit(exp, i)) {
            elem_mul(&tmp, res, &base, a, b, N);
            elem_set(res, &tmp);
        }
        if (i + 1 < bits) {
            elem_mul(&tmp, &base, &base, a, b, N);
            elem_set(&base, &tmp);
        }
    }

    elem_clear(&base); elem_clear(&tmp);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N; mpz_init(N); mpz_set_str(N, argv[1], 10);
    size_t ndig = mpz_sizeinbase(N, 10);
    fprintf(stderr, "Schoof-factor: N = %zu digits\n", ndig);

    struct timespec t0; clock_gettime(CLOCK_MONOTONIC, &t0);

    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, 42);

    int found = 0;
    mpz_t factor, a, b, g;
    mpz_init(factor); mpz_init(a); mpz_init(b); mpz_init(g);

    for (int curve = 0; curve < 10000 && !found; curve++) {
        mpz_urandomm(a, rng, N);
        mpz_urandomm(b, rng, N);

        /* Check discriminant */
        mpz_t disc; mpz_init(disc);
        mpz_mul(disc, a, a); mpz_mul(disc, disc, a); mpz_mul_ui(disc, disc, 4);
        mpz_t t27; mpz_init(t27);
        mpz_mul(t27, b, b); mpz_mul_ui(t27, t27, 27);
        mpz_add(disc, disc, t27);
        mpz_mod(disc, disc, N);
        mpz_clear(t27);
        if (mpz_sgn(disc) == 0) { mpz_clear(disc); continue; }

        /* Check if disc is coprime to N */
        mpz_gcd(g, disc, N);
        mpz_clear(disc);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_set(factor, g); found = 1;
            fprintf(stderr, "Factor from discriminant! curve=%d\n", curve);
            break;
        }

        /* Compute x^N mod (x^3 + ax + b, N) */
        elem_t res; elem_init(&res);
        elem_powx(&res, N, a, b, N);

        /* Check x^N - x: subtract x from c[1] */
        mpz_sub_ui(res.c[1], res.c[1], 1);
        mpz_mod(res.c[1], res.c[1], N);

        /* Check GCD of each coefficient with N */
        for (int i = 0; i < 3; i++) {
            mpz_gcd(g, res.c[i], N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                found = 1;
                fprintf(stderr, "Factor from x^N-x coeff[%d]! curve=%d\n", i, curve);
                break;
            }
        }

        elem_clear(&res);

        struct timespec tn; clock_gettime(CLOCK_MONOTONIC, &tn);
        double el = (tn.tv_sec-t0.tv_sec)+(tn.tv_nsec-t0.tv_nsec)/1e9;
        if (el > 280) { fprintf(stderr, "TIMEOUT at curve %d\n", curve); break; }
        if ((curve+1) % 500 == 0) fprintf(stderr, "  %d curves (%.1fs)\n", curve+1, el);
    }

    struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;

    if (found) {
        mpz_t q; mpz_init(q); mpz_divexact(q, N, factor);
        gmp_printf("%Zd %Zd\n", factor, q);
        fprintf(stderr, "SUCCESS in %.3fs\n", total);
        mpz_clear(q);
    } else {
        fprintf(stderr, "FAIL after %.3fs\n", total);
    }

    mpz_clear(factor); mpz_clear(a); mpz_clear(b); mpz_clear(g);
    gmp_randclear(rng); mpz_clear(N);
    return found ? 0 : 1;
}
