/*
 * Combined factoring: trial division + Pollard rho + ECM
 * Usage: ./factor <N> [deadline_seconds]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

#define TRIAL_LIMIT 100000

static struct timespec t0;
static double deadline = 290.0;

static double now_sec(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (t.tv_sec - t0.tv_sec) + (t.tv_nsec - t0.tv_nsec) / 1e9;
}

/* Trial division up to TRIAL_LIMIT */
static int trial_div(mpz_t n, mpz_t f) {
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; }
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; }
    for (unsigned long d = 5; d <= TRIAL_LIMIT; d += 6) {
        if (mpz_divisible_ui_p(n, d))     { mpz_set_ui(f, d);   return 1; }
        if (mpz_divisible_ui_p(n, d + 2)) { mpz_set_ui(f, d+2); return 1; }
    }
    return 0;
}

/* Pollard rho-Brent, limited iterations based on factor size estimate */
static int rho_brent(mpz_t n, mpz_t factor, unsigned long max_iter) {
    mpz_t x, y, c, d, q, ys, tmp;
    mpz_inits(x, y, c, d, q, ys, tmp, NULL);
    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, (unsigned long)time(NULL) ^ (unsigned long)clock() ^ 12345);

    int found = 0;
    for (int att = 0; att < 10 && !found; att++) {
        mpz_urandomm(c, rng, n); if (mpz_sgn(c)==0) mpz_set_ui(c,1);
        mpz_urandomm(y, rng, n);
        mpz_set(x, y);
        mpz_set_ui(q, 1);
        unsigned long r = 1, m = 128;

        for (unsigned long iter = 0; iter < max_iter && !found;) {
            mpz_set(x, y);
            for (unsigned long i = 0; i < r; i++) {
                mpz_mul(y,y,y); mpz_add(y,y,c); mpz_mod(y,y,n);
            }
            for (unsigned long k = 0; k < r && !found; k += m) {
                mpz_set(ys, y);
                unsigned long bound = (m < r-k) ? m : r-k;
                for (unsigned long i = 0; i < bound; i++) {
                    mpz_mul(y,y,y); mpz_add(y,y,c); mpz_mod(y,y,n);
                    mpz_sub(tmp,x,y); mpz_abs(tmp,tmp);
                    mpz_mul(q,q,tmp); mpz_mod(q,q,n);
                }
                mpz_gcd(d, q, n);
                if (mpz_cmp_ui(d,1)>0 && mpz_cmp(d,n)<0) {
                    mpz_set(factor, d); found = 1;
                } else if (mpz_cmp(d,n)==0) {
                    mpz_set(y, ys);
                    do {
                        mpz_mul(y,y,y); mpz_add(y,y,c); mpz_mod(y,y,n);
                        mpz_sub(tmp,x,y); mpz_abs(tmp,tmp);
                        mpz_gcd(d, tmp, n);
                    } while (mpz_cmp_ui(d,1)==0);
                    if (mpz_cmp(d,n)<0) { mpz_set(factor,d); found=1; }
                    break;
                }
                iter += bound;
            }
            r *= 2;
        }
    }
    mpz_clears(x,y,c,d,q,ys,tmp, NULL);
    gmp_randclear(rng);
    return found;
}

/* ECM with parameters tuned for balanced semiprimes */
static int try_ecm(mpz_t n, mpz_t factor, int ndigits) {
    int factor_digits = (ndigits + 1) / 2;

    /* B1 and expected curves from GMP-ECM documentation */
    struct { int digs; double b1; int curves; } tbl[] = {
        {15, 2000, 25},
        {20, 11000, 90},
        {25, 50000, 300},
        {30, 250000, 700},
        {35, 1000000, 1800},
        {40, 11000000, 5100},
        {45, 43000000, 10600},
        {50, 110000000, 19300},
        {55, 260000000, 42000},
    };

    double b1 = 260000000;
    int curves = 50000;
    for (int i = 0; i < 9; i++) {
        if (factor_digits <= tbl[i].digs) {
            b1 = tbl[i].b1; curves = tbl[i].curves;
            break;
        }
    }

    mpz_t f;
    mpz_init(f);
    int found = 0;

    for (int i = 0; i < curves && !found; i++) {
        ecm_params p;
        ecm_init(p);
        p->B1done = 1.0;
        mpz_set_ui(p->sigma, 0);
        p->method = ECM_ECM;

        mpz_set(f, n);
        int ret = ecm_factor(f, f, b1, p);
        ecm_clear(p);

        if (ret > 0 && mpz_cmp_ui(f,1) > 0 && mpz_cmp(f,n) < 0) {
            mpz_set(factor, f);
            found = 1;
        }
        if (now_sec() > deadline) break;
    }
    mpz_clear(f);
    return found;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N> [deadline]\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &t0);
    if (argc >= 3) deadline = atof(argv[2]);

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);
    mpz_set_str(n, argv[1], 10);
    int ndigits = strlen(argv[1]);

    /* 1. Trial division */
    if (trial_div(n, factor)) {
        gmp_printf("%Zd\n", factor);
        goto done;
    }

    /* 2. Quick Pollard rho — good for factors up to ~15 digits */
    /* For balanced semiprimes, only useful for ~30-digit numbers */
    if (ndigits <= 38) {
        unsigned long iters = (ndigits <= 32) ? 500000 : 2000000;
        if (rho_brent(n, factor, iters)) {
            gmp_printf("%Zd\n", factor);
            goto done;
        }
    }

    /* 3. ECM */
    if (try_ecm(n, factor, ndigits)) {
        gmp_printf("%Zd\n", factor);
        goto done;
    }

    fprintf(stderr, "FAIL: %s (%d digits) %.1fs\n", argv[1], ndigits, now_sec());
    mpz_clears(n, factor, NULL);
    return 1;

done:
    mpz_clears(n, factor, NULL);
    return 0;
}
