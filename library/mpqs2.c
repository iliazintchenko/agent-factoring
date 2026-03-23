/*
 * mpqs2.c - Multiple Polynomial Quadratic Sieve (correct implementation)
 * Usage: ./mpqs2 <number>
 * Single-threaded, seed=42.
 *
 * Uses polynomials: g(x) = (a*x + b)^2 - N where a = q^2 for a prime q.
 * g(x) is divisible by q^2. We check R(x) = g(x)/q^2 for smoothness.
 * Relation: (ax+b)^2 ≡ a * R(x) (mod N).
 * Since a = q^2 is a perfect square, we just need R(x) to factor into
 * perfect-square products over the factor base.
 *
 * Hensel lifting: given r with r^2 ≡ N (mod q), lift to b with b^2 ≡ N (mod q^2).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_FB 25000
#define MAX_RELS 300000

/* Tonelli-Shanks */
static int mpz_sqrtmod(mpz_t r, const mpz_t a, const mpz_t p) {
    if (mpz_sgn(a) == 0) { mpz_set_ui(r, 0); return 1; }
    mpz_t t, q, z, M, c, R, b2;
    mpz_inits(t, q, z, M, c, R, b2, NULL);
    mpz_sub_ui(t, p, 1); mpz_fdiv_q_2exp(t, t, 1);
    mpz_powm(t, a, t, p);
    if (mpz_cmp_ui(t, 1) != 0) { mpz_clears(t,q,z,M,c,R,b2,NULL); return 0; }
    mpz_sub_ui(q, p, 1);
    unsigned long s = 0;
    while (mpz_even_p(q)) { mpz_fdiv_q_2exp(q, q, 1); s++; }
    if (s == 1) {
        mpz_add_ui(t, p, 1); mpz_fdiv_q_2exp(t, t, 2);
        mpz_powm(r, a, t, p);
        mpz_clears(t,q,z,M,c,R,b2,NULL); return 1;
    }
    mpz_set_ui(z, 2);
    while (1) {
        mpz_sub_ui(t, p, 1); mpz_fdiv_q_2exp(t, t, 1);
        mpz_powm(t, z, t, p);
        if (mpz_cmp_ui(t, 1) != 0) break;
        mpz_add_ui(z, z, 1);
    }
    mpz_set_ui(M, s);
    mpz_powm(c, z, q, p);
    mpz_add_ui(t, q, 1); mpz_fdiv_q_2exp(t, t, 1);
    mpz_powm(R, a, t, p);
    mpz_powm(t, a, q, p);
    while (mpz_cmp_ui(t, 1) != 0) {
        unsigned long i = 0;
        mpz_set(b2, t);
        while (mpz_cmp_ui(b2, 1) != 0) { mpz_mul(b2,b2,b2); mpz_mod(b2,b2,p); i++; }
        unsigned long m = mpz_get_ui(M);
        mpz_set(b2, c);
        for (unsigned long j = 0; j < m-i-1; j++) { mpz_mul(b2,b2,b2); mpz_mod(b2,b2,p); }
        mpz_set_ui(M, i);
        mpz_mul(c, b2, b2); mpz_mod(c, c, p);
        mpz_mul(R, R, b2); mpz_mod(R, R, p);
        mpz_mul(t, t, c); mpz_mod(t, t, p);
    }
    mpz_set(r, R);
    mpz_clears(t,q,z,M,c,R,b2,NULL);
    return 1;
}

/* Factor base */
static unsigned int fb_p[MAX_FB];
static int fb_root[MAX_FB];
static double fb_logp[MAX_FB];
static int fb_size;

typedef struct {
    mpz_t x;       /* ax + b where a = q^2 */
    mpz_t q;       /* the polynomial prime q (a = q^2) */
    int *exps;     /* [sign, fb[0], ...], size fb_size+1 */
    unsigned long lp;
} rel_t;

static rel_t rels[MAX_RELS];
static int num_rels, num_smooth, num_partial;
static mpz_t N;
static gmp_randstate_t rstate;

static void build_fb(int target) {
    fb_size = 0;
    fb_p[0] = 2; fb_root[0] = 1; fb_logp[0] = log(2.0); fb_size++;
    unsigned int lim = 2000000;
    char *sv = (char *)calloc(lim, 1);
    for (unsigned int i = 2; (unsigned long)i*i < lim; i++)
        if (!sv[i]) for (unsigned int j = i*i; j < lim; j += i) sv[j] = 1;
    for (unsigned int p = 3; p < lim && fb_size < target; p += 2) {
        if (sv[p]) continue;
        unsigned int nm = mpz_fdiv_ui(N, p);
        mpz_t az, pz, rz;
        mpz_inits(az, pz, rz, NULL);
        mpz_set_ui(az, nm); mpz_set_ui(pz, p);
        if (mpz_sqrtmod(rz, az, pz)) {
            fb_p[fb_size] = p;
            fb_root[fb_size] = (int)mpz_get_ui(rz);
            fb_logp[fb_size] = log((double)p);
            fb_size++;
        }
        mpz_clears(az, pz, rz, NULL);
    }
    free(sv);
}

static void get_params(int digits, int *fb_target, int *half_M) {
    if (digits <= 30)      { *fb_target = 100;   *half_M = 25000; }
    else if (digits <= 35) { *fb_target = 200;   *half_M = 50000; }
    else if (digits <= 40) { *fb_target = 400;   *half_M = 100000; }
    else if (digits <= 45) { *fb_target = 700;   *half_M = 150000; }
    else if (digits <= 50) { *fb_target = 1200;  *half_M = 250000; }
    else if (digits <= 55) { *fb_target = 2000;  *half_M = 400000; }
    else if (digits <= 60) { *fb_target = 3000;  *half_M = 600000; }
    else if (digits <= 65) { *fb_target = 4500;  *half_M = 1000000; }
    else if (digits <= 70) { *fb_target = 6000;  *half_M = 1500000; }
    else if (digits <= 75) { *fb_target = 8000;  *half_M = 2000000; }
    else if (digits <= 80) { *fb_target = 10000; *half_M = 3000000; }
    else if (digits <= 85) { *fb_target = 13000; *half_M = 4000000; }
    else if (digits <= 90) { *fb_target = 16000; *half_M = 5000000; }
    else if (digits <= 95) { *fb_target = 20000; *half_M = 6000000; }
    else                   { *fb_target = 24000; *half_M = 8000000; }
}

int main(int argc, char *argv[]) {
    if (argc != 2) { fprintf(stderr, "Usage: %s <number>\n", argv[0]); return 1; }

    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);

    mpz_init_set_str(N, argv[1], 10);
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    int digits = strlen(argv[1]);
    int fb_target, half_M;
    get_params(digits, &fb_target, &half_M);
    build_fb(fb_target);

    int sieve_len = 2 * half_M;
    int target_rels = fb_size + 30;

    fprintf(stderr, "MPQS2: %d digits, FB=%d (max_p=%u), M=%d, need=%d rels\n",
            digits, fb_size, fb_p[fb_size-1], half_M, target_rels);

    unsigned long lp_bound = (unsigned long)fb_p[fb_size-1] * 300;
    unsigned char *sieve_arr = (unsigned char *)malloc(sieve_len);

    /* target_a = sqrt(2N) / sieve_len. We use a = q^2, so q ≈ sqrt(target_a). */
    mpz_t target_q, tmp, tmp2;
    mpz_inits(target_q, tmp, tmp2, NULL);
    mpz_mul_ui(tmp, N, 2);
    mpz_sqrt(tmp, tmp);
    mpz_fdiv_q_ui(tmp, tmp, sieve_len);  /* tmp = target_a */
    mpz_sqrt(target_q, tmp);             /* target_q = sqrt(target_a) */

    fprintf(stderr, "MPQS2: target_q ≈ %lu bits\n", (unsigned long)mpz_sizeinbase(target_q, 2));

    num_rels = num_smooth = num_partial = 0;
    int polys = 0, found = 0;

    while (num_smooth < target_rels && polys < 2000000) {
        mpz_t q, q2, b, r0, s_val, nmod;
        mpz_inits(q, q2, b, r0, s_val, nmod, NULL);

        /* Pick prime q near target_q */
        if (polys == 0) {
            mpz_set(q, target_q);
        } else {
            unsigned long bits = mpz_sizeinbase(target_q, 2);
            mpz_urandomb(q, rstate, bits > 1 ? bits : 1);
            mpz_fdiv_q_2exp(tmp, target_q, 1);
            mpz_add(q, q, tmp);
        }
        mpz_nextprime(q, q);

        /* Check N is QR mod q and compute r0 = sqrt(N) mod q */
        mpz_mod(nmod, N, q);
        if (!mpz_sqrtmod(r0, nmod, q)) {
            mpz_nextprime(q, q);
            mpz_mod(nmod, N, q);
            if (!mpz_sqrtmod(r0, nmod, q)) {
                polys++;
                mpz_clears(q,q2,b,r0,s_val,nmod,NULL);
                continue;
            }
        }

        /* Hensel lift: b^2 ≡ N (mod q^2) from r0^2 ≡ N (mod q)
         * s = (r0^2 - N) / q
         * t = -s * (2*r0)^{-1} (mod q)
         * b = r0 + t * q (mod q^2)
         */
        mpz_mul(q2, q, q);

        mpz_mul(tmp, r0, r0);
        mpz_sub(tmp, tmp, N);
        /* tmp = r0^2 - N, divisible by q */
        mpz_divexact(s_val, tmp, q);

        /* inv2r0 = (2*r0)^{-1} mod q */
        mpz_mul_ui(tmp, r0, 2);
        mpz_mod(tmp, tmp, q);
        mpz_t inv2r0;
        mpz_init(inv2r0);
        if (!mpz_invert(inv2r0, tmp, q)) {
            polys++;
            mpz_clears(q,q2,b,r0,s_val,nmod,inv2r0,NULL);
            continue;
        }

        /* t = -s * inv2r0 mod q */
        mpz_neg(tmp, s_val);
        mpz_mod(tmp, tmp, q);
        mpz_mul(tmp, tmp, inv2r0);
        mpz_mod(tmp, tmp, q);

        /* b = r0 + t * q */
        mpz_mul(b, tmp, q);
        mpz_add(b, b, r0);
        mpz_mod(b, b, q2);

        mpz_clear(inv2r0);

        /* Verify b^2 ≡ N (mod q^2) */
        mpz_mul(tmp, b, b);
        mpz_sub(tmp, tmp, N);
        if (!mpz_divisible_p(tmp, q2)) {
            polys++;
            mpz_clears(q,q2,b,r0,s_val,nmod,NULL);
            continue;
        }

        polys++;

        /* a = q^2 */
        /* g(x) = (a*(x-M) + b)^2 - N = (q^2*(x-M) + b)^2 - N
         * R(x) = g(x) / q^2 = q^2*(x-M)^2 + 2*b*(x-M) + (b^2-N)/q^2
         *
         * For FB prime p: g(x) ≡ 0 (mod p) when q^2*(x-M) + b ≡ ±r (mod p)
         * where r = sqrt(N) mod p.
         * x ≡ M + (±r - b) * (q^2)^{-1} (mod p)
         */

        memset(sieve_arr, 0, sieve_len);

        unsigned int q2_small = 0; /* q^2 mod p, but might overflow if q > 2^16 */

        for (int i = 1; i < fb_size; i++) {
            unsigned int p = fb_p[i];

            unsigned int qm = mpz_fdiv_ui(q, p);
            unsigned int q2m = (unsigned long long)qm * qm % p;
            unsigned int bm = mpz_fdiv_ui(b, p);

            if (q2m == 0) continue;

            /* (q^2)^{-1} mod p */
            mpz_t qi, pp;
            mpz_inits(qi, pp, NULL);
            mpz_set_ui(qi, q2m); mpz_set_ui(pp, p);
            if (!mpz_invert(qi, qi, pp)) { mpz_clears(qi, pp, NULL); continue; }
            unsigned int q2inv = mpz_get_ui(qi);
            mpz_clears(qi, pp, NULL);

            unsigned int r = fb_root[i];
            unsigned int Mm = half_M % p;

            /* x1 = M + (r - b) * q2_inv (mod p)
             * x2 = M + (-r - b) * q2_inv (mod p) */
            unsigned int x1 = ((unsigned long long)((r + p - bm) % p) * q2inv + Mm) % p;
            unsigned int x2 = ((unsigned long long)((p - r + p - bm) % p) * q2inv + Mm) % p;

            unsigned char lp = (unsigned char)(fb_logp[i] + 0.5);
            for (unsigned int x = x1; x < (unsigned int)sieve_len; x += p) sieve_arr[x] += lp;
            if (x1 != x2)
                for (unsigned int x = x2; x < (unsigned int)sieve_len; x += p) sieve_arr[x] += lp;
        }

        /* Threshold: max |R(x)| at edges ≈ q^2 * M^2 */
        double log_Rmax = 2.0 * log(mpz_get_d(q)) + 2.0 * log((double)half_M);
        unsigned char thresh = (unsigned char)(log_Rmax * 0.55);
        if (thresh < 10) thresh = 10;

        int poly_smooth = 0, poly_cand = 0;

        for (int xi = 0; xi < sieve_len; xi++) {
            if (sieve_arr[xi] < thresh) continue;
            poly_cand++;

            long sx = (long)xi - (long)half_M;
            mpz_t axb, gv, rv;
            mpz_inits(axb, gv, rv, NULL);

            /* axb = q^2 * sx + b */
            mpz_mul_si(axb, q2, sx);
            mpz_add(axb, axb, b);

            /* gv = axb^2 - N */
            mpz_mul(gv, axb, axb);
            mpz_sub(gv, gv, N);

            if (mpz_sgn(gv) == 0) {
                mpz_t g, cof;
                mpz_inits(g, cof, NULL);
                mpz_gcd(g, axb, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_divexact(cof, N, g);
                    gmp_printf("%Zd %Zd\n", g, cof);
                    found = 1;
                }
                mpz_clears(g, cof, axb, gv, rv, NULL);
                if (found) goto done;
                continue;
            }

            int sign = 0;
            mpz_set(rv, gv);
            if (mpz_sgn(rv) < 0) { mpz_neg(rv, rv); sign = 1; }

            /* Divide out q^2 */
            if (!mpz_divisible_p(rv, q2)) {
                mpz_clears(axb, gv, rv, NULL);
                continue;
            }
            mpz_divexact(rv, rv, q2);

            /* Trial divide by FB */
            int *exps = (int *)calloc(fb_size + 1, sizeof(int));
            exps[0] = sign;

            for (int i = 0; i < fb_size; i++) {
                while (mpz_divisible_ui_p(rv, fb_p[i])) {
                    mpz_fdiv_q_ui(rv, rv, fb_p[i]);
                    exps[i+1]++;
                }
            }

            if (mpz_cmp_ui(rv, 1) == 0) {
                if (num_rels < MAX_RELS) {
                    mpz_init_set(rels[num_rels].x, axb);
                    mpz_init_set(rels[num_rels].q, q);
                    rels[num_rels].exps = exps;
                    rels[num_rels].lp = 0;
                    num_rels++; num_smooth++; poly_smooth++;
                } else free(exps);
            } else if (mpz_fits_ulong_p(rv) && mpz_get_ui(rv) <= lp_bound) {
                if (num_rels < MAX_RELS) {
                    mpz_init_set(rels[num_rels].x, axb);
                    mpz_init_set(rels[num_rels].q, q);
                    rels[num_rels].exps = exps;
                    rels[num_rels].lp = mpz_get_ui(rv);
                    num_rels++; num_partial++;
                } else free(exps);
            } else {
                free(exps);
            }

            mpz_clears(axb, gv, rv, NULL);
        }

        if (polys <= 20 || polys % 200 == 0 || num_smooth == target_rels) {
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            double el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)/1e9;
            fprintf(stderr, "poly %d: %d cand, %d sm (total: %d/%d sm, %d part, %.1fs)\n",
                    polys, poly_cand, poly_smooth, num_smooth, target_rels, num_partial, el);
        }

        mpz_clears(q,q2,b,r0,s_val,nmod,NULL);
    }

    {
        clock_gettime(CLOCK_MONOTONIC, &ts1);
        double el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)/1e9;
        fprintf(stderr, "MPQS2: sieve: %d smooth, %d partial, %d polys, %.1fs\n",
                num_smooth, num_partial, polys, el);
    }

    /* Combine large prime partials */
    if (num_smooth < target_rels) {
        for (int i = 0; i < num_rels && num_smooth < target_rels; i++) {
            if (rels[i].lp == 0 || rels[i].lp == (unsigned long)-1) continue;
            for (int j = i+1; j < num_rels; j++) {
                if (rels[j].lp == rels[i].lp) {
                    for (int k = 0; k <= fb_size; k++)
                        rels[i].exps[k] += rels[j].exps[k];
                    mpz_mul(rels[i].x, rels[i].x, rels[j].x);
                    rels[i].lp = 0; rels[j].lp = (unsigned long)-1;
                    num_smooth++; break;
                }
            }
        }
        fprintf(stderr, "MPQS2: after LP combining: %d smooth\n", num_smooth);
    }

    if (num_smooth < fb_size + 1) {
        fprintf(stderr, "MPQS2: FAILED - %d < %d\n", num_smooth, fb_size + 1);
        return 1;
    }

    /* Build GF(2) matrix */
    int mc = fb_size + 1, mr = num_smooth;
    int wm = (mc+63)/64, wh = (mr+63)/64, wt = wm + wh;
    unsigned long **mat = (unsigned long **)malloc(sizeof(unsigned long *) * mr);
    int *smap = (int *)malloc(sizeof(int) * mr);
    int si = 0;

    for (int i = 0; i < num_rels && si < mr; i++) {
        if (rels[i].lp != 0) continue;
        smap[si] = i;
        mat[si] = (unsigned long *)calloc(wt, sizeof(unsigned long));
        mat[si][wm + si/64] |= (1UL << (si%64));
        for (int j = 0; j <= fb_size; j++)
            if (rels[i].exps[j] & 1)
                mat[si][j/64] |= (1UL << (j%64));
        si++;
    }

    /* Gauss */
    int cr = 0;
    for (int col = 0; col < mc && cr < mr; col++) {
        int piv = -1;
        for (int r = cr; r < mr; r++)
            if ((mat[r][col/64] >> (col%64)) & 1) { piv = r; break; }
        if (piv < 0) continue;
        if (piv != cr) { unsigned long *t = mat[cr]; mat[cr] = mat[piv]; mat[piv] = t; }
        for (int r = 0; r < mr; r++)
            if (r != cr && ((mat[r][col/64] >> (col%64)) & 1))
                for (int w = 0; w < wt; w++) mat[r][w] ^= mat[cr][w];
        cr++;
    }

    fprintf(stderr, "MPQS2: rank=%d, nullity=%d\n", cr, mr-cr);

    mpz_t xp, yp, g, d;
    mpz_inits(xp, yp, g, d, NULL);

    for (int r = cr; r < mr && !found; r++) {
        int zr = 1;
        for (int w = 0; w < wm && zr; w++) if (mat[r][w]) zr = 0;
        if (!zr) continue;

        mpz_set_ui(xp, 1);
        mpz_set_ui(yp, 1);  /* y accumulates q_j * product of p^(e/2) */
        int *te = (int *)calloc(fb_size + 1, sizeof(int));

        for (int i = 0; i < mr; i++) {
            if (!((mat[r][wm + i/64] >> (i%64)) & 1)) continue;
            int ri = smap[i];
            mpz_mul(xp, xp, rels[ri].x);
            mpz_mod(xp, xp, N);
            /* Include q in y product: (ax+b)^2 ≡ q^2 * R(x) (mod N) */
            mpz_mul(yp, yp, rels[ri].q);
            mpz_mod(yp, yp, N);
            for (int j = 0; j <= fb_size; j++) te[j] += rels[ri].exps[j];
        }

        /* Multiply in the FB prime contributions: product of p^(e/2) */
        for (int j = 1; j <= fb_size; j++) {
            if (te[j] <= 0) continue;
            mpz_t pw; mpz_init(pw);
            mpz_ui_pow_ui(pw, fb_p[j-1], te[j] / 2);
            mpz_mul(yp, yp, pw);
            mpz_mod(yp, yp, N);
            mpz_clear(pw);
        }

        mpz_sub(d, xp, yp); mpz_gcd(g, d, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t cof; mpz_init(cof);
            mpz_divexact(cof, N, g);
            clock_gettime(CLOCK_MONOTONIC, &ts1);
            double el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)/1e9;
            gmp_printf("%Zd %Zd\n", g, cof);
            fprintf(stderr, "MPQS2: factored in %.3fs\n", el);
            found = 1;
            mpz_clear(cof);
        }
        free(te);
    }

done:
    if (!found) {
        clock_gettime(CLOCK_MONOTONIC, &ts1);
        double el = (ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)/1e9;
        fprintf(stderr, "MPQS2: LA failed, %.3fs\n", el);
    }

    return found ? 0 : 1;
}
