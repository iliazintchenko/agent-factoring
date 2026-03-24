/*
 * qs_engine.c — Clean QS-style factoring engine for scaling measurement
 *
 * Focus: reliable factoring with clean code, supporting novel variations.
 * Uses sieving + trial division + LP matching + GF(2) linear algebra.
 *
 * NOT a production QS — minimal, focused on measuring scaling behavior.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_FB_SIZE 20000
#define MAX_RELS 300000
#define MAX_FACTORS_PER_REL 80

/* Factor base entry */
typedef struct {
    unsigned long p;
    unsigned long r1, r2; /* roots of N mod p */
    double logp;
} fb_t;

/* Relation: (x+m)^2 - N = sign * product(p_i^e_i) * product(LP_j) */
typedef struct {
    int idx;           /* original index */
    mpz_t xm;         /* x + m */
    int sign;          /* -1 flag */
    int nf;            /* number of small factors */
    int primes[MAX_FACTORS_PER_REL]; /* indices into factor base */
    int exps[MAX_FACTORS_PER_REL];
    int nlp;
    unsigned long lp[2]; /* large primes */
} rel_t;

/* --- Tonelli-Shanks --- */
unsigned long modsqrt(unsigned long n, unsigned long p) {
    if (p == 2) return n & 1;
    if (n == 0) return 0;
    n %= p;
    if (n == 0) return 0;

    /* Check QR */
    mpz_t b, m, r;
    mpz_inits(b, m, r, NULL);
    mpz_set_ui(b, n); mpz_set_ui(m, p);
    mpz_powm_ui(r, b, (p-1)/2, m);
    if (mpz_cmp_ui(r, 1) != 0) { mpz_clears(b, m, r, NULL); return 0; }

    if (p % 4 == 3) {
        mpz_powm_ui(r, b, (p+1)/4, m);
        unsigned long res = mpz_get_ui(r);
        mpz_clears(b, m, r, NULL);
        return res;
    }

    /* General Tonelli-Shanks */
    unsigned long Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned long z = 2;
    mpz_t zz, tmp;
    mpz_inits(zz, tmp, NULL);
    while (z < p) {
        mpz_set_ui(zz, z);
        mpz_powm_ui(tmp, zz, (p-1)/2, m);
        if (mpz_cmp_ui(tmp, p-1) == 0) break;
        z++;
    }

    mpz_t c, t, R, bb;
    mpz_inits(c, t, R, bb, NULL);
    mpz_set_ui(zz, z);
    mpz_powm_ui(c, zz, Q, m);
    mpz_set_ui(b, n);
    mpz_powm_ui(t, b, Q, m);
    mpz_powm_ui(R, b, (Q+1)/2, m);

    unsigned long MM = S;
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            unsigned long res = mpz_get_ui(R);
            mpz_clears(b, m, r, zz, tmp, c, t, R, bb, NULL);
            return res;
        }
        unsigned long i = 0;
        mpz_set(tmp, t);
        do { mpz_mul(tmp,tmp,tmp); mpz_mod(tmp,tmp,m); i++; } while (mpz_cmp_ui(tmp,1)!=0);

        mpz_set(bb, c);
        for (unsigned long j = 0; j < MM - i - 1; j++) { mpz_mul(bb,bb,bb); mpz_mod(bb,bb,m); }
        MM = i;
        mpz_mul(c, bb, bb); mpz_mod(c, c, m);
        mpz_mul(t, t, c); mpz_mod(t, t, m);
        mpz_mul(R, R, bb); mpz_mod(R, R, m);
    }
}

/* --- Main factoring --- */
int qs_factor(mpz_t N, mpz_t result) {
    int nbits = mpz_sizeinbase(N, 2);
    int ndig = mpz_sizeinbase(N, 10);
    double lnN = nbits * log(2.0);
    double lnlnN = log(lnN);
    double L = exp(sqrt(lnN * lnlnN));

    /* Parameters — tuned for reliability */
    unsigned long B = (unsigned long)(pow(L, 0.45));
    if (B < 2000) B = 2000;
    if (B > 200000) B = 200000;
    unsigned long B2 = B * 200;
    /* Use chunked sieving: chunk_size fits in memory, sieve multiple chunks */
    int chunk_size = 50000000;
    int n_chunks = 4; /* sieve up to chunk_size * n_chunks total */

    /* Build factor base */
    fb_t *fb = calloc(MAX_FB_SIZE, sizeof(fb_t));
    int fbsz = 0;

    /* -1 pseudo-prime for sign */
    fb[fbsz].p = 2;
    mpz_t nmod;
    mpz_init(nmod);
    mpz_mod_ui(nmod, N, 2);
    fb[fbsz].r1 = fb[fbsz].r2 = mpz_get_ui(nmod) ? 1 : 0;
    fb[fbsz].logp = log(2);
    fbsz++;

    mpz_t pp;
    mpz_init_set_ui(pp, 3);
    while (mpz_cmp_ui(pp, B) <= 0 && fbsz < MAX_FB_SIZE) {
        unsigned long p = mpz_get_ui(pp);
        if (mpz_kronecker_ui(N, p) == 1) {
            mpz_mod_ui(nmod, N, p);
            unsigned long nm = mpz_get_ui(nmod);
            unsigned long r = modsqrt(nm, p);
            if (r > 0) {
                fb[fbsz].p = p;
                fb[fbsz].r1 = r;
                fb[fbsz].r2 = p - r;
                fb[fbsz].logp = log(p);
                fbsz++;
            }
        }
        mpz_nextprime(pp, pp);
    }
    mpz_clears(pp, nmod, NULL);

    int target = fbsz + 20;
    fprintf(stderr, "QS: %d dig, B=%lu, M=%d, FB=%d, target=%d\n",
            ndig, B, M, fbsz, target);

    /* Compute m = ceil(sqrt(N)) */
    mpz_t m, sqN;
    mpz_inits(m, sqN, NULL);
    mpz_sqrt(sqN, N);
    mpz_add_ui(m, sqN, 1);

    /* Sieve */
    float *sv = calloc(M, sizeof(float));
    mpz_t mmod;
    mpz_init(mmod);

    /* Initialize: subtract log of expected Q(x) value */
    double log2sqN = mpz_sizeinbase(sqN, 2) * log(2.0) + log(2.0);
    for (int x = 0; x < M; x++)
        sv[x] = -(float)(log2sqN + log(x + 1.0));

    /* Sieve with factor base */
    for (int i = 0; i < fbsz; i++) {
        unsigned long p = fb[i].p;
        float lp = (float)fb[i].logp;
        mpz_mod_ui(mmod, m, p);
        unsigned long mv = mpz_get_ui(mmod);

        unsigned long s1 = (fb[i].r1 + p - mv) % p;
        unsigned long s2 = (fb[i].r2 + p - mv) % p;

        for (unsigned long x = s1; x < (unsigned long)M; x += p) sv[x] += lp;
        if (s1 != s2)
            for (unsigned long x = s2; x < (unsigned long)M; x += p) sv[x] += lp;

        /* Sieve with p^2 */
        unsigned long p2 = p * p;
        if (p2 < (unsigned long)M && p < 100000) {
            /* Hensel lift: if r is root mod p, lift to mod p^2 */
            /* r' = r - f(r) * (2r)^{-1} mod p^2 */
            /* f(x) = (m+x)^2 - N, f(r) ≡ 0 mod p, f'(r) = 2(m+r) */
            for (int ri = 0; ri < 2; ri++) {
                unsigned long r = ri == 0 ? s1 : s2;
                /* Check: is (m + r)^2 - N ≡ 0 mod p^2? */
                mpz_t test;
                mpz_init(test);
                mpz_add_ui(test, m, r);
                mpz_mul(test, test, test);
                mpz_sub(test, test, N);
                if (mpz_divisible_ui_p(test, p2)) {
                    for (unsigned long x = r; x < (unsigned long)M; x += p2) sv[x] += lp;
                } else {
                    /* Find lifted root */
                    for (unsigned long t = 0; t < p; t++) {
                        unsigned long rr = r + t * p;
                        mpz_add_ui(test, m, rr);
                        mpz_mul(test, test, test);
                        mpz_sub(test, test, N);
                        if (mpz_divisible_ui_p(test, p2)) {
                            for (unsigned long x = rr; x < (unsigned long)M; x += p2) sv[x] += lp;
                            break;
                        }
                    }
                }
                mpz_clear(test);
            }
        }
    }
    mpz_clear(mmod);

    /* Collect candidates */
    rel_t *rels = calloc(MAX_RELS, sizeof(rel_t));
    int nrels = 0, nsmooth = 0, nlp_rels = 0;

    /* LP hash table */
    #define LPH_SIZE (1 << 18)
    int *lp_first = malloc(LPH_SIZE * sizeof(int));
    memset(lp_first, -1, LPH_SIZE * sizeof(int));
    int *lp_next = calloc(MAX_RELS, sizeof(int));
    int nlp_pairs = 0;

    /* Combined LP-pair relations */
    typedef struct {
        int r1, r2; /* indices of the two LP relations */
    } lp_pair_t;
    lp_pair_t *lp_pairs = calloc(MAX_RELS, sizeof(lp_pair_t));

    mpz_t xv, qx, rem;
    mpz_inits(xv, qx, rem, NULL);

    float thresh = -0.35; /* Keep candidates where sieve accounts for >65% */

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int x = 0; x < M && nrels < MAX_RELS - 1; x++) {
        if (sv[x] < thresh * (float)(log2sqN + log(x + 1.0))) continue;

        mpz_add_ui(xv, m, x);
        mpz_mul(qx, xv, xv);
        mpz_sub(qx, qx, N);

        int sign = mpz_sgn(qx) < 0 ? 1 : 0;
        mpz_abs(rem, qx);
        if (mpz_cmp_ui(rem, 0) == 0) {
            mpz_gcd(result, xv, N);
            if (mpz_cmp_ui(result, 1) > 0 && mpz_cmp(result, N) < 0) goto success;
            continue;
        }

        int nf = 0;
        int fp[MAX_FACTORS_PER_REL], fe[MAX_FACTORS_PER_REL];

        for (int i = 0; i < fbsz && mpz_cmp_ui(rem, 1) > 0; i++) {
            unsigned long p = fb[i].p;
            if (mpz_divisible_ui_p(rem, p)) {
                int e = 0;
                while (mpz_divisible_ui_p(rem, p)) { mpz_divexact_ui(rem, rem, p); e++; }
                if (nf < MAX_FACTORS_PER_REL) { fp[nf] = i; fe[nf] = e; nf++; }
            }
        }

        if (mpz_cmp_ui(rem, 1) == 0) {
            /* Fully smooth */
            rels[nrels].idx = nrels;
            mpz_init_set(rels[nrels].xm, xv);
            rels[nrels].sign = sign;
            rels[nrels].nf = nf;
            memcpy(rels[nrels].primes, fp, nf * sizeof(int));
            memcpy(rels[nrels].exps, fe, nf * sizeof(int));
            rels[nrels].nlp = 0;
            nrels++; nsmooth++;
        } else if (mpz_sizeinbase(rem, 2) <= 40 && mpz_probab_prime_p(rem, 10)) {
            unsigned long lp = mpz_get_ui(rem);
            if (lp <= B2) {
                rels[nrels].idx = nrels;
                mpz_init_set(rels[nrels].xm, xv);
                rels[nrels].sign = sign;
                rels[nrels].nf = nf;
                memcpy(rels[nrels].primes, fp, nf * sizeof(int));
                memcpy(rels[nrels].exps, fe, nf * sizeof(int));
                rels[nrels].nlp = 1;
                rels[nrels].lp[0] = lp;

                /* LP matching */
                int h = (int)((lp * 2654435761UL) & (LPH_SIZE - 1));
                int prev = lp_first[h];
                while (prev >= 0) {
                    if (rels[prev].lp[0] == lp) {
                        /* Found match */
                        lp_pairs[nlp_pairs].r1 = prev;
                        lp_pairs[nlp_pairs].r2 = nrels;
                        nlp_pairs++;
                        break;
                    }
                    prev = lp_next[prev];
                }
                lp_next[nrels] = lp_first[h];
                lp_first[h] = nrels;

                nrels++; nlp_rels++;
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    int total_usable = nsmooth + nlp_pairs;
    fprintf(stderr, "QS: %d smooth + %d LP-pairs = %d usable (need %d) [%.1fs]\n",
            nsmooth, nlp_pairs, total_usable, target, sieve_time);

    if (total_usable < target) {
        fprintf(stderr, "QS: Not enough relations\n");
        mpz_set_ui(result, 0);
        goto cleanup;
    }

    /* Build GF(2) matrix */
    /* Rows: smooth relations + LP-pair combined relations */
    /* For LP pair (r1, r2) with shared LP p: combined relation is */
    /* (xm1 * xm2)^2 ≡ product_of_all_factors * p^2 (mod N) */
    /* Since p^2 is a perfect square, the GF(2) row is XOR of the two rows */
    {
        int nrows = nsmooth + nlp_pairs;
        if (nrows > total_usable) nrows = total_usable;
        int ncols = fbsz + 1; /* +1 for sign */
        int words = (ncols + 63) / 64;
        int awords = (nrows + 63) / 64;
        int tw = words + awords;

        unsigned long **mat = malloc(nrows * sizeof(unsigned long *));
        mpz_t *row_xm = malloc(nrows * sizeof(mpz_t)); /* product of xm for each row */
        int **row_fb_exps = malloc(nrows * sizeof(int *)); /* total exps per FB prime */
        /* For LP pairs: track the large prime so Y includes it */
        unsigned long *row_lp = calloc(nrows, sizeof(unsigned long)); /* 0 = smooth, nonzero = LP */

        int ri = 0;

        /* Add smooth relations */
        for (int i = 0; i < nrels && ri < nrows; i++) {
            if (rels[i].nlp != 0) continue;
            mat[ri] = calloc(tw, sizeof(unsigned long));
            mat[ri][words + ri / 64] |= (1UL << (ri % 64));
            mpz_init_set(row_xm[ri], rels[i].xm);
            row_fb_exps[ri] = calloc(fbsz, sizeof(int));

            if (rels[i].sign) mat[ri][0] |= 1UL;
            for (int f = 0; f < rels[i].nf; f++) {
                int col = rels[i].primes[f] + 1;
                row_fb_exps[ri][rels[i].primes[f]] += rels[i].exps[f];
                if (rels[i].exps[f] % 2)
                    mat[ri][col / 64] ^= (1UL << (col % 64));
            }
            ri++;
        }

        /* Add LP-pair combined relations */
        for (int p = 0; p < nlp_pairs && ri < nrows; p++) {
            int i1 = lp_pairs[p].r1, i2 = lp_pairs[p].r2;
            mat[ri] = calloc(tw, sizeof(unsigned long));
            mat[ri][words + ri / 64] |= (1UL << (ri % 64));

            /* xm = xm1 * xm2 mod N */
            mpz_init(row_xm[ri]);
            mpz_mul(row_xm[ri], rels[i1].xm, rels[i2].xm);
            mpz_mod(row_xm[ri], row_xm[ri], N);

            row_fb_exps[ri] = calloc(fbsz, sizeof(int));

            /* Combined sign */
            int csign = rels[i1].sign ^ rels[i2].sign;
            if (csign) mat[ri][0] |= 1UL;

            /* Combined factors (XOR of exponent parities) */
            for (int f = 0; f < rels[i1].nf; f++) {
                int col = rels[i1].primes[f] + 1;
                row_fb_exps[ri][rels[i1].primes[f]] += rels[i1].exps[f];
                if (rels[i1].exps[f] % 2)
                    mat[ri][col / 64] ^= (1UL << (col % 64));
            }
            for (int f = 0; f < rels[i2].nf; f++) {
                int col = rels[i2].primes[f] + 1;
                row_fb_exps[ri][rels[i2].primes[f]] += rels[i2].exps[f];
                if (rels[i2].exps[f] % 2)
                    mat[ri][col / 64] ^= (1UL << (col % 64));
            }
            /* Large prime appears with even total exponent (1+1=2), so no GF(2) contribution */
            /* But we need to track it for Y computation */
            row_lp[ri] = rels[i1].lp[0]; /* = rels[i2].lp[0] since they matched */
            ri++;
        }
        nrows = ri;

        fprintf(stderr, "QS: Matrix %d x %d\n", nrows, ncols);

        /* Gaussian elimination */
        int *pivots = malloc(ncols * sizeof(int));
        memset(pivots, -1, ncols * sizeof(int));

        for (int col = 0; col < ncols; col++) {
            int piv = -1;
            for (int row = 0; row < nrows; row++) {
                if (!((mat[row][col / 64] >> (col % 64)) & 1)) continue;
                int used = 0;
                for (int c = 0; c < col; c++)
                    if (pivots[c] == row) { used = 1; break; }
                if (!used) { piv = row; break; }
            }
            if (piv < 0) continue;
            pivots[col] = piv;
            for (int row = 0; row < nrows; row++) {
                if (row == piv) continue;
                if ((mat[row][col / 64] >> (col % 64)) & 1)
                    for (int w = 0; w < tw; w++) mat[row][w] ^= mat[piv][w];
            }
        }

        /* Try dependencies */
        mpz_t X, Y, g, d;
        mpz_inits(X, Y, g, d, NULL);
        int factored = 0;

        for (int row = 0; row < nrows && !factored; row++) {
            int zero = 1;
            for (int w = 0; w < words && zero; w++)
                if (mat[row][w]) zero = 0;
            if (!zero) continue;

            /* Compute X = product of xm, Y = sqrt(product of values) */
            mpz_set_ui(X, 1);
            int *texp = calloc(fbsz, sizeof(int));
            int dcnt = 0;

            for (int i = 0; i < nrows; i++) {
                if (!((mat[row][words + i / 64] >> (i % 64)) & 1)) continue;
                mpz_mul(X, X, row_xm[i]);
                mpz_mod(X, X, N);
                for (int j = 0; j < fbsz; j++)
                    texp[j] += row_fb_exps[i][j];
                dcnt++;
            }

            if (dcnt < 2) { free(texp); continue; }

            mpz_set_ui(Y, 1);
            for (int j = 0; j < fbsz; j++) {
                if (texp[j] <= 0) continue;
                mpz_t base, pw;
                mpz_inits(base, pw, NULL);
                mpz_set_ui(base, fb[j].p);
                mpz_powm_ui(pw, base, texp[j] / 2, N);
                mpz_mul(Y, Y, pw);
                mpz_mod(Y, Y, N);
                mpz_clears(base, pw, NULL);
            }
            /* Include large primes from LP pairs in the dependency */
            for (int i = 0; i < nrows; i++) {
                if (!((mat[row][words + i / 64] >> (i % 64)) & 1)) continue;
                if (row_lp[i] != 0) {
                    mpz_t lp_mpz;
                    mpz_init_set_ui(lp_mpz, row_lp[i]);
                    mpz_mul(Y, Y, lp_mpz);
                    mpz_mod(Y, Y, N);
                    mpz_clear(lp_mpz);
                }
            }

            mpz_sub(d, X, Y); mpz_gcd(g, d, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(result, g); factored = 1;
            } else {
                mpz_add(d, X, Y); mpz_gcd(g, d, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_set(result, g); factored = 1;
                }
            }
            free(texp);
        }

        mpz_clears(X, Y, g, d, NULL);

        if (factored) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stderr, "QS: Factor found in %.3fs\n",
                    (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9);
        } else {
            fprintf(stderr, "QS: All deps trivial\n");
            mpz_set_ui(result, 0);
        }

        for (int i = 0; i < nrows; i++) {
            free(mat[i]); mpz_clear(row_xm[i]); free(row_fb_exps[i]);
        }
        free(mat); free(row_xm); free(row_fb_exps); free(pivots); free(row_lp);
    }

    goto cleanup;

success:
    clock_gettime(CLOCK_MONOTONIC, &t1);
    fprintf(stderr, "QS: Direct factor in %.3fs\n",
            (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9);

cleanup:
    for (int i = 0; i < nrels; i++) mpz_clear(rels[i].xm);
    free(rels); free(fb); free(sv);
    free(lp_first); free(lp_next); free(lp_pairs);
    mpz_clears(m, sqN, xv, qx, rem, NULL);
    return mpz_cmp_ui(result, 0) > 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N, f;
    mpz_inits(N, f, NULL);
    mpz_set_str(N, argv[1], 10);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    int ok = qs_factor(N, f);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double el = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    if (ok) {
        mpz_t c; mpz_init(c);
        mpz_divexact(c, N, f);
        if (mpz_cmp(f, c) > 0) mpz_swap(f, c);
        gmp_printf("%Zd %Zd\n", f, c);
        fprintf(stderr, "Time: %.3fs\n", el);
        mpz_clear(c);
    } else {
        fprintf(stderr, "FAIL %.3fs\n", el);
    }
    mpz_clears(N, f, NULL);
    return ok ? 0 : 1;
}
