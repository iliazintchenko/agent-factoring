/*
 * sss_factor.c - Smooth Subsum Search for Integer Factoring
 *
 * Based on Hittmeir's "Smooth Subsum Search" (2023).
 * CRT-based candidate generation with collision detection and
 * batch smoothness testing via product/remainder trees.
 *
 * Compile: gcc -O3 -march=native -o sss_factor library/sss_factor.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ========================== PARAMETERS ========================== */

#define MAX_FB       200000
#define MAX_RELS     210000
#define PLIST_DIV    5
#define MAX_SUBSET   8
#define LP_MULT      128

/* ========================== GLOBALS ========================== */

static int fb[MAX_FB], fb_size;
static int plist_idx[MAX_FB]; /* indices into fb[] for CRT primes */
static int plist_size;

static unsigned int roots_arr[MAX_FB][2]; /* sqrt(N) mod fb[i], adjusted for b */

static mpz_t N_val, b_val;
static mpz_t coeffs[MAX_FB]; /* CRT coefficients */
static mpz_t nval;            /* product of plist primes */
static mpz_t mval;            /* powersmooth product for batch test */

typedef struct { mpz_t x, qx; } rel_t;
static rel_t rels[MAX_RELS];
static int nrels;

typedef struct { unsigned long lp; mpz_t x; } partial_t;
static partial_t parts[MAX_RELS * 2];
static int nparts;

/* RNG */
static unsigned long rng_s = 42;
static unsigned long rng(void) {
    rng_s = rng_s * 6364136223846793005ULL + 1442695040888963407ULL;
    return rng_s;
}

/* ========================== PRIMES ========================== */

static int *all_primes;
static int nprimes_total;

static void gen_primes(int limit) {
    char *sieve = calloc(limit + 1, 1);
    all_primes = malloc((limit / 4) * sizeof(int));
    nprimes_total = 0;
    for (int i = 2; i <= limit; i++) {
        if (!sieve[i]) {
            all_primes[nprimes_total++] = i;
            for (long j = (long)i * i; j <= limit; j += i)
                sieve[j] = 1;
        }
    }
    free(sieve);
}

/* ========================== TONELLI-SHANKS ========================== */

static unsigned int tsqrt(unsigned long n_mod, unsigned int p) {
    if (n_mod == 0) return 0;
    if (p == 2) return n_mod & 1;

    /* Legendre symbol check */
    mpz_t nm, pz;
    mpz_init_set_ui(nm, n_mod);
    mpz_init_set_ui(pz, p);
    if (mpz_legendre(nm, pz) != 1) { mpz_clear(nm); mpz_clear(pz); return 0xFFFFFFFF; }

    unsigned int r;
    if (p % 4 == 3) {
        mpz_t e;
        mpz_init_set_ui(e, (p + 1) / 4);
        mpz_powm(nm, nm, e, pz);
        r = mpz_get_ui(nm);
        mpz_clear(e);
    } else {
        /* Full Tonelli-Shanks */
        unsigned long Q = p - 1;
        int S = __builtin_ctzl(Q);
        Q >>= S;

        unsigned long z = 2;
        mpz_t zt;
        mpz_init(zt);
        while (1) {
            mpz_set_ui(zt, z);
            if (mpz_legendre(zt, pz) == -1) break;
            z++;
        }

        mpz_t M_v, c, t, R, e;
        mpz_init_set_ui(M_v, S);
        mpz_init(c); mpz_init(t); mpz_init(R); mpz_init(e);

        mpz_set_ui(e, Q);
        mpz_set_ui(c, z); mpz_powm(c, c, e, pz);
        mpz_set_ui(t, n_mod); mpz_powm(t, t, e, pz);
        mpz_set_ui(e, (Q + 1) / 2);
        mpz_set_ui(R, n_mod); mpz_powm(R, R, e, pz);

        while (1) {
            if (mpz_cmp_ui(t, 1) == 0) break;
            mpz_t tt; mpz_init_set(tt, t);
            int i = 0;
            while (mpz_cmp_ui(tt, 1) != 0) {
                mpz_mul(tt, tt, tt); mpz_mod(tt, tt, pz); i++;
            }
            mpz_clear(tt);

            mpz_set(e, c);
            for (long j = 0; j < (long)mpz_get_ui(M_v) - i - 1; j++) {
                mpz_mul(e, e, e); mpz_mod(e, e, pz);
            }
            mpz_set_ui(M_v, i);
            mpz_mul(c, e, e); mpz_mod(c, c, pz);
            mpz_mul(t, t, c); mpz_mod(t, t, pz);
            mpz_mul(R, R, e); mpz_mod(R, R, pz);
        }
        r = mpz_get_ui(R);
        mpz_clear(M_v); mpz_clear(c); mpz_clear(t); mpz_clear(R); mpz_clear(e);
        mpz_clear(zt);
    }
    mpz_clear(nm); mpz_clear(pz);
    return r;
}

/* ========================== PRODUCT/REMAINDER TREE ========================== */

typedef struct { mpz_t *d; int n; } level_t;

static void prod_tree(mpz_t *xs, int n, mpz_t *result) {
    /* Compute product tree, store root in result */
    if (n == 0) { mpz_set_ui(*result, 1); return; }
    if (n == 1) { mpz_set(*result, xs[0]); return; }

    /* Build tree bottom-up */
    mpz_t *cur = malloc(n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) mpz_init_set(cur[i], xs[i]);
    int cn = n;

    while (cn > 1) {
        int nn = (cn + 1) / 2;
        mpz_t *next = malloc(nn * sizeof(mpz_t));
        for (int i = 0; i < nn; i++) {
            mpz_init(next[i]);
            if (2*i+1 < cn)
                mpz_mul(next[i], cur[2*i], cur[2*i+1]);
            else
                mpz_set(next[i], cur[2*i]);
        }
        for (int i = 0; i < cn; i++) mpz_clear(cur[i]);
        free(cur);
        cur = next;
        cn = nn;
    }
    mpz_set(*result, cur[0]);
    mpz_clear(cur[0]);
    free(cur);
}

/* Compute remainders: result[i] = n mod xs[i] for each i */
static void rem_tree(mpz_t n, mpz_t *xs, int cnt, mpz_t *result) {
    if (cnt == 0) return;
    if (cnt == 1) { mpz_mod(result[0], n, xs[0]); return; }

    /* Build product tree, keeping all levels */
    int max_levels = 0;
    { int t = cnt; while (t > 1) { t = (t+1)/2; max_levels++; } max_levels++; }

    level_t *tree = malloc(max_levels * sizeof(level_t));
    tree[0].n = cnt;
    tree[0].d = malloc(cnt * sizeof(mpz_t));
    for (int i = 0; i < cnt; i++) mpz_init_set(tree[0].d[i], xs[i]);

    int lv = 0;
    while (tree[lv].n > 1) {
        int cn = tree[lv].n;
        int nn = (cn + 1) / 2;
        tree[lv+1].n = nn;
        tree[lv+1].d = malloc(nn * sizeof(mpz_t));
        for (int i = 0; i < nn; i++) {
            mpz_init(tree[lv+1].d[i]);
            if (2*i+1 < cn)
                mpz_mul(tree[lv+1].d[i], tree[lv].d[2*i], tree[lv].d[2*i+1]);
            else
                mpz_set(tree[lv+1].d[i], tree[lv].d[2*i]);
        }
        lv++;
    }
    int nlv = lv + 1;

    /* Top-down remainder */
    mpz_t *rem = malloc(1 * sizeof(mpz_t));
    mpz_init(rem[0]);
    mpz_mod(rem[0], n, tree[nlv-1].d[0]);

    for (int l = nlv - 1; l > 0; l--) {
        int child_n = tree[l-1].n;
        mpz_t *nrem = malloc(child_n * sizeof(mpz_t));
        for (int i = 0; i < child_n; i++) {
            mpz_init(nrem[i]);
            mpz_mod(nrem[i], rem[i/2], tree[l-1].d[i]);
        }
        for (int i = 0; i < tree[l].n; i++) mpz_clear(rem[i]);
        free(rem);
        rem = nrem;
    }

    for (int i = 0; i < cnt; i++) mpz_set(result[i], rem[i]);

    for (int i = 0; i < cnt; i++) mpz_clear(rem[i]);
    free(rem);
    for (int l = 0; l < nlv; l++) {
        for (int i = 0; i < tree[l].n; i++) mpz_clear(tree[l].d[i]);
        free(tree[l].d);
    }
    free(tree);
}

/* ========================== SETUP ========================== */

static int target_fb(int digits) {
    if (digits <= 18) return 120;
    if (digits <= 25) return 300;
    if (digits <= 34) return 400;
    if (digits <= 36) return 600;
    if (digits <= 38) return 800;
    if (digits <= 40) return 1000;
    if (digits <= 42) return 1200;
    if (digits <= 44) return 1400;
    if (digits <= 48) return 2000;
    if (digits <= 52) return 2400;
    if (digits <= 56) return 4000;
    if (digits <= 60) return 8000;
    if (digits <= 66) return 12000;
    if (digits <= 74) return 20000;
    if (digits <= 80) return 60000;
    if (digits <= 88) return 100000;
    if (digits <= 94) return 120000;
    return 200000;
}

static void setup(void) {
    int dig = mpz_sizeinbase(N_val, 10);
    int m = target_fb(dig);
    int plimit = m * 20;
    if (plimit < 3000000) plimit = 3000000;
    gen_primes(plimit);

    fb_size = 0;
    for (int i = 0; i < nprimes_total && fb_size < m; i++) {
        int p = all_primes[i];
        if (mpz_divisible_ui_p(N_val, p)) continue;
        if (p > 2) {
            int leg = mpz_kronecker_ui(N_val, p);
            if (leg != 1) continue;
        }

        unsigned long nmod = mpz_fdiv_ui(N_val, p);
        unsigned int r = tsqrt(nmod, p);
        if (r == 0xFFFFFFFF) continue;

        unsigned int r2 = p - r;
        unsigned long bmod = mpz_fdiv_ui(b_val, p);
        unsigned int adj0 = ((unsigned long)r + p - bmod) % p;
        unsigned int adj1 = ((unsigned long)r2 + p - bmod) % p;
        if (adj0 > adj1) { unsigned int t = adj0; adj0 = adj1; adj1 = t; }

        fb[fb_size] = p;
        roots_arr[fb_size][0] = adj0;
        roots_arr[fb_size][1] = adj1;
        fb_size++;
    }

    plist_size = fb_size / PLIST_DIV;
    if (plist_size < 10) plist_size = 10;
    if (plist_size > fb_size) plist_size = fb_size;

    /* Compute nval = product of plist primes */
    mpz_init_set_ui(nval, 1);
    for (int i = 0; i < plist_size; i++) {
        plist_idx[i] = i;
        mpz_mul_ui(nval, nval, fb[i]);
    }

    /* CRT coefficients */
    for (int i = 0; i < plist_size; i++) {
        mpz_init(coeffs[i]);
        mpz_t ni, pz;
        mpz_init(ni); mpz_init_set_ui(pz, fb[i]);
        mpz_divexact_ui(ni, nval, fb[i]);
        mpz_t inv; mpz_init(inv);
        mpz_invert(inv, ni, pz);
        mpz_mul(coeffs[i], inv, ni);
        mpz_clear(ni); mpz_clear(inv); mpz_clear(pz);
    }

    /* mval = powersmooth product of all FB primes */
    mpz_init_set_ui(mval, 1);
    for (int i = 0; i < fb_size; i++) {
        int p = fb[i];
        int k = (int)(15.0 * log(2.0) / log((double)p));
        if (k < 1) k = 1;
        for (int j = 0; j < k; j++)
            mpz_mul_ui(mval, mval, p);
    }

    fprintf(stderr, "FB: %d primes (largest %d), plist: %d\n",
            fb_size, fb[fb_size-1], plist_size);
}

/* ========================== Q(x) ========================== */

static void qx(mpz_t result, mpz_t x) {
    /* Q(x) = (x + b)^2 - N */
    mpz_add(result, x, b_val);
    mpz_mul(result, result, result);
    mpz_sub(result, result, N_val);
}

/* ========================== CORE SEARCH ========================== */

static void add_relation(mpz_t x_val) {
    if (nrels >= MAX_RELS) return;
    /* Store (x+b) mod N as the "x" so that x^2 ≡ Q(x) (mod N) */
    mpz_init(rels[nrels].x);
    mpz_add(rels[nrels].x, x_val, b_val);
    mpz_mod(rels[nrels].x, rels[nrels].x, N_val);
    mpz_init(rels[nrels].qx);
    qx(rels[nrels].qx, x_val);
    nrels++;
}

static void add_partial(unsigned long lp, mpz_t x_val) {
    if (nparts >= MAX_RELS * 2) return;
    mpz_init_set(parts[nparts].x, x_val);
    parts[nparts].lp = lp;
    nparts++;
}

/*
 * The core SSS search: for a random CRT subset, find x values where Q(x) is
 * divisible by the CRT modulus AND multiple additional large primes (via collision
 * detection), then batch-test the cofactor for smoothness.
 */
static void search_one(int subset_len) {
    /* Pick random subset of plist indices (skip index 0 = prime 2) */
    int indset[MAX_SUBSET];
    int iss = 0;
    char used[MAX_FB];
    memset(used, 0, plist_size);

    while (iss < subset_len && iss < plist_size - 1) {
        int idx = 1 + (rng() % (plist_size - 1)); /* skip 0 */
        if (!used[idx]) { used[idx] = 1; indset[iss++] = idx; }
    }
    int nss = iss;

    /* M = product of subset primes */
    /* primes_sel[0] = 1 (dummy), primes_sel[1..nss] = selected primes */
    unsigned long primes_sel[MAX_SUBSET + 1];
    int ind_sel[MAX_SUBSET + 1];
    primes_sel[0] = 1; ind_sel[0] = -1;
    for (int i = 0; i < nss; i++) {
        primes_sel[i+1] = fb[indset[i]];
        ind_sel[i+1] = indset[i];
    }
    int np = nss + 1;

    mpz_t M;
    mpz_init_set_ui(M, 1);
    for (int i = 0; i < nss; i++)
        mpz_mul_ui(M, M, fb[indset[i]]);

    /* mvals[j] = M / primes_sel[j] */
    mpz_t mv[MAX_SUBSET + 1];
    for (int j = 0; j < np; j++) {
        mpz_init(mv[j]);
        if (primes_sel[j] <= 1) mpz_set(mv[j], M);
        else mpz_divexact_ui(mv[j], M, primes_sel[j]);
    }

    /* xval = CRT solution for roots[][0] */
    mpz_t xval;
    mpz_init_set_ui(xval, 0);
    for (int i = 0; i < nss; i++) {
        mpz_t t; mpz_init(t);
        mpz_mul_ui(t, coeffs[indset[i]], roots_arr[indset[i]][0]);
        mpz_add(xval, xval, t);
        mpz_clear(t);
    }
    mpz_mod(xval, xval, M);

    /* inv[v] = M^(-1) mod fb[v] for large primes */
    unsigned int *inv = calloc(fb_size, sizeof(unsigned int));
    for (int v = plist_size; v < fb_size; v++) {
        unsigned long Mm = mpz_fdiv_ui(M, fb[v]);
        if (Mm == 0) continue;
        mpz_t t, pz;
        mpz_init_set_ui(t, Mm); mpz_init_set_ui(pz, fb[v]);
        if (mpz_invert(t, t, pz)) inv[v] = mpz_get_ui(t);
        mpz_clear(t); mpz_clear(pz);
    }

    /* difflist[i] = coeffs[indset[i]] * (roots[indset[i]][1] - roots[indset[i]][0]) */
    mpz_t diffs[MAX_SUBSET];
    for (int i = 0; i < nss; i++) {
        mpz_init(diffs[i]);
        int idx = indset[i];
        long d = (long)roots_arr[idx][1] - (long)roots_arr[idx][0];
        mpz_mul_si(diffs[i], coeffs[idx], d);
    }

    unsigned long last_p = fb[fb_size - 1];
    unsigned long partial_bound = LP_MULT * last_p;

    /* Buffers for candidate collection */
    int max_cands = 8192;
    mpz_t *cand_args = malloc(max_cands * sizeof(mpz_t));
    mpz_t *cand_vals = malloc(max_cands * sizeof(mpz_t));
    for (int i = 0; i < max_cands; i++) { mpz_init(cand_args[i]); mpz_init(cand_vals[i]); }

    int nlarge = fb_size - plist_size;
    long *kbuf = malloc(nlarge * 4 * sizeof(long));
    unsigned int (*bvs)[2] = calloc(fb_size, sizeof(unsigned int[2]));

    /* Iterate: swap each root from 0→1 */
    for (int ii = 0; ii < nss; ii++) {
        int idx = indset[ii];

        mpz_add(xval, xval, diffs[ii]);
        mpz_mod(xval, xval, M);

        /* Compute base offsets for each large prime */
        unsigned long xvm;
        for (int v = plist_size; v < fb_size; v++) {
            unsigned int p = fb[v];
            xvm = mpz_fdiv_ui(xval, p);
            unsigned long b0 = ((unsigned long)roots_arr[v][0] + p - xvm) % p;
            b0 = (b0 * (unsigned long)inv[v]) % p;
            unsigned long b1 = ((unsigned long)roots_arr[v][1] + p - xvm) % p;
            b1 = (b1 * (unsigned long)inv[v]) % p;
            bvs[v][0] = b0;
            bvs[v][1] = b1;
        }

        /* For each mval (M/prime), find collision offsets */
        for (int j = 0; j < np; j++) {
            if (ind_sel[j] == idx) continue;

            unsigned long pri = primes_sel[j];
            int nk = 0;

            if (j == 0) {
                /* Use base bvs directly */
                for (int v = plist_size; v < fb_size; v++) {
                    unsigned int p = fb[v];
                    long b0 = bvs[v][0], b1 = bvs[v][1];
                    kbuf[nk++] = b0 - p;
                    kbuf[nk++] = b0;
                    kbuf[nk++] = b1 - p;
                    kbuf[nk++] = b1;
                }
            } else {
                /* Multiply bvs by pri mod p */
                for (int v = plist_size; v < fb_size; v++) {
                    unsigned int p = fb[v];
                    long r0 = ((unsigned long)bvs[v][0] * pri) % p;
                    long r1 = ((unsigned long)bvs[v][1] * pri) % p;
                    kbuf[nk++] = r0 - p;
                    kbuf[nk++] = r0;
                    kbuf[nk++] = r1 - p;
                    kbuf[nk++] = r1;
                }
            }

            /* Sort kbuf to find collisions */
            /* Use counting sort or qsort */
            /* For simplicity, use a simple approach: sort then scan */
            /* Actually, for speed, let's use a hash-like approach */
            /* But for correctness first, just sort */

            /* Simple insertion sort for moderate-size arrays would be slow.
             * Use stdlib qsort */
            /* Need a comparator for long */
            /* Actually let's just do a radix-like approach */

            /* Sort */
            {
                /* qsort with long comparator */
                int cmp_long(const void *a, const void *b) {
                    long la = *(const long *)a, lb = *(const long *)b;
                    return (la > lb) - (la < lb);
                }
                qsort(kbuf, nk, sizeof(long), cmp_long);
            }

            /* Find runs of length >= 3 */
            int ncands = 0;
            int i = 0;
            while (i < nk) {
                int run = 1;
                while (i + run < nk && kbuf[i + run] == kbuf[i]) run++;
                if (run >= 3 && ncands < max_cands) {
                    long k = kbuf[i];

                    /* arg = xval + k * mv[j] */
                    mpz_mul_si(cand_args[ncands], mv[j], k);
                    mpz_add(cand_args[ncands], cand_args[ncands], xval);

                    /* val = |Q(arg)| / mv[j] (divide by the submodulus) */
                    mpz_t full_q;
                    mpz_init(full_q);
                    qx(full_q, cand_args[ncands]);
                    mpz_abs(full_q, full_q);

                    /* Divide by mv[j] - should be exact */
                    if (mpz_divisible_p(full_q, mv[j])) {
                        mpz_divexact(cand_vals[ncands], full_q, mv[j]);
                        ncands++;
                    }
                    mpz_clear(full_q);
                }
                i += run;
            }

            /* Batch smoothness test on collected candidates */
            if (ncands > 0) {
                mpz_t *rems = malloc(ncands * sizeof(mpz_t));
                for (int c = 0; c < ncands; c++) mpz_init(rems[c]);

                rem_tree(mval, cand_vals, ncands, rems);

                for (int c = 0; c < ncands; c++) {
                    if (mpz_sgn(rems[c]) == 0) {
                        /* Fully smooth */
                        add_relation(cand_args[c]);
                    } else {
                        mpz_t g, cof;
                        mpz_init(g); mpz_init(cof);
                        mpz_gcd(g, cand_vals[c], rems[c]);
                        mpz_divexact(cof, cand_vals[c], g);

                        if (mpz_cmp_ui(cof, last_p) <= 0) {
                            add_relation(cand_args[c]);
                        } else if (mpz_fits_ulong_p(cof) &&
                                   mpz_get_ui(cof) < partial_bound) {
                            add_partial(mpz_get_ui(cof), cand_args[c]);
                        }
                        mpz_clear(g); mpz_clear(cof);
                    }
                }

                for (int c = 0; c < ncands; c++) mpz_clear(rems[c]);
                free(rems);
            }
        }
    }

    /* Cleanup */
    for (int i = 0; i < max_cands; i++) { mpz_clear(cand_args[i]); mpz_clear(cand_vals[i]); }
    free(cand_args); free(cand_vals);
    free(kbuf); free(bvs);
    for (int i = 0; i < nss; i++) mpz_clear(diffs[i]);
    for (int j = 0; j < np; j++) mpz_clear(mv[j]);
    mpz_clear(M); mpz_clear(xval);
    free(inv);
}

/* ========================== LP MATCHING ========================== */

static int lp_cmp(const void *a, const void *b) {
    unsigned long la = ((const partial_t *)a)->lp;
    unsigned long lb = ((const partial_t *)b)->lp;
    return (la > lb) - (la < lb);
}

static int match_partials(void) {
    if (nparts < 2) return 0;

    /* Sort by large prime */
    qsort(parts, nparts, sizeof(partial_t), lp_cmp);

    int matched = 0;
    int i = 0;
    while (i < nparts - 1) {
        if (parts[i].lp == parts[i+1].lp) {
            /* Combine: x1 and x2 share the same large prime LP */
            /* (x1+b)(x2+b) ≡ ... mod N, Q(x1)*Q(x2)/LP^2 is smooth */
            mpz_t combined_x, combined_q, t1, t2;
            mpz_init(combined_x); mpz_init(combined_q);
            mpz_init(t1); mpz_init(t2);

            /* combined x = (x1+b)*(x2+b) * LP^(-1) mod N */
            /* This is stored directly as x (already includes +b) */
            mpz_add(t1, parts[i].x, b_val);
            mpz_add(t2, parts[i+1].x, b_val);
            mpz_mul(combined_x, t1, t2);
            mpz_t lp_z; mpz_init_set_ui(lp_z, parts[i].lp);
            mpz_t inv_lp; mpz_init(inv_lp);
            if (mpz_invert(inv_lp, lp_z, N_val)) {
                mpz_mul(combined_x, combined_x, inv_lp);
                mpz_mod(combined_x, combined_x, N_val);

                /* combined Q = Q(x1) * Q(x2) / LP^2 */
                mpz_t q1, q2;
                mpz_init(q1); mpz_init(q2);
                qx(q1, parts[i].x);
                qx(q2, parts[i+1].x);
                mpz_mul(combined_q, q1, q2);
                mpz_t lp2; mpz_init(lp2);
                mpz_mul(lp2, lp_z, lp_z);
                mpz_divexact(combined_q, combined_q, lp2);

                if (nrels < MAX_RELS) {
                    /* Store combined_x directly (already has +b baked in) */
                    mpz_init_set(rels[nrels].x, combined_x);
                    mpz_init_set(rels[nrels].qx, combined_q);
                    nrels++;
                    matched++;
                }
                mpz_clear(q1); mpz_clear(q2); mpz_clear(lp2);
            }
            mpz_clear(combined_x); mpz_clear(combined_q);
            mpz_clear(t1); mpz_clear(t2); mpz_clear(lp_z); mpz_clear(inv_lp);
            i += 2;
        } else {
            i++;
        }
    }
    return matched;
}

/* ========================== LINEAR ALGEBRA ========================== */

static int solve(void) {
    fprintf(stderr, "Linear algebra: %d relations, %d FB primes\n", nrels, fb_size);

    int nf = fb_size + 1; /* +1 for sign */

    /* Factorize each relation over FB */
    int **exps = malloc(nrels * sizeof(int *));
    int valid = 0;
    for (int r = 0; r < nrels; r++) {
        exps[r] = calloc(nf, sizeof(int));
        mpz_t tmp; mpz_init_set(tmp, rels[r].qx);
        if (mpz_sgn(tmp) < 0) { exps[r][0] = 1; mpz_neg(tmp, tmp); }
        int ok = 1;
        for (int i = 0; i < fb_size && ok; i++) {
            while (mpz_divisible_ui_p(tmp, fb[i])) {
                mpz_divexact_ui(tmp, tmp, fb[i]);
                exps[r][i + 1]++;
            }
        }
        if (mpz_cmp_ui(tmp, 1) != 0) {
            /* Not smooth - mark invalid */
            exps[r][0] = -99;
        } else {
            valid++;
        }
        mpz_clear(tmp);
    }

    fprintf(stderr, "Valid smooth relations: %d / %d\n", valid, nrels);

    /* Verify all relations: x^2 ≡ qx (mod N) */
    int verified = 0;
    for (int r = 0; r < nrels; r++) {
        if (exps[r][0] == -99) continue;
        mpz_t x2;
        mpz_init(x2);
        mpz_mul(x2, rels[r].x, rels[r].x);
        mpz_mod(x2, x2, N_val);
        mpz_t qmod;
        mpz_init(qmod);
        mpz_mod(qmod, rels[r].qx, N_val);
        if (mpz_cmp(x2, qmod) == 0) {
            verified++;
        } else {
            /* Invalidate this relation */
            exps[r][0] = -99;
            valid--;
        }
        mpz_clear(x2);
        mpz_clear(qmod);
    }
    fprintf(stderr, "Verified x^2 ≡ Q mod N: %d / %d\n", verified, valid + (nrels - verified));

    if (valid < nf / 2) {
        for (int r = 0; r < nrels; r++) free(exps[r]);
        free(exps);
        return 0;
    }

    /* Gaussian elimination mod 2 - bit-packed */
    int nwords = (nrels + 63) / 64;
    unsigned long *cols = calloc((size_t)nf * nwords, sizeof(unsigned long));

    for (int r = 0; r < nrels; r++) {
        if (exps[r][0] == -99) continue;
        for (int c = 0; c < nf; c++) {
            if (exps[r][c] & 1)
                cols[c * nwords + r / 64] ^= (1UL << (r % 64));
        }
    }

    int *pivots = malloc(nf * sizeof(int));
    char *marked = calloc(nrels, 1);

    for (int j = 0; j < nf; j++) {
        pivots[j] = -1;
        for (int w = 0; w < nwords; w++) {
            if (cols[j * nwords + w]) {
                int bit = __builtin_ctzl(cols[j * nwords + w]);
                int i = w * 64 + bit;
                if (i >= nrels) continue;
                pivots[j] = i;
                marked[i] = 1;
                for (int k = 0; k < nf; k++) {
                    if (k != j && (cols[k * nwords + w] & (1UL << bit))) {
                        for (int ww = 0; ww < nwords; ww++)
                            cols[k * nwords + ww] ^= cols[j * nwords + ww];
                    }
                }
                break;
            }
        }
    }

    /* Extract factors from dependencies */
    int found = 0;
    int attempts = 0;
    for (int i = 0; i < nrels && !found; i++) {
        if (marked[i]) continue;
        if (exps[i][0] == -99) continue;

        int dep[MAX_RELS];
        int nd = 0;
        dep[nd++] = i;
        for (int j = 0; j < nf; j++) {
            if (cols[j * nwords + i / 64] & (1UL << (i % 64))) {
                if (pivots[j] >= 0) dep[nd++] = pivots[j];
            }
        }

        attempts++;

        mpz_t s1, s2;
        mpz_init_set_ui(s1, 1);
        mpz_init_set_ui(s2, 1);

        for (int d = 0; d < nd; d++) {
            int idx = dep[d];
            /* s1 *= x mod N */
            mpz_mul(s1, s1, rels[idx].x);
            mpz_mod(s1, s1, N_val);

            /* s2 *= |Q(x)| (accumulate big product, sqrt at end) */
            mpz_t aq; mpz_init(aq);
            mpz_abs(aq, rels[idx].qx);
            mpz_mul(s2, s2, aq);
            mpz_clear(aq);
        }

        /* Integer square root of accumulated product */
        if (!mpz_perfect_square_p(s2)) {
            mpz_clear(s1); mpz_clear(s2);
            continue;
        }
        mpz_sqrt(s2, s2);
        mpz_mod(s2, s2, N_val);

        mpz_t diff, factor;
        mpz_init(diff); mpz_init(factor);

        /* Try s1 - s2 */
        mpz_sub(diff, s1, s2);
        mpz_gcd(factor, diff, N_val);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N_val) < 0) {
            gmp_printf("%Zd\n", factor);
            found = 1;
        }

        /* Try s1 + s2 */
        if (!found) {
            mpz_add(diff, s1, s2);
            mpz_gcd(factor, diff, N_val);
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N_val) < 0) {
                gmp_printf("%Zd\n", factor);
                found = 1;
            }
        }

        mpz_clear(s1); mpz_clear(s2); mpz_clear(diff); mpz_clear(factor);
    }
    fprintf(stderr, "Tried %d dependencies\n", attempts);

    for (int r = 0; r < nrels; r++) free(exps[r]);
    free(exps); free(cols); free(pivots); free(marked);
    return found;
}

/* ========================== MAIN ========================== */

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_init_set_str(N_val, argv[1], 10);
    mpz_init(b_val);
    mpz_sqrt(b_val, N_val);
    mpz_add_ui(b_val, b_val, 1);
    {
        mpz_t t; mpz_init(t);
        mpz_mul(t, b_val, b_val);
        if (mpz_cmp(t, N_val) < 0) mpz_add_ui(b_val, b_val, 1);
        mpz_clear(t);
    }

    int dig = mpz_sizeinbase(N_val, 10);
    fprintf(stderr, "N = %s (%d digits)\n", argv[1], dig);

    setup();

    int target = fb_size + 10;
    int subset_len = 6;
    if (fb_size > 2000) subset_len = 7;
    if (fb_size > 8000) subset_len = 7;

    fprintf(stderr, "Searching (target %d relations, subset_len %d)...\n", target, subset_len);

    int iter = 0;
    while (nrels < target) {
        search_one(subset_len);
        iter++;

        /* Periodically match partials - disabled, combined partials cause
         * degenerate dependencies (structural correlation issue) */
        if (0 && iter % 200 == 0 && nparts > 100) {
            int m = match_partials();
            if (m > 0) {
                fprintf(stderr, "  Matched %d partial relations\n", m);
                /* Reset partials */
                for (int i = 0; i < nparts; i++) mpz_clear(parts[i].x);
                nparts = 0;
            }
        }

        if (iter % 50 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "\rIter %d: %d/%d rels, %d parts, %.1fs   ",
                    iter, nrels, target, nparts, el);
        }

        if (iter > 100000) {
            fprintf(stderr, "\nToo many iterations, giving up.\n");
            break;
        }
    }

    /* Final partial matching disabled */

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double search_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "\nSearch done: %d relations in %.2fs\n", nrels, search_time);

    int found = solve();

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Total time: %.4f seconds\n", total);

    /* Cleanup */
    mpz_clear(N_val); mpz_clear(b_val); mpz_clear(nval); mpz_clear(mval);
    for (int i = 0; i < plist_size; i++) mpz_clear(coeffs[i]);
    for (int i = 0; i < nrels; i++) { mpz_clear(rels[i].x); mpz_clear(rels[i].qx); }
    for (int i = 0; i < nparts; i++) mpz_clear(parts[i].x);
    free(all_primes);

    return found ? 0 : 1;
}
