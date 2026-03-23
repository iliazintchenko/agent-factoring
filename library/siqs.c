/*
 * Quadratic Sieve with Multiple Polynomials (MPQS)
 *
 * For polynomial f(x) = (x + base)^2 - N, sieve f(x) for smooth values.
 * Uses large prime variation with hash table.
 * Uses multiple base points (window shifting) for polynomial variety.
 *
 * Usage: ./siqs <number>
 * Single-threaded, seed=42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_FB 12000
#define MAX_RELS 30000

typedef struct {
    unsigned int p;
    int logp;
    unsigned int sqrt_n_mod_p;
} fb_t;

typedef struct {
    mpz_t x;         /* x + base, the value whose square ≡ Q(x) + N (mod N) */
    int *expo;        /* exponents of Q(x) over the factor base */
    int neg;          /* 1 if Q(x) < 0 */
    unsigned long lp; /* large prime, or 0 if fully smooth */
    unsigned long combined_lp; /* lp value for combined partials, to include in sqrt */
} rel_t;

static fb_t fb[MAX_FB];
static int fb_size;
static rel_t rels[MAX_RELS];
static int nrels = 0;
static mpz_t N_global;

/* ---- primes ---- */
static int *primes_list;
static int nprimes;

static void make_primes(int lim) {
    char *s = calloc(lim+1, 1);
    int c = 0;
    for (int i = 2; i <= lim; i++) if (!s[i]) {
        c++;
        for (long j = (long)i*i; j <= lim; j += i) s[j] = 1;
    }
    primes_list = malloc(c * sizeof(int));
    nprimes = 0;
    memset(s, 0, lim+1);
    for (int i = 2; i <= lim; i++) if (!s[i]) {
        primes_list[nprimes++] = i;
        for (long j = (long)i*i; j <= lim; j += i) s[j] = 1;
    }
    free(s);
}

/* Tonelli-Shanks: returns sqrt(a) mod p, or -1 if not a QR */
static unsigned long tsqrt(unsigned long a, unsigned long p) {
    if (p == 2) return a & 1;
    if (a == 0) return 0;
    mpz_t az, pz, rz;
    mpz_init_set_ui(az, a);
    mpz_init_set_ui(pz, p);
    mpz_init(rz);
    mpz_powm_ui(rz, az, (p-1)/2, pz);
    if (mpz_cmp_ui(rz, 1) != 0) { mpz_clear(az); mpz_clear(pz); mpz_clear(rz); return (unsigned long)-1; }
    /* p ≡ 3 mod 4 fast path */
    if ((p & 3) == 3) {
        mpz_powm_ui(rz, az, (p+1)/4, pz);
        unsigned long r = mpz_get_ui(rz);
        mpz_clear(az); mpz_clear(pz); mpz_clear(rz);
        return r;
    }
    /* General case */
    unsigned long Q = p-1, S = 0;
    while (!(Q&1)) { Q >>= 1; S++; }
    unsigned long z = 2;
    mpz_t zz;
    mpz_init(zz);
    while (1) { mpz_set_ui(zz,z); mpz_powm_ui(zz,zz,(p-1)/2,pz); if (mpz_cmp_ui(zz,p-1)==0) break; z++; }
    mpz_t M,c,t,R,b;
    mpz_init_set_ui(M,S); mpz_init(c); mpz_init(t); mpz_init(R); mpz_init(b);
    mpz_set_ui(c,z); mpz_powm_ui(c,c,Q,pz);
    mpz_set_ui(t,a); mpz_powm_ui(t,t,Q,pz);
    mpz_set_ui(R,a); mpz_powm_ui(R,R,(Q+1)/2,pz);
    while (1) {
        if (mpz_cmp_ui(t,1)==0) { unsigned long r=mpz_get_ui(R); mpz_clear(M);mpz_clear(c);mpz_clear(t);mpz_clear(R);mpz_clear(b);mpz_clear(zz);mpz_clear(az);mpz_clear(pz);mpz_clear(rz); return r; }
        mpz_set(zz,t); unsigned long i=0;
        while (mpz_cmp_ui(zz,1)!=0) { mpz_mul(zz,zz,zz); mpz_mod(zz,zz,pz); i++; }
        unsigned long m=mpz_get_ui(M);
        mpz_set(b,c);
        for (unsigned long j=0;j<m-i-1;j++) { mpz_mul(b,b,b); mpz_mod(b,b,pz); }
        mpz_set_ui(M,i);
        mpz_mul(c,b,b); mpz_mod(c,c,pz);
        mpz_mul(t,t,c); mpz_mod(t,t,pz);
        mpz_mul(R,R,b); mpz_mod(R,R,pz);
    }
}

/* ---- parameters by digit count ---- */
typedef struct { int fb_sz; int M; double thresh_sub; int lp_mult; } params_t;

static params_t get_params(int d) {
    /* M = sieve half-width, thresh_sub = how much below log2(Qmax) to set threshold */
    if (d<=30)  return (params_t){100,  32768,  25.0, 30};
    if (d<=35)  return (params_t){200,  65536,  27.0, 40};
    if (d<=40)  return (params_t){350,  65536,  27.0, 50};
    if (d<=45)  return (params_t){550,  65536,  27.0, 50};
    if (d<=50)  return (params_t){1000, 131072, 27.0, 60};
    if (d<=55)  return (params_t){1500, 131072, 27.0, 60};
    if (d<=60)  return (params_t){2200, 196608, 26.0, 60};
    if (d<=65)  return (params_t){3200, 262144, 26.0, 60};
    if (d<=70)  return (params_t){4200, 262144, 25.0, 60};
    if (d<=75)  return (params_t){5500, 327680, 25.0, 60};
    if (d<=80)  return (params_t){7000, 393216, 24.0, 60};
    if (d<=85)  return (params_t){8500, 393216, 23.0, 60};
    if (d<=90)  return (params_t){10000,524288, 23.0, 60};
    if (d<=95)  return (params_t){11000,524288, 22.0, 60};
    return (params_t){12000,524288, 22.0, 60};
}

/* ---- factor base construction ---- */
static int build_fb(int target) {
    int cnt = 0;
    for (int i = 0; i < nprimes && cnt < target; i++) {
        int p = primes_list[i];
        unsigned long nm = mpz_fdiv_ui(N_global, p);
        if (p == 2) {
            fb[cnt].p = 2; fb[cnt].logp = 1; fb[cnt].sqrt_n_mod_p = nm; cnt++;
            continue;
        }
        unsigned long sr = tsqrt(nm, p);
        if (sr == (unsigned long)-1) continue;
        fb[cnt].p = p;
        fb[cnt].logp = (int)(log2((double)p) + 0.5);
        if (fb[cnt].logp < 1) fb[cnt].logp = 1;
        fb[cnt].sqrt_n_mod_p = (unsigned int)sr;
        cnt++;
    }
    return cnt;
}

/* ---- large prime hash ---- */
#define HTSIZE (1<<20)
typedef struct ht_entry { unsigned long lp; int idx; struct ht_entry *next; } ht_entry;
static ht_entry *htable[HTSIZE];

static int ht_lookup(unsigned long lp) {
    int h = (int)((lp * 2654435761UL) >> 12) & (HTSIZE-1);
    for (ht_entry *e = htable[h]; e; e = e->next)
        if (e->lp == lp) return e->idx;
    return -1;
}
static void ht_insert(unsigned long lp, int idx) {
    int h = (int)((lp * 2654435761UL) >> 12) & (HTSIZE-1);
    ht_entry *e = malloc(sizeof(ht_entry));
    e->lp = lp; e->idx = idx; e->next = htable[h]; htable[h] = e;
}

/* ---- main QS logic ---- */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_init_set_str(N_global, argv[1], 10);
    int digits = strlen(argv[1]);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Trial division up to 1M */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N_global, p)) {
            mpz_t q; mpz_init(q); mpz_divexact_ui(q, N_global, p);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            gmp_printf("%lu %Zd\n", p, q);
            fprintf(stderr, "Trial division: time=%.3fs\n",
                (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9);
            return 0;
        }
    }

    params_t par = get_params(digits);
    make_primes(par.fb_sz * 50);
    fb_size = build_fb(par.fb_sz);
    int M = par.M;
    int sieve_len = 2 * M;

    fprintf(stderr, "QS: %d digits, FB=%d (max=%u), M=%d\n", digits, fb_size, fb[fb_size-1].p, M);

    int needed = fb_size + 50;
    unsigned long lp_bound = (unsigned long)fb[fb_size-1].p * par.lp_mult;

    /* Compute sqrt(N) */
    mpz_t sqrt_n, base, Q, xv, cofactor, tmp;
    mpz_init(sqrt_n); mpz_init(base); mpz_init(Q); mpz_init(xv);
    mpz_init(cofactor); mpz_init(tmp);
    mpz_sqrt(sqrt_n, N_global);
    mpz_add_ui(sqrt_n, sqrt_n, 1); /* ceil */

    unsigned char *sieve = malloc(sieve_len);
    int *soln1 = malloc(fb_size * sizeof(int));
    int *soln2 = malloc(fb_size * sizeof(int));
    int *expo = malloc(fb_size * sizeof(int));

    /* log2 of max |Q(x)|: Q(x) = (base+x)^2 - N ≈ 2*sqrt(N)*M when x≈M */
    double log2q = log2((double)M) + mpz_sizeinbase(N_global, 2)/2.0 + 1.0;
    double log2_lpb = log2((double)lp_bound);
    /* Threshold = log2(Q) - log2(lp_bound) - slack (for accumulated rounding errors in sieve) */
    int thresh = (int)(log2q - log2_lpb - 4.0);
    if (thresh < 15) thresh = 15;
    fprintf(stderr, "Threshold: %d (log2q=%.1f)\n", thresh, log2q);

    int smooth_found = 0, partials_combined = 0, partials_stored = 0;
    int window = 0, max_windows = 2000000;

    memset(htable, 0, sizeof(htable));

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 42);

    while (smooth_found + partials_combined < needed && window < max_windows) {
        /* base = ceil(sqrt(N)) + window * sieve_len/2 (overlapping windows) */
        /* Use alternating: sqrt(N) + k, sqrt(N) - k to stay close */
        mpz_set(base, sqrt_n);
        if (window & 1)
            mpz_sub_ui(base, base, (unsigned long)(window/2 + 1) * M);
        else
            mpz_add_ui(base, base, (unsigned long)(window/2) * M);

        /* Compute sieve offsets:
         * Q(x) = (base + x)^2 - N, x in [0, sieve_len)
         * We want base + x ≡ ± sqrt(N) (mod p)
         * x ≡ sqrt_n_mod_p - base (mod p)  or  x ≡ -sqrt_n_mod_p - base (mod p)
         */
        for (int i = 0; i < fb_size; i++) {
            unsigned int p = fb[i].p;
            unsigned long bmod = mpz_fdiv_ui(base, p);
            unsigned long sr = fb[i].sqrt_n_mod_p;
            if (p == 2) {
                /* (base+x)^2 - N mod 2: base is near sqrt(N), N is odd */
                /* (base+x) must be chosen so (base+x)^2 ≡ N mod 2 */
                /* N odd means N mod 2 = 1, so (base+x) must be odd */
                soln1[i] = (bmod & 1) ? 0 : 1;
                soln2[i] = -1;
                continue;
            }
            long s1 = ((long)sr - (long)bmod) % (long)p;
            if (s1 < 0) s1 += p;
            long s2 = ((long)(p - sr) - (long)bmod) % (long)p;
            if (s2 < 0) s2 += p;
            soln1[i] = (int)s1;
            soln2[i] = (int)s2;
        }

        /* Sieve */
        memset(sieve, 0, sieve_len);
        for (int i = 0; i < fb_size; i++) {
            unsigned int p = fb[i].p;
            int lp = fb[i].logp;
            if (soln1[i] >= 0)
                for (int j = soln1[i]; j < sieve_len; j += p) sieve[j] += lp;
            if (soln2[i] >= 0 && soln2[i] != soln1[i])
                for (int j = soln2[i]; j < sieve_len; j += p) sieve[j] += lp;
        }

        /* Scan candidates */
        for (int i = 0; i < sieve_len && nrels < MAX_RELS - 10; i++) {
            if (sieve[i] < thresh) continue;

            /* Q = (base + i)^2 - N */
            mpz_set(xv, base);
            mpz_add_ui(xv, xv, i);
            mpz_mul(Q, xv, xv);
            mpz_sub(Q, Q, N_global);

            int neg = (mpz_sgn(Q) < 0);
            if (neg) mpz_neg(Q, Q);

            /* Trial divide */
            mpz_set(cofactor, Q);
            memset(expo, 0, fb_size * sizeof(int));
            for (int j = 0; j < fb_size; j++) {
                unsigned int p = fb[j].p;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    expo[j]++;
                }
            }

            int smooth = 0;
            unsigned long lp = 0;
            if (mpz_cmp_ui(cofactor, 1) == 0) {
                smooth = 1;
            } else if (mpz_fits_ulong_p(cofactor)) {
                unsigned long r = mpz_get_ui(cofactor);
                if (r <= lp_bound && r > 1) {
                    smooth = 2;
                    lp = r;
                }
            }

            if (!smooth) continue;

            /* Store relation */
            rel_t *r = &rels[nrels];
            mpz_init_set(r->x, xv);
            r->expo = malloc(fb_size * sizeof(int));
            memcpy(r->expo, expo, fb_size * sizeof(int));
            r->neg = neg;
            r->lp = lp;
            r->combined_lp = 0;
            int this_idx = nrels;
            nrels++;

            if (smooth == 1) {
                smooth_found++;
            } else {
                /* Partial: check if we can combine */
                int match = ht_lookup(lp);
                if (match >= 0) {
                    /* Combine: create a new relation that is the "product" */
                    /* (x1 * x2)^2 ≡ Q1 * Q2 (mod N), and Q1*Q2/lp^2 is smooth */
                    rel_t *combined = &rels[nrels];
                    mpz_init(combined->x);
                    mpz_mul(combined->x, rels[match].x, rels[this_idx].x);
                    mpz_mod(combined->x, combined->x, N_global);
                    combined->expo = malloc(fb_size * sizeof(int));
                    for (int j = 0; j < fb_size; j++)
                        combined->expo[j] = rels[match].expo[j] + rels[this_idx].expo[j];
                    combined->neg = rels[match].neg ^ rels[this_idx].neg;
                    combined->lp = 0; /* fully smooth now (lp^2 divides product) */
                    combined->combined_lp = lp; /* need this for computing Y */
                    nrels++;
                    partials_combined++;
                } else {
                    ht_insert(lp, this_idx);
                    partials_stored++;
                }
            }
        }

        window++;
        if (window % 200 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
            fprintf(stderr, "Win %d: smooth=%d combined=%d stored=%d total=%d/%d (%.1fs)\n",
                    window, smooth_found, partials_combined, partials_stored,
                    smooth_found+partials_combined, needed, el);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_time = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
    fprintf(stderr, "Sieve: smooth=%d combined=%d, total=%d/%d (%.1fs, %d windows)\n",
            smooth_found, partials_combined, smooth_found+partials_combined, needed, sieve_time, window);

    /* Collect usable relations (lp == 0) */
    int *usable = malloc(nrels * sizeof(int));
    int nusable = 0;
    for (int i = 0; i < nrels; i++)
        if (rels[i].lp == 0)
            usable[nusable++] = i;

    fprintf(stderr, "Usable: %d (need %d)\n", nusable, fb_size + 1);

    if (nusable < fb_size + 1) {
        fprintf(stderr, "Not enough relations\n");
        printf("FAIL\n");
        return 1;
    }

    /* ---- Gaussian elimination mod 2 ---- */
    int nrows = nusable;
    int ncols = fb_size + 1; /* +1 for sign */
    int total_cols = ncols + nrows;

    unsigned char **mat = malloc(nrows * sizeof(unsigned char *));
    for (int i = 0; i < nrows; i++) {
        mat[i] = calloc(total_cols, 1);
        int ri = usable[i];
        mat[i][0] = rels[ri].neg & 1;
        for (int j = 0; j < fb_size; j++)
            mat[i][j+1] = rels[ri].expo[j] & 1;
        mat[i][ncols + i] = 1; /* identity for tracking */
    }

    int rank = 0;
    for (int col = 0; col < ncols && rank < nrows; col++) {
        int pr = -1;
        for (int row = rank; row < nrows; row++)
            if (mat[row][col]) { pr = row; break; }
        if (pr < 0) continue;
        if (pr != rank) { unsigned char *t = mat[pr]; mat[pr] = mat[rank]; mat[rank] = t; }
        for (int row = 0; row < nrows; row++) {
            if (row != rank && mat[row][col]) {
                for (int c = 0; c < total_cols; c++)
                    mat[row][c] ^= mat[rank][c];
            }
        }
        rank++;
    }

    fprintf(stderr, "Matrix %dx%d, rank=%d, deps=%d\n", nrows, ncols, rank, nrows-rank);

    /* Try null-space vectors */
    int found = 0;
    mpz_t X, Y, g, p_mpz;
    mpz_init(X); mpz_init(Y); mpz_init(g); mpz_init(p_mpz);

    for (int row = rank; row < nrows && !found; row++) {
        mpz_set_ui(X, 1);
        int *total = calloc(fb_size, sizeof(int));
        int neg_cnt = 0;

        for (int i = 0; i < nrows; i++) {
            if (!mat[row][ncols + i]) continue;
            int ri = usable[i];
            mpz_mul(X, X, rels[ri].x);
            mpz_mod(X, X, N_global);
            for (int j = 0; j < fb_size; j++)
                total[j] += rels[ri].expo[j];
            if (rels[ri].neg) neg_cnt++;
        }

        /* All exponents should be even (including neg_cnt) */
        mpz_set_ui(Y, 1);
        for (int j = 0; j < fb_size; j++) {
            if (total[j] <= 0) continue;
            mpz_set_ui(p_mpz, fb[j].p);
            mpz_powm_ui(p_mpz, p_mpz, total[j] / 2, N_global);
            mpz_mul(Y, Y, p_mpz);
            mpz_mod(Y, Y, N_global);
        }

        /* Include combined_lp factors: each combined partial contributes lp to Y */
        for (int i = 0; i < nrows; i++) {
            if (!mat[row][ncols + i]) continue;
            int ri = usable[i];
            if (rels[ri].combined_lp > 0) {
                mpz_set_ui(p_mpz, rels[ri].combined_lp);
                mpz_mul(Y, Y, p_mpz);
                mpz_mod(Y, Y, N_global);
            }
        }

        mpz_sub(g, X, Y); mpz_gcd(g, g, N_global);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) { found = 1; }
        else { mpz_add(g, X, Y); mpz_gcd(g, g, N_global);
               if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) found = 1; }

        free(total);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;

    if (found) {
        mpz_t cof; mpz_init(cof); mpz_divexact(cof, N_global, g);
        gmp_printf("%Zd %Zd\n", g, cof);
        fprintf(stderr, "QS found factor, time=%.3fs\n", elapsed);
        mpz_clear(cof);
        return 0;
    } else {
        fprintf(stderr, "QS failed in linear algebra, time=%.3fs\n", elapsed);
        printf("FAIL\n");
        return 1;
    }
}
