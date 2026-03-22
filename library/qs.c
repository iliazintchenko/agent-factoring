/*
 * Basic Quadratic Sieve (single polynomial)
 * f(x) = (x + m)^2 - N where m = ceil(sqrt(N))
 * Relation: (x + m)^2 ≡ f(x) (mod N)
 *
 * Usage: ./qs <N> [deadline_seconds]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

static struct timespec t0;
static double deadline = 290.0;
static double now_sec(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (t.tv_sec - t0.tv_sec) + (t.tv_nsec - t0.tv_nsec) / 1e9;
}

/* ---- Parameters ---- */
typedef struct {
    int fb_size;
    int sieve_radius; /* sieve x from -M to M */
    int lp_mult;
} qs_params_t;

static qs_params_t get_params(int ndigits) {
    if (ndigits <= 35) return (qs_params_t){150, 50000, 30};
    if (ndigits <= 40) return (qs_params_t){300, 100000, 40};
    if (ndigits <= 45) return (qs_params_t){500, 200000, 50};
    if (ndigits <= 50) return (qs_params_t){900, 300000, 60};
    if (ndigits <= 55) return (qs_params_t){1500, 500000, 60};
    if (ndigits <= 60) return (qs_params_t){2500, 700000, 70};
    if (ndigits <= 65) return (qs_params_t){3500, 1000000, 70};
    if (ndigits <= 70) return (qs_params_t){5000, 1500000, 80};
    if (ndigits <= 75) return (qs_params_t){7000, 2000000, 90};
    if (ndigits <= 80) return (qs_params_t){9500, 3000000, 100};
    if (ndigits <= 85) return (qs_params_t){13000, 4000000, 100};
    if (ndigits <= 90) return (qs_params_t){18000, 5000000, 120};
    if (ndigits <= 95) return (qs_params_t){25000, 7000000, 120};
    return (qs_params_t){35000, 10000000, 150};
}

/* ---- Tonelli-Shanks ---- */
static unsigned long mod_sqrt(mpz_t n_big, unsigned long p) {
    unsigned long n = mpz_fdiv_ui(n_big, p);
    if (p == 2) return n & 1;
    if (n == 0) return 0;

    /* p ≡ 3 mod 4: easy case */
    if (p % 4 == 3) {
        unsigned long exp = (p + 1) / 4;
        /* compute n^exp mod p */
        unsigned long long result = 1, base = n;
        unsigned long e = exp;
        while (e > 0) {
            if (e & 1) result = (result * base) % p;
            base = (base * base) % p;
            e >>= 1;
        }
        return (unsigned long)result;
    }

    /* General Tonelli-Shanks */
    unsigned long q = p - 1, s = 0;
    while (q % 2 == 0) { q /= 2; s++; }

    /* Find quadratic non-residue */
    unsigned long z = 2;
    while (1) {
        unsigned long long tmp = 1, b2 = z;
        unsigned long e2 = (p-1)/2;
        while (e2 > 0) {
            if (e2 & 1) tmp = (tmp * b2) % p;
            b2 = (b2 * b2) % p;
            e2 >>= 1;
        }
        if (tmp == p - 1) break;
        z++;
    }

    unsigned long long M = s;
    unsigned long long c = 1, b2 = z;
    unsigned long eq = q;
    while (eq > 0) { if (eq & 1) c = (c * b2) % p; b2 = (b2 * b2) % p; eq >>= 1; }

    unsigned long long t = 1; b2 = n;
    eq = q;
    while (eq > 0) { if (eq & 1) t = (t * b2) % p; b2 = (b2 * b2) % p; eq >>= 1; }

    unsigned long long R = 1; b2 = n;
    eq = (q + 1) / 2;
    while (eq > 0) { if (eq & 1) R = (R * b2) % p; b2 = (b2 * b2) % p; eq >>= 1; }

    while (t != 1) {
        unsigned long long t2 = t;
        unsigned long i;
        for (i = 1; i < M; i++) { t2 = (t2 * t2) % p; if (t2 == 1) break; }
        unsigned long long bb = c;
        for (unsigned long j = 0; j < M - i - 1; j++) bb = (bb * bb) % p;
        R = (R * bb) % p;
        c = (bb * bb) % p;
        t = (t * c) % p;
        M = i;
    }
    return (unsigned long)R;
}

/* ---- Multiplier selection ---- */
static int select_multiplier(mpz_t n) {
    static int mults[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,30,31,33,37,41,43};
    int nm = sizeof(mults)/sizeof(mults[0]);
    double best = -1e30;
    int best_k = 1;
    mpz_t kn;
    mpz_init(kn);
    for (int mi = 0; mi < nm; mi++) {
        int k = mults[mi];
        mpz_mul_ui(kn, n, k);
        double score = -0.5 * log(k);
        for (unsigned long p = 2; p < 200; p++) {
            /* Check if p is prime */
            int is_p = 1;
            for (unsigned long d = 2; d * d <= p; d++) if (p % d == 0) { is_p = 0; break; }
            if (!is_p) continue;
            if (p == 2) {
                int r = mpz_fdiv_ui(kn, 8);
                if (r == 1) score += 2*log(2);
                else if (r == 5) score += 1.5*log(2);
                else if (r % 2 == 1) score += log(2);
            } else {
                int kr = mpz_kronecker_ui(kn, p);
                if (kr == 1) score += 2.0 * log(p) / (p - 1);
                else if (kr == 0) score += log(p) / p;
            }
        }
        if (score > best) { best = score; best_k = k; }
    }
    mpz_clear(kn);
    return best_k;
}

/* ---- Factor Base ---- */
typedef struct {
    int *p;           /* primes */
    unsigned long *r;  /* sqrt(kN) mod p */
    unsigned char *logp; /* scaled log(p) */
    int size;
} fb_t;

static fb_t* build_fb(mpz_t kn, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    fb->p = malloc((target + 10) * sizeof(int));
    fb->r = malloc((target + 10) * sizeof(unsigned long));
    fb->logp = malloc((target + 10) * sizeof(unsigned char));
    fb->size = 0;

    fb->p[0] = 2;
    fb->r[0] = 1;
    fb->logp[0] = 1; /* log2(2) = 1 */
    fb->size = 1;

    int limit = target * 30;
    if (limit < 10000) limit = 10000;

    /* Simple sieve for primes */
    char *sieve = calloc(limit + 1, 1);
    memset(sieve, 1, limit + 1);
    for (int i = 2; (long)i*i <= limit; i++)
        if (sieve[i]) for (int j = i*i; j <= limit; j += i) sieve[j] = 0;

    for (int i = 3; i <= limit && fb->size < target; i++) {
        if (!sieve[i]) continue;
        if (mpz_kronecker_ui(kn, i) != 1) continue;
        fb->p[fb->size] = i;
        fb->r[fb->size] = mod_sqrt(kn, i);
        fb->logp[fb->size] = (unsigned char)(log(i) / log(2) * 1.5 + 0.5);
        if (fb->logp[fb->size] < 1) fb->logp[fb->size] = 1;
        fb->size++;
    }
    free(sieve);
    return fb;
}

/* ---- Relation ---- */
typedef struct {
    int *exp;       /* exponent vector (fb_size) */
    mpz_t x;        /* x + m */
    int sign;       /* 1 if f(x) < 0 */
    unsigned long lp; /* large prime (0 if full) */
} rel_t;

/* ---- Partial relation hash table ---- */
#define PHASH_SIZE (1 << 18)
typedef struct pnode {
    unsigned long lp;
    int *exp;
    mpz_t x;
    int sign;
    struct pnode *next;
} pnode_t;

static pnode_t **phash;
static int n_partials = 0;

static void phash_init(void) {
    phash = calloc(PHASH_SIZE, sizeof(pnode_t*));
}

/* ---- GF(2) matrix & Gaussian elimination ---- */
static int gf2_solve(int nrels, int ncols, unsigned long long **mat_rows,
                     int ***out_deps, int **out_sizes) {
    int words = (ncols + 63) / 64;
    int aug_words = (nrels + 63) / 64;

    /* Augmentation: identity */
    unsigned long long **aug = malloc(nrels * sizeof(unsigned long long*));
    for (int i = 0; i < nrels; i++) {
        aug[i] = calloc(aug_words, sizeof(unsigned long long));
        aug[i][i / 64] |= (1ULL << (i % 64));
    }

    int *pivot_for_col = malloc(ncols * sizeof(int));
    memset(pivot_for_col, -1, ncols * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        /* Find pivot row */
        int pr = -1;
        for (int r = 0; r < nrels; r++) {
            if ((mat_rows[r][col/64] >> (col%64)) & 1) {
                /* Check not already a pivot */
                int used = 0;
                for (int c2 = 0; c2 < col; c2++)
                    if (pivot_for_col[c2] == r) { used = 1; break; }
                if (!used) { pr = r; break; }
            }
        }
        if (pr == -1) continue;
        pivot_for_col[col] = pr;

        for (int r = 0; r < nrels; r++) {
            if (r != pr && ((mat_rows[r][col/64] >> (col%64)) & 1)) {
                for (int w = 0; w < words; w++) mat_rows[r][w] ^= mat_rows[pr][w];
                for (int w = 0; w < aug_words; w++) aug[r][w] ^= aug[pr][w];
            }
        }
    }

    /* Extract null space */
    int ndeps = 0;
    for (int r = 0; r < nrels; r++) {
        int zero = 1;
        for (int w = 0; w < words && zero; w++) if (mat_rows[r][w]) zero = 0;
        if (zero) ndeps++;
    }

    *out_deps = malloc(ndeps * sizeof(int*));
    *out_sizes = malloc(ndeps * sizeof(int));
    int di = 0;
    for (int r = 0; r < nrels; r++) {
        int zero = 1;
        for (int w = 0; w < words && zero; w++) if (mat_rows[r][w]) zero = 0;
        if (!zero) continue;

        int cnt = 0;
        for (int i = 0; i < nrels; i++)
            if ((aug[r][i/64] >> (i%64)) & 1) cnt++;
        (*out_deps)[di] = malloc(cnt * sizeof(int));
        (*out_sizes)[di] = cnt;
        int idx = 0;
        for (int i = 0; i < nrels; i++)
            if ((aug[r][i/64] >> (i%64)) & 1)
                (*out_deps)[di][idx++] = i;
        di++;
    }

    for (int i = 0; i < nrels; i++) free(aug[i]);
    free(aug);
    free(pivot_for_col);
    return ndeps;
}

/* ---- Main QS ---- */
static int qs_factor(mpz_t n, mpz_t factor) {
    int ndigits = mpz_sizeinbase(n, 10);
    qs_params_t par = get_params(ndigits);

    fprintf(stderr, "QS: %d digits, fb=%d, M=%d\n", ndigits, par.fb_size, par.sieve_radius);

    /* Multiplier */
    int k = select_multiplier(n);
    mpz_t kn;
    mpz_init(kn);
    mpz_mul_ui(kn, n, k);
    fprintf(stderr, "QS: multiplier k=%d\n", k);

    /* Factor base */
    fb_t *fb = build_fb(kn, par.fb_size);
    fprintf(stderr, "QS: fb built, %d primes, largest=%d\n", fb->size, fb->p[fb->size-1]);
    if (fb->size < par.fb_size) {
        fprintf(stderr, "QS: WARNING: only got %d/%d fb primes\n", fb->size, par.fb_size);
    }

    /* m = ceil(sqrt(kN)) */
    mpz_t m, m2;
    mpz_inits(m, m2, NULL);
    mpz_sqrt(m, kn);
    mpz_mul(m2, m, m);
    if (mpz_cmp(m2, kn) < 0) mpz_add_ui(m, m, 1);

    /* Large prime bound */
    unsigned long lp_bound = (unsigned long)fb->p[fb->size-1] * par.lp_mult;

    /* Sieve threshold: approximate log2(f(M)) */
    /* f(x) ≈ 2*m*x for small x, f(M) ≈ 2*m*M */
    mpz_t fM;
    mpz_init(fM);
    mpz_mul_ui(fM, m, 2);
    mpz_mul_ui(fM, fM, par.sieve_radius);
    double log2_fM = mpz_sizeinbase(fM, 2);
    double log2_lp = log(lp_bound) / log(2);
    /* Threshold: sieve only counts each prime once (not higher powers),
     * so we need a lower threshold. Use about 60% of optimal. */
    double log2_sqrtN = mpz_sizeinbase(m, 2);
    int threshold = (int)((log2_sqrtN - log2_lp) * 1.5 * 0.65);
    if (threshold < 15) threshold = 15;
    fprintf(stderr, "QS: lp_bound=%lu, log2(fM)=%.0f, threshold=%d\n",
            lp_bound, log2_fM, threshold);
    mpz_clear(fM);

    /* Relations */
    int need_rels = fb->size + 50;
    rel_t *rels = calloc(need_rels * 3, sizeof(rel_t));
    int nrels = 0;

    /* Partial table */
    phash_init();

    /* Sieve array (byte sieve) */
    int sieve_size = 2 * par.sieve_radius + 1;
    unsigned char *sieve_arr = malloc(sieve_size);

    /* Compute sieve start positions */
    /* f(x) = (x+m)^2 - kN ≡ (x+m)^2 (mod p) */
    /* For each p, (x+m)^2 ≡ 0 mod p when x ≡ r - m or x ≡ -r - m (mod p) */
    /* where r = sqrt(kN) mod p */

    /* Sieve */
    memset(sieve_arr, 0, sieve_size);

    for (int fi = 0; fi < fb->size; fi++) {
        int p = fb->p[fi];
        unsigned char lp_val = fb->logp[fi];

        if (p == 2) {
            /* Handle p=2 specially */
            unsigned long m_mod2 = mpz_fdiv_ui(m, 2);
            long start = (m_mod2 == 0) ? 0 : 1;
            start += par.sieve_radius; /* offset for sieve array */
            start %= 2;
            for (long i = start; i < sieve_size; i += 2)
                sieve_arr[i] += lp_val;
            continue;
        }

        unsigned long r = fb->r[fi];
        unsigned long m_mod = mpz_fdiv_ui(m, p);

        /* soln1: x ≡ r - m (mod p) */
        long s1 = ((long)r - (long)m_mod) % p;
        if (s1 < 0) s1 += p;
        /* soln2: x ≡ -r - m ≡ p - r - m (mod p) */
        long s2 = ((long)(p - r) - (long)m_mod) % p;
        if (s2 < 0) s2 += p;

        /* Offset for sieve array: index i corresponds to x = i - M */
        long start1 = (s1 + par.sieve_radius) % p;
        long start2 = (s2 + par.sieve_radius) % p;

        for (long i = start1; i < sieve_size; i += p)
            sieve_arr[i] += lp_val;
        if (s1 != s2) {
            for (long i = start2; i < sieve_size; i += p)
                sieve_arr[i] += lp_val;
        }

        /* Also sieve with p^2, p^3 etc for small primes */
        if (p < 50) {
            long pp = (long)p * p;
            while (pp < sieve_size) {
                /* For p^k, need to find start position */
                /* Skip for now - the basic sieve is enough */
                break;
            }
        }
    }

    /* Debug: print sieve statistics */
    int max_sieve = 0, sum_over_thresh = 0;
    for (long i = 0; i < sieve_size; i++) {
        if (sieve_arr[i] > max_sieve) max_sieve = sieve_arr[i];
        if (sieve_arr[i] >= threshold) sum_over_thresh++;
    }
    fprintf(stderr, "QS: max_sieve=%d, #over_thresh=%d (thresh=%d)\n",
            max_sieve, sum_over_thresh, threshold);

    /* Scan for smooth candidates */
    mpz_t fx, xm, tmp, residue;
    mpz_inits(fx, xm, tmp, residue, NULL);

    for (long i = 0; i < sieve_size && nrels < need_rels; i++) {
        if (sieve_arr[i] < threshold) continue;

        long x = i - par.sieve_radius;

        /* f(x) = (x+m)^2 - kN */
        mpz_set(xm, m);
        if (x >= 0) mpz_add_ui(xm, xm, x);
        else mpz_sub_ui(xm, xm, -x);

        mpz_mul(fx, xm, xm);
        mpz_sub(fx, fx, kn);

        /* Trial divide |f(x)| by factor base */
        int *exp = calloc(fb->size, sizeof(int));
        int sign = (mpz_sgn(fx) < 0) ? 1 : 0;
        mpz_abs(residue, fx);

        if (mpz_sgn(residue) == 0) { free(exp); continue; }

        for (int fi = 0; fi < fb->size; fi++) {
            while (mpz_divisible_ui_p(residue, fb->p[fi])) {
                exp[fi]++;
                mpz_divexact_ui(residue, residue, fb->p[fi]);
            }
        }

        if (mpz_cmp_ui(residue, 1) == 0) {
            /* Full relation */
            rels[nrels].exp = exp;
            mpz_init_set(rels[nrels].x, xm);
            rels[nrels].sign = sign;
            rels[nrels].lp = 0;
            nrels++;
        } else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) <= lp_bound) {
            /* Partial: single large prime */
            unsigned long lp = mpz_get_ui(residue);
            unsigned long hash = lp % PHASH_SIZE;

            /* Look for match */
            pnode_t *pe = phash[hash];
            int matched = 0;
            while (pe) {
                if (pe->lp == lp) {
                    /* Combine */
                    int *cexp = calloc(fb->size, sizeof(int));
                    for (int fi = 0; fi < fb->size; fi++)
                        cexp[fi] = exp[fi] + pe->exp[fi];
                    rels[nrels].exp = cexp;
                    mpz_init(rels[nrels].x);
                    mpz_mul(rels[nrels].x, xm, pe->x);
                    rels[nrels].sign = sign ^ pe->sign;
                    rels[nrels].lp = lp;
                    nrels++;
                    matched = 1;
                    break;
                }
                pe = pe->next;
            }
            if (!matched) {
                pnode_t *nn = malloc(sizeof(pnode_t));
                nn->lp = lp;
                nn->exp = exp;
                exp = NULL;
                mpz_init_set(nn->x, xm);
                nn->sign = sign;
                nn->next = phash[hash];
                phash[hash] = nn;
                n_partials++;
            } else {
                free(exp);
                exp = NULL;
            }
        }

        if (exp) free(exp);
    }

    fprintf(stderr, "QS: %d relations (%d partials stored), %.1fs\n",
            nrels, n_partials, now_sec());

    if (nrels < fb->size + 1) {
        fprintf(stderr, "QS: not enough relations (need %d, got %d)\n", fb->size + 1, nrels);
        goto fail;
    }

    /* Build GF(2) matrix */
    int ncols = fb->size + 1; /* +1 for sign */
    int words = (ncols + 63) / 64;
    unsigned long long **mat = malloc(nrels * sizeof(unsigned long long*));
    for (int i = 0; i < nrels; i++) {
        mat[i] = calloc(words, sizeof(unsigned long long));
        if (rels[i].sign) mat[i][0] |= 1ULL;
        for (int fi = 0; fi < fb->size; fi++)
            if (rels[i].exp[fi] & 1)
                mat[i][(fi+1)/64] |= (1ULL << ((fi+1)%64));
    }

    fprintf(stderr, "QS: starting linear algebra (%d rels, %d cols)\n", nrels, ncols);

    int **deps;
    int *dep_sizes;
    int ndeps = gf2_solve(nrels, ncols, mat, &deps, &dep_sizes);
    fprintf(stderr, "QS: %d dependencies found\n", ndeps);

    /* Try dependencies */
    int found = 0;
    mpz_t X, Y, g;
    mpz_inits(X, Y, g, NULL);

    for (int di = 0; di < ndeps && !found; di++) {
        /* X = product of x values mod n */
        mpz_set_ui(X, 1);
        int *total_exp = calloc(fb->size, sizeof(int));

        for (int j = 0; j < dep_sizes[di]; j++) {
            int ri = deps[di][j];
            mpz_mul(X, X, rels[ri].x);
            mpz_mod(X, X, n);
            for (int fi = 0; fi < fb->size; fi++)
                total_exp[fi] += rels[ri].exp[fi];
        }

        /* Y = product of p^(exp/2) mod n * product of large_primes */
        mpz_set_ui(Y, 1);
        for (int fi = 0; fi < fb->size; fi++) {
            if (total_exp[fi] > 0) {
                mpz_t pe;
                mpz_init(pe);
                mpz_set_ui(pe, fb->p[fi]);
                mpz_powm_ui(pe, pe, total_exp[fi] / 2, n);
                mpz_mul(Y, Y, pe);
                mpz_mod(Y, Y, n);
                mpz_clear(pe);
            }
        }
        /* Include large primes from combined partials */
        for (int j = 0; j < dep_sizes[di]; j++) {
            int ri = deps[di][j];
            if (rels[ri].lp > 0) {
                mpz_mul_ui(Y, Y, rels[ri].lp);
                mpz_mod(Y, Y, n);
            }
        }

        /* gcd(X ± Y, n) */
        mpz_sub(g, X, Y);
        mpz_gcd(g, g, n);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
            mpz_set(factor, g);
            found = 1;
        } else {
            mpz_add(g, X, Y);
            mpz_gcd(g, g, n);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                mpz_set(factor, g);
                found = 1;
            }
        }
        free(total_exp);
    }

    /* Cleanup */
    mpz_clears(X, Y, g, NULL);
    for (int i = 0; i < ndeps; i++) free(deps[i]);
    free(deps); free(dep_sizes);
    for (int i = 0; i < nrels; i++) free(mat[i]);
    free(mat);

    if (found) {
        fprintf(stderr, "QS: found factor in %.1fs\n", now_sec());
        /* Cleanup */
        for (int i = 0; i < nrels; i++) { free(rels[i].exp); mpz_clear(rels[i].x); }
        free(rels);
        free(fb->p); free(fb->r); free(fb->logp); free(fb);
        mpz_clears(kn, m, m2, fx, xm, tmp, residue, NULL);
        free(sieve_arr);
        return 1;
    }

fail:
    /* Simple cleanup - just free what we can safely */
    free(sieve_arr);
    mpz_clears(kn, m, m2, fx, xm, tmp, residue, NULL);
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N> [deadline]\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &t0);
    if (argc >= 3) deadline = atof(argv[2]);

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);
    mpz_set_str(n, argv[1], 10);

    if (qs_factor(n, factor)) {
        gmp_printf("%Zd\n", factor);
        mpz_clears(n, factor, NULL);
        return 0;
    }

    fprintf(stderr, "FAIL\n");
    mpz_clears(n, factor, NULL);
    return 1;
}
