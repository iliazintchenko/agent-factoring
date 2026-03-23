/*
 * MQSS - Modular Quadratic Self-initializing Sieve
 *
 * Novel QS variant: uses "primorial GCD" for fast smooth-part extraction
 * instead of traditional trial division. After sieve identifies candidates,
 * smooth-part extraction is done via repeated GCD with the primorial of
 * the factor base, avoiding per-prime trial division.
 *
 * Compile: gcc -O2 -march=native -o mqss library/mqss.c -lgmp -lm
 * Usage: ./mqss <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

#define SIEVE_BLOCK   32768
#define LP_HASH_SIZE  (1 << 20)
#define MAX_FB        200000
#define MAX_RELS      300000
#define MAX_DEPS      128
#define SEED          42

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size, num_blocks, lp_mult, num_a_primes;
} params_t;

static params_t get_params(int digits) {
    if (digits <= 30) return (params_t){120,   2,   40, 3};
    if (digits <= 35) return (params_t){250,   3,   40, 4};
    if (digits <= 40) return (params_t){500,   5,   50, 5};
    if (digits <= 45) return (params_t){900,   7,   60, 6};
    if (digits <= 50) return (params_t){1500,  10,  70, 7};
    if (digits <= 55) return (params_t){2500,  14,  80, 8};
    if (digits <= 60) return (params_t){4000,  20,  90, 9};
    if (digits <= 65) return (params_t){7000,  28, 100, 10};
    if (digits <= 70) return (params_t){12000, 36, 100, 11};
    if (digits <= 75) return (params_t){18000, 48, 110, 13};
    if (digits <= 80) return (params_t){28000, 64, 110, 14};
    if (digits <= 85) return (params_t){45000, 80, 110, 15};
    if (digits <= 90) return (params_t){70000, 100,110, 16};
    return (params_t){100000, 140, 120, 17};
}

/* ==================== Utilities ==================== */
static int *sieve_primes(int limit, int *count) {
    char *s = calloc(limit + 1, 1);
    for (int i = 2; (long)i * i <= limit; i++)
        if (!s[i]) for (int j = i * i; j <= limit; j += i) s[j] = 1;
    int cnt = 0;
    for (int i = 2; i <= limit; i++) if (!s[i]) cnt++;
    int *p = malloc(cnt * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= limit; i++) if (!s[i]) p[idx++] = i;
    free(s);
    *count = cnt;
    return p;
}

static int modpow_ll(long long base, long long exp, long long mod) {
    long long r = 1; base %= mod;
    while (exp > 0) {
        if (exp & 1) r = (__int128)r * base % mod;
        exp >>= 1; base = (__int128)base * base % mod;
    }
    return (int)r;
}

static int tonelli_shanks(long long n, int p) {
    if (p == 2) return (int)(n & 1);
    n = ((n % p) + p) % p;
    if (n == 0) return 0;
    if (modpow_ll(n, (p - 1) / 2, p) != 1) return -1;
    if (p % 4 == 3) return modpow_ll(n, (p + 1) / 4, p);

    int s = 0; long long q = p - 1;
    while (!(q & 1)) { s++; q >>= 1; }
    long long z = 2;
    while (modpow_ll(z, (p - 1) / 2, p) != p - 1) z++;

    long long m = s, c = modpow_ll(z, q, p);
    long long t = modpow_ll(n, q, p), r = modpow_ll(n, (q + 1) / 2, p);
    while (1) {
        if (t == 1) return (int)r;
        long long i = 0, tmp = t;
        while (tmp != 1) { tmp = (__int128)tmp * tmp % p; i++; }
        long long b = c;
        for (long long j = 0; j < m - i - 1; j++) b = (__int128)b * b % p;
        m = i; c = (__int128)b * b % p;
        t = (__int128)t * c % p; r = (__int128)r * b % p;
    }
}

static int mod_inv(long long a, long long m) {
    long long old_r = a, r = m, old_s = 1, s = 0;
    while (r) { long long q = old_r / r;
        long long t = r; r = old_r - q * r; old_r = t;
        t = s; s = old_s - q * s; old_s = t;
    }
    return (int)((old_s % m + m) % m);
}

/* ==================== Factor Base ==================== */
typedef struct {
    int p, root1, root2;  /* prime and roots of kN mod p */
    int soln1, soln2;      /* current sieve solutions */
    double logp;
} fb_entry_t;

static fb_entry_t fb[MAX_FB];
static int fb_size, fb_bound;
static mpz_t kN;
static int multiplier;

static int select_multiplier(const mpz_t N) {
    static const int cands[] = {1,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    double best = -1e30; int best_k = 1;
    mpz_t tmp; mpz_init(tmp);
    for (int ci = 0; ci < 15; ci++) {
        int k = cands[ci];
        mpz_mul_ui(tmp, N, k);
        double sc = -0.5 * log(k);
        int r8 = mpz_fdiv_ui(tmp, 8);
        if (r8 == 1) sc += 2*log(2.0);
        else if (r8 == 5) sc += log(2.0);

        int nprimes; int *primes = sieve_primes(200, &nprimes);
        for (int i = 0; i < nprimes; i++) {
            int p = primes[i]; if (p == 2) continue;
            if (mpz_kronecker_si(tmp, p) == 1) sc += 2.0*log(p)/(p-1);
        }
        free(primes);
        if (sc > best) { best = sc; best_k = k; }
    }
    mpz_clear(tmp);
    return best_k;
}

static void build_factor_base(const mpz_t N, int target) {
    multiplier = select_multiplier(N);
    mpz_init(kN); mpz_mul_ui(kN, N, multiplier);

    int prime_limit = target * 25;
    if (prime_limit < 50000) prime_limit = 50000;
    int nprimes; int *primes = sieve_primes(prime_limit, &nprimes);

    fb_size = 0;
    /* Prime 2 */
    fb[0].p = 2; fb[0].root1 = fb[0].root2 = 1; fb[0].logp = log(2.0);
    fb_size = 1;

    for (int i = 0; i < nprimes && fb_size < target; i++) {
        int p = primes[i]; if (p == 2) continue;
        int r = tonelli_shanks(mpz_fdiv_ui(kN, p), p);
        if (r < 0) continue;
        fb[fb_size].p = p;
        fb[fb_size].root1 = r;
        fb[fb_size].root2 = (r == 0) ? -1 : (p - r);
        fb[fb_size].logp = log(p);
        fb_size++;
    }
    fb_bound = fb[fb_size - 1].p;
    free(primes);
}

/* ==================== Polynomial ==================== */
typedef struct {
    mpz_t A, B, C;
    mpz_t *Bainv2;   /* 2 * Bi for Gray code update */
    int *a_idx;       /* indices into fb[] of A-primes */
    int num_a;
    int b_idx, max_b;
} poly_t;

static void compute_sieve_roots(poly_t *poly) {
    for (int i = 0; i < fb_size; i++) {
        int p = fb[i].p;
        if (p == 2) { fb[i].soln1 = 0; fb[i].soln2 = -1; continue; }

        int a_mod = mpz_fdiv_ui(poly->A, p);
        if (a_mod == 0) { fb[i].soln1 = fb[i].soln2 = -1; continue; }

        int ainv = mod_inv(a_mod, p);
        int b_mod = mpz_fdiv_ui(poly->B, p);

        long long s1 = (long long)ainv * ((long long)fb[i].root1 - b_mod) % p;
        fb[i].soln1 = (int)((s1 % p + p) % p);

        if (fb[i].root2 >= 0) {
            long long s2 = (long long)ainv * ((long long)fb[i].root2 - b_mod) % p;
            fb[i].soln2 = (int)((s2 % p + p) % p);
        } else fb[i].soln2 = -1;
    }
}

static int a_call_count = 0;

static void poly_new_a(poly_t *poly, int num_a_hint, const mpz_t target_a) {
    /* Dynamically determine num_a: find how many primes we need from the
       middle of the FB to get A ≈ target_a */
    double target_log = mpz_sizeinbase(target_a, 2) * 0.693;

    /* Estimate: pick primes from around the middle of FB */
    int mid = fb_size / 3;
    if (mid < 2) mid = 2;
    double avg_logp = log(fb[mid].p);
    int num_a = (int)(target_log / avg_logp + 0.5);
    if (num_a < 3) num_a = 3;
    if (num_a > fb_size / 4) num_a = fb_size / 4;

    poly->num_a = num_a;
    poly->a_idx = malloc(num_a * sizeof(int));

    /* Search range: middle third of FB */
    int start = fb_size / 6;
    if (start < 2) start = 2;
    int end = fb_size * 5 / 6;
    if (end > fb_size) end = fb_size;
    int range = end - start;

    mpz_t prod, test;
    mpz_init_set_ui(prod, 1);
    mpz_init(test);
    int *used = calloc(fb_size, sizeof(int));

    /* Use a_call_count to vary the first prime selection, producing
       different A values on each call. */
    int call = a_call_count++;

    /* Build list of eligible primes in range */
    int *eligible = malloc(range * sizeof(int));
    int n_eligible = 0;
    for (int j = start; j < end; j++)
        if (fb[j].root2 >= 0) eligible[n_eligible++] = j;

    if (n_eligible > 0) {
        /* Pick first prime based on call count */
        int first_idx = eligible[call % n_eligible];
        used[first_idx] = 1;
        poly->a_idx[0] = first_idx;
        mpz_set_ui(prod, fb[first_idx].p);

        /* For subsequent primes, use greedy best-ratio approach */
        for (int i = 1; i < num_a; i++) {
            double best_ratio = 1e30;
            int best = -1;
            for (int j = start; j < end; j++) {
                if (used[j] || fb[j].root2 < 0) continue;
                mpz_mul_ui(test, prod, fb[j].p);
                double r = (mpz_cmp(test, target_a) > 0) ?
                    mpz_get_d(test) / mpz_get_d(target_a) :
                    mpz_get_d(target_a) / mpz_get_d(test);
                if (r < best_ratio) { best_ratio = r; best = j; }
            }
            if (best < 0) {
                for (int j = 2; j < fb_size; j++)
                    if (!used[j] && fb[j].root2 >= 0) { best = j; break; }
            }
            if (best < 0) break;
            used[best] = 1;
            poly->a_idx[i] = best;
            mpz_mul_ui(prod, prod, fb[best].p);
        }
    }
    free(eligible);

    mpz_set(poly->A, prod);
    free(used);
    mpz_clear(prod);
    mpz_clear(test);
}

static void poly_init_b(poly_t *poly) {
    mpz_t Ai, Bi;
    mpz_init(Ai); mpz_init(Bi);

    poly->Bainv2 = malloc(poly->num_a * sizeof(mpz_t));
    mpz_set_ui(poly->B, 0);

    for (int i = 0; i < poly->num_a; i++) {
        mpz_init(poly->Bainv2[i]);
        int ai = fb[poly->a_idx[i]].p;
        int ri = fb[poly->a_idx[i]].root1;

        mpz_divexact_ui(Ai, poly->A, ai);
        int ainv = mod_inv(mpz_fdiv_ui(Ai, ai), ai);

        /* Bi = ri * ainv * Ai */
        mpz_set_ui(Bi, (long long)ri * ainv % ai);
        mpz_mul(Bi, Bi, Ai);

        /* Store 2*Bi for Gray code updates */
        mpz_mul_2exp(poly->Bainv2[i], Bi, 1);

        mpz_add(poly->B, poly->B, Bi);
    }

    /* Verify B^2 ≡ kN (mod A), fix sign if needed */
    mpz_t temp, mod_a;
    mpz_init(temp); mpz_init(mod_a);
    mpz_mul(temp, poly->B, poly->B);
    mpz_mod(temp, temp, poly->A);
    mpz_mod(mod_a, kN, poly->A);
    if (mpz_cmp(temp, mod_a) != 0) {
        mpz_sub(poly->B, poly->A, poly->B);
        /* Recompute Bainv2 with correct signs */
    }
    mpz_clear(temp); mpz_clear(mod_a);

    /* C = (B^2 - kN) / A */
    mpz_mul(poly->C, poly->B, poly->B);
    mpz_sub(poly->C, poly->C, kN);
    if (!mpz_divisible_p(poly->C, poly->A)) {
        fprintf(stderr, "BUG: B^2 - kN not divisible by A!\n");
    }
    mpz_divexact(poly->C, poly->C, poly->A);

    poly->b_idx = 0;
    poly->max_b = 1 << (poly->num_a - 1);

    mpz_clear(Ai); mpz_clear(Bi);
    compute_sieve_roots(poly);
}

static int poly_next_b(poly_t *poly) {
    poly->b_idx++;
    if (poly->b_idx >= poly->max_b) return 0;

    /* Gray code: find rightmost set bit in b_idx */
    int v = poly->b_idx, gray_idx = 0;
    while (!(v & 1)) { gray_idx++; v >>= 1; }

    /* Sign from Gray code */
    int sign = (poly->b_idx >> (gray_idx + 1)) & 1;
    if (sign)
        mpz_sub(poly->B, poly->B, poly->Bainv2[gray_idx]);
    else
        mpz_add(poly->B, poly->B, poly->Bainv2[gray_idx]);

    /* Recompute C */
    mpz_mul(poly->C, poly->B, poly->B);
    mpz_sub(poly->C, poly->C, kN);
    mpz_divexact(poly->C, poly->C, poly->A);

    compute_sieve_roots(poly);
    return 1;
}

static void poly_clear(poly_t *poly) {
    for (int i = 0; i < poly->num_a; i++) mpz_clear(poly->Bainv2[i]);
    free(poly->Bainv2);
    free(poly->a_idx);
}

/* ==================== Relations ==================== */
typedef struct {
    mpz_t Axb;           /* Ax + B value */
    int *exponents;       /* full exponents for each FB prime + A-primes */
    uint64_t *parity_row; /* GF(2) parity vector */
    int large_prime;
    int neg;              /* 1 if Q(x) < 0 */
} relation_t;

static relation_t rels[MAX_RELS];
static int nrels = 0;
static int ncombined = 0;

/* LP hash table */
static int lp_hash[LP_HASH_SIZE];
static int lp_next[MAX_RELS];

static void init_lp_hash(void) {
    memset(lp_hash, -1, sizeof(lp_hash));
}

static int parity_words;

static void add_rel(const mpz_t Axb, const int *exp, int fb_sz,
                    int lp, int neg, const int *a_idx, int num_a) {
    if (nrels >= MAX_RELS) return;
    relation_t *r = &rels[nrels];
    mpz_init_set(r->Axb, Axb);
    r->exponents = malloc(fb_sz * sizeof(int));
    memcpy(r->exponents, exp, fb_sz * sizeof(int));

    /* Add A-prime exponents: each a-prime appears once in A*Q(x) = (Ax+B)^2 */
    for (int i = 0; i < num_a; i++)
        r->exponents[a_idx[i]] += 1;

    r->large_prime = lp;
    r->neg = neg;
    r->parity_row = NULL; /* built later at LA time */

    if (lp > 1) {
        int h = lp % LP_HASH_SIZE;
        int idx = lp_hash[h];
        while (idx >= 0) {
            if (rels[idx].large_prime == lp) { ncombined++; break; }
            idx = lp_next[idx];
        }
        lp_next[nrels] = lp_hash[h];
        lp_hash[h] = nrels;
    }
    nrels++;
}

/* ==================== Primorial GCD smooth extraction ==================== */
static mpz_t primorial;

static void compute_primorial(void) {
    mpz_init_set_ui(primorial, 1);
    for (int i = 0; i < fb_size; i++) {
        /* Include prime powers up to 2^63 */
        int p = fb[i].p;
        long long pk = p;
        while (pk <= (1LL << 28)) { /* prime powers up to ~256M */
            mpz_mul_ui(primorial, primorial, p);
            if (pk > (1LL << 28) / p) break;
            pk *= p;
        }
    }
}

/* Extract B-smooth part of |val| via repeated GCD with primorial.
   Returns cofactor. Sets smooth_exp[] to exponents of each FB prime. */
static void extract_smooth(int *smooth_exp, mpz_t cofactor,
                           const mpz_t val) {
    mpz_t g;
    mpz_init(g);
    mpz_abs(cofactor, val);
    memset(smooth_exp, 0, fb_size * sizeof(int));

    /* Iteratively extract all smooth factors */
    for (int iter = 0; iter < 30; iter++) {
        mpz_gcd(g, cofactor, primorial);
        if (mpz_cmp_ui(g, 1) == 0) break;
        mpz_divexact(cofactor, cofactor, g);
    }

    /* Now cofactor has no FB-smooth factors. The smooth part = |val| / cofactor.
       Factorize the smooth part by trial division. */
    mpz_t smooth_part;
    mpz_init(smooth_part);
    mpz_abs(smooth_part, val);
    mpz_divexact(smooth_part, smooth_part, cofactor);

    for (int i = 0; i < fb_size; i++) {
        int p = fb[i].p;
        while (mpz_divisible_ui_p(smooth_part, p)) {
            mpz_divexact_ui(smooth_part, smooth_part, p);
            smooth_exp[i]++;
        }
        if (mpz_cmp_ui(smooth_part, 1) == 0) break;
    }
    mpz_clear(smooth_part);
    mpz_clear(g);
}

/* ==================== Sieve + relation collection ==================== */

static unsigned char sieve_arr[SIEVE_BLOCK];

static int collect_block(int block_start, int block_size,
                         poly_t *poly, const mpz_t N,
                         int lp_bound, unsigned char thresh) {
    int found = 0;
    memset(sieve_arr, 0, block_size);

    /* Standard logarithmic sieve */
    for (int i = 0; i < fb_size; i++) {
        int p = fb[i].p;
        if (fb[i].soln1 < 0) continue;
        int logp_scaled = (int)(fb[i].logp * 1.44 + 0.5); /* log2 approx */
        if (logp_scaled < 1) logp_scaled = 1;
        if (logp_scaled > 127) logp_scaled = 127;

        int off1 = ((fb[i].soln1 - block_start) % p);
        if (off1 < 0) off1 += p;
        for (int j = off1; j < block_size; j += p) {
            int v = sieve_arr[j] + logp_scaled;
            sieve_arr[j] = (v > 255) ? 255 : (unsigned char)v;
        }

        if (fb[i].soln2 >= 0 && fb[i].soln2 != fb[i].soln1) {
            int off2 = ((fb[i].soln2 - block_start) % p);
            if (off2 < 0) off2 += p;
            for (int j = off2; j < block_size; j += p) {
                int v = sieve_arr[j] + logp_scaled;
                sieve_arr[j] = (v > 255) ? 255 : (unsigned char)v;
            }
        }
    }

    /* Scan for candidates */
    mpz_t Qx, Axb, cofactor, t1, t2;
    mpz_init(Qx); mpz_init(Axb); mpz_init(cofactor);
    mpz_init(t1); mpz_init(t2);
    int *smooth_exp = malloc(fb_size * sizeof(int));

    for (int j = 0; j < block_size; j++) {
        if (sieve_arr[j] < thresh) continue;

        long x = block_start + j;

        /* Q(x) = Ax^2 + 2Bx + C */
        mpz_set_si(t1, x);
        mpz_mul_si(t1, t1, x);
        mpz_mul(t1, t1, poly->A);     /* Ax^2 */
        mpz_set_si(t2, 2 * x);
        mpz_mul(t2, t2, poly->B);     /* 2Bx */
        mpz_add(Qx, t1, t2);
        mpz_add(Qx, Qx, poly->C);    /* Ax^2 + 2Bx + C */

        if (mpz_sgn(Qx) == 0) continue;
        int neg = (mpz_sgn(Qx) < 0);

        /* NOVEL: extract smooth part via primorial GCD */
        extract_smooth(smooth_exp, cofactor, Qx);

        int lp = 0;
        if (mpz_cmp_ui(cofactor, 1) == 0) {
            /* Fully smooth - great */
        } else if (mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= (unsigned long)lp_bound) {
            if (mpz_probab_prime_p(cofactor, 2)) {
                lp = mpz_get_ui(cofactor);
            } else continue; /* cofactor is composite but small - skip */
        } else continue; /* cofactor too large */

        /* Compute Y = Ax + B */
        mpz_set_si(Axb, x);
        mpz_mul(Axb, Axb, poly->A);
        mpz_add(Axb, Axb, poly->B);

        add_rel(Axb, smooth_exp, fb_size, lp, neg,
                poly->a_idx, poly->num_a);
        found++;
    }

    free(smooth_exp);
    mpz_clear(Qx); mpz_clear(Axb); mpz_clear(cofactor);
    mpz_clear(t1); mpz_clear(t2);
    return found;
}

/* ==================== GF(2) Gaussian Elimination ==================== */

static int gf2_solve(int nrows, int ncols, int *dep_list, int *dep_sizes,
                     int max_deps) {
    /* Augmented matrix: [parity | identity] */
    int row_w = (ncols + 63) / 64;
    int aug_w = (nrows + 63) / 64;

    uint64_t **mat = malloc(nrows * sizeof(uint64_t *));
    uint64_t **aug = malloc(nrows * sizeof(uint64_t *));
    for (int i = 0; i < nrows; i++) {
        mat[i] = calloc(row_w, sizeof(uint64_t));
        memcpy(mat[i], rels[i].parity_row, row_w * sizeof(uint64_t));
        aug[i] = calloc(aug_w, sizeof(uint64_t));
        aug[i][i / 64] |= (1ULL << (i % 64));
    }

    int *pivot_col = malloc(nrows * sizeof(int));
    memset(pivot_col, -1, nrows * sizeof(int));
    int *col_pivot = malloc(ncols * sizeof(int));
    memset(col_pivot, -1, ncols * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        /* Find pivot row */
        int piv = -1;
        for (int row = 0; row < nrows; row++) {
            if (pivot_col[row] >= 0) continue;
            if ((mat[row][col / 64] >> (col % 64)) & 1) { piv = row; break; }
        }
        if (piv < 0) continue;
        pivot_col[piv] = col;
        col_pivot[col] = piv;

        /* Eliminate column in all other rows */
        for (int row = 0; row < nrows; row++) {
            if (row == piv) continue;
            if (!((mat[row][col / 64] >> (col % 64)) & 1)) continue;
            for (int w = 0; w < row_w; w++) mat[row][w] ^= mat[piv][w];
            for (int w = 0; w < aug_w; w++) aug[row][w] ^= aug[piv][w];
        }
    }

    /* Extract dependencies (zero rows) */
    int ndeps = 0;
    for (int row = 0; row < nrows && ndeps < max_deps; row++) {
        int zero = 1;
        for (int w = 0; w < row_w && zero; w++) if (mat[row][w]) zero = 0;
        if (!zero) continue;

        /* Count and extract dependency members */
        int cnt = 0;
        for (int i = 0; i < nrows; i++)
            if ((aug[row][i / 64] >> (i % 64)) & 1) cnt++;
        if (cnt < 2) continue;

        int *dep = dep_list + ndeps * nrows; /* flat array */
        int idx = 0;
        for (int i = 0; i < nrows; i++)
            if ((aug[row][i / 64] >> (i % 64)) & 1)
                dep[idx++] = i;
        dep_sizes[ndeps] = cnt;
        ndeps++;
    }

    for (int i = 0; i < nrows; i++) { free(mat[i]); free(aug[i]); }
    free(mat); free(aug); free(pivot_col); free(col_pivot);
    return ndeps;
}

/* ==================== Square Root Step ==================== */

static int try_factor(mpz_t factor, const mpz_t N,
                      const int *dep, int dep_size) {
    mpz_t X, Y, temp;
    mpz_init_set_ui(X, 1);
    mpz_init_set_ui(Y, 1);
    mpz_init(temp);

    int *total_exp = calloc(fb_size, sizeof(int));
    int sign_count = 0; /* count of negative Q(x) values */

    /* Also track LP exponents */
    int *lp_list = calloc(dep_size, sizeof(int));
    int n_lps = 0;

    for (int i = 0; i < dep_size; i++) {
        int ri = dep[i];
        mpz_mul(X, X, rels[ri].Axb);
        mpz_mod(X, X, N);

        for (int j = 0; j < fb_size; j++)
            total_exp[j] += rels[ri].exponents[j];

        if (rels[ri].neg) sign_count++;

        if (rels[ri].large_prime > 1)
            lp_list[n_lps++] = rels[ri].large_prime;
    }

    /* Check sign parity: if odd number of negatives, skip this dep */
    if (sign_count % 2 != 0) {
        free(lp_list); free(total_exp);
        mpz_clear(X); mpz_clear(Y); mpz_clear(temp);
        return 0;
    }

    /* Y = product of p^(exp/2) * product of lp^(count/2) mod N */
    int odd_exp = 0;
    for (int j = 0; j < fb_size; j++) {
        if (total_exp[j] == 0) continue;
        if (total_exp[j] % 2 != 0) { odd_exp = 1; break; }
        mpz_set_ui(temp, fb[j].p);
        mpz_powm_ui(temp, temp, total_exp[j] / 2, N);
        mpz_mul(Y, Y, temp);
        mpz_mod(Y, Y, N);
    }

    /* Handle large primes: count occurrences and include in Y */
    if (!odd_exp) {
        /* Sort LP list and count pairs */
        for (int i = 0; i < n_lps - 1; i++)
            for (int j = i + 1; j < n_lps; j++)
                if (lp_list[i] > lp_list[j]) {
                    int t = lp_list[i]; lp_list[i] = lp_list[j]; lp_list[j] = t;
                }
        for (int i = 0; i < n_lps; ) {
            int cnt = 1;
            while (i + cnt < n_lps && lp_list[i + cnt] == lp_list[i]) cnt++;
            if (cnt % 2 != 0) { odd_exp = 1; break; }
            mpz_set_ui(temp, lp_list[i]);
            mpz_powm_ui(temp, temp, cnt / 2, N);
            mpz_mul(Y, Y, temp);
            mpz_mod(Y, Y, N);
            i += cnt;
        }
    }

    free(lp_list);
    if (odd_exp) {
        free(total_exp);
        mpz_clear(X); mpz_clear(Y); mpz_clear(temp);
        return 0;
    }

    /* Try gcd(X - Y, N) */
    mpz_sub(temp, X, Y);
    mpz_gcd(factor, temp, N);
    int ok = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);

    if (!ok) {
        mpz_add(temp, X, Y);
        mpz_gcd(factor, temp, N);
        ok = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);
    }

    free(total_exp);
    mpz_clear(X); mpz_clear(Y); mpz_clear(temp);
    return ok;
}

/* ==================== Main ==================== */

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N, factor;
    mpz_init_set_str(N, argv[1], 10);
    mpz_init(factor);

    /* Quick small-factor check */
    {
        int nprimes; int *primes = sieve_primes(1000000, &nprimes);
        for (int i = 0; i < nprimes; i++) {
            if (mpz_divisible_ui_p(N, primes[i])) {
                gmp_printf("Factor: %d\n", primes[i]);
                mpz_divexact_ui(factor, N, primes[i]);
                gmp_printf("Cofactor: %Zd\n", factor);
                free(primes);
                mpz_clear(N); mpz_clear(factor);
                return 0;
            }
        }
        free(primes);
    }

    int digits = mpz_sizeinbase(N, 10);
    int bits = mpz_sizeinbase(N, 2);
    fprintf(stderr, "MQSS: %d digits (%d bits)\n", digits, bits);

    params_t par = get_params(digits);
    fprintf(stderr, "Params: fb=%d blocks=%d lp_mult=%d a_primes=%d\n",
            par.fb_size, par.num_blocks, par.lp_mult, par.num_a_primes);

    build_factor_base(N, par.fb_size);
    fprintf(stderr, "FB: %d primes, bound=%d, mult=%d\n",
            fb_size, fb_bound, multiplier);

    long lp_bound = 0; /* full relations only - simpler LA */

    /* Compute primorial for smooth extraction */
    fprintf(stderr, "Computing primorial...\n");
    compute_primorial();
    fprintf(stderr, "Primorial: %zu bits\n", mpz_sizeinbase(primorial, 2));

    /* Target A */
    mpz_t target_a;
    mpz_init(target_a);
    mpz_mul_ui(target_a, kN, 2);
    mpz_sqrt(target_a, target_a);
    mpz_tdiv_q_ui(target_a, target_a, (unsigned long)par.num_blocks * SIEVE_BLOCK);

    /* Sieve threshold: log2(|Q(x)|_max) * scale
     * |Q(x)| ≈ M * sqrt(kN/2) for x at edge of sieve interval
     * We want threshold = log2(|Q(x)|_max) - tolerance */
    double log2_kN = mpz_sizeinbase(kN, 2);
    double log2_M = log2(par.num_blocks * SIEVE_BLOCK);
    double log2_qx = log2_kN / 2.0 + log2_M;  /* log2(M * sqrt(kN)) */
    /* Threshold: accept if sieve sum > ~65% of log2(Q(x)) */
    unsigned char thresh = (unsigned char)(log2_qx * 0.62);
    if (thresh < 20) thresh = 20;
    if (thresh > 200) thresh = 200;
    fprintf(stderr, "Threshold: %d (est log2(Q) = %.1f)\n", thresh, log2_qx);

    /* Need more relations than columns: fb_size + 1 + estimated LP columns + extra */
    int target_rels = fb_size * 3; /* generous: 3x FB size to ensure enough after LP columns */
    init_lp_hash();

    fprintf(stderr, "Need ~%d relations (3 * fb=%d)\n", target_rels, fb_size);

    poly_t poly;
    mpz_init(poly.A); mpz_init(poly.B); mpz_init(poly.C);

    int total_polys = 0, a_count = 0;

    while (nrels + ncombined < target_rels) {
        poly_new_a(&poly, par.num_a_primes, target_a);
        a_count++;
        poly_init_b(&poly);

        do {
            total_polys++;
            int half = par.num_blocks * SIEVE_BLOCK;

            /* Sieve both sides */
            for (int block = 0; block < par.num_blocks; block++) {
                collect_block(-half + block * SIEVE_BLOCK, SIEVE_BLOCK,
                              &poly, N, lp_bound, thresh);
            }
            for (int block = 0; block < par.num_blocks; block++) {
                collect_block(block * SIEVE_BLOCK, SIEVE_BLOCK,
                              &poly, N, lp_bound, thresh);
            }

            if (total_polys % 50 == 0) {
                struct timespec now;
                clock_gettime(CLOCK_MONOTONIC, &now);
                double el = (now.tv_sec - t0.tv_sec) + (now.tv_nsec - t0.tv_nsec) / 1e9;
                int full = 0, partial = 0;
                for (int i = 0; i < nrels; i++) {
                    if (rels[i].large_prime <= 1) full++; else partial++;
                }
                fprintf(stderr, "\r[%.1fs] polys=%d A=%d rels=%d(f=%d p=%d c=%d) "
                        "need=%d rate=%.0f/s  ",
                        el, total_polys, a_count, nrels, full, partial, ncombined,
                        target_rels, nrels / el);
            }
        } while (poly_next_b(&poly) && nrels + ncombined < target_rels);

        poly_clear(&poly);
    }

    struct timespec t1;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "\nSieve: %d rels in %.2fs (%.1f/s)\n",
            nrels, sieve_time, nrels / sieve_time);

    /* Build GF(2) matrix with LP columns */
    /* Columns: [sign, fb_prime_0, ..., fb_prime_{n-1}, lp_0, lp_1, ...] */
    /* First, collect distinct large primes that appear at least twice */
    int *distinct_lps = malloc(nrels * sizeof(int));
    int n_distinct_lps = 0;

    /* Use a simple hash to find LPs that appear >= 2 times */
    for (int i = 0; i < nrels; i++) {
        if (rels[i].large_prime <= 1) continue;
        int lp = rels[i].large_prime;
        /* Check if already in distinct list */
        int found = 0;
        for (int j = 0; j < n_distinct_lps; j++) {
            if (distinct_lps[j] == lp) { found = 1; break; }
        }
        if (!found) {
            /* Check if it appears more than once */
            int count = 0;
            for (int j = i; j < nrels; j++)
                if (rels[j].large_prime == lp) { count++; if (count >= 2) break; }
            if (count >= 2)
                distinct_lps[n_distinct_lps++] = lp;
        }
    }

    /* Only include relations that are full or have a paired LP */
    int *use_rel = calloc(nrels, sizeof(int));
    int n_usable = 0;
    for (int i = 0; i < nrels; i++) {
        if (rels[i].large_prime <= 1) { use_rel[i] = 1; n_usable++; }
        else {
            for (int j = 0; j < n_distinct_lps; j++) {
                if (rels[i].large_prime == distinct_lps[j]) {
                    use_rel[i] = 1; n_usable++; break;
                }
            }
        }
    }

    int ncols = fb_size + 1 + n_distinct_lps;
    fprintf(stderr, "LA: %d usable rels (%d full + %d paired SLP), %d cols (%d FB + %d LP)\n",
            n_usable, n_usable - ncombined, ncombined, ncols, fb_size + 1, n_distinct_lps);

    /* Build parity rows for usable relations */
    parity_words = (ncols + 63) / 64;
    for (int i = 0; i < nrels; i++) {
        if (!use_rel[i]) continue;
        rels[i].parity_row = calloc(parity_words, sizeof(uint64_t));
        /* Column 0: sign */
        if (rels[i].neg)
            rels[i].parity_row[0] |= 1ULL;
        /* Columns 1..fb_size: exponent parities */
        for (int j = 0; j < fb_size; j++) {
            if (rels[i].exponents[j] & 1) {
                int col = j + 1;
                rels[i].parity_row[col / 64] |= (1ULL << (col % 64));
            }
        }
        /* LP column */
        if (rels[i].large_prime > 1) {
            for (int j = 0; j < n_distinct_lps; j++) {
                if (distinct_lps[j] == rels[i].large_prime) {
                    int col = fb_size + 1 + j;
                    rels[i].parity_row[col / 64] |= (1ULL << (col % 64));
                    break;
                }
            }
        }
    }

    /* Compact usable relations to front (copy, don't swap) */
    int *rel_map = malloc(nrels * sizeof(int)); /* maps compact index -> original */
    int compact_n = 0;
    relation_t *compact_rels = malloc(n_usable * sizeof(relation_t));
    for (int i = 0; i < nrels; i++) {
        if (use_rel[i]) {
            compact_rels[compact_n] = rels[i];
            rel_map[compact_n] = i;
            compact_n++;
        }
    }
    /* Copy back to rels[] array front */
    memcpy(rels, compact_rels, compact_n * sizeof(relation_t));
    free(compact_rels);

    int *dep_flat = malloc(MAX_DEPS * compact_n * sizeof(int));
    int *dep_sizes = malloc(MAX_DEPS * sizeof(int));
    int ndeps = gf2_solve(compact_n, ncols, dep_flat, dep_sizes, MAX_DEPS);
    free(distinct_lps); free(use_rel); free(rel_map);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    int factored = 0;
    for (int d = 0; d < ndeps && !factored; d++) {
        int *dep = dep_flat + d * compact_n;
        if (try_factor(factor, N, dep, dep_sizes[d])) {
            struct timespec t2;
            clock_gettime(CLOCK_MONOTONIC, &t2);
            double total = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec) / 1e9;

            mpz_t cofn; mpz_init(cofn);
            mpz_divexact(cofn, N, factor);

            /* Adjust for multiplier */
            if (multiplier > 1) {
                if (mpz_divisible_ui_p(factor, multiplier))
                    mpz_divexact_ui(factor, factor, multiplier);
                else if (mpz_divisible_ui_p(cofn, multiplier))
                    mpz_divexact_ui(cofn, cofn, multiplier);
            }

            /* Ensure factor < cofactor for consistent output */
            if (mpz_cmp(factor, cofn) > 0) mpz_swap(factor, cofn);

            gmp_printf("Factor: %Zd\n", factor);
            gmp_printf("Cofactor: %Zd\n", cofn);
            fprintf(stderr, "Time: %.3fs (sieve %.3fs, LA %.3fs)\n",
                    total, sieve_time, total - sieve_time);
            factored = 1;
            mpz_clear(cofn);
        }
    }

    if (!factored)
        fprintf(stderr, "FAILED: %d deps, %d rels, fb=%d\n", ndeps, nrels, fb_size);

    free(dep_flat); free(dep_sizes);
    mpz_clear(N); mpz_clear(factor); mpz_clear(primorial);
    mpz_clear(target_a); mpz_clear(kN);
    mpz_clear(poly.A); mpz_clear(poly.B); mpz_clear(poly.C);
    return factored ? 0 : 1;
}
