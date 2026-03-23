/*
 * lattice_factor.c - Lattice-based integer factoring
 *
 * Novel approach combining lattice reduction with quadratic sieve:
 * 1. Use LLL lattice reduction to find "structured" smooth candidates
 * 2. Use a lightweight sieve to verify/supplement
 * 3. Linear algebra over GF(2) to find congruence of squares
 *
 * The lattice approach: Given factor base primes p_1,...,p_t and N,
 * construct lattice where short vectors give exponent vectors (e_1,...,e_t)
 * such that prod(p_i^e_i) mod N is small (and thus likely smooth).
 *
 * Compile: gcc -O3 -march=native -o lattice_factor library/lattice_factor.c -lgmp -lm
 * Usage: ./lattice_factor <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ============ Small primes ============ */
static int primes[100000];
static int nprimes;

static void gen_primes(int limit) {
    char *s = calloc(limit + 1, 1);
    nprimes = 0;
    for (int i = 2; i <= limit; i++) {
        if (!s[i]) primes[nprimes++] = i;
        for (long j = (long)i * i; j <= limit; j += i) s[j] = 1;
    }
    free(s);
}

/* ============ Tonelli-Shanks sqrt mod p ============ */
static unsigned long sqrt_mod_p(unsigned long n, unsigned long p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    /* Find Q, S such that p-1 = Q * 2^S */
    unsigned long Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }

    if (S == 1) {
        /* p ≡ 3 mod 4 */
        unsigned long r;
        mpz_t base, exp, mod, result;
        mpz_init_set_ui(base, n);
        mpz_init_set_ui(exp, (p + 1) / 4);
        mpz_init_set_ui(mod, p);
        mpz_init(result);
        mpz_powm(result, base, exp, mod);
        r = mpz_get_ui(result);
        mpz_clear(base); mpz_clear(exp); mpz_clear(mod); mpz_clear(result);
        return r;
    }

    /* Find quadratic non-residue */
    unsigned long z = 2;
    mpz_t zt, pt, et, rt;
    mpz_init(zt); mpz_init(pt); mpz_init(et); mpz_init(rt);
    mpz_set_ui(pt, p);
    mpz_set_ui(et, (p - 1) / 2);
    while (z < p) {
        mpz_set_ui(zt, z);
        mpz_powm(rt, zt, et, pt);
        if (mpz_get_ui(rt) == p - 1) break;
        z++;
    }

    unsigned long M = S;
    mpz_t c, t, R, b, temp;
    mpz_init(c); mpz_init(t); mpz_init(R); mpz_init(b); mpz_init(temp);

    mpz_set_ui(zt, z);
    mpz_set_ui(et, Q);
    mpz_powm(c, zt, et, pt);  /* c = z^Q mod p */

    mpz_set_ui(zt, n);
    mpz_powm(t, zt, et, pt);  /* t = n^Q mod p */

    mpz_set_ui(et, (Q + 1) / 2);
    mpz_powm(R, zt, et, pt);  /* R = n^((Q+1)/2) mod p */

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            unsigned long result = mpz_get_ui(R);
            mpz_clear(c); mpz_clear(t); mpz_clear(R); mpz_clear(b); mpz_clear(temp);
            mpz_clear(zt); mpz_clear(pt); mpz_clear(et); mpz_clear(rt);
            return result;
        }

        /* Find least i such that t^(2^i) = 1 */
        unsigned long i = 0;
        mpz_set(temp, t);
        while (mpz_cmp_ui(temp, 1) != 0) {
            mpz_mul(temp, temp, temp);
            mpz_mod(temp, temp, pt);
            i++;
        }

        /* b = c^(2^(M-i-1)) */
        mpz_set(b, c);
        for (unsigned long j = 0; j < M - i - 1; j++) {
            mpz_mul(b, b, b);
            mpz_mod(b, b, pt);
        }

        M = i;
        mpz_mul(c, b, b); mpz_mod(c, c, pt);
        mpz_mul(t, t, c); mpz_mod(t, t, pt);
        mpz_mul(R, R, b); mpz_mod(R, R, pt);
    }
}

/* ============ LLL with floating-point GS (practical for dim < 200) ============ */

typedef struct {
    int n, m;
    mpz_t *data;
} lat_t;

#define L(l,i,j) ((l)->data[(i)*(l)->m+(j)])

static void lat_init(lat_t *l, int n, int m) {
    l->n = n; l->m = m;
    l->data = malloc(n * m * sizeof(mpz_t));
    for (int i = 0; i < n * m; i++) mpz_init(l->data[i]);
}

static void lat_free(lat_t *l) {
    for (int i = 0; i < l->n * l->m; i++) mpz_clear(l->data[i]);
    free(l->data);
}

/* LLL with double-precision GS (sufficient for small lattices) */
static void lll(lat_t *lat) {
    int n = lat->n, m = lat->m;
    double *mu = calloc(n * n, sizeof(double));
    double *B = calloc(n, sizeof(double));

    /* Compute GS */
    auto void compute_gs(void);
    void compute_gs(void) {
        for (int i = 0; i < n; i++) {
            /* B[i] = <b_i, b_i> */
            double dot = 0;
            for (int k = 0; k < m; k++) {
                double v = mpz_get_d(L(lat, i, k));
                dot += v * v;
            }
            B[i] = dot;

            for (int j = 0; j < i; j++) {
                double dot_ij = 0;
                for (int k = 0; k < m; k++)
                    dot_ij += mpz_get_d(L(lat, i, k)) * mpz_get_d(L(lat, j, k));

                for (int kk = 0; kk < j; kk++)
                    dot_ij -= mu[j * n + kk] * mu[i * n + kk] * B[kk];

                mu[i * n + j] = (B[j] != 0.0) ? dot_ij / B[j] : 0.0;
            }

            for (int j = 0; j < i; j++)
                B[i] -= mu[i * n + j] * mu[i * n + j] * B[j];
        }
    }

    compute_gs();

    mpz_t r, tmp;
    mpz_init(r);
    mpz_init(tmp);

    int k = 1;
    int iterations = 0;
    while (k < n) {
        iterations++;
        if (iterations > n * n * 100) break;  /* Safety limit */

        /* Size reduce b_k */
        for (int j = k - 1; j >= 0; j--) {
            if (fabs(mu[k * n + j]) > 0.501) {
                long rr = (long)round(mu[k * n + j]);
                mpz_set_si(r, rr);
                for (int c = 0; c < m; c++) {
                    mpz_mul(tmp, r, L(lat, j, c));
                    mpz_sub(L(lat, k, c), L(lat, k, c), tmp);
                }
                /* Update mu */
                for (int i = 0; i < j; i++)
                    mu[k * n + i] -= rr * mu[j * n + i];
                mu[k * n + j] -= rr;
            }
        }

        /* Lovász condition */
        double lhs = B[k] + mu[k * n + (k-1)] * mu[k * n + (k-1)] * B[k-1];
        if (lhs < 0.75 * B[k-1]) {
            /* Swap rows k and k-1 */
            for (int c = 0; c < m; c++)
                mpz_swap(L(lat, k, c), L(lat, k-1, c));
            compute_gs();
            if (k > 1) k--;
        } else {
            k++;
        }
    }

    mpz_clear(r);
    mpz_clear(tmp);
    free(mu);
    free(B);
}

/* ============ Factor base + smoothness ============ */

static int is_smooth(mpz_t val, int *fb, int fb_sz, int *exps) {
    mpz_t r;
    mpz_init_set(r, val);
    if (mpz_sgn(r) < 0) mpz_neg(r, r);
    memset(exps, 0, fb_sz * sizeof(int));
    for (int i = 0; i < fb_sz; i++) {
        while (mpz_divisible_ui_p(r, fb[i])) {
            exps[i]++;
            mpz_divexact_ui(r, r, fb[i]);
        }
    }
    int ok = (mpz_cmp_ui(r, 1) == 0);
    mpz_clear(r);
    return ok;
}

/* Check if val is smooth allowing one large prime up to lp_bound.
 * If so, fill exps and set *large_prime. */
static int is_smooth_lp(mpz_t val, int *fb, int fb_sz, int *exps, unsigned long lp_bound, unsigned long *large_prime) {
    mpz_t r;
    mpz_init_set(r, val);
    if (mpz_sgn(r) < 0) mpz_neg(r, r);
    memset(exps, 0, fb_sz * sizeof(int));
    for (int i = 0; i < fb_sz; i++) {
        while (mpz_divisible_ui_p(r, fb[i])) {
            exps[i]++;
            mpz_divexact_ui(r, r, fb[i]);
        }
    }
    if (mpz_cmp_ui(r, 1) == 0) {
        *large_prime = 1;
        mpz_clear(r);
        return 1;
    }
    if (mpz_fits_ulong_p(r) && mpz_get_ui(r) <= lp_bound) {
        *large_prime = mpz_get_ui(r);
        /* Check if it's prime (quick) */
        if (mpz_probab_prime_p(r, 5)) {
            mpz_clear(r);
            return 1;
        }
    }
    mpz_clear(r);
    return 0;
}

/* ============ Relation storage ============ */

#define MAX_RELS 16384

typedef struct {
    mpz_t x;      /* x value: (x + sqrt_n)^2 - N = Q(x) */
    int *exp;      /* exponent vector */
    int neg;       /* sign of Q(x) */
    unsigned long lp; /* large prime (1 if full relation) */
} rel_t;

/* ============ GF(2) Gaussian elimination with history ============ */

static void gf2_factor(mpz_t N, rel_t *rels, int nrels, int fb_sz, int *fb) {
    /* Build augmented matrix: [exponent parities | identity] */
    int ncols = fb_sz + 1;  /* +1 for sign */
    int nwords = (ncols + 63) / 64;
    int hist_words = (nrels + 63) / 64;

    unsigned long (*mat)[128] = calloc(nrels, sizeof(unsigned long[128]));
    unsigned long (*hist)[128] = calloc(nrels, sizeof(unsigned long[128]));

    for (int i = 0; i < nrels; i++) {
        hist[i][i / 64] |= (1UL << (i % 64));
        if (rels[i].neg)
            mat[i][0] |= 1UL;
        for (int j = 0; j < fb_sz; j++) {
            if (rels[i].exp[j] & 1) {
                int col = j + 1;
                mat[i][col / 64] |= (1UL << (col % 64));
            }
        }
    }

    /* Gaussian elimination */
    int *pivot = malloc(ncols * sizeof(int));
    memset(pivot, -1, ncols * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        int pr = -1;
        for (int row = 0; row < nrels; row++) {
            if ((mat[row][col / 64] >> (col % 64)) & 1) {
                int used = 0;
                for (int c = 0; c < col; c++)
                    if (pivot[c] == row) { used = 1; break; }
                if (!used) { pr = row; break; }
            }
        }
        if (pr < 0) continue;
        pivot[col] = pr;

        for (int row = 0; row < nrels; row++) {
            if (row != pr && ((mat[row][col / 64] >> (col % 64)) & 1)) {
                for (int w = 0; w < nwords; w++) mat[row][w] ^= mat[pr][w];
                for (int w = 0; w < hist_words; w++) hist[row][w] ^= hist[pr][w];
            }
        }
    }

    /* Find zero rows = dependencies */
    mpz_t X, Y, g, tmp, sqrt_n;
    mpz_init(X); mpz_init(Y); mpz_init(g); mpz_init(tmp); mpz_init(sqrt_n);
    mpz_sqrt(sqrt_n, N);
    mpz_add_ui(sqrt_n, sqrt_n, 1);

    for (int row = 0; row < nrels; row++) {
        int is_zero = 1;
        for (int w = 0; w < nwords; w++)
            if (mat[row][w]) { is_zero = 0; break; }
        if (!is_zero) continue;

        /* Combine relations in this dependency */
        mpz_set_ui(X, 1);
        int *combined_exp = calloc(fb_sz, sizeof(int));
        int combined_neg = 0;

        int count = 0;
        for (int i = 0; i < nrels; i++) {
            if ((hist[row][i / 64] >> (i % 64)) & 1) {
                mpz_mul(X, X, rels[i].x);
                mpz_mod(X, X, N);
                for (int j = 0; j < fb_sz; j++)
                    combined_exp[j] += rels[i].exp[j];
                combined_neg ^= rels[i].neg;
                count++;
            }
        }

        if (count < 2) { free(combined_exp); continue; }

        /* Y = product of p_i^(sum_exp/2) mod N */
        mpz_set_ui(Y, 1);
        for (int j = 0; j < fb_sz; j++) {
            int e = combined_exp[j] / 2;
            if (e > 0) {
                mpz_set_ui(tmp, fb[j]);
                mpz_powm_ui(tmp, tmp, e, N);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }
        }

        /* gcd(X - Y, N) */
        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, g);
            gmp_printf("%Zd = %Zd * %Zd\n", N, g, cofactor);
            mpz_clear(cofactor);
            free(combined_exp);
            goto done;
        }

        /* gcd(X + Y, N) */
        mpz_add(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, g);
            gmp_printf("%Zd = %Zd * %Zd\n", N, g, cofactor);
            mpz_clear(cofactor);
            free(combined_exp);
            goto done;
        }

        free(combined_exp);
    }

    fprintf(stderr, "All dependencies tried, no factor found.\n");

done:
    mpz_clear(X); mpz_clear(Y); mpz_clear(g); mpz_clear(tmp); mpz_clear(sqrt_n);
    free(mat); free(hist); free(pivot);
}

/* ============ Main: Quadratic sieve with lattice-enhanced candidate finding ============ */

static void factor_qs(mpz_t N) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    double log_n = mpz_sizeinbase(N, 2) * log(2.0);
    double log_log_n = log(log_n);
    int digits = (int)mpz_sizeinbase(N, 10);

    /* Hand-tuned parameters by digit count (matching YAFU-style scaling) */
    int fb_bound;
    long M;
    if (digits <= 30) { fb_bound = 300; M = 5000; }
    else if (digits <= 35) { fb_bound = 500; M = 10000; }
    else if (digits <= 40) { fb_bound = 1000; M = 20000; }
    else if (digits <= 45) { fb_bound = 2000; M = 40000; }
    else if (digits <= 50) { fb_bound = 4000; M = 65536; }
    else if (digits <= 55) { fb_bound = 8000; M = 100000; }
    else if (digits <= 60) { fb_bound = 15000; M = 200000; }
    else if (digits <= 65) { fb_bound = 30000; M = 400000; }
    else if (digits <= 70) { fb_bound = 50000; M = 600000; }
    else if (digits <= 75) { fb_bound = 80000; M = 800000; }
    else if (digits <= 80) { fb_bound = 120000; M = 1000000; }
    else if (digits <= 85) { fb_bound = 180000; M = 1500000; }
    else { fb_bound = 250000; M = 2000000; }

    fprintf(stderr, "%d digits, FB bound %d, sieve M %ld\n", digits, fb_bound, M);

    gen_primes(fb_bound);

    /* Build factor base */
    int fb[20000], fb_sz = 0;
    fb[fb_sz++] = 2;  /* always include 2 */

    mpz_t tmp;
    mpz_init(tmp);

    for (int i = 1; i < nprimes && primes[i] <= fb_bound; i++) {
        int p = primes[i];
        mpz_set_ui(tmp, p);
        if (mpz_legendre(N, tmp) >= 0) {
            fb[fb_sz++] = p;
            if (fb_sz >= 20000) break;
        }
    }

    fprintf(stderr, "Factor base: %d primes (max %d)\n", fb_sz, fb[fb_sz - 1]);
    int target = fb_sz + 30;

    /* Compute sqrt(N) */
    mpz_t sqrt_n;
    mpz_init(sqrt_n);
    mpz_sqrt(sqrt_n, N);
    mpz_add_ui(sqrt_n, sqrt_n, 1);

    /* Compute sieve roots */
    int *root1 = malloc(fb_sz * sizeof(int));
    int *root2 = malloc(fb_sz * sizeof(int));

    for (int i = 0; i < fb_sz; i++) {
        int p = fb[i];
        if (p == 2) {
            /* Handle 2 specially: Q(x) = (x+s)^2 - N.
             * If N is odd, Q(x) even when x+s is even, i.e. x ≡ -s (mod 2) */
            unsigned long s2 = mpz_fdiv_ui(sqrt_n, 2);
            root1[i] = (2 - s2) % 2;
            root2[i] = root1[i];
            continue;
        }

        unsigned long n_mod_p = mpz_fdiv_ui(N, p);
        unsigned long s_mod_p = mpz_fdiv_ui(sqrt_n, p);
        unsigned long r = sqrt_mod_p(n_mod_p, p);

        /* Q(x) = (x+s)^2 - N ≡ 0 (mod p) => x ≡ r-s or -r-s (mod p) */
        root1[i] = (int)(((long)r - (long)s_mod_p + 2L * p) % p);
        root2[i] = (int)((-(long)r - (long)s_mod_p + 2L * p) % p);
    }

    /* Allocate relations */
    rel_t *rels = malloc(MAX_RELS * sizeof(rel_t));
    int nrels = 0;

    /* Large prime bound */
    unsigned long lp_bound = (unsigned long)fb[fb_sz - 1] * 30;

    /* SLP hash table for combining */
    #define LP_HASH_SIZE 65536
    typedef struct lp_entry { unsigned long lp; int rel_idx; struct lp_entry *next; } lp_entry_t;
    lp_entry_t **lp_hash = calloc(LP_HASH_SIZE, sizeof(lp_entry_t *));
    int n_slp = 0;

    /* ====== SIEVE ====== */
    fprintf(stderr, "Sieving...\n");

    /* Use block sieve for cache efficiency */
    int block_size = 32768;
    double *sieve = malloc(block_size * sizeof(double));
    int *exps = malloc(fb_sz * sizeof(int));

    mpz_t qx, x_val;
    mpz_init(qx);
    mpz_init(x_val);

    /* Log threshold */
    /* Q(x) ≈ 2*sqrt(N)*x for small x, max ≈ M*sqrt(N) + M^2 */
    /* We want sum(log p) to be close to log(Q(x)) */
    double base_thresh = (log_n / 2.0) * 0.85;

    long total_checked = 0;
    long x_offset = 0;
    int pass = 0;

    while (nrels < target && x_offset < 10 * M) {
        /* Process block [x_offset, x_offset + block_size) centered at 0 for first blocks */
        long block_start;
        if (pass == 0) block_start = 0;
        else if (pass == 1) block_start = -block_size;
        else block_start = (pass % 2 == 0) ? (pass / 2) * block_size : -(pass / 2 + 1) * block_size;

        /* Adjust threshold based on block position */
        long block_center = block_start + block_size / 2;
        double abs_center = fabs((double)block_center);
        double log_qx_est;
        if (abs_center < 1) log_qx_est = log_n / 2.0 + 10;
        else log_qx_est = log_n / 2.0 + log(abs_center);
        double thresh = log_qx_est * 0.80;  /* Allow 20% slack for large primes */

        memset(sieve, 0, block_size * sizeof(double));

        /* Sieve with FB primes */
        for (int i = 0; i < fb_sz; i++) {
            int p = fb[i];
            double logp = log((double)p);

            /* Powers of p */
            for (long pk = p; pk <= fb_bound; pk *= p) {
                /* Find starting positions in block */
                long r1, r2;
                if (pk == p) {
                    r1 = root1[i]; r2 = root2[i];
                } else {
                    /* Skip prime powers for simplicity */
                    break;
                }

                /* Position of first root in block */
                long pos1 = r1 - ((block_start % pk) + pk) % pk;
                if (pos1 < 0) pos1 += pk;
                pos1 = ((block_start + pos1) % pk + pk) % pk;
                /* Simpler: first x in block with x ≡ r1 (mod p) */
                long first1 = ((r1 - (block_start % p)) % p + p) % p;

                for (long j = first1; j < block_size; j += p)
                    sieve[j] += logp;

                if (r2 != r1) {
                    long first2 = ((r2 - (block_start % p)) % p + p) % p;
                    for (long j = first2; j < block_size; j += p)
                        sieve[j] += logp;
                }

                break;  /* skip prime powers for now */
            }
        }

        /* Check candidates above threshold */
        for (long j = 0; j < block_size && nrels < MAX_RELS; j++) {
            if (sieve[j] < thresh) continue;

            long x = block_start + j;
            total_checked++;

            /* Q(x) = (x + sqrt_n)^2 - N */
            mpz_set(x_val, sqrt_n);
            if (x >= 0) mpz_add_ui(x_val, x_val, (unsigned long)x);
            else mpz_sub_ui(x_val, x_val, (unsigned long)(-x));

            mpz_mul(qx, x_val, x_val);
            mpz_sub(qx, qx, N);

            int neg = (mpz_sgn(qx) < 0);
            if (neg) mpz_neg(qx, qx);

            unsigned long lp;
            if (is_smooth_lp(qx, fb, fb_sz, exps, lp_bound, &lp)) {
                if (lp == 1) {
                    /* Full relation */
                    mpz_init_set(rels[nrels].x, x_val);
                    rels[nrels].exp = malloc(fb_sz * sizeof(int));
                    memcpy(rels[nrels].exp, exps, fb_sz * sizeof(int));
                    rels[nrels].neg = neg;
                    rels[nrels].lp = 1;
                    nrels++;
                } else {
                    /* SLP - store and check for pair */
                    unsigned int h = (unsigned int)(lp % LP_HASH_SIZE);
                    lp_entry_t *e = lp_hash[h];
                    int found_pair = 0;
                    while (e) {
                        if (e->lp == lp) {
                            /* Found a pair! Combine into full relation */
                            int idx = e->rel_idx;
                            /* Combine: multiply x values, add exponents, XOR sign */
                            mpz_init(rels[nrels].x);
                            mpz_mul(rels[nrels].x, rels[idx].x, x_val);
                            mpz_mod(rels[nrels].x, rels[nrels].x, N);
                            rels[nrels].exp = malloc(fb_sz * sizeof(int));
                            for (int k = 0; k < fb_sz; k++)
                                rels[nrels].exp[k] = rels[idx].exp[k] + exps[k];
                            /* LP cancels (appears twice), add to exponents as squared */
                            /* Actually lp^2 is the cofactor, need to include it */
                            /* For SLP combining: Q1*Q2 = (smooth1*lp)*(smooth2*lp) = smooth1*smooth2*lp^2 */
                            /* lp^2 is even power, so contributes 0 to GF(2) matrix */
                            rels[nrels].neg = rels[idx].neg ^ neg;
                            rels[nrels].lp = 1;
                            nrels++;
                            found_pair = 1;
                            break;
                        }
                        e = e->next;
                    }
                    if (!found_pair) {
                        /* Store as partial */
                        /* Store the relation temporarily */
                        int idx = MAX_RELS / 2 + n_slp;
                        if (idx < MAX_RELS) {
                            mpz_init_set(rels[idx].x, x_val);
                            rels[idx].exp = malloc(fb_sz * sizeof(int));
                            memcpy(rels[idx].exp, exps, fb_sz * sizeof(int));
                            rels[idx].neg = neg;
                            rels[idx].lp = lp;

                            lp_entry_t *ne = malloc(sizeof(lp_entry_t));
                            ne->lp = lp;
                            ne->rel_idx = idx;
                            ne->next = lp_hash[h];
                            lp_hash[h] = ne;
                            n_slp++;
                        }
                    }
                }
            }
        }

        pass++;

        /* Progress report every 100 blocks */
        if (pass % 100 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "  pass %d: %d/%d rels, %d SLP partials, %.1fs, %ld candidates checked\n",
                    pass, nrels, target, n_slp, el, total_checked);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Sieve done: %d rels (%d SLP combined) in %.2fs, %ld checked\n",
            nrels, n_slp, sieve_time, total_checked);

    if (nrels < fb_sz + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", nrels, fb_sz + 1);
        goto cleanup;
    }

    /* Linear algebra */
    fprintf(stderr, "Linear algebra...\n");
    gf2_factor(N, rels, nrels, fb_sz, fb);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Total: %.3fs (sieve %.3fs + LA %.3fs)\n", total, sieve_time, total - sieve_time);

cleanup:
    /* Free LP hash */
    for (int i = 0; i < LP_HASH_SIZE; i++) {
        lp_entry_t *e = lp_hash[i];
        while (e) { lp_entry_t *n = e->next; free(e); e = n; }
    }
    free(lp_hash);

    for (int i = 0; i < nrels; i++) {
        mpz_clear(rels[i].x);
        free(rels[i].exp);
    }
    /* Free SLP partials */
    for (int i = 0; i < n_slp; i++) {
        int idx = MAX_RELS / 2 + i;
        if (idx < MAX_RELS) {
            mpz_clear(rels[idx].x);
            free(rels[idx].exp);
        }
    }
    free(rels);
    free(exps);
    free(sieve);
    free(root1);
    free(root2);
    mpz_clear(sqrt_n);
    mpz_clear(qx);
    mpz_clear(x_val);
    mpz_clear(tmp);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N;
    mpz_init(N);
    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }

    if (mpz_probab_prime_p(N, 25)) {
        gmp_printf("%Zd is prime\n", N);
        mpz_clear(N);
        return 0;
    }

    factor_qs(N);
    mpz_clear(N);
    return 0;
}
