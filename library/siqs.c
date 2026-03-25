/*
 * Basic Quadratic Sieve (QS) for factoring semiprimes.
 * NOT a true SIQS — uses single polynomial with shifted starting points.
 * L[1/2] baseline implementation. Uses GMP.
 * Factors 30-46 digit semiprimes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ---------- tunables ---------- */
#define MAX_FB       4096
#define MAX_SIEVE    131072
#define EXTRA_RELS   20      /* extra relations beyond FB size */
#define MAX_RELS     (MAX_FB + EXTRA_RELS + 200)
#define BLOCK_SIZE   65536

/* Parameters indexed by digit count */
typedef struct {
    int fb_size;       /* factor base size */
    int sieve_radius;  /* M: sieve from -M to M */
    int large_prime;   /* large prime bound (multiple of largest FB prime) */
    double thresh_adj; /* threshold adjustment */
} SieveParams;

static SieveParams get_params(int digits) {
    if (digits <= 30) return (SieveParams){200,  32768, 0, 0.0};
    if (digits <= 32) return (SieveParams){300,  50000, 0, 0.0};
    if (digits <= 34) return (SieveParams){400,  65536, 0, 0.0};
    if (digits <= 36) return (SieveParams){600,  65536, 0, 0.0};
    if (digits <= 38) return (SieveParams){900,  80000, 0, 0.0};
    if (digits <= 40) return (SieveParams){1200, 100000, 0, 0.0};
    if (digits <= 42) return (SieveParams){1600, 100000, 0, 0.0};
    if (digits <= 44) return (SieveParams){2000, 120000, 0, 0.0};
    return (SieveParams){2500, 131072, 0, 0.0};
}

/* ---------- factor base ---------- */
typedef struct {
    unsigned int p;    /* prime */
    unsigned int r1;   /* sqrt(n) mod p, root 1 */
    unsigned int r2;   /* sqrt(n) mod p, root 2 */
    double logp;       /* log2(p) */
    unsigned int ainv; /* used in SIQS poly switching */
} FBEntry;

static int primes[MAX_FB + 100];
static int nprimes;

static void sieve_primes(int limit) {
    char *is_composite = calloc(limit + 1, 1);
    nprimes = 0;
    for (int i = 2; i <= limit; i++) {
        if (!is_composite[i]) {
            primes[nprimes++] = i;
            if (nprimes >= MAX_FB + 100) break;
            for (long long j = (long long)i * i; j <= limit; j += i)
                is_composite[j] = 1;
        }
    }
    free(is_composite);
}

/* Tonelli-Shanks for sqrt(n) mod p */
static unsigned int mod_sqrt(mpz_t n, unsigned int p) {
    unsigned long nmod = mpz_fdiv_ui(n, p);
    if (p == 2) return nmod & 1;
    if (nmod == 0) return 0;

    /* Check if QR */
    mpz_t tmp;
    mpz_init(tmp);

    /* Use GMP's mpz_legendre or manual */
    /* Simple: a^((p-1)/2) mod p */
    unsigned long exp = (p - 1) / 2;
    unsigned long result = 1, base = nmod;
    unsigned long mod = p;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    if (result != 1) { mpz_clear(tmp); return 0; } /* not QR */

    /* Tonelli-Shanks */
    unsigned long q = p - 1, s = 0;
    while ((q & 1) == 0) { q >>= 1; s++; }

    if (s == 1) {
        /* p ≡ 3 (mod 4) */
        unsigned long r = 1;
        base = nmod;
        exp = (p + 1) / 4;
        while (exp > 0) {
            if (exp & 1) r = (r * base) % mod;
            base = (base * base) % mod;
            exp >>= 1;
        }
        mpz_clear(tmp);
        return (unsigned int)r;
    }

    /* Find non-residue */
    unsigned long z = 2;
    while (1) {
        unsigned long t = 1;
        base = z; exp = (p - 1) / 2;
        while (exp > 0) {
            if (exp & 1) t = (t * base) % mod;
            base = (base * base) % mod;
            exp >>= 1;
        }
        if (t == p - 1) break;
        z++;
    }

    unsigned long M = s;
    /* c = z^q mod p */
    unsigned long c = 1;
    base = z; exp = q;
    while (exp > 0) {
        if (exp & 1) c = (c * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    /* t = n^q mod p */
    unsigned long t = 1;
    base = nmod; exp = q;
    while (exp > 0) {
        if (exp & 1) t = (t * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    /* r = n^((q+1)/2) mod p */
    unsigned long r = 1;
    base = nmod; exp = (q + 1) / 2;
    while (exp > 0) {
        if (exp & 1) r = (r * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }

    while (1) {
        if (t == 1) { mpz_clear(tmp); return (unsigned int)r; }
        unsigned long i = 0, tt = t;
        while (tt != 1) { tt = (tt * tt) % mod; i++; }
        unsigned long b = c;
        for (unsigned long j = 0; j < M - i - 1; j++)
            b = (b * b) % mod;
        M = i;
        c = (b * b) % mod;
        t = (t * c) % mod;
        r = (r * b) % mod;
    }
}

/* ---------- relation storage ---------- */
typedef struct {
    mpz_t x;           /* the x value (or Q(x) after adjustment) */
    mpz_t y;           /* the Y value = (ax+b) for sqrt extraction */
    int *exponents;    /* exponent vector over factor base (mod 2 stored in matrix) */
    int smooth;        /* 1 if smooth */
} Relation;

/* ---------- matrix for Gaussian elimination mod 2 ---------- */
/* Store as array of bit-vectors using unsigned long */
#define BITS_PER_WORD 64
typedef unsigned long BitWord;

typedef struct {
    int rows, cols;
    int words_per_row;
    BitWord *data;  /* rows * words_per_row */
} BitMatrix;

static BitMatrix *bmat_alloc(int rows, int cols) {
    BitMatrix *m = malloc(sizeof(BitMatrix));
    m->rows = rows;
    m->cols = cols;
    m->words_per_row = (cols + BITS_PER_WORD - 1) / BITS_PER_WORD;
    m->data = calloc((size_t)rows * m->words_per_row, sizeof(BitWord));
    return m;
}

static void bmat_free(BitMatrix *m) {
    free(m->data);
    free(m);
}

static inline void bmat_set(BitMatrix *m, int r, int c) {
    m->data[(size_t)r * m->words_per_row + c / BITS_PER_WORD] |= (1UL << (c % BITS_PER_WORD));
}

static inline int bmat_get(BitMatrix *m, int r, int c) {
    return (m->data[(size_t)r * m->words_per_row + c / BITS_PER_WORD] >> (c % BITS_PER_WORD)) & 1;
}

static inline void bmat_xor_row(BitMatrix *m, int dst, int src) {
    BitWord *d = m->data + (size_t)dst * m->words_per_row;
    BitWord *s = m->data + (size_t)src * m->words_per_row;
    for (int i = 0; i < m->words_per_row; i++)
        d[i] ^= s[i];
}

/* Gaussian elimination mod 2 on transpose approach:
 * We have nrels relations (rows) and mat_cols columns.
 * We want to find linear dependencies among the rows.
 * Standard approach: work column by column, use row swaps.
 */
static int gaussian_elim(BitMatrix *mat, BitWord **deps, int *ndeps) {
    int nrows = mat->rows;
    int ncols = mat->cols;

    /* Create augmented matrix: [M | I] where M is the exponent matrix
     * and I is the identity. After elimination, zero rows in M part
     * give us the dependency info from I part. */
    int aug_cols = ncols + nrows;
    int aug_wpr = (aug_cols + BITS_PER_WORD - 1) / BITS_PER_WORD;
    BitWord *aug = calloc((size_t)nrows * aug_wpr, sizeof(BitWord));

    for (int r = 0; r < nrows; r++) {
        /* Copy left part */
        for (int c = 0; c < ncols; c++) {
            if (bmat_get(mat, r, c))
                aug[(size_t)r * aug_wpr + c / BITS_PER_WORD] |= (1UL << (c % BITS_PER_WORD));
        }
        /* Set identity bit */
        int ic = ncols + r;
        aug[(size_t)r * aug_wpr + ic / BITS_PER_WORD] |= (1UL << (ic % BITS_PER_WORD));
    }

    /* Row-reduce column by column */
    int cur_row = 0;
    for (int col = 0; col < ncols && cur_row < nrows; col++) {
        /* Find pivot in this column at or below cur_row */
        int piv = -1;
        for (int r = cur_row; r < nrows; r++) {
            if ((aug[(size_t)r * aug_wpr + col / BITS_PER_WORD] >> (col % BITS_PER_WORD)) & 1) {
                piv = r;
                break;
            }
        }
        if (piv == -1) continue;

        /* Swap piv and cur_row */
        if (piv != cur_row) {
            for (int w = 0; w < aug_wpr; w++) {
                BitWord tmp = aug[(size_t)piv * aug_wpr + w];
                aug[(size_t)piv * aug_wpr + w] = aug[(size_t)cur_row * aug_wpr + w];
                aug[(size_t)cur_row * aug_wpr + w] = tmp;
            }
        }

        /* Eliminate all other rows with a 1 in this column */
        for (int r = 0; r < nrows; r++) {
            if (r == cur_row) continue;
            if ((aug[(size_t)r * aug_wpr + col / BITS_PER_WORD] >> (col % BITS_PER_WORD)) & 1) {
                for (int w = 0; w < aug_wpr; w++)
                    aug[(size_t)r * aug_wpr + w] ^= aug[(size_t)cur_row * aug_wpr + w];
            }
        }
        cur_row++;
    }

    /* Rows from cur_row..nrows-1 have zero left part => dependencies */
    *ndeps = 0;
    int dep_words = (nrows + BITS_PER_WORD - 1) / BITS_PER_WORD;
    *deps = NULL;

    int left_full_words = ncols / BITS_PER_WORD;
    int left_remaining_bits = ncols % BITS_PER_WORD;
    BitWord left_last_mask = left_remaining_bits ? ((1UL << left_remaining_bits) - 1) : 0;

    for (int r = 0; r < nrows; r++) {
        /* Check if left part (columns 0..ncols-1) is zero */
        int is_zero = 1;
        for (int w = 0; w < left_full_words && is_zero; w++) {
            if (aug[(size_t)r * aug_wpr + w] != 0) is_zero = 0;
        }
        if (is_zero && left_remaining_bits) {
            if (aug[(size_t)r * aug_wpr + left_full_words] & left_last_mask)
                is_zero = 0;
        }

        if (is_zero) {
            *deps = realloc(*deps, (size_t)(*ndeps + 1) * dep_words * sizeof(BitWord));
            BitWord *dep = *deps + (size_t)(*ndeps) * dep_words;
            memset(dep, 0, dep_words * sizeof(BitWord));
            for (int i = 0; i < nrows; i++) {
                int ic = ncols + i;
                if ((aug[(size_t)r * aug_wpr + ic / BITS_PER_WORD] >> (ic % BITS_PER_WORD)) & 1)
                    dep[i / BITS_PER_WORD] |= (1UL << (i % BITS_PER_WORD));
            }
            (*ndeps)++;
        }
    }

    free(aug);
    return *ndeps;
}

/* ---------- SIQS core ---------- */

/* Compute modular inverse using GMP */
static unsigned long mod_inverse(unsigned long a, unsigned long m) {
    mpz_t ga, gm, gi;
    mpz_inits(ga, gm, gi, NULL);
    mpz_set_ui(ga, a);
    mpz_set_ui(gm, m);
    mpz_invert(gi, ga, gm);
    unsigned long result = mpz_get_ui(gi);
    mpz_clears(ga, gm, gi, NULL);
    return result;
}

static int factor_siqs(mpz_t n, mpz_t factor1, mpz_t factor2) {
    /* Quick trial division first */
    for (int i = 0; i < 1000 && i < nprimes; i++) {
        if (mpz_divisible_ui_p(n, primes[i])) {
            mpz_set_ui(factor1, primes[i]);
            mpz_divexact(factor2, n, factor1);
            return 1;
        }
    }

    int digits = mpz_sizeinbase(n, 10);
    SieveParams params = get_params(digits);
    int fb_size = params.fb_size;
    int M = params.sieve_radius;

    /* Build factor base */
    FBEntry *fb = calloc(fb_size + 1, sizeof(FBEntry));
    int fb_len = 0;

    /* Add -1 as index 0 (for sign) */
    /* Actually, let's handle sign separately */

    /* Find primes where n is a QR */
    for (int i = 0; i < nprimes && fb_len < fb_size; i++) {
        int p = primes[i];
        unsigned int r = mod_sqrt(n, p);
        if (r == 0 && mpz_divisible_ui_p(n, p)) continue; /* p divides n */
        if (r == 0 && p != 2) continue; /* not QR */
        if (p == 2) {
            if (mpz_fdiv_ui(n, 8) != 1) {
                if (mpz_odd_p(n)) {
                    fb[fb_len].p = 2;
                    fb[fb_len].r1 = 1;
                    fb[fb_len].r2 = 1;
                    fb[fb_len].logp = 1.0;
                    fb_len++;
                }
                continue;
            }
        }
        /* Verify: r^2 ≡ n (mod p) */
        if (p > 2) {
            unsigned long check = ((unsigned long)r * r) % p;
            if (check != mpz_fdiv_ui(n, p)) continue;
        }
        fb[fb_len].p = p;
        fb[fb_len].r1 = r;
        fb[fb_len].r2 = (p - r) % p;
        fb[fb_len].logp = log2(p);
        fb_len++;
    }

    if (fb_len < 20) {
        free(fb);
        return 0;
    }

    int target_rels = fb_len + EXTRA_RELS;

    /* Relation storage */
    int nrels = 0;
    mpz_t *rel_y = malloc(MAX_RELS * sizeof(mpz_t));
    mpz_t *rel_qx = malloc(MAX_RELS * sizeof(mpz_t));
    int **rel_exp = malloc(MAX_RELS * sizeof(int *));
    for (int i = 0; i < MAX_RELS; i++) {
        mpz_init(rel_y[i]);
        mpz_init(rel_qx[i]);
        rel_exp[i] = calloc(fb_len + 1, sizeof(int)); /* +1 for sign */
    }

    /* Sieve array */
    double *sieve = malloc(2 * M * sizeof(double));

    /* Threshold for smoothness */
    mpz_t sqrtn, temp, temp2, ax, bpoly, qx;
    mpz_inits(sqrtn, temp, temp2, ax, bpoly, qx, NULL);
    mpz_sqrt(sqrtn, n);

    /* For SIQS, we generate polynomials g(x) = (ax+b)^2 - n
     * where a = q^2 for some prime q, b chosen so b^2 ≡ n (mod a)
     *
     * Simplified: use single-polynomial QS for now with polynomial
     * Q(x) = (x + ceil(sqrt(n)))^2 - n
     */

    /* Actually, let's do proper MPQS with multiple polynomials */
    /* For simplicity and correctness, use basic QS first:
     * Q(x) = (x + s)^2 - n, where s = ceil(sqrt(n))
     * Sieve x from -M to M-1
     */

    mpz_add_ui(sqrtn, sqrtn, 1); /* s = ceil(sqrt(n)) */

    /* Compute approximate log2 of Q(x) at edges for threshold */
    /* Q(M) ≈ 2*M*sqrt(n), so log2(Q) ≈ log2(2*M) + bits(n)/2 */
    double log2_target = log2(2.0 * M) + mpz_sizeinbase(n, 2) / 2.0;
    double threshold = log2_target - params.thresh_adj;
    /* We want to detect values that are smooth over FB.
     * A value is "probably smooth" if the sum of logs of FB primes
     * that divide it is close to log2(Q(x)).
     * Use a threshold somewhat below the target. */
    double thresh_cutoff = threshold * 0.75;  /* Generous cutoff */

    int poly_count = 0;
    int max_polys = 10000;

    /* MPQS: choose 'a' values from products of FB primes */
    /* For simplicity and speed, let's do basic QS with multiple starting points */

    /* Basic QS sieving */
    while (nrels < target_rels && poly_count < max_polys) {
        /* Polynomial: Q(x) = (x + s + poly_count*2*M)^2 - n */
        /* Actually shift s by poly_count * 2M */
        mpz_t s;
        mpz_init(s);
        mpz_sqrt(s, n);
        mpz_add_ui(s, s, 1);
        if (poly_count > 0) {
            mpz_add_ui(s, s, (unsigned long)poly_count * 2 * M);
        }

        /* Initialize sieve to 0 */
        memset(sieve, 0, 2 * M * sizeof(double));

        /* For each prime in FB, find starting positions */
        for (int i = 0; i < fb_len; i++) {
            unsigned int p = fb[i].p;
            if (p == 0) continue;

            /* We need (x + s)^2 ≡ 0 (mod p), so x + s ≡ ±r (mod p)
             * x ≡ r - s (mod p) and x ≡ -r - s (mod p)
             * Map x in [-M, M-1] to index x + M in [0, 2M-1]
             */
            unsigned long s_mod = mpz_fdiv_ui(s, p);
            int r1 = ((int)fb[i].r1 - (int)s_mod % (int)p + 2 * (int)p) % (int)p;
            int r2 = ((int)fb[i].r2 - (int)s_mod % (int)p + 2 * (int)p) % (int)p;

            /* Adjust for sieve offset (-M maps to index 0) */
            int start1 = (r1 + M % (int)p) % (int)p;  /* position of first hit in [-M..] */
            if (start1 < 0) start1 += p;

            /* Sieve with this prime */
            double lp = fb[i].logp;

            /* For powers of 2, handle separately */
            if (p == 2) {
                for (int j = start1; j < 2 * M; j += p)
                    sieve[j] += lp;
                continue;
            }

            /* Two roots */
            /* x ≡ r1 - s (mod p): starting from -M, first x is: */
            /* We want the smallest non-negative idx s.t. (idx - M + s) ≡ r1 (mod p) */
            /* idx ≡ r1 - s + M (mod p) */
            long long idx1 = (((long long)fb[i].r1 - (long long)s_mod + (long long)M) % (long long)p + p) % p;
            long long idx2 = (((long long)fb[i].r2 - (long long)s_mod + (long long)M) % (long long)p + p) % p;

            for (long long j = idx1; j < 2 * M; j += p)
                sieve[j] += lp;
            if (idx1 != idx2) {
                for (long long j = idx2; j < 2 * M; j += p)
                    sieve[j] += lp;
            }
        }

        /* Check sieve for smooth values */
        for (int idx = 0; idx < 2 * M && nrels < target_rels; idx++) {
            if (sieve[idx] < thresh_cutoff) continue;

            /* Compute Q(x) = (x + s)^2 - n, where x = idx - M */
            long x = idx - M;
            mpz_set_si(temp, x);
            mpz_add(temp, temp, s);   /* temp = x + s */
            mpz_set(temp2, temp);     /* save x + s for Y value */
            mpz_mul(temp, temp, temp); /* temp = (x+s)^2 */
            mpz_sub(qx, temp, n);     /* qx = (x+s)^2 - n */

            if (mpz_sgn(qx) == 0) continue;

            /* Trial divide by factor base */
            mpz_t rem;
            mpz_init(rem);
            mpz_abs(rem, qx);

            int *exps = rel_exp[nrels];
            memset(exps, 0, (fb_len + 1) * sizeof(int));

            /* Handle sign */
            if (mpz_sgn(qx) < 0) {
                exps[0] = 1; /* sign bit */
            }

            int fully_factored = 1;
            for (int i = 0; i < fb_len; i++) {
                unsigned int p = fb[i].p;
                while (mpz_divisible_ui_p(rem, p)) {
                    mpz_divexact_ui(rem, rem, p);
                    exps[i + 1]++;
                }
            }

            if (mpz_cmp_ui(rem, 1) == 0) {
                /* Smooth! */
                mpz_set(rel_y[nrels], temp2);  /* Y = x + s */
                mpz_set(rel_qx[nrels], qx);
                nrels++;
            }
            mpz_clear(rem);
        }

        mpz_clear(s);
        poly_count++;

        /* Print progress every 100 polys */
        if (poly_count % 200 == 0) {
            fprintf(stderr, "  poly %d, rels %d/%d\n", poly_count, nrels, target_rels);
        }
    }

    fprintf(stderr, "  Collected %d relations with %d factor base primes\n", nrels, fb_len);

    if (nrels < fb_len + 2) {
        fprintf(stderr, "  Not enough relations\n");
        /* cleanup */
        for (int i = 0; i < MAX_RELS; i++) {
            mpz_clear(rel_y[i]);
            mpz_clear(rel_qx[i]);
            free(rel_exp[i]);
        }
        free(rel_y); free(rel_qx); free(rel_exp);
        free(fb); free(sieve);
        mpz_clears(sqrtn, temp, temp2, ax, bpoly, qx, NULL);
        return 0;
    }

    /* Build exponent matrix mod 2 */
    /* Rows = relations, Cols = fb_len + 1 (sign + primes) */
    int mat_cols = fb_len + 1;
    BitMatrix *mat = bmat_alloc(nrels, mat_cols);

    for (int r = 0; r < nrels; r++) {
        for (int c = 0; c < mat_cols; c++) {
            if (rel_exp[r][c] & 1)
                bmat_set(mat, r, c);
        }
    }

    /* Gaussian elimination */
    BitWord *deps = NULL;
    int ndeps = 0;
    gaussian_elim(mat, &deps, &ndeps);

    fprintf(stderr, "  Found %d dependencies\n", ndeps);

    int dep_words = (nrels + BITS_PER_WORD - 1) / BITS_PER_WORD;
    int found = 0;

    for (int d = 0; d < ndeps && !found; d++) {
        BitWord *dep = deps + (size_t)d * dep_words;

        /* Compute X = product of Y values mod n */
        /* Compute total exponents, then Y = product of primes^(exp/2) mod n */
        mpz_t X, Y;
        mpz_init_set_ui(X, 1);
        mpz_init_set_ui(Y, 1);

        int *total_exp = calloc(mat_cols, sizeof(int));

        for (int r = 0; r < nrels; r++) {
            if (!((dep[r / BITS_PER_WORD] >> (r % BITS_PER_WORD)) & 1)) continue;

            /* X *= Y_r (mod n) */
            mpz_mul(X, X, rel_y[r]);
            mpz_mod(X, X, n);

            /* Accumulate exponents */
            for (int c = 0; c < mat_cols; c++)
                total_exp[c] += rel_exp[r][c];
        }

        /* Check all exponents are even */
        int all_even = 1;
        for (int c = 0; c < mat_cols; c++) {
            if (total_exp[c] & 1) { all_even = 0; break; }
        }

        if (!all_even) {
            free(total_exp);
            mpz_clear(X); mpz_clear(Y);
            continue;
        }

        /* Y = product of p_i^(exp_i/2) mod n */
        /* Skip sign (index 0) */
        for (int c = 1; c < mat_cols; c++) {
            if (total_exp[c] == 0) continue;
            mpz_t pp;
            mpz_init(pp);
            mpz_set_ui(pp, fb[c - 1].p);
            mpz_powm_ui(pp, pp, total_exp[c] / 2, n);
            mpz_mul(Y, Y, pp);
            mpz_mod(Y, Y, n);
            mpz_clear(pp);
        }

        /* gcd(X - Y, n) and gcd(X + Y, n) */
        mpz_t g;
        mpz_init(g);

        mpz_sub(temp, X, Y);
        mpz_gcd(g, temp, n);

        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
            mpz_set(factor1, g);
            mpz_divexact(factor2, n, g);
            found = 1;
        } else {
            mpz_add(temp, X, Y);
            mpz_gcd(g, temp, n);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                mpz_set(factor1, g);
                mpz_divexact(factor2, n, g);
                found = 1;
            }
        }

        free(total_exp);
        mpz_clear(X); mpz_clear(Y); mpz_clear(g);
    }

    /* Cleanup */
    if (deps) free(deps);
    bmat_free(mat);
    for (int i = 0; i < MAX_RELS; i++) {
        mpz_clear(rel_y[i]);
        mpz_clear(rel_qx[i]);
        free(rel_exp[i]);
    }
    free(rel_y); free(rel_qx); free(rel_exp);
    free(fb); free(sieve);
    mpz_clears(sqrtn, temp, temp2, ax, bpoly, qx, NULL);

    return found;
}

/* ---------- main ---------- */

int main(int argc, char *argv[]) {
    /* Generate primes up to 1,000,000 */
    sieve_primes(1000000);

    /* Read semiprimes from stdin or command line */
    /* Usage: ./siqs "number" OR ./siqs < file_with_numbers */

    if (argc > 1) {
        /* Factor each argument */
        for (int i = 1; i < argc; i++) {
            mpz_t n, f1, f2;
            mpz_inits(n, f1, f2, NULL);
            mpz_set_str(n, argv[i], 10);

            int digits = mpz_sizeinbase(n, 10);
            fprintf(stderr, "Factoring %s (%d digits)...\n", argv[i], digits);

            struct timespec start, end;
            clock_gettime(CLOCK_MONOTONIC, &start);

            int ok = factor_siqs(n, f1, f2);

            clock_gettime(CLOCK_MONOTONIC, &end);
            double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

            if (ok) {
                gmp_printf("FACTORED %Zd = %Zd * %Zd  (%.3fs)\n", n, f1, f2, elapsed);
            } else {
                gmp_printf("FAILED %Zd  (%.3fs)\n", n, elapsed);
            }

            mpz_clears(n, f1, f2, NULL);
        }
    } else {
        /* Read numbers from stdin, one per line */
        char line[1024];
        while (fgets(line, sizeof(line), stdin)) {
            /* Strip whitespace */
            char *p = line;
            while (*p == ' ' || *p == '\t') p++;
            char *end = p + strlen(p) - 1;
            while (end > p && (*end == '\n' || *end == '\r' || *end == ' ')) *end-- = '\0';
            if (*p == '\0' || *p == '#') continue;

            mpz_t n, f1, f2;
            mpz_inits(n, f1, f2, NULL);
            if (mpz_set_str(n, p, 10) != 0) {
                fprintf(stderr, "Invalid number: %s\n", p);
                mpz_clears(n, f1, f2, NULL);
                continue;
            }

            int digits = mpz_sizeinbase(n, 10);
            fprintf(stderr, "Factoring %s (%d digits)...\n", p, digits);

            struct timespec start_t, end_t;
            clock_gettime(CLOCK_MONOTONIC, &start_t);

            int ok = factor_siqs(n, f1, f2);

            clock_gettime(CLOCK_MONOTONIC, &end_t);
            double elapsed = (end_t.tv_sec - start_t.tv_sec) + (end_t.tv_nsec - start_t.tv_nsec) / 1e9;

            if (ok) {
                gmp_printf("FACTORED %Zd = %Zd * %Zd  (%.3fs)\n", n, f1, f2, elapsed);
            } else {
                gmp_printf("FAILED %Zd  (%.3fs)\n", n, elapsed);
            }

            mpz_clears(n, f1, f2, NULL);
        }
    }

    return 0;
}
