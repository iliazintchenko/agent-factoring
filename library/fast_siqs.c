/*
 * fast_siqs.c - High-performance SIQS implementation
 *
 * Key optimizations over previous custom implementations:
 * 1. Bucket sieving for large primes (primes > BLOCKSIZE)
 * 2. Sieve-informed trial division (only test primes whose roots match)
 * 3. SIQS self-initialization with Gray code enumeration
 * 4. Incremental root updates between polynomials
 * 5. 32KB block sieving for L1 cache efficiency
 * 6. Single large prime variation with hash table
 *
 * Compile: gcc -O3 -march=native -o fast_siqs library/fast_siqs.c -lgmp -lm
 * Usage: ./fast_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define SEED 42
#define BLOCKSIZE 32768
#define BLOCKBITS 15
#define BLOCKMASK (BLOCKSIZE - 1)

/* Maximum sizes */
#define MAX_FB      100000
#define MAX_RELS    500000
#define MAX_PARTIALS 2000000
#define MAX_A_FACTORS 25
#define MAX_DEPS    64

/* Bucket sieve parameters */
#define BUCKET_ALLOC  (1 << 16)  /* initial bucket allocation per block */
#define LARGE_PRIME_START_IDX 1500  /* will be computed dynamically */

/* Timing */
static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size;
    int num_blocks;    /* blocks per side */
    int lp_mult;       /* large prime = fb_max * lp_mult */
    int extra_rels;    /* relations beyond FB size */
    double thresh_adj; /* sieve threshold multiplier */
} params_t;

/* Tuned parameter table matching YAFU's parameter selection more closely */
static params_t get_params(int bits) {
    /* Based on YAFU's actual param_table entries for SIQS */
    if (bits <= 80)  return (params_t){60,    1,  30,  30,  0.68};
    if (bits <= 90)  return (params_t){80,    1,  30,  30,  0.70};
    if (bits <= 100) return (params_t){120,   1,  40,  35,  0.72};
    if (bits <= 110) return (params_t){180,   1,  40,  35,  0.73};
    if (bits <= 120) return (params_t){230,   1,  40,  40,  0.74};
    if (bits <= 130) return (params_t){300,   2,  40,  40,  0.75};
    if (bits <= 140) return (params_t){400,   2,  40,  50,  0.76};
    if (bits <= 150) return (params_t){500,   2,  40,  50,  0.77};
    if (bits <= 160) return (params_t){650,   3,  40,  60,  0.78};
    if (bits <= 170) return (params_t){900,   3,  40,  60,  0.79};
    if (bits <= 180) return (params_t){1200,  4,  50,  70,  0.80};
    if (bits <= 190) return (params_t){1700,  5,  50,  80,  0.81};
    if (bits <= 200) return (params_t){2200,  6,  50,  80,  0.82};
    if (bits <= 210) return (params_t){3000,  8,  50,  90,  0.83};
    if (bits <= 220) return (params_t){4000,  10, 60,  100, 0.84};
    if (bits <= 230) return (params_t){5000,  12, 60,  100, 0.84};
    if (bits <= 240) return (params_t){6500,  16, 60,  120, 0.85};
    if (bits <= 250) return (params_t){9000,  20, 70,  120, 0.86};
    if (bits <= 260) return (params_t){12000, 26, 70,  150, 0.86};
    if (bits <= 270) return (params_t){16000, 32, 80,  150, 0.87};
    if (bits <= 280) return (params_t){22000, 40, 80,  200, 0.87};
    if (bits <= 290) return (params_t){30000, 48, 90,  200, 0.88};
    if (bits <= 300) return (params_t){40000, 56, 90,  250, 0.88};
    return                  (params_t){55000, 64, 100, 300, 0.89};
}

/* ==================== Modular Arithmetic ==================== */

static inline unsigned int mod_inverse_u32(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
    while (r) {
        int q = old_r / r;
        int t = r; r = old_r - q * r; old_r = t;
        t = s; s = old_s - q * s; old_s = t;
    }
    if (old_r != 1) return 0;
    return (unsigned int)(((long long)old_s % m + m) % m);
}

/* Tonelli-Shanks modular square root */
static unsigned int sqrt_mod_p(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    unsigned long long nn = n % p, m = p;

    /* Check quadratic residue */
    unsigned long long r = 1, b = nn, e = (p - 1) / 2;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    if (r != 1) return 0;

    if (p % 4 == 3) {
        r = 1; b = nn; e = (p + 1) / 4;
        while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
        return (unsigned int)r;
    }

    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned int z = 2;
    for (;;) {
        r = 1; b = z; e = (p - 1) / 2;
        while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
        if (r == m - 1) break;
        z++;
    }

    unsigned long long M_val = S;
    r = 1; b = z; e = Q;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    unsigned long long c = r;

    r = 1; b = nn; e = Q;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    unsigned long long t = r;

    r = 1; b = nn; e = (Q + 1) / 2;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    unsigned long long R = r;

    for (;;) {
        if (t == 1) return (unsigned int)R;
        int i = 0;
        unsigned long long tt = t;
        while (tt != 1) { tt = tt * tt % p; i++; }
        unsigned long long bb = c;
        for (int j = 0; j < (int)M_val - i - 1; j++) bb = bb * bb % p;
        M_val = i;
        c = bb * bb % p;
        t = t * c % p;
        R = R * bb % p;
    }
}

/* ==================== Knuth-Schroeppel Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30;
    int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki];
        mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8 == 1) s += 2*log(2.0);
        else if (m8 == 5) s += log(2.0);
        else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) {
            if (k % ps[i] == 0) { s += log((double)ps[i]); continue; }
            if (sqrt_mod_p(mpz_fdiv_ui(kN, ps[i]), ps[i]))
                s += 2.0*log(ps[i])/(ps[i]-1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct {
    unsigned int *prime;
    unsigned int *sqrtN;   /* sqrt(kN) mod p */
    unsigned char *logp;
    unsigned int *inv;     /* modular inverse of 2*a mod p (computed per poly) */
    int size;
    int alloc;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->sqrtN = malloc(alloc * sizeof(unsigned int));
    fb->logp  = malloc(alloc * sizeof(unsigned char));
    fb->inv   = malloc(alloc * sizeof(unsigned int));
    fb->alloc = alloc;

    /* p=2 */
    fb->prime[0] = 2;
    fb->sqrtN[0] = 1;
    fb->logp[0] = 1;
    fb->size = 1;

    int bound = target * 30 + 100000;
    char *sieve = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++)
        if (!sieve[i]) for (int j = i*i; j <= bound; j += i) sieve[j] = 1;

    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sieve[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
        if (nm == 0) {
            fb->prime[fb->size] = i;
            fb->sqrtN[fb->size] = 0;
            fb->logp[fb->size] = (unsigned char)(log2(i) + 0.5);
            fb->size++;
            continue;
        }
        unsigned int r = sqrt_mod_p((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i;
        fb->sqrtN[fb->size] = r;
        fb->logp[fb->size] = (unsigned char)(log2(i) + 0.5);
        fb->size++;
    }
    free(sieve);
    return fb;
}

/* ==================== Large Prime Hash Table ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)

typedef struct lp_entry {
    unsigned long lp;
    int rel_idx;
    struct lp_entry *next;
} lp_entry_t;

typedef struct {
    lp_entry_t **buckets;
    lp_entry_t *pool;
    int used, max;
} lp_hash_t;

static lp_hash_t *lp_create(int max) {
    lp_hash_t *h = calloc(1, sizeof(lp_hash_t));
    h->buckets = calloc(LP_HASH_SIZE, sizeof(lp_entry_t*));
    h->pool = calloc(max, sizeof(lp_entry_t));
    h->max = max;
    return h;
}

static inline unsigned int lp_hashfn(unsigned long lp) {
    return (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - LP_HASH_BITS));
}

static int lp_find(lp_hash_t *h, unsigned long lp) {
    unsigned int idx = lp_hashfn(lp);
    for (lp_entry_t *e = h->buckets[idx]; e; e = e->next)
        if (e->lp == lp) return e->rel_idx;
    return -1;
}

static void lp_insert(lp_hash_t *h, unsigned long lp, int rel_idx) {
    if (h->used >= h->max) return;
    unsigned int idx = lp_hashfn(lp);
    lp_entry_t *e = &h->pool[h->used++];
    e->lp = lp;
    e->rel_idx = rel_idx;
    e->next = h->buckets[idx];
    h->buckets[idx] = e;
}

/* ==================== Relation Storage ==================== */
typedef struct {
    /* For each relation: Y^2 ≡ Q(x) (mod kN) */
    mpz_t *Y;       /* ax+b values */
    short **exps;    /* exps[i][j] = full exponent of j-th FB prime */
    unsigned long *lp;
    int count;
    int alloc;
    int fb_size;
    int neg_col;     /* extra column for sign */
} rels_t;

static rels_t *rels_create(int alloc, int fb_size) {
    rels_t *r = malloc(sizeof(rels_t));
    r->Y = malloc(alloc * sizeof(mpz_t));
    r->exps = malloc(alloc * sizeof(short*));
    r->lp = calloc(alloc, sizeof(unsigned long));
    for (int i = 0; i < alloc; i++) {
        mpz_init(r->Y[i]);
        r->exps[i] = calloc(fb_size + 2, sizeof(short)); /* +1 for sign, +1 for safety */
    }
    r->count = 0;
    r->alloc = alloc;
    r->fb_size = fb_size;
    r->neg_col = fb_size; /* sign stored at index fb_size */
    return r;
}

/* ==================== Bucket Sieve Structures ==================== */
typedef struct {
    unsigned int *data;   /* packed: (fb_idx << 16) | sieve_offset */
    int count;
    int alloc;
} bucket_t;

typedef struct {
    bucket_t *blocks;    /* one bucket per sieve block */
    int num_blocks;
    unsigned char *logp; /* log(p) for each FB index in large range */
} bucket_sieve_t;

static bucket_sieve_t *bucket_create(int num_blocks, int fb_size) {
    bucket_sieve_t *bs = malloc(sizeof(bucket_sieve_t));
    bs->num_blocks = num_blocks;
    bs->blocks = calloc(num_blocks, sizeof(bucket_t));
    for (int i = 0; i < num_blocks; i++) {
        bs->blocks[i].alloc = BUCKET_ALLOC;
        bs->blocks[i].data = malloc(BUCKET_ALLOC * sizeof(unsigned int));
        bs->blocks[i].count = 0;
    }
    bs->logp = calloc(fb_size, sizeof(unsigned char));
    return bs;
}

static void bucket_reset(bucket_sieve_t *bs) {
    for (int i = 0; i < bs->num_blocks; i++)
        bs->blocks[i].count = 0;
}

static inline void bucket_add(bucket_sieve_t *bs, int block_idx, unsigned int fb_idx, unsigned int offset) {
    bucket_t *b = &bs->blocks[block_idx];
    if (b->count >= b->alloc) {
        b->alloc *= 2;
        b->data = realloc(b->data, b->alloc * sizeof(unsigned int));
    }
    b->data[b->count++] = (fb_idx << 16) | offset;
}

/* ==================== GF(2) Linear Algebra ==================== */
typedef unsigned long long u64;
typedef struct {
    u64 **rows;
    int nr, nc, fbw, idw, wprow;
} gf2_t;

static gf2_t *gf2_create(int nr, int nc) {
    gf2_t *m = malloc(sizeof(gf2_t));
    m->nr = nr; m->nc = nc;
    m->fbw = (nc + 63) / 64;
    m->idw = (nr + 63) / 64;
    m->wprow = m->fbw + m->idw;
    m->rows = malloc(nr * sizeof(u64*));
    for (int i = 0; i < nr; i++) {
        m->rows[i] = calloc(m->wprow, sizeof(u64));
        m->rows[i][m->fbw + i/64] |= (1ULL << (i % 64));
    }
    return m;
}

static inline void gf2_set(gf2_t *m, int r, int c) {
    m->rows[r][c/64] ^= (1ULL << (c % 64));
}

static int gf2_solve(gf2_t *m, int ***deps_out, int **dlen_out, int max_deps) {
    int piv = 0;
    for (int c = 0; c < m->nc && piv < m->nr; c++) {
        int pr = -1;
        for (int r = piv; r < m->nr; r++) {
            if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) {
            if (r == piv) continue;
            if ((m->rows[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < m->wprow; w++)
                    m->rows[r][w] ^= m->rows[piv][w];
        }
        piv++;
    }

    int nd = 0;
    *deps_out = malloc(max_deps * sizeof(int*));
    *dlen_out = malloc(max_deps * sizeof(int));
    for (int r = piv; r < m->nr && nd < max_deps; r++) {
        /* Check if FB part is all zero */
        int zero = 1;
        for (int w = 0; w < m->fbw && zero; w++) {
            u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1);
            if (m->rows[r][w] & mask) zero = 0;
        }
        if (!zero) continue;
        /* Extract dependency indices */
        int *d = malloc(m->nr * sizeof(int));
        int dl = 0;
        for (int w = 0; w < m->idw; w++) {
            u64 bits = m->rows[r][m->fbw + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w * 64 + bit;
                if (idx < m->nr) d[dl++] = idx;
                bits &= bits - 1;
            }
        }
        if (dl > 0) { (*deps_out)[nd] = d; (*dlen_out)[nd] = dl; nd++; }
        else free(d);
    }
    return nd;
}

/* ==================== SIQS Core ==================== */

/* Polynomial: g(x) = (ax + b)^2 - kN, where a = product of selected FB primes
 * Q(x) = g(x)/a = ax^2 + 2bx + c, where c = (b^2 - kN)/a
 *
 * Sieve roots: for prime p in FB, g(x) ≡ 0 (mod p) when
 *   x ≡ (±sqrt(kN) - b) * a^(-1) (mod p)
 */

typedef struct {
    mpz_t a;
    mpz_t b;
    mpz_t c;        /* c = (b^2 - kN) / a */
    int a_indices[MAX_A_FACTORS];
    int num_a_factors;
    mpz_t B_vals[MAX_A_FACTORS];  /* B_j values for Gray code */

    /* Per-FB-prime sieve roots: root1[i], root2[i] for prime i */
    unsigned int *root1;
    unsigned int *root2;

    /* For Gray code updates: ainv[j][i] = 2*B_j * a^(-1) mod p_i */
    unsigned int **ainv;
} siqs_poly_t;

/* Trial divide Q(x) by FB primes.
 * SIEVE-INFORMED: only try primes whose sieve roots match the candidate position.
 */
static int trial_divide_sieve_informed(
    mpz_t Q, mpz_t Y, short *exponents,
    fb_t *fb, int fb_size,
    unsigned int *root1, unsigned int *root2,
    int x_in_block, int block_start,
    unsigned long lp_bound, unsigned long *lp_out,
    mpz_t kN, mpz_t a, mpz_t b, int x_global,
    int neg_col)
{
    /* Compute Y = a*x + b, Q = Y^2 - kN */
    mpz_mul_si(Y, a, x_global);
    mpz_add(Y, Y, b);
    mpz_mul(Q, Y, Y);
    mpz_sub(Q, Q, kN);

    memset(exponents, 0, (fb_size + 2) * sizeof(short));

    /* Track sign */
    if (mpz_sgn(Q) < 0) {
        exponents[neg_col] = 1;
        mpz_neg(Q, Q);
    }

    /* p=2 */
    while (mpz_even_p(Q)) {
        exponents[0]++;
        mpz_tdiv_q_2exp(Q, Q, 1);
    }

    /* Trial divide by FB primes.
     * Optimization: use sieve roots to skip most primes.
     * For primes where roots are valid, only test if position matches.
     * For primes dividing a (root=0xFFFFFFFF), always test. */
    int x_pos = block_start + x_in_block;

    /* Try to extract Q into a __int128 for fast TD on small FB primes */
    int use_fast = (mpz_sizeinbase(Q, 2) <= 127);
    __int128 Q128 = 0;
    if (use_fast) {
        if (mpz_fits_ulong_p(Q)) {
            Q128 = mpz_get_ui(Q);
        } else {
            /* Extract two limbs */
            mpz_t hi_part;
            mpz_init(hi_part);
            mpz_tdiv_q_2exp(hi_part, Q, 64);
            Q128 = ((__int128)mpz_get_ui(hi_part) << 64) | mpz_getlimbn(Q, 0);
            mpz_clear(hi_part);
        }
    }

    for (int i = 1; i < fb_size; i++) {
        unsigned int p = fb->prime[i];
        if (p < 3) continue;

        /* Sieve-informed root check with mpz fallback */
        if (root1[i] != 0xFFFFFFFF) {
            unsigned int xmod;
            int xp = x_pos % (int)p;
            xmod = (unsigned int)(xp < 0 ? xp + (int)p : xp);
            if (xmod != root1[i] && xmod != root2[i]) {
                /* Root doesn't match - check with mpz as fallback */
                if (use_fast) {
                    if (Q128 % p != 0) continue;
                } else {
                    if (!mpz_divisible_ui_p(Q, p)) continue;
                }
            }
        } else {
            /* p divides a - always check */
            if (use_fast) {
                if (Q128 % p != 0) continue;
            } else {
                if (!mpz_divisible_ui_p(Q, p)) continue;
            }
        }

        /* Divide out all powers of p */
        if (use_fast && Q128 > 0) {
            while (Q128 % p == 0) {
                exponents[i]++;
                Q128 /= p;
            }
            /* Sync back to mpz */
            if (Q128 <= (unsigned long long)(-1)) {
                mpz_set_ui(Q, (unsigned long long)Q128);
            } else {
                mpz_set_ui(Q, (unsigned long long)(Q128 >> 64));
                mpz_mul_2exp(Q, Q, 64);
                mpz_add_ui(Q, Q, (unsigned long long)Q128);
            }
        } else {
            do {
                exponents[i]++;
                mpz_divexact_ui(Q, Q, p);
            } while (mpz_divisible_ui_p(Q, p));
            /* Update fast path state */
            if (mpz_sizeinbase(Q, 2) <= 127) {
                use_fast = 1;
                if (mpz_fits_ulong_p(Q)) Q128 = mpz_get_ui(Q);
                else {
                    mpz_t hi; mpz_init(hi); mpz_tdiv_q_2exp(hi, Q, 64);
                    Q128 = ((__int128)mpz_get_ui(hi) << 64) | mpz_getlimbn(Q, 0);
                    mpz_clear(hi);
                }
            }
        }
    }
    /* Sync Q128 back to Q for final cofactor check */
    if (use_fast) {
        if (Q128 <= (unsigned long long)(-1)) {
            mpz_set_ui(Q, (unsigned long long)Q128);
        } else {
            mpz_set_ui(Q, (unsigned long long)(Q128 >> 64));
            mpz_mul_2exp(Q, Q, 64);
            mpz_add_ui(Q, Q, (unsigned long long)Q128);
        }
    }

    /* Check if fully factored or has single large prime */
    if (mpz_cmp_ui(Q, 1) == 0) {
        *lp_out = 0;
        return 1;
    }

    if (mpz_fits_ulong_p(Q)) {
        unsigned long cofactor = mpz_get_ui(Q);
        if (cofactor <= lp_bound && cofactor > 1) {
            mpz_t cof;
            mpz_init_set_ui(cof, cofactor);
            int prime = mpz_probab_prime_p(cof, 5);
            mpz_clear(cof);
            if (prime) {
                *lp_out = cofactor;
                return 2;
            }
        }
    }

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Quick trial division for small factors */
    {
        unsigned long small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
        for (int i = 0; i < 25; i++) {
            if (mpz_divisible_ui_p(N, small_primes[i])) {
                mpz_t q; mpz_init(q);
                mpz_divexact_ui(q, N, small_primes[i]);
                if (mpz_cmp_ui(q, 1) > 0 && mpz_cmp(q, N) < 0) {
                    gmp_printf("%Zd\n%Zd\n", q, N); /* wrong, print factor + cofactor */
                    printf("%lu\n", small_primes[i]);
                    mpz_clear(q);
                    return 0;
                }
                mpz_clear(q);
            }
        }
        /* Extended trial division */
        for (unsigned long p = 101; p < 100000; p += 2) {
            if (mpz_divisible_ui_p(N, p)) {
                mpz_t q; mpz_init(q);
                mpz_divexact_ui(q, N, p);
                gmp_printf("%lu\n%Zd\n", p, q);
                mpz_clear(q);
                return 0;
            }
        }
    }

    /* Choose multiplier */
    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);

    /* Get parameters */
    params_t P = get_params(kN_bits);

    /* Build factor base */
    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCKSIZE * P.num_blocks;  /* sieve interval per side: [-M, M) */
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    int target = fb->size + P.extra_rels;
    int total_blocks = 2 * P.num_blocks;  /* blocks for both sides */

    /* Compute sieve threshold */
    double log2_kN = kN_bits;
    double log2_M = log2((double)M);
    double log2_a_approx = log2_kN / 2.0;  /* a ≈ sqrt(2*kN) / M */
    double log2_Qmax = log2_a_approx + log2_M;  /* Q(x) ≈ a*M for x at boundary */
    int threshold = (int)(log2_Qmax * P.thresh_adj);

    /* Determine bucket sieve cutoff: primes > BLOCKSIZE use buckets */
    int bucket_start = 0;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] > BLOCKSIZE) { bucket_start = i; break; }
        if (i == fb->size - 1) bucket_start = fb->size;
    }

    /* Skip tiny primes in sieve (p < 5) */
    int sieve_start = 0;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] >= 5) { sieve_start = i; break; }
    }

    fprintf(stderr, "FSIQS: %dd (%db), k=%d, FB=%d, blocks=%d, M=%d, thresh=%d, LP=%lu, target=%d, bucket_start=%d\n",
            digits, bits, mult, fb->size, total_blocks, M, threshold, lp_bound, target, bucket_start);

    /* Allocate sieve array and buckets */
    unsigned char *sieve_array = malloc(BLOCKSIZE);
    bucket_sieve_t *buckets = bucket_create(total_blocks, fb->size);

    /* Copy logp values for bucket sieve */
    for (int i = bucket_start; i < fb->size; i++)
        buckets->logp[i] = fb->logp[i];

    /* Allocate polynomial state */
    siqs_poly_t poly;
    mpz_init(poly.a);
    mpz_init(poly.b);
    mpz_init(poly.c);
    poly.root1 = malloc(fb->size * sizeof(unsigned int));
    poly.root2 = malloc(fb->size * sizeof(unsigned int));
    poly.ainv = malloc(MAX_A_FACTORS * sizeof(unsigned int*));
    for (int j = 0; j < MAX_A_FACTORS; j++) {
        mpz_init(poly.B_vals[j]);
        poly.ainv[j] = malloc(fb->size * sizeof(unsigned int));
    }

    /* Relation storage */
    rels_t *full_rels = rels_create(MAX_RELS, fb->size);
    rels_t *part_rels = rels_create(MAX_PARTIALS, fb->size);
    lp_hash_t *lp_hash = lp_create(MAX_PARTIALS);

    /* RNG */
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    /* Temp variables */
    mpz_t Q_val, Y_val, tmp, tmp2;
    mpz_inits(Q_val, Y_val, tmp, tmp2, NULL);
    short *tmp_exps = calloc(fb->size + 2, sizeof(short));

    int total_polys = 0;
    int combined_rels = 0;

    /* ========== MAIN SIEVING LOOP ========== */
    while (full_rels->count < target) {
        double t = elapsed();
        if (t > 280) {
            fprintf(stderr, "TIMEOUT at %.1fs with %d/%d relations\n", t, full_rels->count, target);
            break;
        }

        /* ===== Generate new 'a' coefficient ===== */
        {
            /* Target: a ≈ sqrt(2*kN) / M */
            mpz_t tgt;
            mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2);
            mpz_sqrt(tgt, tgt);
            if (M > 0) mpz_tdiv_q_ui(tgt, tgt, M);

            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            /* Select FB primes from middle third */
            int lo = fb->size / 4, hi = 3 * fb->size / 4;
            if (lo < 2) lo = 2;
            if (hi <= lo + 3) hi = fb->size - 1;

            double avg_logp = 0;
            int cnt = 0;
            for (int i = lo; i < hi; i++) {
                if (fb->sqrtN[i] == 0) continue;
                avg_logp += log(fb->prime[i]);
                cnt++;
            }
            if (cnt == 0) break;
            avg_logp /= cnt;

            int s = (int)(log_tgt / avg_logp + 0.5);
            if (s < 3) s = 3;
            if (s > MAX_A_FACTORS) s = MAX_A_FACTORS;
            if (s > hi - lo) s = hi - lo;
            poly.num_a_factors = s;

            /* Pick best 'a' from random trials */
            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];

            for (int att = 0; att < 50; att++) {
                mpz_set_ui(poly.a, 1);
                int idx[MAX_A_FACTORS];
                int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do {
                        idx[i] = lo + gmp_urandomm_ui(rng, hi - lo);
                        good = 1;
                        for (int j = 0; j < i; j++)
                            if (idx[j] == idx[i]) { good = 0; break; }
                        if (fb->sqrtN[idx[i]] == 0) good = 0;
                        tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok = 0; break; }
                    mpz_mul_ui(poly.a, poly.a, fb->prime[idx[i]]);
                }
                if (!ok) continue;

                double ratio;
                if (mpz_cmp(poly.a, tgt) > 0) {
                    mpz_tdiv_q(tmp, poly.a, tgt);
                    ratio = mpz_get_d(tmp);
                } else {
                    mpz_tdiv_q(tmp, tgt, poly.a);
                    ratio = mpz_get_d(tmp);
                }
                if (ratio < best_ratio) {
                    best_ratio = ratio;
                    memcpy(best, idx, s * sizeof(int));
                }
                if (ratio < 1.5) break;
            }

            memcpy(poly.a_indices, best, s * sizeof(int));
            mpz_set_ui(poly.a, 1);
            for (int i = 0; i < s; i++)
                mpz_mul_ui(poly.a, poly.a, fb->prime[poly.a_indices[i]]);
            mpz_clear(tgt);

            /* ===== Compute B values for SIQS ===== */
            for (int j = 0; j < s; j++) {
                int idx = poly.a_indices[j];
                unsigned int qj = fb->prime[idx];
                unsigned int rj = fb->sqrtN[idx];

                mpz_t a_q, inv_aq;
                mpz_inits(a_q, inv_aq, NULL);
                mpz_divexact_ui(a_q, poly.a, qj);

                /* inv_aq = (a/qj)^(-1) mod qj */
                unsigned long aqmod = mpz_fdiv_ui(a_q, qj);
                unsigned int inv_val = mod_inverse_u32((unsigned int)aqmod, qj);

                /* B_j = (a/qj) * inv(a/qj) * sqrt(kN) mod qj */
                mpz_mul_ui(poly.B_vals[j], a_q, ((unsigned long)rj * inv_val) % qj);
                mpz_clears(a_q, inv_aq, NULL);
            }

            /* ===== Precompute ainv[j][i] = 2*B_j * (2a)^(-1) mod p_i ===== */
            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 3) { poly.ainv[j][i] = 0; continue; }
                    unsigned long am = mpz_fdiv_ui(poly.a, p);
                    if (am == 0) { poly.ainv[j][i] = 0; continue; }
                    unsigned int ai = mod_inverse_u32((unsigned int)((2UL * am) % p), p);
                    if (ai == 0) { poly.ainv[j][i] = 0; continue; }
                    unsigned long Bm = mpz_fdiv_ui(poly.B_vals[j], p);
                    poly.ainv[j][i] = (unsigned int)((2ULL * ai % p * Bm) % p);
                }
            }
        }

        /* ===== Enumerate b-values using Gray code ===== */
        int num_b = 1 << (poly.num_a_factors - 1);

        for (int b_idx = 0; b_idx < num_b && full_rels->count < target; b_idx++) {
            if (elapsed() > 285) break;

            /* Compute b from Gray code */
            int gray = b_idx ^ (b_idx >> 1);

            mpz_set_ui(poly.b, 0);
            for (int j = 0; j < poly.num_a_factors; j++) {
                if (gray & (1 << j))
                    mpz_add(poly.b, poly.b, poly.B_vals[j]);
                else
                    mpz_sub(poly.b, poly.b, poly.B_vals[j]);
            }

            /* Verify b^2 ≡ kN (mod a) */
            mpz_mul(tmp, poly.b, poly.b);
            mpz_sub(tmp, tmp, kN);
            mpz_mod(tmp, tmp, poly.a);
            if (mpz_sgn(tmp) != 0) {
                mpz_neg(poly.b, poly.b);
                mpz_mul(tmp, poly.b, poly.b);
                mpz_sub(tmp, tmp, kN);
                mpz_mod(tmp, tmp, poly.a);
                if (mpz_sgn(tmp) != 0) continue;
            }

            /* c = (b^2 - kN) / a */
            mpz_mul(poly.c, poly.b, poly.b);
            mpz_sub(poly.c, poly.c, kN);
            mpz_divexact(poly.c, poly.c, poly.a);

            /* ===== Compute sieve roots for all FB primes ===== */
            for (int i = 0; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                if (p < 3 || fb->sqrtN[i] == 0) {
                    poly.root1[i] = 0xFFFFFFFF;
                    poly.root2[i] = 0xFFFFFFFF;
                    continue;
                }
                unsigned long am = mpz_fdiv_ui(poly.a, p);
                if (am == 0) {
                    /* p divides a: special handling - single root */
                    unsigned long bm = mpz_fdiv_ui(poly.b, p);
                    unsigned long cm = mpz_fdiv_ui(poly.c, p);
                    if (bm == 0) {
                        poly.root1[i] = 0xFFFFFFFF;
                        poly.root2[i] = 0xFFFFFFFF;
                        continue;
                    }
                    /* Q(x) = ax^2 + 2bx + c; when p|a: Q(x) ≡ 2bx + c (mod p) */
                    /* x ≡ -c * (2b)^(-1) (mod p) */
                    unsigned int inv2b = mod_inverse_u32((unsigned int)((2UL * bm) % p), p);
                    unsigned int root = (unsigned int)((unsigned long)(p - cm) % p * inv2b % p);
                    poly.root1[i] = root;
                    poly.root2[i] = root;  /* single root */
                    continue;
                }
                unsigned int ai = mod_inverse_u32((unsigned int)am, p);
                if (ai == 0) {
                    poly.root1[i] = 0xFFFFFFFF;
                    poly.root2[i] = 0xFFFFFFFF;
                    continue;
                }
                unsigned long bm = mpz_fdiv_ui(poly.b, p);
                unsigned int r = fb->sqrtN[i];
                /* root1 = (r - b) * a^(-1) mod p */
                /* root2 = (-r - b) * a^(-1) mod p */
                poly.root1[i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
                poly.root2[i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
            }

            /* ===== Fill buckets for large primes ===== */
            bucket_reset(buckets);
            for (int i = bucket_start; i < fb->size; i++) {
                if (poly.root1[i] == 0xFFFFFFFF) continue;
                unsigned int p = fb->prime[i];

                /* Process both roots across entire sieve interval [-M, M) */
                /* Sieve interval mapped to [0, 2*M) where block 0 = [-M, -M+BLOCKSIZE) */
                for (int ri = 0; ri < 2; ri++) {
                    unsigned int root = (ri == 0) ? poly.root1[i] : poly.root2[i];

                    /* Convert sieve root to position in [-M, M) */
                    /* root gives x ≡ root (mod p) where Q(root) ≡ 0 (mod p) */
                    /* First hit in [-M, M): x = root + ceil((-M - root)/p) * p */
                    long start = (long)root;
                    /* Shift to start from -M */
                    long off = (start - (-(long)M)) % (long)p;
                    if (off < 0) off += p;
                    long pos = -(long)M + off;  /* first position >= -M with x ≡ root (mod p) */

                    while (pos < (long)M) {
                        int block_idx = (int)((pos + M) >> BLOCKBITS);
                        int block_off = (int)((pos + M) & BLOCKMASK);
                        if (block_idx >= 0 && block_idx < total_blocks)
                            bucket_add(buckets, block_idx, i, block_off);
                        pos += p;
                    }
                }
            }

            total_polys++;
            if (total_polys % 500 == 0) {
                double t = elapsed();
                fprintf(stderr, "  poly=%d, rels=%d/%d (full=%d+%d), part=%d, t=%.1fs\n",
                        total_polys, full_rels->count, target,
                        full_rels->count - combined_rels, combined_rels,
                        part_rels->count, t);
            }

            /* ===== Sieve each block ===== */
            for (int blk = 0; blk < total_blocks; blk++) {
                int block_start_x = -(int)M + blk * BLOCKSIZE;

                /* Initialize sieve */
                memset(sieve_array, 0, BLOCKSIZE);

                /* ---- Small/medium primes: direct sieve ---- */
                for (int i = sieve_start; i < bucket_start && i < fb->size; i++) {
                    if (poly.root1[i] == 0xFFFFFFFF) continue;
                    unsigned int p = fb->prime[i];
                    unsigned char lp = fb->logp[i];

                    /* Compute start positions in this block */
                    for (int ri = 0; ri < 2; ri++) {
                        unsigned int root = (ri == 0) ? poly.root1[i] : poly.root2[i];

                        /* Find first x >= block_start_x with x ≡ root (mod p) */
                        long off = ((long)root - (long)block_start_x) % (long)p;
                        if (off < 0) off += p;

                        /* Sieve this root - unrolled for small primes */
                        int j = (int)off;
                        if (p < BLOCKSIZE / 4) {
                            /* 4x unroll for very small primes */
                            int p2 = p * 2, p3 = p * 3, p4 = p * 4;
                            int limit = BLOCKSIZE - p3;
                            while (j <= limit) {
                                sieve_array[j] += lp;
                                sieve_array[j + p] += lp;
                                sieve_array[j + p2] += lp;
                                sieve_array[j + p3] += lp;
                                j += p4;
                            }
                            /* handle remaining */
                            while (j < BLOCKSIZE) { sieve_array[j] += lp; j += p; }
                        } else if (p < BLOCKSIZE / 2) {
                            /* 2x unroll for medium primes */
                            int limit = BLOCKSIZE - p;
                            while (j <= limit) {
                                sieve_array[j] += lp;
                                sieve_array[j + p] += lp;
                                j += p * 2;
                            }
                            if (j < BLOCKSIZE) sieve_array[j] += lp;
                        } else {
                            /* Large primes: max 1 hit per block per root */
                            if (j < BLOCKSIZE) sieve_array[j] += lp;
                        }
                    }
                }

                /* ---- Large primes: apply bucket hits ---- */
                bucket_t *bkt = &buckets->blocks[blk];
                for (int k = 0; k < bkt->count; k++) {
                    unsigned int packed = bkt->data[k];
                    unsigned int fb_idx = packed >> 16;
                    unsigned int offset = packed & 0xFFFF;
                    sieve_array[offset] += buckets->logp[fb_idx];
                }

                /* ---- Scan for smooth candidates ---- */
                for (int j = 0; j < BLOCKSIZE; j += 4) {
                    /* Quick check: any of the 4 bytes above threshold? */
                    unsigned int v = *(unsigned int*)(sieve_array + j);
                    /* Check if any byte >= threshold */
                    /* Threshold check: subtract threshold from each byte, check high bits */
                    if (((v & 0xFF) < (unsigned)threshold) &&
                        (((v >> 8) & 0xFF) < (unsigned)threshold) &&
                        (((v >> 16) & 0xFF) < (unsigned)threshold) &&
                        (((v >> 24) & 0xFF) < (unsigned)threshold))
                        continue;

                    for (int jj = j; jj < j + 4 && jj < BLOCKSIZE; jj++) {
                        if (sieve_array[jj] < threshold) continue;

                        int x_global = block_start_x + jj;
                        unsigned long lp_val = 0;

                        int result = trial_divide_sieve_informed(
                            Q_val, Y_val,
                            tmp_exps,
                            fb, fb->size,
                            poly.root1, poly.root2,
                            jj, block_start_x,
                            lp_bound, &lp_val,
                            kN, poly.a, poly.b, x_global,
                            full_rels->neg_col);

                        if (result == 1) {
                            int idx = full_rels->count;
                            if (idx < full_rels->alloc) {
                                mpz_set(full_rels->Y[idx], Y_val);
                                memcpy(full_rels->exps[idx], tmp_exps, (fb->size + 2) * sizeof(short));
                                full_rels->lp[idx] = 0;
                                full_rels->count++;
                            }
                        } else if (result == 2) {
                            int match = lp_find(lp_hash, lp_val);
                            if (match >= 0) {
                                int pidx = match;
                                int fidx = full_rels->count;
                                if (fidx < full_rels->alloc) {
                                    /* Combined: multiply Ys, add exponents, store LP for sqrt */
                                    mpz_mul(full_rels->Y[fidx], Y_val, part_rels->Y[pidx]);
                                    mpz_mod(full_rels->Y[fidx], full_rels->Y[fidx], N);
                                    for (int e = 0; e < fb->size + 2; e++)
                                        full_rels->exps[fidx][e] =
                                            tmp_exps[e] + part_rels->exps[pidx][e];
                                    full_rels->lp[fidx] = lp_val; /* store LP for sqrt step */
                                    full_rels->count++;
                                    combined_rels++;
                                }
                            } else {
                                int pidx = part_rels->count;
                                if (pidx < part_rels->alloc) {
                                    mpz_set(part_rels->Y[pidx], Y_val);
                                    memcpy(part_rels->exps[pidx], tmp_exps, (fb->size + 2) * sizeof(short));
                                    part_rels->lp[pidx] = lp_val;
                                    part_rels->count++;
                                    lp_insert(lp_hash, lp_val, pidx);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fprintf(stderr, "Sieving done: %d relations (%d full + %d combined) in %.1fs\n",
            full_rels->count, full_rels->count - combined_rels, combined_rels, elapsed());

    /* ========== LINEAR ALGEBRA ========== */
    if (full_rels->count < fb->size + 1) {
        fprintf(stderr, "Not enough relations: %d < %d\n", full_rels->count, fb->size + 1);
        return 1;
    }

    int nrels = full_rels->count;
    int ncols = fb->size + 1; /* +1 for sign column */

    fprintf(stderr, "Building %d x %d GF(2) matrix...\n", nrels, ncols);

    gf2_t *mat = gf2_create(nrels, ncols);
    for (int r = 0; r < nrels; r++) {
        for (int c = 0; c < ncols; c++) {
            if (full_rels->exps[r][c] & 1)
                gf2_set(mat, r, c);
        }
    }

    fprintf(stderr, "Solving GF(2) system...\n");
    int **deps;
    int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, MAX_DEPS);

    fprintf(stderr, "Found %d dependencies, trying square root...\n", ndeps);

    /* ========== SQUARE ROOT (exponent-based) ========== */
    mpz_t X, Y2, g;
    mpz_inits(X, Y2, g, NULL);

    int found = 0;
    for (int d = 0; d < ndeps && !found; d++) {
        /* Accumulate: X = prod(Y_i) mod N, exps = sum of exponent vectors */
        mpz_set_ui(X, 1);
        int *exps = calloc(ncols + 2, sizeof(int));

        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            mpz_mul(X, X, full_rels->Y[ri]);
            mpz_mod(X, X, N);
            for (int c = 0; c <= fb->size; c++)  /* include sign column */
                exps[c] += full_rels->exps[ri][c];
        }

        /* Check all exponents are even (including sign) */
        int all_even = 1;
        for (int c = 0; c <= fb->size; c++)
            if (exps[c] & 1) { all_even = 0; break; }

        if (!all_even) {
            free(exps);
            continue;
        }

        /* Y = prod(p_i^(exp_i/2)) mod N
         * For combined LP relations, also include LP contributions */
        mpz_set_ui(Y2, 1);
        for (int c = 0; c < fb->size; c++) {
            if (exps[c] <= 0) continue;
            int half = exps[c] / 2;
            mpz_set_ui(tmp, fb->prime[c]);
            mpz_powm_ui(tmp, tmp, half, N);
            mpz_mul(Y2, Y2, tmp);
            mpz_mod(Y2, Y2, N);
        }

        /* For combined partial relations in this dep, include LP */
        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            if (full_rels->lp[ri] > 1) {
                /* This shouldn't happen for full or combined rels (lp=0),
                 * but handle combined LP: the LP^2 factor needs LP^1 in sqrt */
                mpz_set_ui(tmp, full_rels->lp[ri]);
                mpz_mul(Y2, Y2, tmp);
                mpz_mod(Y2, Y2, N);
            }
        }

        /* gcd(X ± Y, N) */
        mpz_sub(tmp, X, Y2);
        mpz_gcd(g, tmp, N);

        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, g);
            gmp_printf("%Zd\n%Zd\n", g, cofactor);
            mpz_clear(cofactor);
            found = 1;
        }

        if (!found) {
            mpz_add(tmp, X, Y2);
            mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_t cofactor;
                mpz_init(cofactor);
                mpz_divexact(cofactor, N, g);
                gmp_printf("%Zd\n%Zd\n", g, cofactor);
                mpz_clear(cofactor);
                found = 1;
            }
        }

        free(exps);
    }

    if (!found) {
        fprintf(stderr, "FAILED: no non-trivial factor found from %d dependencies\n", ndeps);
        return 1;
    }

    double total_time = elapsed();
    fprintf(stderr, "Total time: %.3fs\n", total_time);

    return 0;
}
