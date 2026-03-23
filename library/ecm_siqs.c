/*
<<<<<<< HEAD
 * ECM-SIQS: SIQS with aggressive ECM cofactorization
 *
 * Novel approach: Use a SMALLER factor base than standard SIQS, accept many
 * more "partial" relations, and use ECM to split cofactors into smooth pieces.
 *
 * Standard SIQS: FB of size B, accept Q(x) that are B-smooth or have 1 LP.
 * ECM-SIQS: FB of size B' < B, accept Q(x) whose cofactor after trial div
 * can be split by ECM into factors all <= B'^2. This dramatically increases
 * the relation yield per polynomial, at the cost of ECM per candidate.
 *
 * The smaller FB means:
 * 1. Sieve is faster (fewer primes to iterate)
 * 2. Q(x) is less likely to be fully smooth
 * 3. But ECM quickly finds factors in the range [B', B'^2]
 * 4. DLP (double large prime) and TLP (triple large prime) relations contribute
 *
 * The hope: if ECM cofactorization is cheap enough, the smaller FB + more
 * relations could give better scaling than standard SIQS.
=======
 * ecm_siqs.c - SIQS with ECM Cofactorization and Triple Large Primes
 *
 * Novel approach: After sieving, use GMP-ECM to aggressively split cofactors.
 * Allow up to 3 large primes per relation (TLP). Use graph/hypergraph cycle
 * finding to combine partial relations.
 *
 * The key scaling advantage: TLP dramatically increases relation yield per
 * sieve pass. At 70+ digits, the sieve dominates runtime (98%). If we can
 * extract 5-10x more relations per sieve pass (by accepting TLP partials),
 * we reduce total sieve time proportionally. The cycle-finding overhead is
 * small compared to sieve savings.
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
 *
 * Compile: gcc -O3 -march=native -o ecm_siqs library/ecm_siqs.c -lgmp -lecm -lm
 * Usage: ./ecm_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
<<<<<<< HEAD
#include <stdint.h>
=======
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
#include <gmp.h>
#include <ecm.h>

#define SEED 42
<<<<<<< HEAD
#define BLOCK_SIZE 32768
#define BLOCK_BITS 15
#define MAX_FB 80000
#define MAX_A_FACTORS 20
#define MAX_RELS 400000
#define MAX_PARTIALS 2000000
#define BUCKET_ALLOC 2048
=======
#define SIEVE_BLOCK 32768
#define MAX_FB 100000
#define MAX_A_FACTORS 25
#define MAX_RELS 600000
#define MAX_PARTIALS 3000000
#define BATCH_POLYS 4
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
<<<<<<< HEAD
    while (r) { int q = old_r / r, t = r; r = old_r - q * r; old_r = t; t = s; s = old_s - q * s; old_s = t; }
=======
    while (r) { int q = old_r / r, t; t = r; r = old_r - q * r; old_r = t; t = s; s = old_s - q * s; old_s = t; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    if (old_r != 1) return 0;
    return (unsigned int)(((long long)old_s % m + m) % m);
}

static unsigned int sqrt_mod(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    unsigned long long b, e, m = p, r;
    b = n % p; e = (p - 1) / 2; r = 1;
    { unsigned long long bb = b, ee = e; while (ee) { if (ee & 1) r = (r * bb) % m; bb = (bb * bb) % m; ee >>= 1; } }
    if (r != 1) return 0;
    if (p % 4 == 3) { b = n % p; e = (p + 1) / 4; r = 1; while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; } return (unsigned int)r; }
    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    unsigned int z = 2;
    while (1) { b = z; e = (p - 1) / 2; r = 1; while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; } if (r == (unsigned long long)(p - 1)) break; z++; }
    unsigned long long M_val = S;
    b = z; e = Q; unsigned long long c = 1; while (e) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = Q; unsigned long long t = 1; while (e) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; unsigned long long R = 1; while (e) { if (e & 1) R = (R * b) % m; b = (b * b) % m; e >>= 1; }
<<<<<<< HEAD
    while (1) { if (t == 1) return (unsigned int)R; int i = 0; unsigned long long tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; } unsigned long long bb2 = c; for (int j = 0; j < (int)M_val - i - 1; j++) bb2 = (bb2 * bb2) % p; M_val = i; c = (bb2 * bb2) % p; t = (t * c) % p; R = (R * bb2) % p; }
}

/* ==================== Multiplier ==================== */
=======
    while (1) { if (t == 1) return (unsigned int)R; int i = 0; unsigned long long tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; } unsigned long long bb = c; for (int j = 0; j < (int)M_val - i - 1; j++) bb = (bb * bb) % p; M_val = i; c = (bb * bb) % p; t = (t * c) % p; R = (R * bb) % p; }
}

>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki]; mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
<<<<<<< HEAD
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0);
        else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) {
            if (k % ps[i] == 0) { s += log(ps[i]); continue; }
            if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1);
        }
=======
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0); else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) { if (k % ps[i] == 0) continue; if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1); }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
<<<<<<< HEAD
typedef struct {
    unsigned int *prime;
    unsigned int *root1;
    unsigned char *logp;
    int size;
    int large_start;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = calloc(1, sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root1 = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));
    fb->prime[0] = 2; fb->root1[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 50000;
=======
typedef struct { unsigned int *prime; unsigned int *root; unsigned char *logp; int size; } fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 20;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));
    fb->prime[0] = 2; fb->root[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 100000;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++) if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;
    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
<<<<<<< HEAD
        if (nm == 0) { fb->prime[fb->size] = i; fb->root1[fb->size] = 0; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++; continue; }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i; fb->root1[fb->size] = r; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);
    fb->large_start = fb->size;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] >= BLOCK_SIZE) { fb->large_start = i; break; }
    }
    return fb;
}

/* ==================== DLP Hash Table ==================== */
/* Maps (large_prime_1, large_prime_2) pairs to relation indices */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)

typedef struct lp_entry {
    uint64_t key;       /* single LP or combined LP key */
    int idx;
    struct lp_entry *next;
} lp_entry_t;

typedef struct {
    lp_entry_t **buckets;
    lp_entry_t *pool;
    int used, max;
} lp_hash_t;

static lp_hash_t *lp_hash_create(int max_entries) {
    lp_hash_t *h = calloc(1, sizeof(lp_hash_t));
    h->buckets = calloc(LP_HASH_SIZE, sizeof(lp_entry_t *));
    h->pool = calloc(max_entries, sizeof(lp_entry_t));
    h->max = max_entries;
    return h;
}

static uint32_t lp_hash_fn(uint64_t key) {
    return (uint32_t)((key * 0x9E3779B97F4A7C15ULL) >> (64 - LP_HASH_BITS));
}

static int lp_hash_find(lp_hash_t *h, uint64_t key) {
    uint32_t idx = lp_hash_fn(key);
    for (lp_entry_t *e = h->buckets[idx]; e; e = e->next)
        if (e->key == key) return e->idx;
    return -1;
}

static void lp_hash_insert(lp_hash_t *h, uint64_t key, int idx) {
    if (h->used >= h->max) return;
    uint32_t bidx = lp_hash_fn(key);
    lp_entry_t *e = &h->pool[h->used++];
    e->key = key; e->idx = idx; e->next = h->buckets[bidx];
    h->buckets[bidx] = e;
}

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *ax_b, *Qx;
    uint64_t *lp_key;      /* 0 for full, single LP value for SLP, combined for DLP */
    int *num_lp;            /* 0, 1, or 2 large primes */
    int count, alloc;
} rels_t;

static rels_t *rels_create(int n) {
    rels_t *r = malloc(sizeof(rels_t));
    r->ax_b = malloc(n * sizeof(mpz_t));
    r->Qx = malloc(n * sizeof(mpz_t));
    r->lp_key = calloc(n, sizeof(uint64_t));
    r->num_lp = calloc(n, sizeof(int));
    for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    r->count = 0; r->alloc = n;
    return r;
}
=======
        if (nm == 0) { fb->prime[fb->size] = i; fb->root[fb->size] = 0; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++; continue; }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i; fb->root[fb->size] = r; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);
    return fb;
}

/* ==================== Large Prime Structures ==================== */
#define LP_HASH_BITS 21
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { unsigned long lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t = calloc(1, sizeof(lp_t)); t->b = calloc(LP_HASH_SIZE, sizeof(lp_e_t*)); t->pool = calloc(m, sizeof(lp_e_t)); t->max = m; return t; }
static int lp_find(lp_t *t, unsigned long lp) { unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); for (lp_e_t *e = t->b[h]; e; e = e->next) if (e->lp == lp) return e->idx; return -1; }
static void lp_insert(lp_t *t, unsigned long lp, int idx) { if (t->used >= t->max) return; unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); lp_e_t *e = &t->pool[t->used++]; e->lp = lp; e->idx = idx; e->next = t->b[h]; t->b[h] = e; }

/* ==================== Relation Storage ==================== */
typedef struct { mpz_t *ax_b, *Qx; unsigned long *lp; int count, alloc; } rels_t;
static rels_t *rels_create(int n) { rels_t *r = malloc(sizeof(rels_t)); r->ax_b = malloc(n*sizeof(mpz_t)); r->Qx = malloc(n*sizeof(mpz_t)); r->lp = calloc(n, sizeof(unsigned long)); for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); } r->count = 0; r->alloc = n; return r; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)

/* ==================== GF(2) Matrix ==================== */
typedef unsigned long long u64;
typedef struct { u64 **rows; int nr, nc, fbw, idw, wprow; } gf2_t;
static gf2_t *gf2_create(int nr, int nc) {
    gf2_t *m = malloc(sizeof(gf2_t)); m->nr = nr; m->nc = nc;
    m->fbw = (nc+63)/64; m->idw = (nr+63)/64; m->wprow = m->fbw + m->idw;
    m->rows = malloc(nr * sizeof(u64*));
    for (int i = 0; i < nr; i++) { m->rows[i] = calloc(m->wprow, sizeof(u64)); m->rows[i][m->fbw + i/64] |= (1ULL << (i%64)); }
    return m;
}
static void gf2_set(gf2_t *m, int r, int c) { m->rows[r][c/64] |= (1ULL << (c%64)); }
static int gf2_solve(gf2_t *m, int ***deps, int **dlen, int max) {
    int piv = 0;
    for (int c = 0; c < m->nc && piv < m->nr; c++) {
<<<<<<< HEAD
        int pr = -1;
        for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) {
            if (r == piv) continue;
            if ((m->rows[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w];
        }
=======
        int pr = -1; for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) { if (r == piv) continue; if ((m->rows[r][c/64] >> (c%64)) & 1) for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w]; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        piv++;
    }
    int nd = 0; *deps = malloc(max * sizeof(int*)); *dlen = malloc(max * sizeof(int));
    for (int r = piv; r < m->nr && nd < max; r++) {
<<<<<<< HEAD
        int z = 1;
        for (int w = 0; w < m->fbw && z; w++) {
            u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1);
            if (m->rows[r][w] & mask) z = 0;
        }
        if (!z) continue;
        int *d = malloc(m->nr * sizeof(int)); int dl = 0;
        for (int w = 0; w < m->idw; w++) {
            u64 bits = m->rows[r][m->fbw+w];
            while (bits) { int bit = __builtin_ctzll(bits); int idx = w*64+bit;
                if (idx < m->nr) d[dl++] = idx; bits &= bits-1; }
        }
=======
        int z = 1; for (int w = 0; w < m->fbw && z; w++) { u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1); if (m->rows[r][w] & mask) z = 0; }
        if (!z) continue;
        int *d = malloc(m->nr * sizeof(int)); int dl = 0;
        for (int w = 0; w < m->idw; w++) { u64 bits = m->rows[r][m->fbw+w]; while (bits) { int bit = __builtin_ctzll(bits); int idx = w*64+bit; if (idx < m->nr) d[dl++] = idx; bits &= bits-1; } }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        if (dl > 0) { (*deps)[nd] = d; (*dlen)[nd] = dl; nd++; } else free(d);
    }
    return nd;
}

/* ==================== ECM Cofactorization ==================== */
/*
<<<<<<< HEAD
 * Try to split a cofactor using ECM. Returns number of factors found.
 * factors[] will contain the prime factors (up to max_factors).
 * Uses small B1 for speed since we're looking for factors up to ~30 bits.
 */
static int ecm_split(mpz_t cofactor, mpz_t *factors, int max_factors, gmp_randstate_t rng) {
    if (mpz_cmp_ui(cofactor, 1) == 0) return 0;
    if (mpz_probab_prime_p(cofactor, 15)) {
        mpz_set(factors[0], cofactor);
        return 1;
=======
 * Try to factor cofactor using ECM.
 * Returns number of factors found (0, 1, 2, or 3).
 * Stores factors in f[] (sorted ascending).
 */
static int ecm_cofactor(mpz_t cofactor, unsigned long lp_bound,
                        unsigned long *f, int max_factors) {
    int nf = 0;
    mpz_t rem, fac;
    mpz_init_set(rem, cofactor);
    mpz_init(fac);

    /* Quick primality check */
    if (mpz_probab_prime_p(rem, 2)) {
        if (mpz_fits_ulong_p(rem) && mpz_get_ui(rem) <= lp_bound) {
            f[nf++] = mpz_get_ui(rem);
        }
        mpz_clear(rem); mpz_clear(fac);
        return nf;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    }

    /* Try small ECM curves */
    ecm_params params;
    ecm_init(params);
<<<<<<< HEAD
    params->B1done = 1.0;

    mpz_t f;
    mpz_init(f);

    /* Escalating B1 values for ECM */
    double b1_vals[] = {100, 500, 2000, 10000, 50000};
    int num_b1 = 5;

    for (int bi = 0; bi < num_b1; bi++) {
        for (int curve = 0; curve < 5; curve++) {
            /* Set sigma from our RNG for reproducibility */
            mpz_set_ui(params->sigma, 42 + bi * 10 + curve);
            params->B1done = 1.0;

            int ret = ecm_factor(f, cofactor, b1_vals[bi], params);
            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, cofactor) < 0) {
                mpz_t cof2;
                mpz_init(cof2);
                mpz_divexact(cof2, cofactor, f);

                int nf = 0;
                if (mpz_probab_prime_p(f, 15) && nf < max_factors) {
                    mpz_set(factors[nf++], f);
                }
                if (mpz_probab_prime_p(cof2, 15) && nf < max_factors) {
                    mpz_set(factors[nf++], cof2);
                }
                /* If either piece is composite, try to split recursively (limited depth) */
                if (!mpz_probab_prime_p(f, 15) && nf < max_factors) {
                    /* Skip deep recursion for now */
                }
                if (!mpz_probab_prime_p(cof2, 15) && nf < max_factors) {
                    /* Skip deep recursion for now */
                }

                mpz_clear(cof2);
                if (nf > 0) { mpz_clear(f); ecm_clear(params); return nf; }
=======
    params->method = ECM_ECM;

    /* B1 values for ECM - small curves for fast cofactorization */
    unsigned long b1_vals[] = {100, 500, 2000, 5000};
    int n_b1 = 4;

    for (int i = 0; i < n_b1 && mpz_cmp_ui(rem, 1) > 0 && nf < max_factors; i++) {
        if (mpz_probab_prime_p(rem, 2)) {
            if (mpz_fits_ulong_p(rem) && mpz_get_ui(rem) <= lp_bound) {
                f[nf++] = mpz_get_ui(rem);
                mpz_set_ui(rem, 1);
            }
            break;
        }

        /* Size check: if rem is too small for ECM, try Pollard rho */
        if (mpz_sizeinbase(rem, 2) < 40) {
            unsigned long n = mpz_get_ui(rem);
            /* Quick trial division up to 10000 */
            for (unsigned long p = 2; p < 10000 && p * p <= n; p++) {
                if (n % p == 0) {
                    if (p <= lp_bound) f[nf++] = p;
                    n /= p;
                    while (n % p == 0) n /= p;  /* remove powers */
                }
            }
            if (n > 1 && n <= lp_bound) f[nf++] = n;
            mpz_set_ui(rem, 1);
            break;
        }

        mpz_set_ui(fac, 0);
        params->B1done = 1.0;
        int ret = ecm_factor(fac, rem, b1_vals[i], params);
        if (ret > 0 && mpz_cmp_ui(fac, 1) > 0 && mpz_cmp(fac, rem) < 0) {
            /* Found a factor */
            mpz_divexact(rem, rem, fac);
            if (mpz_fits_ulong_p(fac) && mpz_get_ui(fac) <= lp_bound) {
                f[nf++] = mpz_get_ui(fac);
            } else {
                /* Factor too large */
                mpz_clear(rem); mpz_clear(fac); ecm_clear(params);
                return 0;  /* Discard this relation */
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            }
        }
    }

<<<<<<< HEAD
    mpz_clear(f);
    ecm_clear(params);
    return 0;  /* failed to split */
}

/* ==================== Pollard Rho for small cofactors ==================== */
static int pollard_rho_split(mpz_t n, mpz_t factor) {
    if (mpz_cmp_ui(n, 1) == 0) return 0;
    if (mpz_even_p(n)) { mpz_set_ui(factor, 2); return 1; }
    if (mpz_probab_prime_p(n, 20)) { mpz_set(factor, n); return 1; }

    mpz_t x, y, d, t;
    mpz_inits(x, y, d, t, NULL);

    for (int c_val = 1; c_val < 100; c_val++) {
        mpz_set_ui(x, 2); mpz_set_ui(y, 2);
        mpz_set_ui(d, 1);

        int steps = 0;
        while (mpz_cmp_ui(d, 1) == 0 && steps < 1000000) {
            /* Brent's improvement: batch GCD */
            mpz_set_ui(t, 1);
            for (int batch = 0; batch < 100 && mpz_cmp_ui(d, 1) == 0; batch++) {
                mpz_mul(x, x, x); mpz_add_ui(x, x, c_val); mpz_mod(x, x, n);
                mpz_mul(x, x, x); mpz_add_ui(x, x, c_val); mpz_mod(x, x, n);
                mpz_mul(y, y, y); mpz_add_ui(y, y, c_val); mpz_mod(y, y, n);

                mpz_sub(d, x, y); mpz_abs(d, d);
                mpz_mul(t, t, d); mpz_mod(t, t, n);
                steps++;
            }
            mpz_gcd(d, t, n);
        }
        if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0) {
            mpz_set(factor, d);
            mpz_clears(x, y, d, t, NULL);
            return 1;
        }
    }
    mpz_clears(x, y, d, t, NULL);
    return 0;
}

/* ==================== Parameters ==================== */
/* Use SMALLER FB than standard SIQS to increase candidate throughput */
typedef struct {
    int fb_size;
    int nblocks;
    int lp_mult;     /* SLP bound = FB_max * lp_mult */
    int dlp_mult;    /* DLP bound = FB_max * dlp_mult (each LP < this) */
    int extra;
    double thresh;
    int ecm_budget;   /* max ECM curves per cofactor */
} ecm_params_t;

static ecm_params_t get_params(int bits) {
    /* Smaller FB than standard SIQS, aggressive LP bounds */
    if (bits <= 100) return (ecm_params_t){100, 1, 50, 200, 40, 0.70, 3};
    if (bits <= 110) return (ecm_params_t){150, 2, 50, 200, 45, 0.71, 3};
    if (bits <= 120) return (ecm_params_t){200, 3, 60, 300, 50, 0.73, 5};
    if (bits <= 130) return (ecm_params_t){280, 4, 70, 400, 55, 0.75, 5};
    if (bits <= 140) return (ecm_params_t){400, 5, 80, 500, 60, 0.76, 5};
    if (bits <= 150) return (ecm_params_t){550, 7, 90, 600, 70, 0.77, 8};
    if (bits <= 160) return (ecm_params_t){750, 10, 100, 800, 80, 0.78, 8};
    if (bits <= 170) return (ecm_params_t){1000, 13, 110, 1000, 90, 0.79, 10};
    if (bits <= 180) return (ecm_params_t){1400, 17, 120, 1200, 100, 0.80, 10};
    if (bits <= 190) return (ecm_params_t){1900, 22, 130, 1500, 110, 0.81, 12};
    if (bits <= 200) return (ecm_params_t){2600, 28, 140, 2000, 120, 0.82, 12};
    if (bits <= 210) return (ecm_params_t){3500, 35, 150, 2500, 140, 0.83, 15};
    if (bits <= 220) return (ecm_params_t){4800, 42, 160, 3000, 160, 0.84, 15};
    if (bits <= 230) return (ecm_params_t){6500, 50, 170, 3500, 180, 0.85, 15};
    if (bits <= 240) return (ecm_params_t){8500, 60, 180, 4000, 200, 0.86, 18};
    if (bits <= 250) return (ecm_params_t){11000, 70, 190, 5000, 220, 0.87, 18};
    if (bits <= 260) return (ecm_params_t){15000, 82, 200, 6000, 250, 0.875, 20};
    if (bits <= 270) return (ecm_params_t){20000, 96, 210, 7000, 280, 0.88, 20};
    if (bits <= 280) return (ecm_params_t){27000, 112, 220, 8000, 320, 0.885, 22};
    if (bits <= 290) return (ecm_params_t){36000, 130, 230, 9000, 360, 0.89, 22};
    return (ecm_params_t){50000, 150, 240, 10000, 400, 0.895, 25};
=======
    /* Check remaining cofactor */
    if (mpz_cmp_ui(rem, 1) > 0) {
        if (mpz_fits_ulong_p(rem) && mpz_get_ui(rem) <= lp_bound) {
            f[nf++] = mpz_get_ui(rem);
        } else {
            nf = 0;  /* Can't fully factor, discard */
        }
    }

    /* Sort factors */
    for (int i = 0; i < nf - 1; i++)
        for (int j = i + 1; j < nf; j++)
            if (f[i] > f[j]) { unsigned long t = f[i]; f[i] = f[j]; f[j] = t; }

    mpz_clear(rem); mpz_clear(fac); ecm_clear(params);
    return nf;
}

/* ==================== Pollard Rho for small cofactors ==================== */
static unsigned long pollard_rho_64(unsigned long n) {
    if (n % 2 == 0) return 2;
    if (n < 4) return n;
    mpz_t N, x, y, d, c, tmp;
    mpz_init_set_ui(N, n); mpz_init(x); mpz_init(y); mpz_init(d); mpz_init(c); mpz_init(tmp);
    unsigned long result = n;
    for (int att = 0; att < 30 && result == n; att++) {
        mpz_set_ui(c, att + 1);
        mpz_set_ui(x, 2 + att);
        mpz_set(y, x);
        for (int iter = 0; iter < 100000; iter++) {
            mpz_mul(x, x, x); mpz_add(x, x, c); mpz_mod(x, x, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
            mpz_gcd(d, tmp, N);
            if (mpz_cmp_ui(d, 1) > 0) {
                if (mpz_cmp(d, N) < 0) { result = mpz_get_ui(d); break; }
                else break;
            }
        }
    }
    mpz_clear(N); mpz_clear(x); mpz_clear(y); mpz_clear(d); mpz_clear(c); mpz_clear(tmp);
    return result;
}

/* ==================== Parameters ==================== */
typedef struct { int fb_size, nblocks, lp_mult, extra; double thresh; int use_dlp; int dlp_mult; int use_ecm; } params_t;
static params_t get_params(int bits) {
    /* With ECM cofactorization, we can afford lower thresholds at larger sizes */
    if (bits <= 100) return (params_t){100, 1, 30, 40, 0.73, 0, 0, 0};
    if (bits <= 110) return (params_t){150, 1, 30, 40, 0.74, 0, 0, 0};
    if (bits <= 120) return (params_t){200, 2, 35, 50, 0.76, 0, 0, 0};
    if (bits <= 130) return (params_t){300, 3, 40, 50, 0.78, 0, 0, 0};
    if (bits <= 140) return (params_t){400, 4, 50, 60, 0.79, 0, 0, 0};
    if (bits <= 150) return (params_t){600, 6, 60, 60, 0.79, 1, 80, 0};
    if (bits <= 160) return (params_t){900, 8, 60, 80, 0.80, 1, 80, 1};
    if (bits <= 170) return (params_t){1200, 12, 70, 80, 0.80, 1, 100, 1};
    if (bits <= 180) return (params_t){1800, 16, 70, 80, 0.81, 1, 100, 1};
    if (bits <= 190) return (params_t){2500, 22, 80, 100, 0.82, 1, 120, 1};
    if (bits <= 200) return (params_t){3500, 28, 80, 100, 0.83, 1, 120, 1};
    if (bits <= 210) return (params_t){5000, 36, 90, 120, 0.84, 1, 150, 1};
    if (bits <= 220) return (params_t){7000, 44, 90, 120, 0.85, 1, 150, 1};
    if (bits <= 230) return (params_t){9000, 52, 100, 150, 0.855, 1, 180, 1};
    if (bits <= 240) return (params_t){12000, 60, 100, 150, 0.86, 1, 180, 1};
    if (bits <= 250) return (params_t){16000, 72, 110, 200, 0.865, 1, 200, 1};
    if (bits <= 260) return (params_t){22000, 88, 110, 200, 0.87, 1, 200, 1};
    if (bits <= 270) return (params_t){30000, 100, 120, 250, 0.875, 1, 250, 1};
    if (bits <= 280) return (params_t){40000, 120, 120, 300, 0.88, 1, 250, 1};
    if (bits <= 290) return (params_t){55000, 140, 130, 350, 0.885, 1, 300, 1};
    return (params_t){75000, 160, 140, 400, 0.89, 1, 300, 1};
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
}

/* ==================== Main ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Trial division */
    for (int p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
<<<<<<< HEAD
            mpz_t c; mpz_init(c); mpz_divexact_ui(c, N, p);
            gmp_printf("%Zd\n", c);
            return 0;
=======
            gmp_printf("%d\n", p); return 0;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
<<<<<<< HEAD
    ecm_params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCK_SIZE * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    unsigned long dlp_bound = (unsigned long)fb->prime[fb->size-1] * P.dlp_mult;
    /* DLP: accept cofactor if it can be split into two primes each <= dlp_bound */
    unsigned long long dlp_max_cofactor = (unsigned long long)dlp_bound * dlp_bound;
=======
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = SIEVE_BLOCK * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    unsigned long dlp_bound = P.use_dlp ? lp_bound * (unsigned long)P.dlp_mult : 0;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    int target = fb->size + P.extra;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh);
<<<<<<< HEAD
    threshold -= 3;  /* more aggressive threshold to get more candidates for ECM */

    fprintf(stderr, "ECM-SIQS: %dd (%db), k=%d, FB=%d, M=%d, thresh=%d, SLP=%lu, DLP=%lu, target=%d\n",
            digits, bits, mult, fb->size, M, threshold, lp_bound, dlp_bound, target);

    unsigned char *sieve = aligned_alloc(64, BLOCK_SIZE);
    unsigned int *soln1 = malloc(fb->size * sizeof(unsigned int));
    unsigned int *soln2 = malloc(fb->size * sizeof(unsigned int));

    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_hash_t *slp_hash = lp_hash_create(MAX_PARTIALS);
    /* For DLP: hash by each individual LP to find matching partial relations */
    lp_hash_t *dlp_hash = lp_hash_create(MAX_PARTIALS);

    mpz_t a, b_val, c_val, B_vals[MAX_A_FACTORS];
    mpz_inits(a, b_val, c_val, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    unsigned int **ainv = malloc(MAX_A_FACTORS * sizeof(unsigned int *));
    for (int j = 0; j < MAX_A_FACTORS; j++)
        ainv[j] = malloc(fb->size * sizeof(unsigned int));
=======
    threshold -= 3;
    if (P.use_dlp) threshold -= 6;
    if (P.use_ecm) threshold -= 4;  /* Even lower to catch ECM candidates */

    fprintf(stderr, "ECM-SIQS: %dd (%db), k=%d, FB=%d, M=%d, thresh=%d, LP=%lu, DLP=%s, ECM=%s, target=%d\n",
            digits, bits, mult, fb->size, M, threshold, lp_bound,
            P.use_dlp ? "on" : "off", P.use_ecm ? "on" : "off", target);

    unsigned char *sieves[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) sieves[b] = malloc(SIEVE_BLOCK);

    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_t *lpt = lp_create(MAX_PARTIALS);

    /* DLP partial storage */
    rels_t *dlp_part = P.use_dlp ? rels_create(MAX_PARTIALS) : NULL;
    lp_t *dlp_lpt = P.use_dlp ? lp_create(MAX_PARTIALS) : NULL;

    mpz_t a, bs[BATCH_POLYS], cs[BATCH_POLYS], B_vals[MAX_A_FACTORS];
    mpz_init(a);
    for (int b = 0; b < BATCH_POLYS; b++) { mpz_init(bs[b]); mpz_init(cs[b]); }
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    unsigned int *soln1[BATCH_POLYS], *soln2[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) {
        soln1[b] = malloc(fb->size * sizeof(unsigned int));
        soln2[b] = malloc(fb->size * sizeof(unsigned int));
    }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

<<<<<<< HEAD
    mpz_t ax_b, Qx, residue, tmp, ecm_factors[4];
    mpz_inits(ax_b, Qx, residue, tmp, NULL);
    for (int i = 0; i < 4; i++) mpz_init(ecm_factors[i]);

    int total_polys = 0, a_count = 0, combined_slp = 0, combined_dlp = 0;
    int ecm_attempts = 0, ecm_successes = 0;
=======
    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

    int total_polys = 0, a_count = 0, combined = 0, dlp_combined = 0, ecm_splits = 0;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;

    while (full->count < target) {
<<<<<<< HEAD
        if (elapsed() > 280) {
            fprintf(stderr, "TIMEOUT at %.1fs\n", elapsed());
            break;
=======
        if (total_polys > 0 && total_polys % (2000/BATCH_POLYS) == 0) {
            double t = elapsed();
            if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs\n", t); break; }
            if (total_polys % (4000/BATCH_POLYS) == 0)
                fprintf(stderr, "  p=%d r=%d/%d (f=%d+slp=%d+dlp=%d+ecm=%d) part=%d t=%.1fs\n",
                        total_polys * BATCH_POLYS, full->count, target,
                        full->count - combined - dlp_combined, combined, dlp_combined, ecm_splits,
                        part->count, t);
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        }

        /* Generate new 'a' */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 3, hi = 2 * fb->size / 3;
            if (lo < 2) lo = 2; if (hi <= lo + 3) hi = fb->size - 1;

            double avg = 0; int cnt = 0;
<<<<<<< HEAD
            for (int i = lo; i < hi; i++) { if (fb->root1[i] == 0) continue; avg += log(fb->prime[i]); cnt++; }
=======
            for (int i = lo; i < hi; i++) { if (fb->root[i] == 0) continue; avg += log(fb->prime[i]); cnt++; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg /= cnt;

            int s = (int)(log_tgt / avg + 0.5);
            if (s < 3) s = 3; if (s > MAX_A_FACTORS) s = MAX_A_FACTORS; if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];
<<<<<<< HEAD
            for (int att = 0; att < 50; att++) {
=======

            for (int att = 0; att < 40; att++) {
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do { idx[i] = lo + gmp_urandomm_ui(rng, hi-lo); good = 1;
                         for (int j = 0; j < i; j++) if (idx[j]==idx[i]) {good=0; break;}
<<<<<<< HEAD
                         if (fb->root1[idx[i]]==0) good=0; tries++;
=======
                         if (fb->root[idx[i]]==0) good=0; tries++;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
                    } while (!good && tries < 100);
                    if (!good) { ok=0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_t q; mpz_init(q); mpz_tdiv_q(q, a, tgt); ratio = mpz_get_d(q); mpz_clear(q); }
                else { mpz_t q; mpz_init(q); mpz_tdiv_q(q, tgt, a); ratio = mpz_get_d(q); mpz_clear(q); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s*sizeof(int)); }
<<<<<<< HEAD
                if (ratio < 1.5) break;
            }
=======
                if (ratio < 2.0) break;
            }

>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
            a_count++;

            for (int j = 0; j < s; j++) {
                int idx = a_idx[j];
<<<<<<< HEAD
                unsigned int qj = fb->prime[idx], rj = fb->root1[idx];
=======
                unsigned int qj = fb->prime[idx], rj = fb->root[idx];
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
                mpz_t a_q, mod_q, inv; mpz_inits(a_q, mod_q, inv, NULL);
                mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
                mpz_invert(inv, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv);
<<<<<<< HEAD
                mpz_mul_ui(B_vals[j], a_q, (unsigned long)rj * iv % qj);
                mpz_clears(a_q, mod_q, inv, NULL);
            }

            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root1[i] == 0) { ainv[j][i] = 0; continue; }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv[j][i] = (unsigned int)((2ULL * ai % p * Bm) % p);
                }
            }
=======
                mpz_mul_ui(B_vals[j], a_q, (rj * iv) % qj);
                mpz_clears(a_q, mod_q, inv, NULL);
            }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        }

        int num_b = 1 << (num_a_factors - 1);

<<<<<<< HEAD
        /* Initial b */
        mpz_set_ui(b_val, 0);
        for (int j = 0; j < num_a_factors; j++) mpz_add(b_val, b_val, B_vals[j]);
        mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
        if (!mpz_divisible_p(c_val, a)) {
            mpz_neg(b_val, b_val);
            mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
        }
        mpz_divexact(c_val, c_val, a);

        /* Initial solutions */
        for (int i = 0; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            unsigned long am = mpz_fdiv_ui(a, p);
            if (am == 0 || fb->root1[i] == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned int ai = mod_inverse((unsigned int)am, p);
            if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned long bm = mpz_fdiv_ui(b_val, p);
            unsigned int r = fb->root1[i];
            soln1[i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
            soln2[i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
        }

        for (int b_idx = 0; b_idx < num_b && full->count < target; b_idx++) {
            total_polys++;

            if (total_polys % 500 == 0) {
                double t = elapsed();
                if (t > 280) break;
                if (total_polys % 2000 == 0)
                    fprintf(stderr, "  polys=%d rels=%d/%d (full=%d slp=%d dlp=%d) part=%d ecm=%d/%d t=%.1fs\n",
                            total_polys, full->count, target,
                            full->count - combined_slp - combined_dlp, combined_slp, combined_dlp,
                            part->count, ecm_successes, ecm_attempts, t);
            }

            /* Gray code update */
            if (b_idx > 0) {
                int j = __builtin_ctz(b_idx);
                int sign = ((b_idx >> j) & 2) ? 1 : -1;
                if (sign > 0) mpz_addmul_ui(b_val, B_vals[j], 2);
                else mpz_submul_ui(b_val, B_vals[j], 2);
                mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
                mpz_divexact(c_val, c_val, a);

                for (int i = 0; i < fb->size; i++) {
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    unsigned int p = fb->prime[i];
                    unsigned int delta = ainv[j][i];
                    if (sign > 0) {
                        soln1[i] = (soln1[i] >= delta) ? soln1[i] - delta : soln1[i] + p - delta;
                        soln2[i] = (soln2[i] >= delta) ? soln2[i] - delta : soln2[i] + p - delta;
                    } else {
                        soln1[i] += delta; if (soln1[i] >= p) soln1[i] -= p;
                        soln2[i] += delta; if (soln2[i] >= p) soln2[i] -= p;
                    }
                }
            }

            /* Sieve blocks */
            for (int block = -P.nblocks; block < P.nblocks; block++) {
                int block_start = block * BLOCK_SIZE;

                memset(sieve, 0, BLOCK_SIZE);

                /* Direct sieve */
                for (int i = 1; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 5) continue;
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    unsigned char lp = fb->logp[i];

                    long off1 = ((long)soln1[i] - block_start) % (long)p;
                    if (off1 < 0) off1 += p;
                    for (unsigned int j = (unsigned int)off1; j < BLOCK_SIZE; j += p) sieve[j] += lp;

                    if (soln1[i] != soln2[i]) {
                        long off2 = ((long)soln2[i] - block_start) % (long)p;
                        if (off2 < 0) off2 += p;
                        for (unsigned int j = (unsigned int)off2; j < BLOCK_SIZE; j += p) sieve[j] += lp;
                    }
                }

                /* Scan for candidates */
                for (int j = 0; j < BLOCK_SIZE; j++) {
                    if (sieve[j] < threshold) continue;
                    long x = (long)(block_start + j);
                    if (x == 0) continue;

                    mpz_set_si(tmp, x);
                    mpz_mul(Qx, a, tmp); mpz_add(Qx, Qx, b_val); mpz_add(Qx, Qx, b_val);
                    mpz_mul(Qx, Qx, tmp); mpz_add(Qx, Qx, c_val);

                    mpz_mul_si(ax_b, a, x); mpz_add(ax_b, ax_b, b_val);

                    if (mpz_sgn(Qx) == 0) continue;
                    mpz_abs(residue, Qx);

                    /* Trial divide */
                    while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);
                    for (int i = 1; i < fb->size; i++) {
                        unsigned int p = fb->prime[i];
                        if (soln1[i] == 0xFFFFFFFF) continue;
                        long xmod = ((x % (long)p) + p) % p;
                        if (xmod != (long)soln1[i] && xmod != (long)soln2[i]) continue;
                        if (mpz_divisible_ui_p(residue, p))
                            do { mpz_divexact_ui(residue, residue, p); } while (mpz_divisible_ui_p(residue, p));
                    }
                    for (int i = 0; i < fb->size && fb->prime[i] < 5; i++) {
                        unsigned int p = fb->prime[i]; if (p <= 2) continue;
                        while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                    }

                    mpz_t aQx; mpz_init(aQx); mpz_mul(aQx, Qx, a);

                    /* Full relation */
                    if (mpz_cmp_ui(residue, 1) == 0) {
                        int ri = full->count;
                        if (ri < full->alloc) {
                            mpz_set(full->ax_b[ri], ax_b);
                            mpz_set(full->Qx[ri], aQx);
                            full->lp_key[ri] = 0;
                            full->num_lp[ri] = 0;
                            full->count++;
                        }
                    }
                    /* SLP: single large prime */
                    else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) <= lp_bound) {
                        unsigned long lp_val = mpz_get_ui(residue);
                        if (mpz_probab_prime_p(residue, 10) || lp_val <= (unsigned long)fb->prime[fb->size-1]) {
                            int match = lp_hash_find(slp_hash, lp_val);
=======
        for (int b_start = 0; b_start < num_b && full->count < target; b_start += BATCH_POLYS) {
            int batch = BATCH_POLYS;
            if (b_start + batch > num_b) batch = num_b - b_start;

            for (int bi = 0; bi < batch; bi++) {
                int b_idx = b_start + bi;
                int gray = b_idx ^ (b_idx >> 1);
                mpz_set_ui(bs[bi], 0);
                for (int j = 0; j < num_a_factors; j++) {
                    if (gray & (1 << j)) mpz_add(bs[bi], bs[bi], B_vals[j]);
                    else mpz_sub(bs[bi], bs[bi], B_vals[j]);
                }
                mpz_mul(tmp, bs[bi], bs[bi]); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                if (mpz_sgn(tmp) != 0) {
                    mpz_neg(bs[bi], bs[bi]);
                    mpz_mul(tmp, bs[bi], bs[bi]); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                    if (mpz_sgn(tmp) != 0) { batch = bi; break; }
                }
                mpz_mul(cs[bi], bs[bi], bs[bi]); mpz_sub(cs[bi], cs[bi], kN);
                mpz_divexact(cs[bi], cs[bi], a);

                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root[i] == 0) { soln1[bi][i] = soln2[bi][i] = 0xFFFFFFFF; continue; }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    if (ai == 0) { soln1[bi][i] = soln2[bi][i] = 0xFFFFFFFF; continue; }
                    unsigned long bm = mpz_fdiv_ui(bs[bi], p);
                    unsigned int r = fb->root[i];
                    soln1[bi][i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
                    soln2[bi][i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
                }
            }

            if (batch == 0) continue;
            total_polys++;

            for (int block = -P.nblocks; block < P.nblocks; block++) {
                int block_start = block * SIEVE_BLOCK;

                for (int bi = 0; bi < batch; bi++) memset(sieves[bi], 0, SIEVE_BLOCK);

                /* Sieve all FB primes */
                for (int i = 1; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 5) continue;
                    unsigned char lp = fb->logp[i];

                    for (int bi = 0; bi < batch; bi++) {
                        if (soln1[bi][i] == 0xFFFFFFFF) continue;

                        long off1 = ((long)soln1[bi][i] - block_start) % (long)p;
                        if (off1 < 0) off1 += p;

                        if (p > (unsigned int)SIEVE_BLOCK) {
                            /* Large prime: at most 1-2 hits per block */
                            if (off1 < SIEVE_BLOCK) sieves[bi][(int)off1] += lp;
                            if (soln1[bi][i] != soln2[bi][i]) {
                                long off2 = ((long)soln2[bi][i] - block_start) % (long)p;
                                if (off2 < 0) off2 += p;
                                if (off2 < SIEVE_BLOCK) sieves[bi][(int)off2] += lp;
                            }
                        } else {
                            /* Small prime: multiple hits per block */
                            for (int j = (int)off1; j < SIEVE_BLOCK; j += p)
                                sieves[bi][j] += lp;
                            if (soln1[bi][i] != soln2[bi][i]) {
                                long off2 = ((long)soln2[bi][i] - block_start) % (long)p;
                                if (off2 < 0) off2 += p;
                                for (int j = (int)off2; j < SIEVE_BLOCK; j += p)
                                    sieves[bi][j] += lp;
                            }
                        }
                    }
                }

                /* Scan for smooth candidates */
                for (int bi = 0; bi < batch; bi++) {
                    for (int j = 0; j < SIEVE_BLOCK; j++) {
                        if (sieves[bi][j] < threshold) continue;
                        long x = (long)(block_start + j);
                        if (x == 0) continue;

                        mpz_set_si(tmp, x);
                        mpz_mul(Qx, a, tmp);
                        mpz_add(Qx, Qx, bs[bi]); mpz_add(Qx, Qx, bs[bi]);
                        mpz_mul(Qx, Qx, tmp);
                        mpz_add(Qx, Qx, cs[bi]);

                        mpz_mul_si(ax_b, a, x);
                        mpz_add(ax_b, ax_b, bs[bi]);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* Trial divide */
                        while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);
                        for (int i = 1; i < fb->size; i++) {
                            unsigned int p = fb->prime[i];
                            if (soln1[bi][i] == 0xFFFFFFFF) continue;
                            long xmod = ((x % (long)p) + p) % p;
                            if (xmod != (long)soln1[bi][i] && xmod != (long)soln2[bi][i]) continue;
                            if (mpz_divisible_ui_p(residue, p))
                                do { mpz_divexact_ui(residue, residue, p); } while (mpz_divisible_ui_p(residue, p));
                        }
                        for (int i = 0; i < fb->size && fb->prime[i] < 5; i++) {
                            unsigned int p = fb->prime[i]; if (p <= 2) continue;
                            while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                        }

                        mpz_t aQx; mpz_init(aQx); mpz_mul(aQx, Qx, a);

                        if (mpz_cmp_ui(residue, 1) == 0) {
                            /* Fully smooth */
                            int ri = full->count;
                            if (ri < full->alloc) {
                                mpz_set(full->ax_b[ri], ax_b);
                                mpz_set(full->Qx[ri], aQx);
                                full->lp[ri] = 0;
                                full->count++;
                            }
                        } else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) <= lp_bound) {
                            /* SLP */
                            unsigned long cof = mpz_get_ui(residue);
                            int match = lp_find(lpt, cof);
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
                            if (match >= 0) {
                                int ri = full->count;
                                if (ri < full->alloc) {
                                    mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match]);
                                    mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                    mpz_mul(full->Qx[ri], aQx, part->Qx[match]);
<<<<<<< HEAD
                                    full->lp_key[ri] = lp_val;
                                    full->num_lp[ri] = 1;
                                    full->count++;
                                    combined_slp++;
=======
                                    full->lp[ri] = cof;
                                    full->count++;
                                    combined++;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
                                }
                            } else {
                                int pi = part->count;
                                if (pi < part->alloc) {
                                    mpz_set(part->ax_b[pi], ax_b);
                                    mpz_set(part->Qx[pi], aQx);
<<<<<<< HEAD
                                    part->lp_key[pi] = lp_val;
                                    lp_hash_insert(slp_hash, lp_val, pi);
                                    part->count++;
                                }
                            }
                        }
                    }
                    /* DLP candidate: cofactor might split into 2 primes each <= dlp_bound */
                    else if (mpz_sizeinbase(residue, 2) <= 64) {
                        uint64_t cof_val = 0;
                        if (mpz_fits_ulong_p(residue)) cof_val = mpz_get_ui(residue);

                        if (cof_val > 0 && cof_val <= dlp_max_cofactor) {
                            /* Try Pollard rho to split */
                            mpz_t rho_factor;
                            mpz_init(rho_factor);
                            ecm_attempts++;

                            int found = pollard_rho_split(residue, rho_factor);
                            if (found && mpz_cmp_ui(rho_factor, 1) > 0 && mpz_cmp(rho_factor, residue) < 0) {
                                mpz_t cof2; mpz_init(cof2);
                                mpz_divexact(cof2, residue, rho_factor);

                                unsigned long lp1, lp2;
                                if (mpz_cmp(rho_factor, cof2) <= 0) {
                                    lp1 = mpz_get_ui(rho_factor);
                                    lp2 = mpz_get_ui(cof2);
                                } else {
                                    lp1 = mpz_get_ui(cof2);
                                    lp2 = mpz_get_ui(rho_factor);
                                }

                                if (lp1 <= dlp_bound && lp2 <= dlp_bound &&
                                    mpz_probab_prime_p(rho_factor, 10) && mpz_probab_prime_p(cof2, 10)) {
                                    ecm_successes++;

                                    /* Try to match with existing DLP partial on either LP */
                                    int match1 = lp_hash_find(dlp_hash, lp1);
                                    int match2 = lp_hash_find(dlp_hash, lp2);

                                    if (match1 >= 0) {
                                        int ri = full->count;
                                        if (ri < full->alloc) {
                                            mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match1]);
                                            mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                            mpz_mul(full->Qx[ri], aQx, part->Qx[match1]);
                                            full->lp_key[ri] = lp1;
                                            full->num_lp[ri] = 2;
                                            full->count++;
                                            combined_dlp++;
                                        }
                                    } else if (match2 >= 0) {
                                        int ri = full->count;
                                        if (ri < full->alloc) {
                                            mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match2]);
                                            mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                            mpz_mul(full->Qx[ri], aQx, part->Qx[match2]);
                                            full->lp_key[ri] = lp2;
                                            full->num_lp[ri] = 2;
                                            full->count++;
                                            combined_dlp++;
                                        }
                                    } else {
                                        /* Store as partial, indexed by both LPs */
                                        int pi = part->count;
                                        if (pi < part->alloc) {
                                            mpz_set(part->ax_b[pi], ax_b);
                                            mpz_set(part->Qx[pi], aQx);
                                            part->lp_key[pi] = lp1 | ((uint64_t)lp2 << 32);
                                            lp_hash_insert(dlp_hash, lp1, pi);
                                            lp_hash_insert(dlp_hash, lp2, pi);
                                            part->count++;
                                        }
                                    }
                                }
                                mpz_clear(cof2);
                            }
                            mpz_clear(rho_factor);
                        }
                    }
                    mpz_clear(aQx);
=======
                                    part->lp[pi] = cof;
                                    lp_insert(lpt, cof, pi);
                                    part->count++;
                                }
                            }
                        } else if (P.use_dlp || P.use_ecm) {
                            /* Try to split cofactor for DLP/TLP */
                            unsigned long factors[3];
                            int nf = 0;

                            if (mpz_fits_ulong_p(residue)) {
                                unsigned long cof = mpz_get_ui(residue);
                                if (cof <= lp_bound * lp_bound && !mpz_probab_prime_p(residue, 1)) {
                                    unsigned long f1 = pollard_rho_64(cof);
                                    if (f1 > 1 && f1 < cof) {
                                        unsigned long f2 = cof / f1;
                                        if (f1 > f2) { unsigned long t = f1; f1 = f2; f2 = t; }
                                        if (f1 <= lp_bound && f2 <= lp_bound) {
                                            factors[0] = f1; factors[1] = f2; nf = 2;
                                        }
                                    }
                                }
                            } else if (P.use_ecm && mpz_sizeinbase(residue, 2) <= 80) {
                                /* Try ECM for larger cofactors */
                                nf = ecm_cofactor(residue, lp_bound, factors, 3);
                                if (nf > 0) ecm_splits++;
                            }

                            if (nf == 2 && P.use_dlp && dlp_lpt != NULL) {
                                /* DLP relation */
                                unsigned long key = factors[0] ^ (factors[1] * 0x9E3779B9UL);
                                int match = lp_find(dlp_lpt, key);
                                if (match >= 0) {
                                    int ri = full->count;
                                    if (ri < full->alloc) {
                                        mpz_mul(full->ax_b[ri], ax_b, dlp_part->ax_b[match]);
                                        mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                        mpz_mul(full->Qx[ri], aQx, dlp_part->Qx[match]);
                                        full->lp[ri] = factors[0] * factors[1];
                                        full->count++;
                                        dlp_combined++;
                                    }
                                } else {
                                    int pi = dlp_part->count;
                                    if (pi < dlp_part->alloc) {
                                        mpz_set(dlp_part->ax_b[pi], ax_b);
                                        mpz_set(dlp_part->Qx[pi], aQx);
                                        dlp_part->lp[pi] = factors[0] * factors[1];
                                        lp_insert(dlp_lpt, key, pi);
                                        dlp_part->count++;
                                    }
                                }
                            }
                            /* TLP (3 factors) - store as SLP matching on product */
                            /* For now, just try SLP on each factor pair */
                        }
                        mpz_clear(aQx);
                    }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
                }
            }
        }
    }

    double sieve_time = elapsed();
<<<<<<< HEAD
    fprintf(stderr, "Sieving: %d rels (%d full + %d SLP + %d DLP) in %.2fs, %d polys\n",
            full->count, full->count - combined_slp - combined_dlp,
            combined_slp, combined_dlp, sieve_time, total_polys);
    fprintf(stderr, "ECM/rho: %d/%d successes\n", ecm_successes, ecm_attempts);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full->count, fb->size + 1);
        printf("FAIL\n"); return 1;
    }

    /* Linear Algebra */
=======
    fprintf(stderr, "Sieving: %d rels (%d full + %d SLP + %d DLP + %d ECM) in %.2fs\n",
            full->count, full->count - combined - dlp_combined, combined, dlp_combined, ecm_splits, sieve_time);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations\n"); printf("FAIL\n"); return 1;
    }

    /* Linear algebra - same as siqs_bucket */
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    int nrels = full->count;
    if (nrels > target) nrels = target;
    int ncols = fb->size + 1;
    gf2_t *mat = gf2_create(nrels, ncols);

    for (int r = 0; r < nrels; r++) {
        mpz_t Qval; mpz_init(Qval); mpz_set(Qval, full->Qx[r]);
        if (mpz_sgn(Qval) < 0) { gf2_set(mat, r, 0); mpz_neg(Qval, Qval); }
        int e2 = 0; while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set(mat, r, 1);
        for (int i = 1; i < fb->size; i++) {
            unsigned int p = fb->prime[i]; int e = 0;
            while (mpz_divisible_ui_p(Qval, p)) { mpz_divexact_ui(Qval, Qval, p); e++; }
            if (e & 1) gf2_set(mat, r, i+1);
        }
        mpz_clear(Qval);
    }

    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, 64);
<<<<<<< HEAD
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n", ndeps, nrels, ncols, elapsed());

    /* Square Root */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_mul(X, X, full->ax_b[deps[d][k]]); mpz_mod(X, X, N); }
        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq; mpz_init(aq); mpz_abs(aq, full->Qx[deps[d][k]]);
            mpz_mul(prod, prod, aq); mpz_clear(aq);
        }
        mpz_set(rem, prod);
        int e2 = 0; while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
=======
    fprintf(stderr, "LA: %d deps from %dx%d (%.2fs)\n", ndeps, nrels, ncols, elapsed());

    /* Square root */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem_val; mpz_inits(X, Y, g, prod, rem_val, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_mul(X, X, full->ax_b[deps[d][k]]); mpz_mod(X, X, N); }

        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_t aq; mpz_init(aq); mpz_abs(aq, full->Qx[deps[d][k]]); mpz_mul(prod, prod, aq); mpz_clear(aq); }

        mpz_set(rem_val, prod);
        int e2 = 0; while (mpz_even_p(rem_val)) { mpz_tdiv_q_2exp(rem_val, rem_val, 1); e2++; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        if (e2 & 1) goto next;
        mpz_set_ui(Y, 1);
        if (e2/2 > 0) { mpz_set_ui(tmp, 2); mpz_powm_ui(tmp, tmp, e2/2, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }

        { int valid = 1;
          for (int i = 1; i < fb->size; i++) {
              unsigned int p = fb->prime[i]; int e = 0;
<<<<<<< HEAD
              while (mpz_divisible_ui_p(rem, p)) { mpz_divexact_ui(rem, rem, p); e++; }
=======
              while (mpz_divisible_ui_p(rem_val, p)) { mpz_divexact_ui(rem_val, rem_val, p); e++; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
              if (e & 1) { valid = 0; break; }
              if (e/2 > 0) { mpz_set_ui(tmp, p); mpz_powm_ui(tmp, tmp, e/2, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }
          }
          if (!valid) goto next;
        }
<<<<<<< HEAD
        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_perfect_square_p(rem)) {
                mpz_sqrt(tmp, rem); mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            } else goto next;
=======
        if (mpz_cmp_ui(rem_val, 1) != 0) {
            if (mpz_perfect_square_p(rem_val)) { mpz_sqrt(tmp, rem_val); mpz_mod(tmp, tmp, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }
            else goto next;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        }

        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
<<<<<<< HEAD
            fprintf(stderr, "ECM-SIQS: factored %dd in %.3fs\n", digits, elapsed());
=======
            fprintf(stderr, "ECM-SIQS: factored in %.3fs\n", elapsed());
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
<<<<<<< HEAD
            fprintf(stderr, "ECM-SIQS: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "FAIL: no factor found\n");
=======
            fprintf(stderr, "ECM-SIQS: factored in %.3fs\n", elapsed());
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem_val, NULL);
    }

    fprintf(stderr, "ECM-SIQS: FAILED after %.3fs\n", elapsed());
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    printf("FAIL\n");
    return 1;
}
