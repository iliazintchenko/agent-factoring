/*
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
 *
 * Compile: gcc -O3 -march=native -o ecm_siqs library/ecm_siqs.c -lgmp -lecm -lm
 * Usage: ./ecm_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>
#include <ecm.h>

#define SEED 42
#define BLOCK_SIZE 32768
#define BLOCK_BITS 15
#define MAX_FB 80000
#define MAX_A_FACTORS 20
#define MAX_RELS 400000
#define MAX_PARTIALS 2000000
#define BUCKET_ALLOC 2048

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
    while (r) { int q = old_r / r, t = r; r = old_r - q * r; old_r = t; t = s; s = old_s - q * s; old_s = t; }
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
    while (1) { if (t == 1) return (unsigned int)R; int i = 0; unsigned long long tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; } unsigned long long bb2 = c; for (int j = 0; j < (int)M_val - i - 1; j++) bb2 = (bb2 * bb2) % p; M_val = i; c = (bb2 * bb2) % p; t = (t * c) % p; R = (R * bb2) % p; }
}

/* ==================== Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki]; mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0);
        else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) {
            if (k % ps[i] == 0) { s += log(ps[i]); continue; }
            if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
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
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++) if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;
    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
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
        int pr = -1;
        for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) {
            if (r == piv) continue;
            if ((m->rows[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w];
        }
        piv++;
    }
    int nd = 0; *deps = malloc(max * sizeof(int*)); *dlen = malloc(max * sizeof(int));
    for (int r = piv; r < m->nr && nd < max; r++) {
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
        if (dl > 0) { (*deps)[nd] = d; (*dlen)[nd] = dl; nd++; } else free(d);
    }
    return nd;
}

/* ==================== ECM Cofactorization ==================== */
/*
 * Try to split a cofactor using ECM. Returns number of factors found.
 * factors[] will contain the prime factors (up to max_factors).
 * Uses small B1 for speed since we're looking for factors up to ~30 bits.
 */
static int ecm_split(mpz_t cofactor, mpz_t *factors, int max_factors, gmp_randstate_t rng) {
    if (mpz_cmp_ui(cofactor, 1) == 0) return 0;
    if (mpz_probab_prime_p(cofactor, 15)) {
        mpz_set(factors[0], cofactor);
        return 1;
    }

    /* Try small ECM curves */
    ecm_params params;
    ecm_init(params);
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
            }
        }
    }

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
            mpz_t c; mpz_init(c); mpz_divexact_ui(c, N, p);
            gmp_printf("%Zd\n", c);
            return 0;
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    ecm_params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCK_SIZE * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    unsigned long dlp_bound = (unsigned long)fb->prime[fb->size-1] * P.dlp_mult;
    /* DLP: accept cofactor if it can be split into two primes each <= dlp_bound */
    unsigned long long dlp_max_cofactor = (unsigned long long)dlp_bound * dlp_bound;
    int target = fb->size + P.extra;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh);
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

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp, ecm_factors[4];
    mpz_inits(ax_b, Qx, residue, tmp, NULL);
    for (int i = 0; i < 4; i++) mpz_init(ecm_factors[i]);

    int total_polys = 0, a_count = 0, combined_slp = 0, combined_dlp = 0;
    int ecm_attempts = 0, ecm_successes = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;

    while (full->count < target) {
        if (elapsed() > 280) {
            fprintf(stderr, "TIMEOUT at %.1fs\n", elapsed());
            break;
        }

        /* Generate new 'a' */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 3, hi = 2 * fb->size / 3;
            if (lo < 2) lo = 2; if (hi <= lo + 3) hi = fb->size - 1;

            double avg = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) { if (fb->root1[i] == 0) continue; avg += log(fb->prime[i]); cnt++; }
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg /= cnt;

            int s = (int)(log_tgt / avg + 0.5);
            if (s < 3) s = 3; if (s > MAX_A_FACTORS) s = MAX_A_FACTORS; if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];
            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do { idx[i] = lo + gmp_urandomm_ui(rng, hi-lo); good = 1;
                         for (int j = 0; j < i; j++) if (idx[j]==idx[i]) {good=0; break;}
                         if (fb->root1[idx[i]]==0) good=0; tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok=0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_t q; mpz_init(q); mpz_tdiv_q(q, a, tgt); ratio = mpz_get_d(q); mpz_clear(q); }
                else { mpz_t q; mpz_init(q); mpz_tdiv_q(q, tgt, a); ratio = mpz_get_d(q); mpz_clear(q); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s*sizeof(int)); }
                if (ratio < 1.5) break;
            }
            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
            a_count++;

            for (int j = 0; j < s; j++) {
                int idx = a_idx[j];
                unsigned int qj = fb->prime[idx], rj = fb->root1[idx];
                mpz_t a_q, mod_q, inv; mpz_inits(a_q, mod_q, inv, NULL);
                mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
                mpz_invert(inv, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv);
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
        }

        int num_b = 1 << (num_a_factors - 1);

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
                            if (match >= 0) {
                                int ri = full->count;
                                if (ri < full->alloc) {
                                    mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match]);
                                    mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                    mpz_mul(full->Qx[ri], aQx, part->Qx[match]);
                                    full->lp_key[ri] = lp_val;
                                    full->num_lp[ri] = 1;
                                    full->count++;
                                    combined_slp++;
                                }
                            } else {
                                int pi = part->count;
                                if (pi < part->alloc) {
                                    mpz_set(part->ax_b[pi], ax_b);
                                    mpz_set(part->Qx[pi], aQx);
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
                }
            }
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving: %d rels (%d full + %d SLP + %d DLP) in %.2fs, %d polys\n",
            full->count, full->count - combined_slp - combined_dlp,
            combined_slp, combined_dlp, sieve_time, total_polys);
    fprintf(stderr, "ECM/rho: %d/%d successes\n", ecm_successes, ecm_attempts);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full->count, fb->size + 1);
        printf("FAIL\n"); return 1;
    }

    /* Linear Algebra */
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
        if (e2 & 1) goto next;
        mpz_set_ui(Y, 1);
        if (e2/2 > 0) { mpz_set_ui(tmp, 2); mpz_powm_ui(tmp, tmp, e2/2, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }

        { int valid = 1;
          for (int i = 1; i < fb->size; i++) {
              unsigned int p = fb->prime[i]; int e = 0;
              while (mpz_divisible_ui_p(rem, p)) { mpz_divexact_ui(rem, rem, p); e++; }
              if (e & 1) { valid = 0; break; }
              if (e/2 > 0) { mpz_set_ui(tmp, p); mpz_powm_ui(tmp, tmp, e/2, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }
          }
          if (!valid) goto next;
        }
        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_perfect_square_p(rem)) {
                mpz_sqrt(tmp, rem); mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            } else goto next;
        }

        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "ECM-SIQS: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "ECM-SIQS: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "FAIL: no factor found\n");
    printf("FAIL\n");
    return 1;
}
