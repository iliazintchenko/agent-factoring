/*
 * Turbo TLP - Triple Large Prime SIQS
 *
 * Key features (extends turbo_siqs with TLP):
 * 1. 48KB L1-cache-optimized sieve blocks (AMD EPYC 9R45)
 * 2. Bucket sieving for large primes
 * 3. Gray code self-initialization with incremental solution updates
 * 4. Double Large Primes (DLP) with graph-based cycle finding (union-find)
 * 5. Triple Large Primes (TLP) - stored as DLP with lp1=product of two smallest
 * 6. Two-threshold sieve scanning (tight for full/SLP, loose for DLP)
 * 7. 64-bit fast path trial division
 * 8. Sieve-informed trial division (only test matching roots)
 * 9. Lower sieve threshold when DLP is enabled (dlp_bonus)
 * 10. factor64: complete factorization of 64-bit numbers for TLP detection
 *
 * Compile: gcc -O3 -march=native -o turbo_tlp library/turbo_tlp.c -lgmp -lm
 * Usage: ./turbo_tlp <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>
#ifdef __AVX512BW__
#include <immintrin.h>
#endif

/* ==================== Constants ==================== */
#define SEED 42
#define BLOCK_SIZE 49152       /* 48KB = L1d cache */
#define MAX_FB 50000
#define MAX_A_FACTORS 25
#define MAX_FULL_RELS 500000
#define MAX_PARTIAL_RELS 4000000
#define BUCKET_ALLOC 65536

/* ==================== Timing ==================== */
static struct timespec g_start;
static double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static inline uint32_t mod_inv32(uint32_t a, uint32_t m) {
    int64_t old_r = a, r = m, old_s = 1, s = 0;
    while (r) {
        int64_t q = old_r / r;
        int64_t t = r; r = old_r - q * r; old_r = t;
        t = s; s = old_s - q * s; old_s = t;
    }
    return old_r == 1 ? (uint32_t)(((old_s % (int64_t)m) + m) % m) : 0;
}

static uint32_t sqrt_mod_p(uint32_t n, uint32_t p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    uint64_t mod = p;

    /* Euler criterion */
    uint64_t b = n % p, e = (p - 1) / 2, r = 1;
    while (e) { if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
    if (r != 1) return 0;

    if ((p & 3) == 3) {
        b = n % p; e = (p + 1) / 4; r = 1;
        while (e) { if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
        return (uint32_t)r;
    }

    /* Tonelli-Shanks */
    uint32_t Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }
    uint32_t z = 2;
    while (1) {
        b = z; e = (p - 1) / 2; r = 1;
        while (e) { if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
        if (r == (uint64_t)(p - 1)) break;
        z++;
    }
    uint64_t M = S;
    b = z; e = Q; uint64_t c = 1;
    while (e) { if (e & 1) c = c * b % mod; b = b * b % mod; e >>= 1; }
    b = n % p; e = Q; uint64_t t = 1;
    while (e) { if (e & 1) t = t * b % mod; b = b * b % mod; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; uint64_t R = 1;
    while (e) { if (e & 1) R = R * b % mod; b = b * b % mod; e >>= 1; }
    while (1) {
        if (t == 1) return (uint32_t)R;
        uint64_t i = 0, tt = t;
        while (tt != 1) { tt = tt * tt % p; i++; }
        uint64_t bb = c;
        for (uint64_t j = 0; j < M - i - 1; j++) bb = bb * bb % p;
        M = i; c = bb * bb % p; t = t * c % p; R = R * bb % p;
    }
}

/* ==================== Factor Base ==================== */
typedef struct {
    uint32_t *prime, *root;
    uint8_t *logp;
    int size;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = calloc(1, sizeof(fb_t));
    fb->prime = malloc((target + 100) * sizeof(uint32_t));
    fb->root  = malloc((target + 100) * sizeof(uint32_t));
    fb->logp  = malloc((target + 100) * sizeof(uint8_t));
    fb->prime[0] = 2; fb->root[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 100000;
    char *sieve = calloc(bound + 1, 1);
    for (int i = 2; (long)i * i <= bound; i++)
        if (!sieve[i]) for (int j = i * i; j <= bound; j += i) sieve[j] = 1;
    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sieve[i]) continue;
        uint32_t p = (uint32_t)i;
        unsigned long nm = mpz_fdiv_ui(kN, p);
        if (nm == 0) { fb->prime[fb->size] = p; fb->root[fb->size] = 0; fb->logp[fb->size] = (uint8_t)(log2(p)+0.5); fb->size++; continue; }
        uint32_t r = sqrt_mod_p((uint32_t)nm, p);
        if (!r) continue;
        fb->prime[fb->size] = p; fb->root[fb->size] = r; fb->logp[fb->size] = (uint8_t)(log2(p)+0.5); fb->size++;
    }
    free(sieve);
    return fb;
}

/* ==================== Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki]; mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0); else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
        for (int i = 0; i < 18; i++) {
            if (k % ps[i] == 0) { s += log((double)ps[i])/ps[i]; continue; }
            if (sqrt_mod_p(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1);
        }
        if (s > best) { best = s; best_k = k; } mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Bucket Sieve ==================== */
typedef struct { uint32_t pos; uint8_t logp; } bucket_entry_t;
typedef struct { bucket_entry_t *entries; int count, alloc; } bucket_t;

static void bucket_add(bucket_t *b, uint32_t pos, uint8_t logp) {
    if (b->count >= b->alloc) { b->alloc *= 2; b->entries = realloc(b->entries, b->alloc * sizeof(bucket_entry_t)); }
    b->entries[b->count++] = (bucket_entry_t){pos, logp};
}

/* ==================== DLP Graph (Union-Find) ==================== */
#define DLP_HASH_BITS 22
#define DLP_HASH_SIZE (1 << DLP_HASH_BITS)

typedef struct dlp_node {
    uint64_t lp;          /* large prime value */
    int parent;           /* union-find parent (index into node array) */
    int rank;
    int edge_count;       /* number of edges incident to this node */
    struct dlp_node *hash_next;
} dlp_node_t;

typedef struct {
    dlp_node_t *nodes;
    int num_nodes;
    int max_nodes;
    dlp_node_t **hash_table;
    int num_edges;        /* total edges (relations) */
    int num_trees;        /* number of connected components */
    int num_cycles;       /* num_edges - (num_nodes - num_trees) = free edges = cycles */
} dlp_graph_t;

static dlp_graph_t *dlp_create(int max_nodes) {
    dlp_graph_t *g = calloc(1, sizeof(dlp_graph_t));
    g->nodes = calloc(max_nodes, sizeof(dlp_node_t));
    g->max_nodes = max_nodes;
    g->hash_table = calloc(DLP_HASH_SIZE, sizeof(dlp_node_t*));
    return g;
}

static int dlp_find_or_create(dlp_graph_t *g, uint64_t lp) {
    uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - DLP_HASH_BITS));
    for (dlp_node_t *n = g->hash_table[h]; n; n = n->hash_next)
        if (n->lp == lp) return (int)(n - g->nodes);

    if (g->num_nodes >= g->max_nodes) return -1;
    int idx = g->num_nodes++;
    g->nodes[idx].lp = lp;
    g->nodes[idx].parent = idx;
    g->nodes[idx].rank = 0;
    g->nodes[idx].edge_count = 0;
    g->nodes[idx].hash_next = g->hash_table[h];
    g->hash_table[h] = &g->nodes[idx];
    g->num_trees++;
    return idx;
}

static int dlp_find_root(dlp_graph_t *g, int x) {
    while (g->nodes[x].parent != x) {
        g->nodes[x].parent = g->nodes[g->nodes[x].parent].parent; /* path compression */
        x = g->nodes[x].parent;
    }
    return x;
}

/* Returns 1 if this edge creates a cycle (= usable relation), 0 if it's a tree edge */
static int dlp_add_edge(dlp_graph_t *g, uint64_t lp1, uint64_t lp2) {
    int n1 = dlp_find_or_create(g, lp1);
    int n2 = dlp_find_or_create(g, lp2);
    if (n1 < 0 || n2 < 0) return 0;

    g->nodes[n1].edge_count++;
    g->nodes[n2].edge_count++;
    g->num_edges++;

    int r1 = dlp_find_root(g, n1);
    int r2 = dlp_find_root(g, n2);
    if (r1 == r2) {
        g->num_cycles++;
        return 1; /* cycle! */
    }
    /* Union by rank */
    if (g->nodes[r1].rank < g->nodes[r2].rank) { int t = r1; r1 = r2; r2 = t; }
    g->nodes[r2].parent = r1;
    if (g->nodes[r1].rank == g->nodes[r2].rank) g->nodes[r1].rank++;
    g->num_trees--;
    return 0;
}

/* ==================== LP Hash for SLP matching ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE2 (1 << LP_HASH_BITS)
typedef struct lp_e { uint64_t lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_hash_t;
static lp_hash_t *lp_create(int m) { lp_hash_t *t = calloc(1, sizeof(lp_hash_t)); t->b = calloc(LP_HASH_SIZE2, sizeof(lp_e_t*)); t->pool = calloc(m, sizeof(lp_e_t)); t->max = m; return t; }
static int lp_find(lp_hash_t *t, uint64_t lp) { uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); for (lp_e_t *e = t->b[h]; e; e = e->next) if (e->lp == lp) return e->idx; return -1; }
static void lp_insert(lp_hash_t *t, uint64_t lp, int idx) { if (t->used >= t->max) return; uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); lp_e_t *e = &t->pool[t->used++]; e->lp = lp; e->idx = idx; e->next = t->b[h]; t->b[h] = e; }

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *ax_b, *Qx;
    uint64_t *lp1, *lp2;  /* 0 = full, nonzero = partial */
    int count, alloc;
} rels_t;

static rels_t *rels_create(int n) {
    rels_t *r = calloc(1, sizeof(rels_t));
    r->alloc = n; r->ax_b = malloc(n*sizeof(mpz_t)); r->Qx = malloc(n*sizeof(mpz_t));
    r->lp1 = calloc(n, sizeof(uint64_t)); r->lp2 = calloc(n, sizeof(uint64_t));
    for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    return r;
}

static int rels_add(rels_t *r, mpz_t ax_b, mpz_t Qx, uint64_t lp1, uint64_t lp2) {
    if (r->count >= r->alloc) return -1;
    int i = r->count++;
    mpz_set(r->ax_b[i], ax_b); mpz_set(r->Qx[i], Qx);
    r->lp1[i] = lp1; r->lp2[i] = lp2;
    return i;
}

/* ==================== Linear Algebra Headers ==================== */
/* Save our BLOCK_SIZE and restore it after headers (headers redefine BLOCK_SIZE=64) */
#pragma push_macro("BLOCK_SIZE")
#include "block_lanczos.h"
#include "structured_gauss.h"
#pragma pop_macro("BLOCK_SIZE")

typedef uint64_t gf2w;
/* sparse_matrix_t is provided by block_lanczos.h; use it as sparse_t alias */
typedef sparse_matrix_t sparse_t;

/* Inline LA functions removed: using block_lanczos.h / structured_gauss.h instead */



/* ==================== Pollard Rho (64-bit cofactor splitting) ==================== */
static int is_prime64(uint64_t n) {
    if (n < 2) return 0; if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    uint64_t d = n - 1; int r = 0;
    while (!(d & 1)) { d >>= 1; r++; }
    uint64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int wi = 0; wi < 12; wi++) {
        uint64_t a = witnesses[wi];
        if (a >= n) continue;
        __uint128_t x = 1, base = a; uint64_t exp = d;
        while (exp) { if (exp & 1) x = x * base % n; base = base * base % n; exp >>= 1; }
        uint64_t xv = (uint64_t)x;
        if (xv == 1 || xv == n - 1) continue;
        int found = 0;
        for (int i = 0; i < r - 1; i++) { xv = (uint64_t)((__uint128_t)xv * xv % n); if (xv == n - 1) { found = 1; break; } }
        if (!found) return 0;
    }
    return 1;
}

static int split64(uint64_t n, uint64_t *f1, uint64_t *f2) {
    if (n < 4) return 0;
    if (n % 2 == 0) { *f1 = 2; *f2 = n/2; return 1; }
    if (is_prime64(n)) return 0;
    for (uint64_t c = 1; c < 200; c++) {
        uint64_t x = 2, y = 2, p = 1;
        /* Brent's improved rho with product accumulation */
        uint64_t ys = 2, q = 1;
        int m = 128;
        do {
            x = y;
            for (int i = 0; i < m; i++) y = (uint64_t)((__uint128_t)y * y % n) + c;
            int k = 0;
            do {
                ys = y;
                int lim = m - k < 32 ? m - k : 32;
                for (int i = 0; i < lim; i++) {
                    y = (uint64_t)((__uint128_t)y * y % n) + c;
                    uint64_t diff = x > y ? x - y : y - x;
                    q = (uint64_t)((__uint128_t)q * diff % n);
                }
                /* GCD(q, n) */
                { uint64_t a = q, b = n; while (b) { uint64_t t = b; b = a % b; a = t; } p = a; }
                k += 32;
            } while (p == 1 && k < m);
            m *= 2;
        } while (p == 1);
        if (p == n) {
            /* Backtrack */
            do {
                ys = (uint64_t)((__uint128_t)ys * ys % n) + c;
                uint64_t diff = x > ys ? x - ys : ys - x;
                uint64_t a = diff, b = n;
                while (b) { uint64_t t = b; b = a % b; a = t; }
                p = a;
            } while (p == 1);
        }
        if (p != n && p != 1) { *f1 = p; *f2 = n / p; return 1; }
    }
    return 0;
}

/* ==================== factor64: complete factorization of 64-bit number ==================== */
/* Factors n into prime factors, returns count. factors[] must have room for at least 64 entries. */
static int factor64(uint64_t n, uint64_t *factors) {
    if (n <= 1) return 0;
    int count = 0;
    /* Small prime trial division up to 1000 */
    static const uint32_t small_primes[] = {
        2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
        73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,
        179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,
        283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,
        419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,
        547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,
        661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
        811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,
        947,953,967,971,977,983,991,997,0
    };
    for (int i = 0; small_primes[i] && n > 1; i++) {
        uint32_t p = small_primes[i];
        if ((uint64_t)p * p > n) break;
        while (n % p == 0) { factors[count++] = p; n /= p; }
    }
    if (n == 1) return count;

    /* n is now > 1000 and has no small prime factors */
    /* Use a stack-based approach with split64 for composite detection */
    uint64_t stack[64];
    int ssize = 0;
    stack[ssize++] = n;

    while (ssize > 0) {
        uint64_t v = stack[--ssize];
        if (v == 1) continue;
        if (is_prime64(v)) {
            factors[count++] = v;
            continue;
        }
        uint64_t f1, f2;
        if (!split64(v, &f1, &f2)) {
            /* Failed to factor - treat as prime (shouldn't happen for small enough inputs) */
            factors[count++] = v;
            continue;
        }
        stack[ssize++] = f1;
        stack[ssize++] = f2;
    }
    return count;
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size, nblocks, lp_mult, extra;
    double thresh_adj;
    int dlp;
} params_t;

static params_t get_params(int bits) {
    /* Tuned for single-core 300s with DLP graph matching + TLP support.
     * For 235-258 bits: smaller FB + higher LP multiplier to increase DLP/TLP yield. */
    if (bits <= 100) return (params_t){120,   1,  40,  30, 0.73, 0};
    if (bits <= 110) return (params_t){170,   1,  50,  35, 0.74, 0};
    if (bits <= 120) return (params_t){240,   2,  60,  45, 0.76, 1};
    if (bits <= 130) return (params_t){340,   3,  70,  50, 0.77, 1};
    if (bits <= 140) return (params_t){480,   4,  80,  60, 0.78, 1};
    if (bits <= 150) return (params_t){650,   5,  90,  70, 0.79, 1};
    if (bits <= 160) return (params_t){900,   7, 100,  80, 0.80, 1};
    if (bits <= 170) return (params_t){1200, 10, 110,  90, 0.81, 1};
    if (bits <= 180) return (params_t){1700, 14, 120, 100, 0.82, 1};
    if (bits <= 190) return (params_t){2300, 18, 140, 120, 0.825, 1};
    if (bits <= 200) return (params_t){3000, 24, 160, 140, 0.83, 1};
    if (bits <= 210) return (params_t){4000, 30, 180, 160, 0.835, 1};
    if (bits <= 220) return (params_t){5500,  30, 200, 180, 0.84, 1};
    if (bits <= 230) return (params_t){7000,  38, 200, 200, 0.845, 1};
    /* TLP-tuned params: smaller FB + higher LP multiplier for 235-258 bits */
    if (bits <= 235) return (params_t){10000, 46, 200, 180, 0.855, 1};
    if (bits <= 245) return (params_t){10000, 56, 300, 200, 0.86,  1};
    if (bits <= 253) return (params_t){14000, 60, 250, 250, 0.86,  1};
    if (bits <= 256) return (params_t){16000, 64, 200, 280, 0.858, 1};
    if (bits <= 260) return (params_t){18000, 68, 180, 300, 0.858, 1};
    if (bits <= 265) return (params_t){16000, 80, 250, 300, 0.87,  1};
    if (bits <= 270) return (params_t){22000, 80, 200, 350, 0.865, 1};
    if (bits <= 280) return (params_t){30000, 96, 200, 400, 0.87, 1};
    return (params_t){40000, 120, 200, 500, 0.875, 1};
}

/* ==================== MAIN ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N, 10);
    int bits   = (int)mpz_sizeinbase(N, 2);

    /* Trial division */
    for (unsigned long p = 2; p < 100000; p++) {
        if (p > 2 && p % 2 == 0) continue;
        if (mpz_divisible_ui_p(N, p)) { printf("%lu\n", p); return 0; }
    }
    { mpz_t sq; mpz_init(sq); if (mpz_perfect_square_p(N)) { mpz_sqrt(sq, N); gmp_printf("%Zd\n", sq); return 0; } mpz_clear(sq); }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);
    fb_t *fb = fb_create(kN, P.fb_size);

    int M = BLOCK_SIZE * P.nblocks;
    uint64_t lp_bound = (uint64_t)fb->prime[fb->size-1] * P.lp_mult;
    uint64_t dlp_max = (uint64_t)lp_bound * lp_bound; /* max cofactor for DLP */
    /* TLP max: lp_bound^3, but we cap to avoid overflow. Use __uint128_t check. */
    /* We'll check tlp inline using 128-bit arithmetic */
    int target = fb->size + P.extra;

    /* Threshold: lower when DLP is enabled to catch more LP candidates */
    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2((double)M);
    int dlp_bonus = P.dlp ? 5 : 0; /* moderate DLP threshold bonus */
    int threshold = (int)(log2_Qmax * P.thresh_adj) - 3 - dlp_bonus;
    if (threshold < 20) threshold = 20;

    /* Bucket sieve setup */
    int bucket_thresh = 0;
    for (int i = 0; i < fb->size; i++) { if (fb->prime[i] > BLOCK_SIZE) { bucket_thresh = i; break; } }
    if (bucket_thresh == 0) bucket_thresh = fb->size;

    int total_blocks = 2 * P.nblocks;
    bucket_t *buckets = calloc(total_blocks, sizeof(bucket_t));
    for (int i = 0; i < total_blocks; i++) { buckets[i].alloc = BUCKET_ALLOC; buckets[i].entries = malloc(BUCKET_ALLOC * sizeof(bucket_entry_t)); }

    fprintf(stderr, "TurboTLP: %dd (%db), k=%d, FB=%d, M=%d, thresh=%d (dlp_bonus=%d), LP=%lu%s, target=%d\n",
            digits, bits, mult, fb->size, M, threshold, dlp_bonus, lp_bound, P.dlp?" [DLP+TLP]":"", target);

    uint8_t *sieve = malloc(BLOCK_SIZE);
    rels_t *full = rels_create(MAX_FULL_RELS);
    rels_t *part = rels_create(MAX_PARTIAL_RELS);
    lp_hash_t *slp = lp_create(MAX_PARTIAL_RELS);
    dlp_graph_t *dlp_g = P.dlp ? dlp_create(2000000) : NULL;

    /* DLP relation tracking: for each partial with 2 LPs, store index */
    #define DLP_PAIR_HASH_BITS 22
    #define DLP_PAIR_HASH_SIZE (1 << DLP_PAIR_HASH_BITS)
    typedef struct dp_e { uint64_t lp; int rel_idx; struct dp_e *next; } dp_e_t;
    dp_e_t **dp_hash = calloc(DLP_PAIR_HASH_SIZE, sizeof(dp_e_t*));
    dp_e_t *dp_pool = calloc(MAX_PARTIAL_RELS * 2, sizeof(dp_e_t));
    int dp_pool_used = 0;

    mpz_t a, b_val, c_val, B_vals[MAX_A_FACTORS];
    mpz_inits(a, b_val, c_val, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    uint32_t *soln1 = malloc(fb->size * sizeof(uint32_t));
    uint32_t *soln2 = malloc(fb->size * sizeof(uint32_t));
    uint32_t *ainv_data = malloc((size_t)MAX_A_FACTORS * fb->size * sizeof(uint32_t));
    #define AINV(j,i) ainv_data[(j)*fb->size+(i)]

    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng, SEED);
    mpz_t ax_b, Qx, aQx, residue, tmp, tmp2;
    mpz_inits(ax_b, Qx, aQx, residue, tmp, tmp2, NULL);

    int total_polys = 0, combined_slp = 0, combined_dlp = 0, combined_tlp = 0;
    int dlp_found = 0, dlp_cycles = 0, tlp_found = 0;
    int a_idx[MAX_A_FACTORS], num_a_factors = 0;

    /* ==================== Main Sieving Loop ==================== */
    while (full->count < target) {
        double t = elapsed_sec();
        if (t > 288.0) { fprintf(stderr, "TIMEOUT at %.1fs with %d/%d rels\n", t, full->count, target); break; }

        /* === Generate new 'a' === */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 4, hi = 3 * fb->size / 4;
            if (lo < 2) lo = 2; if (hi <= lo + 5) hi = fb->size - 1;

            double avg_logp = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) { if (fb->root[i] == 0) continue; avg_logp += log(fb->prime[i]); cnt++; }
            if (cnt == 0) break;
            avg_logp /= cnt;

            int s = (int)(log_tgt / avg_logp + 0.5);
            if (s < 3) s = 3; if (s > MAX_A_FACTORS) s = MAX_A_FACTORS; if (s > hi-lo) s = hi-lo;
            num_a_factors = s;

            double best_ratio = 1e30; int best_idx[MAX_A_FACTORS];
            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1); int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do { idx[i] = lo + gmp_urandomm_ui(rng, hi-lo); good = 1;
                         for (int j = 0; j < i; j++) if (idx[j]==idx[i]) {good=0;break;}
                         if (fb->root[idx[i]]==0) good=0; tries++;
                    } while (!good && tries < 100);
                    if (!good) {ok=0;break;} mpz_mul_ui(a, a, fb->prime[idx[i]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_tdiv_q(tmp, a, tgt); ratio = mpz_get_d(tmp); }
                else { mpz_tdiv_q(tmp, tgt, a); ratio = mpz_get_d(tmp); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best_idx, idx, s*sizeof(int)); }
                if (ratio < 1.5) break;
            }
            memcpy(a_idx, best_idx, s*sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
        }

        /* === B values === */
        for (int j = 0; j < num_a_factors; j++) {
            int idx = a_idx[j]; uint32_t qj = fb->prime[idx], rj = fb->root[idx];
            mpz_t a_q, mod_q, inv; mpz_inits(a_q, mod_q, inv, NULL);
            mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
            mpz_invert(inv, a_q, mod_q);
            mpz_mul_ui(B_vals[j], a_q, ((unsigned long)rj * mpz_get_ui(inv)) % qj);
            mpz_clears(a_q, mod_q, inv, NULL);
        }

        /* === Precompute ainv for Gray code updates === */
        for (int j = 0; j < num_a_factors; j++) {
            for (int i = 0; i < fb->size; i++) {
                uint32_t p = fb->prime[i];
                unsigned long am = mpz_fdiv_ui(a, p);
                if (am == 0 || fb->root[i] == 0) { AINV(j,i) = 0; continue; }
                uint32_t ai = mod_inv32((uint32_t)am, p);
                unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                AINV(j,i) = (uint32_t)((uint64_t)2 * ai % p * Bm % p);
            }
        }

        /* === Gray code enumeration === */
        int num_b = 1 << (num_a_factors - 1);
        int prev_gray = 0;

        /* Initialize b for gray=0 */
        mpz_set_ui(b_val, 0);
        for (int j = 0; j < num_a_factors; j++) mpz_sub(b_val, b_val, B_vals[j]);

        /* Initial solutions */
        for (int i = 0; i < fb->size; i++) {
            uint32_t p = fb->prime[i];
            unsigned long am = mpz_fdiv_ui(a, p);
            if (am == 0 || fb->root[i] == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            uint32_t ai = mod_inv32((uint32_t)am, p);
            if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned long bm = mpz_fdiv_ui(b_val, p);
            uint32_t r = fb->root[i];
            soln1[i] = (uint32_t)((uint64_t)ai * ((r + p - (uint32_t)bm) % p) % p);
            soln2[i] = (uint32_t)((uint64_t)ai * ((p - r + p - (uint32_t)bm) % p) % p);
        }

        for (int b_idx = 0; b_idx < num_b && full->count < target; b_idx++) {
            int gray = b_idx ^ (b_idx >> 1);

            if (b_idx > 0) {
                int changed = gray ^ prev_gray;
                int j = __builtin_ctz(changed);
                int sign = (gray >> j) & 1;

                if (sign) mpz_addmul_ui(b_val, B_vals[j], 2);
                else      mpz_submul_ui(b_val, B_vals[j], 2);

                for (int i = 0; i < fb->size; i++) {
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    uint32_t p = fb->prime[i], delta = AINV(j,i);
                    if (delta == 0) continue;
                    if (sign) {
                        soln1[i] = soln1[i] >= delta ? soln1[i] - delta : soln1[i] + p - delta;
                        soln2[i] = soln2[i] >= delta ? soln2[i] - delta : soln2[i] + p - delta;
                    } else {
                        soln1[i] += delta; if (soln1[i] >= p) soln1[i] -= p;
                        soln2[i] += delta; if (soln2[i] >= p) soln2[i] -= p;
                    }
                }
            }
            prev_gray = gray;

            /* c = (b^2 - kN) / a */
            mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN); mpz_divexact(c_val, c_val, a);
            total_polys++;

            /* === Fill buckets for large primes === */
            for (int i = 0; i < total_blocks; i++) buckets[i].count = 0;
            for (int i = bucket_thresh; i < fb->size; i++) {
                if (soln1[i] == 0xFFFFFFFF) continue;
                uint32_t p = fb->prime[i]; uint8_t lp = fb->logp[i];
                for (int root = 0; root < 2; root++) {
                    uint32_t s = root == 0 ? soln1[i] : soln2[i];
                    if (root == 1 && soln1[i] == soln2[i]) continue;
                    int64_t pos = ((int64_t)s - (-M)) % (int64_t)p;
                    if (pos < 0) pos += p;
                    int64_t x = -M + pos;
                    while (x < M) {
                        int bi = (int)((x + M) / BLOCK_SIZE);
                        if (bi >= 0 && bi < total_blocks) {
                            uint32_t bp = (uint32_t)((x + M) - (int64_t)bi * BLOCK_SIZE);
                            if (bp < BLOCK_SIZE) bucket_add(&buckets[bi], bp, lp);
                        }
                        x += p;
                    }
                }
            }

            /* === Sieve each block === */
            for (int blk = 0; blk < total_blocks; blk++) {
                int64_t bbase = (int64_t)blk * BLOCK_SIZE - M;
                memset(sieve, 0, BLOCK_SIZE);

                /* Small/medium primes */
                for (int i = 1; i < bucket_thresh; i++) {
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    uint32_t p = fb->prime[i];
                    if (p < 4) continue;
                    uint8_t lp = fb->logp[i];
                    int64_t off1 = ((int64_t)soln1[i] - bbase) % (int64_t)p; if (off1 < 0) off1 += p;
                    int64_t off2 = ((int64_t)soln2[i] - bbase) % (int64_t)p; if (off2 < 0) off2 += p;
                    if (soln1[i] == soln2[i]) {
                        for (int64_t j = off1; j < BLOCK_SIZE; j += p) sieve[j] += lp;
                    } else {
                        for (int64_t j = off1; j < BLOCK_SIZE; j += p) sieve[j] += lp;
                        for (int64_t j = off2; j < BLOCK_SIZE; j += p) sieve[j] += lp;
                    }
                }

                /* Bucket entries */
                for (int e = 0; e < buckets[blk].count; e++)
                    sieve[buckets[blk].entries[e].pos] += buckets[blk].entries[e].logp;

                /* === Scan for candidates (AVX512BW accelerated) === */
#ifdef __AVX512BW__
                __m512i thresh_vec = _mm512_set1_epi8((char)threshold);
                for (int j512 = 0; j512 < BLOCK_SIZE; j512 += 64) {
                    __m512i sv = _mm512_loadu_si512((__m512i*)(sieve + j512));
                    uint64_t mask = _mm512_cmpge_epu8_mask(sv, thresh_vec);
                    while (mask) {
                        int bit = __builtin_ctzll(mask);
                        mask &= mask - 1;
                        int j = j512 + bit;
#else
                for (int j = 0; j < BLOCK_SIZE; j++) {
                    if (sieve[j] < threshold) continue;
                    {
#endif
                    int64_t x = bbase + j;
                    if (x == 0) continue;

                    /* Compute Q(x) = a*x^2 + 2*b*x + c and ax+b */
                    mpz_set_si(tmp, (long)x);
                    mpz_mul_si(ax_b, a, (long)x); mpz_add(ax_b, ax_b, b_val);
                    /* Q(x) = (ax+b)^2/a - kN/a = ax^2 + 2bx + c */
                    mpz_mul(Qx, tmp, tmp); mpz_mul(Qx, Qx, a);
                    mpz_mul(tmp2, b_val, tmp); mpz_addmul_ui(Qx, tmp2, 2);
                    mpz_add(Qx, Qx, c_val);

                    if (mpz_sgn(Qx) == 0) continue;
                    mpz_abs(residue, Qx);

                    /* Trial division - divide by 2 first */
                    while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);

                    /* Sieve-informed TD: only divide by primes matching sieve roots */
                    for (int i = 1; i < fb->size; i++) {
                        uint32_t p = fb->prime[i];
                        if (soln1[i] == 0xFFFFFFFF) continue;
                        int64_t xmod = ((x % (int64_t)p) + p) % p;
                        if (xmod != (int64_t)soln1[i] && xmod != (int64_t)soln2[i]) continue;
                        while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                    }
                    /* a-factor primes */
                    for (int i = 0; i < num_a_factors; i++) {
                        uint32_t p = fb->prime[a_idx[i]];
                        while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                    }

                    mpz_mul(aQx, Qx, a);

                    if (mpz_cmp_ui(residue, 1) == 0) {
                        /* Full relation */
                        rels_add(full, ax_b, aQx, 0, 0);
                    } else if (mpz_fits_ulong_p(residue)) {
                        uint64_t cof = mpz_get_ui(residue);

                        /* Determine what kind of partial this cofactor represents.
                         * We try: SLP, DLP, TLP in order of increasing cost. */
                        int handled = 0;

                        if (cof <= lp_bound) {
                            /* SLP: single large prime */
                            handled = 1;
                            int match = lp_find(slp, cof);
                            if (match >= 0) {
                                /* Combine: product Q has lp^2 as factor; LA handles it
                                 * since lp^2 contributes even exponents (zero mod 2). */
                                mpz_mul(tmp, ax_b, part->ax_b[match]); mpz_mod(tmp, tmp, N);
                                mpz_mul(tmp2, aQx, part->Qx[match]);
                                rels_add(full, tmp, tmp2, 0, 0);
                                combined_slp++;
                            } else {
                                int pi = rels_add(part, ax_b, aQx, cof, 0);
                                if (pi >= 0) lp_insert(slp, cof, pi);
                            }
                        } else if (P.dlp && cof <= dlp_max) {
                            /* Extended SLP: if cofactor is prime but > lp_bound,
                             * still store in SLP hash for matching (like YAFU) */
                            /* Extended SLP disabled (causes correctness issues with LA) */
                            {
                                /* Composite cofactor: factor completely */
                                uint64_t primes[10]; int nprimes = factor64(cof, primes);
                                for (int a2=0;a2<nprimes-1;a2++) for (int b2=a2+1;b2<nprimes;b2++)
                                    if (primes[a2]>primes[b2]) {uint64_t t=primes[a2];primes[a2]=primes[b2];primes[b2]=t;}
                                int all_ok = 1;
                                for (int pi2=0;pi2<nprimes;pi2++) if (primes[pi2] > lp_bound) all_ok=0;

                                if (nprimes == 2 && all_ok) {
                                    uint64_t f1 = primes[0], f2 = primes[1];
                                    if (1) {
                                    /* Valid DLP: two primes both within bound */
                                    handled = 1;
                                    dlp_found++;

                                    /* DLP→SLP pipeline: check if either LP matches an SLP partial */
                                    int m1 = lp_find(slp, f1);
                                    int m2 = lp_find(slp, f2);
                                    if (m1 >= 0 && m2 >= 0 && m1 != m2) {
                                        /* Both LPs match SLP partials! DLP+SLP+SLP → full */
                                        mpz_mul(tmp, ax_b, part->ax_b[m1]);
                                        mpz_mul(tmp, tmp, part->ax_b[m2]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m1]);
                                        mpz_mul(tmp2, tmp2, part->Qx[m2]);
                                        rels_add(full, tmp, tmp2, 0, 0);
                                        combined_dlp++;
                                    } else if (m1 >= 0) {
                                        /* f1 matched SLP: DLP+SLP → new SLP partial with LP=f2 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m1]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m1]);
                                        int pi2 = rels_add(part, tmp, tmp2, f2, 0);
                                        if (pi2 >= 0) lp_insert(slp, f2, pi2);
                                        combined_dlp++;
                                    } else if (m2 >= 0) {
                                        /* f2 matched SLP: DLP+SLP → new SLP partial with LP=f1 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m2]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m2]);
                                        int pi2 = rels_add(part, tmp, tmp2, f1, 0);
                                        if (pi2 >= 0) lp_insert(slp, f1, pi2);
                                        combined_dlp++;
                                    } else {
                                        /* No SLP match: store as DLP partial with both LPs */
                                        int pi = rels_add(part, ax_b, aQx, f1, f2);
                                        if (pi >= 0) {
                                            lp_insert(slp, f1, pi);
                                            lp_insert(slp, f2, pi);
                                        }

                                        /* DLP graph and exact pair matching */
                                        int cycle = dlp_add_edge(dlp_g, f1, f2);
                                        if (cycle) dlp_cycles++;

                                        uint32_t h1 = (uint32_t)((f1 * 0x9E3779B97F4A7C15ULL) >> (64-DLP_PAIR_HASH_BITS));
                                        uint32_t h2 = (uint32_t)((f2 * 0x9E3779B97F4A7C15ULL) >> (64-DLP_PAIR_HASH_BITS));

                                        for (dp_e_t *e = dp_hash[h1]; e; e = e->next) {
                                            if (e->lp == f1) {
                                                int oi = e->rel_idx;
                                                if (part->lp2[oi] == f2 || part->lp1[oi] == f2) {
                                                    mpz_mul(tmp, ax_b, part->ax_b[oi]); mpz_mod(tmp, tmp, N);
                                                    mpz_mul(tmp2, aQx, part->Qx[oi]);
                                                    rels_add(full, tmp, tmp2, 0, 0);
                                                    combined_dlp++;
                                                    break;
                                                }
                                            }
                                        }

                                        if (dp_pool_used + 2 <= MAX_PARTIAL_RELS * 2) {
                                            dp_e_t *e1 = &dp_pool[dp_pool_used++];
                                            e1->lp = f1; e1->rel_idx = pi; e1->next = dp_hash[h1]; dp_hash[h1] = e1;
                                            dp_e_t *e2 = &dp_pool[dp_pool_used++];
                                            e2->lp = f2; e2->rel_idx = pi; e2->next = dp_hash[h2]; dp_hash[h2] = e2;
                                        }
                                    }
                                }
                                /* If f1 or f2 > lp_bound, fall through to TLP handling below */
                            }
                            } /* end else (composite cofactor) */
                            /* If factoring failed or primes out of range, falls through */
                        }

                        /* TLP handling: try to factor cof into 3 primes all <= lp_bound.
                         * This covers:
                         * (a) cof > lp_bound^2 (direct TLP range)
                         * (b) cof <= lp_bound^2 but split64 gave a composite or out-of-range factor */
                        if (!handled && P.dlp) {
                            __uint128_t lp3 = (__uint128_t)lp_bound * lp_bound * lp_bound;
                            if ((__uint128_t)cof <= lp3) {
                                uint64_t factors[64];
                                int nf = factor64(cof, factors);
                                /* Sort factors ascending */
                                for (int fi = 0; fi < nf; fi++)
                                    for (int fj = fi+1; fj < nf; fj++)
                                        if (factors[fi] > factors[fj]) {
                                            uint64_t tmp_f = factors[fi];
                                            factors[fi] = factors[fj];
                                            factors[fj] = tmp_f;
                                        }

                                if (nf == 3 &&
                                    factors[0] <= lp_bound &&
                                    factors[1] <= lp_bound &&
                                    factors[2] <= lp_bound) {
                                    /* Valid TLP! */
                                    tlp_found++;
                                    /* Store as DLP partial: use the two smallest primes as lp1,lp2.
                                     * Also register under the largest prime for cross-matching.
                                     * When another relation matches on one of these primes and the
                                     * DLP matching pipeline resolves the other two, we get a full rel. */
                                    int m0 = lp_find(slp, factors[0]);
                                    int m1 = lp_find(slp, factors[1]);
                                    int m2 = lp_find(slp, factors[2]);

                                    if (m2 >= 0 && m1 >= 0 && m0 >= 0) {
                                        /* All 3 primes match existing partials - rare but possible */
                                        /* Skip complex combination, just store as DLP */
                                    }
                                    if (m2 >= 0 && m1 >= 0 && m0 < 0) {
                                        /* Two largest match: TLP + SLP(f2) + SLP(f1) → new SLP with f0 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m2]);
                                        mpz_mul(tmp, tmp, part->ax_b[m1]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m2]);
                                        mpz_mul(tmp2, tmp2, part->Qx[m1]);
                                        int pi2 = rels_add(part, tmp, tmp2, factors[0], 0);
                                        if (pi2 >= 0) lp_insert(slp, factors[0], pi2);
                                        combined_tlp++;
                                    } else if (m2 >= 0 && m0 >= 0 && m1 < 0) {
                                        /* f2 and f0 match: TLP + SLP(f2) + SLP(f0) → new SLP with f1 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m2]);
                                        mpz_mul(tmp, tmp, part->ax_b[m0]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m2]);
                                        mpz_mul(tmp2, tmp2, part->Qx[m0]);
                                        int pi2 = rels_add(part, tmp, tmp2, factors[1], 0);
                                        if (pi2 >= 0) lp_insert(slp, factors[1], pi2);
                                        combined_tlp++;
                                    } else if (m1 >= 0 && m0 >= 0 && m2 < 0) {
                                        /* f1 and f0 match: TLP + SLP(f1) + SLP(f0) → new SLP with f2 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m1]);
                                        mpz_mul(tmp, tmp, part->ax_b[m0]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m1]);
                                        mpz_mul(tmp2, tmp2, part->Qx[m0]);
                                        int pi2 = rels_add(part, tmp, tmp2, factors[2], 0);
                                        if (pi2 >= 0) lp_insert(slp, factors[2], pi2);
                                        combined_tlp++;
                                    } else if (m2 >= 0) {
                                        /* TLP + SLP(f2) → DLP partial with f0,f1 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m2]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m2]);
                                        int pi2 = rels_add(part, tmp, tmp2, factors[0], factors[1]);
                                        if (pi2 >= 0) {
                                            lp_insert(slp, factors[0], pi2);
                                            lp_insert(slp, factors[1], pi2);
                                        }
                                        combined_tlp++;
                                    } else if (m1 >= 0) {
                                        /* TLP + SLP(f1) → DLP partial with f0,f2 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m1]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m1]);
                                        int pi2 = rels_add(part, tmp, tmp2, factors[0], factors[2]);
                                        if (pi2 >= 0) {
                                            lp_insert(slp, factors[0], pi2);
                                            lp_insert(slp, factors[2], pi2);
                                        }
                                        combined_tlp++;
                                    } else if (m0 >= 0) {
                                        /* TLP + SLP(f0) → DLP partial with f1,f2 */
                                        mpz_mul(tmp, ax_b, part->ax_b[m0]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m0]);
                                        int pi2 = rels_add(part, tmp, tmp2, factors[1], factors[2]);
                                        if (pi2 >= 0) {
                                            lp_insert(slp, factors[1], pi2);
                                            lp_insert(slp, factors[2], pi2);
                                        }
                                        combined_tlp++;
                                    } else {
                                        /* No existing matches: store TLP as DLP partial (two smallest primes) */
                                        int pi = rels_add(part, ax_b, aQx, factors[0], factors[1]);
                                        if (pi >= 0) {
                                            lp_insert(slp, factors[0], pi);
                                            lp_insert(slp, factors[1], pi);
                                            lp_insert(slp, factors[2], pi);  /* register 3rd factor too */
                                        }
                                    }
                                } else if (nf == 2 &&
                                           factors[0] <= lp_bound &&
                                           factors[1] <= lp_bound) {
                                    /* Two primes within bound - treat as DLP */
                                    uint64_t f1 = factors[0], f2 = factors[1];
                                    dlp_found++;
                                    int m1 = lp_find(slp, f1);
                                    int m2_dlp = lp_find(slp, f2);
                                    if (m1 >= 0 && m2_dlp >= 0 && m1 != m2_dlp) {
                                        mpz_mul(tmp, ax_b, part->ax_b[m1]);
                                        mpz_mul(tmp, tmp, part->ax_b[m2_dlp]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m1]);
                                        mpz_mul(tmp2, tmp2, part->Qx[m2_dlp]);
                                        rels_add(full, tmp, tmp2, 0, 0);
                                        combined_dlp++;
                                    } else if (m1 >= 0) {
                                        mpz_mul(tmp, ax_b, part->ax_b[m1]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m1]);
                                        int pi2 = rels_add(part, tmp, tmp2, f2, 0);
                                        if (pi2 >= 0) lp_insert(slp, f2, pi2);
                                        combined_dlp++;
                                    } else if (m2_dlp >= 0) {
                                        mpz_mul(tmp, ax_b, part->ax_b[m2_dlp]);
                                        mpz_mod(tmp, tmp, N);
                                        mpz_mul(tmp2, aQx, part->Qx[m2_dlp]);
                                        int pi2 = rels_add(part, tmp, tmp2, f1, 0);
                                        if (pi2 >= 0) lp_insert(slp, f1, pi2);
                                        combined_dlp++;
                                    } else {
                                        int pi = rels_add(part, ax_b, aQx, f1, f2);
                                        if (pi >= 0) {
                                            lp_insert(slp, f1, pi);
                                            lp_insert(slp, f2, pi);
                                        }
                                    }
                                }
                            }
                        }
                    } else if (P.dlp) {
                        /* Residue > 64 bits. Skip */
                    }
#ifdef __AVX512BW__
                    }  /* end while(mask) */
                }  /* end j512 loop */
#else
                }  /* end if threshold check */
                }  /* end j loop */
#endif
            }

            /* Progress */
            if (total_polys % 500 == 0) {
                double t = elapsed_sec();
                if (t > 275.0) break;
                if (total_polys % 2000 == 0)
                    fprintf(stderr, "  poly=%d rels=%d/%d (full=%d slp=%d dlp=%d tlp=%d) part=%d dlp_found=%d tlp_found=%d t=%.1fs\n",
                            total_polys, full->count, target,
                            full->count-combined_slp-combined_dlp-combined_tlp, combined_slp, combined_dlp, combined_tlp,
                            part->count, dlp_found, tlp_found, t);
            }
        }
    }

    double sieve_time = elapsed_sec();
    fprintf(stderr, "Sieve: %d rels in %.2fs (%d polys, %d SLP, %d DLP, %d TLP, %d dlp_found, %d tlp_found)\n",
            full->count, sieve_time, total_polys, combined_slp, combined_dlp, combined_tlp, dlp_found, tlp_found);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full->count, fb->size + 1);
        return 1;
    }

    /* ==================== Linear Algebra ==================== */
    int ncols = fb->size + 1;
    fprintf(stderr, "Building sparse matrix: %d x %d\n", full->count, ncols);

    /* Build sparse matrix in CRS format */
    /* First pass: count nonzeros per row */
    int *nnz_per_row = calloc(full->count, sizeof(int));
    int total_nnz = 0;

    for (int r = 0; r < full->count; r++) {
        mpz_abs(residue, full->Qx[r]);
        if (mpz_sgn(full->Qx[r]) < 0) nnz_per_row[r]++; /* sign col */

        int exp2 = 0;
        while (mpz_even_p(residue)) { mpz_tdiv_q_2exp(residue, residue, 1); exp2++; }
        if (exp2 & 1) nnz_per_row[r]++;

        for (int i = 1; i < fb->size; i++) {
            uint32_t p = fb->prime[i]; int exp = 0;
            while (mpz_divisible_ui_p(residue, p)) { mpz_divexact_ui(residue, residue, p); exp++; }
            if (exp & 1) nnz_per_row[r]++;
        }
        total_nnz += nnz_per_row[r];
    }

    sparse_t *smat = calloc(1, sizeof(sparse_t));
    smat->nrows = full->count;
    smat->ncols = ncols;
    smat->nnz = total_nnz;
    smat->row_start = malloc((full->count + 1) * sizeof(int));
    smat->col_idx = malloc(total_nnz * sizeof(int));

    /* Second pass: fill sparse matrix */
    smat->row_start[0] = 0;
    for (int r = 0; r < full->count; r++) {
        int pos = smat->row_start[r];
        mpz_abs(residue, full->Qx[r]);

        if (mpz_sgn(full->Qx[r]) < 0) smat->col_idx[pos++] = 0;

        int exp2 = 0;
        while (mpz_even_p(residue)) { mpz_tdiv_q_2exp(residue, residue, 1); exp2++; }
        if (exp2 & 1) smat->col_idx[pos++] = 1;

        for (int i = 1; i < fb->size; i++) {
            uint32_t p = fb->prime[i]; int exp = 0;
            while (mpz_divisible_ui_p(residue, p)) { mpz_divexact_ui(residue, residue, p); exp++; }
            if (exp & 1) smat->col_idx[pos++] = i + 1;
        }
        smat->row_start[r + 1] = pos;
    }

    fprintf(stderr, "Sparse matrix: %d nnz (avg %.1f per row), solving...\n",
            total_nnz, (double)total_nnz / full->count);

    int **deps; int *dlen;
    int ndeps = block_lanczos_solve(smat, &deps, &dlen, 64);
    fprintf(stderr, "BL found %d deps in %.2fs\n", ndeps, elapsed_sec());

    if (ndeps < 5) {
        fprintf(stderr, "BL insufficient deps, falling back to structured Gauss\n");
        sg_mat_t *sg = sg_create(full->count, ncols);
        for (int r = 0; r < full->count; r++) {
            for (int j = smat->row_start[r]; j < smat->row_start[r+1]; j++)
                sg_set(sg, r, smat->col_idx[j]);
        }
        if (deps) { for (int i = 0; i < ndeps; i++) free(deps[i]); free(deps); free(dlen); }
        ndeps = sg_solve(sg, &deps, &dlen, 64);
        fprintf(stderr, "SG found %d deps in %.2fs\n", ndeps, elapsed_sec());
        free(sg);
    }

    free(nnz_per_row);

    if (ndeps == 0) { fprintf(stderr, "FAIL: no dependencies\n"); return 1; }

    /* ==================== Square Root ==================== */
    mpz_t X, Y, g;
    mpz_inits(X, Y, g, NULL);
    int factored = 0;

    for (int d = 0; d < ndeps && !factored; d++) {
        mpz_set_ui(X, 1); mpz_set_ui(Y, 1);
        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            mpz_mul(X, X, full->ax_b[ri]); mpz_mod(X, X, N);
            mpz_mul(Y, Y, full->Qx[ri]);
        }
        mpz_abs(Y, Y);
        if (!mpz_perfect_square_p(Y)) continue;
        mpz_sqrt(Y, Y);
        mpz_mod(Y, Y, N);

        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) { gmp_printf("%Zd\n", g); factored = 1; }
        else { mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
               if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) { gmp_printf("%Zd\n", g); factored = 1; } }
    }

    if (!factored) { fprintf(stderr, "FAIL: no factor from %d deps\n", ndeps); return 1; }
    fprintf(stderr, "Total time: %.3fs\n", elapsed_sec());

    /* Cleanup */
    free(smat->row_start); free(smat->col_idx); free(smat);
    mpz_clears(N, kN, a, b_val, c_val, ax_b, Qx, aQx, residue, tmp, tmp2, X, Y, g, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_clear(B_vals[j]);
    gmp_randclear(rng);
    return 0;
}
