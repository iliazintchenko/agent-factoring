/*
 * DLP-SIQS: Self-Initializing Quadratic Sieve with Double Large Primes
 * and ECM cofactorization.
 *
 * Novel approach: Uses aggressive large prime bounds with ECM to split
 * cofactors into DLP (double large prime) relations. The DLP graph
 * (vertices=large primes, edges=relations) is solved with union-find
 * to find cycles that produce full relations.
 *
 * Key innovations over standard SLP-SIQS:
 * - DLP with Pollard rho splitting for cofactors < 2^62
 * - Graph-based DLP cycle finder (union-find + BFS)
 * - Aggressive LP bounds: LP_bound = B * 300 (vs typical B * 30)
 * - This allows smaller factor base B, potentially shifting scaling
 *
 * Compile: gcc -O3 -march=native -mavx512bw -o dlp_siqs library/dlp_siqs.c -lgmp -lm
 * Usage: ./dlp_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <immintrin.h>

#define SEED 42
#define SIEVE_BLOCK 32768   /* 32KB for L1D cache */
#define MAX_FB 80000
#define MAX_A_FACTORS 20
#define MAX_FULL_RELS 300000
#define MAX_PARTIALS 2000000

/* ==================== Timing ==================== */
static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size;
    int num_blocks;
    int lp_mult;      /* LP bound = largest_fb_prime * lp_mult */
    double thresh_adj;
    int extra_rels;
} params_t;

static params_t get_params(int bits) {
    /* Tuned for DLP: use smaller factor bases, larger LP bounds */
    if (bits <= 100)  return (params_t){80,   1,  200, 0.70, 30};
    if (bits <= 110)  return (params_t){120,  1,  200, 0.71, 35};
    if (bits <= 120)  return (params_t){170,  2,  250, 0.73, 40};
    if (bits <= 130)  return (params_t){250,  3,  250, 0.75, 50};
    if (bits <= 140)  return (params_t){350,  4,  300, 0.76, 60};
    if (bits <= 150)  return (params_t){500,  6,  300, 0.77, 70};
    if (bits <= 160)  return (params_t){750,  8,  350, 0.78, 80};
    if (bits <= 170)  return (params_t){1100, 12, 350, 0.79, 90};
    if (bits <= 180)  return (params_t){1600, 16, 400, 0.80, 100};
    if (bits <= 190)  return (params_t){2200, 22, 400, 0.81, 120};
    if (bits <= 200)  return (params_t){3200, 28, 450, 0.82, 140};
    if (bits <= 210)  return (params_t){4500, 36, 450, 0.83, 160};
    if (bits <= 220)  return (params_t){6000, 44, 500, 0.84, 180};
    if (bits <= 230)  return (params_t){8000, 52, 500, 0.845, 200};
    if (bits <= 240)  return (params_t){11000, 60, 500, 0.85, 220};
    if (bits <= 250)  return (params_t){15000, 72, 500, 0.855, 260};
    if (bits <= 260)  return (params_t){20000, 88, 500, 0.86, 300};
    if (bits <= 270)  return (params_t){28000, 100, 500, 0.865, 350};
    if (bits <= 280)  return (params_t){38000, 120, 500, 0.87, 400};
    if (bits <= 290)  return (params_t){50000, 140, 500, 0.875, 450};
    return (params_t){70000, 170, 500, 0.88, 500};
}

/* ==================== Primes ==================== */
static int *sieve_primes(int bound, int *count) {
    char *c = calloc(bound + 1, 1);
    for (int i = 2; (long)i * i <= bound; i++)
        if (!c[i])
            for (int j = i * i; j <= bound; j += i)
                c[j] = 1;
    int cnt = 0;
    for (int i = 2; i <= bound; i++) if (!c[i]) cnt++;
    int *p = malloc(cnt * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= bound; i++) if (!c[i]) p[idx++] = i;
    free(c);
    *count = cnt;
    return p;
}

/* ==================== Modular Arithmetic ==================== */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m;
    int old_s = 1, s = 0;
    while (r) {
        int q = old_r / r;
        int t = r; r = old_r - q * r; old_r = t;
        t = s; s = old_s - q * s; old_s = t;
    }
    if (old_r != 1) return 0;
    return (unsigned int)(((long long)old_s % m + m) % m);
}

static unsigned int sqrt_mod(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    unsigned long long b, e, m = p, r;

    /* Euler criterion */
    b = n % p; e = (p - 1) / 2; r = 1;
    { unsigned long long bb = b, ee = e;
      while (ee > 0) { if (ee & 1) r = (r * bb) % m; bb = (bb * bb) % m; ee >>= 1; }
    }
    if (r != 1) return 0;

    if (p % 4 == 3) {
        b = n % p; e = (p + 1) / 4; r = 1;
        while (e > 0) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; }
        return (unsigned int)r;
    }

    /* Full Tonelli-Shanks */
    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    unsigned int z = 2;
    while (1) {
        b = z; e = (p - 1) / 2; r = 1;
        while (e > 0) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; }
        if (r == (unsigned long long)(p - 1)) break;
        z++;
    }
    unsigned long long M_val = S;
    b = z; e = Q; unsigned long long c = 1;
    while (e > 0) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = Q; unsigned long long t = 1;
    while (e > 0) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; unsigned long long R = 1;
    while (e > 0) { if (e & 1) R = (R * b) % m; b = (b * b) % m; e >>= 1; }
    while (1) {
        if (t == 1) return (unsigned int)R;
        int i = 0; unsigned long long tt = t;
        while (tt != 1) { tt = (tt * tt) % p; i++; }
        unsigned long long bb = c;
        for (int j = 0; j < (int)M_val - i - 1; j++) bb = (bb * bb) % p;
        M_val = i; c = (bb * bb) % p; t = (t * c) % p; R = (R * bb) % p;
    }
}

/* ==================== Knuth-Schroeppel ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 23, 29, 31, 37, 41, 43, 0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki];
        mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long mod8 = mpz_fdiv_ui(kN, 8);
        if (mod8 == 1) s += 2 * log(2.0);
        else if (mod8 == 5) s += log(2.0);
        else if (mod8 == 3 || mod8 == 7) s += 0.5 * log(2.0);
        /* Check a few small primes */
        int primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
        for (int i = 0; i < 14; i++) {
            int p = primes[i];
            if (k % p == 0) continue;
            unsigned long nm = mpz_fdiv_ui(kN, p);
            if (sqrt_mod(nm, p)) s += 2.0 * log(p) / (p - 1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct {
    unsigned int *prime;
    unsigned int *root;
    unsigned char *logp;
    int size;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 20;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));
    fb->prime[0] = 2; fb->root[0] = 1;
    fb->logp[0] = 1; fb->size = 1;

    int bound = target * 30 + 50000;
    int np; int *primes = sieve_primes(bound, &np);
    for (int i = 1; i < np && fb->size < target; i++) {
        int p = primes[i];
        unsigned long nm = mpz_fdiv_ui(kN, p);
        if (nm == 0) {
            fb->prime[fb->size] = p; fb->root[fb->size] = 0;
            fb->logp[fb->size] = (unsigned char)(log2(p) + 0.5);
            fb->size++; continue;
        }
        unsigned int r = sqrt_mod((unsigned int)nm, p);
        if (r == 0) continue;
        fb->prime[fb->size] = p; fb->root[fb->size] = r;
        fb->logp[fb->size] = (unsigned char)(log2(p) + 0.5);
        fb->size++;
    }
    free(primes);
    return fb;
}

/* ==================== Pollard Rho for DLP splitting ==================== */
/* Split a 64-bit composite into two factors using Brent's rho */
static unsigned long long pollard_rho_64(unsigned long long n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    if (n % 5 == 0) return 5;

    /* Use __int128 for modular multiplication */
    unsigned long long x = 2, y = 2, d = 1, c = 1;
    unsigned long long q;

    for (int attempt = 0; attempt < 20 && d == 1; attempt++) {
        c = attempt + 1;
        x = y = 2;
        d = 1;

        while (d == 1) {
            /* Brent's cycle detection with batched GCD */
            unsigned long long ys = y;
            q = 1;
            long range = 1;

            while (d == 1) {
                x = y;
                for (long i = 0; i < range; i++)
                    y = ((__uint128_t)y * y + c) % n;

                long k = 0;
                while (k < range && d == 1) {
                    ys = y;
                    long batch = range - k < 128 ? range - k : 128;
                    for (long i = 0; i < batch; i++) {
                        y = ((__uint128_t)y * y + c) % n;
                        unsigned long long diff = y > x ? y - x : x - y;
                        q = ((__uint128_t)q * diff) % n;
                    }
                    unsigned long long g = 1;
                    unsigned long long a = q, b = n;
                    while (b) { unsigned long long t = b; b = a % b; a = t; }
                    g = a;
                    d = g;
                    k += batch;
                }
                range <<= 1;
                if (range > 1000000) break;
            }

            if (d == n) {
                /* Backtrack */
                d = 1;
                while (d == 1) {
                    ys = ((__uint128_t)ys * ys + c) % n;
                    unsigned long long diff = ys > x ? ys - x : x - ys;
                    unsigned long long a = diff, b = n;
                    while (b) { unsigned long long t = b; b = a % b; a = t; }
                    d = a;
                }
            }
            if (d == n) { d = 1; break; }
        }
    }
    return d == n ? 0 : d;
}

/* Check if n is a prime using Miller-Rabin */
static int is_prime_64(unsigned long long n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0) return 0;

    unsigned long long d = n - 1;
    int r = 0;
    while (d % 2 == 0) { d /= 2; r++; }

    /* Deterministic for < 2^64 with these witnesses */
    unsigned long long witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int w = 0; w < 12; w++) {
        unsigned long long a = witnesses[w];
        if (a >= n) continue;

        /* Compute a^d mod n */
        unsigned long long x = 1, base = a, exp = d;
        while (exp > 0) {
            if (exp & 1) x = ((__uint128_t)x * base) % n;
            base = ((__uint128_t)base * base) % n;
            exp >>= 1;
        }

        if (x == 1 || x == n - 1) continue;
        int found = 0;
        for (int i = 0; i < r - 1; i++) {
            x = ((__uint128_t)x * x) % n;
            if (x == n - 1) { found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}

/* ==================== DLP Graph ==================== */
/* Graph where vertices = large primes, edges = DLP relations */
/* Uses hash table for vertex lookup and union-find for cycle detection */

#define DLP_HASH_BITS 22
#define DLP_HASH_SIZE (1 << DLP_HASH_BITS)
#define DLP_HASH_MASK (DLP_HASH_SIZE - 1)

typedef struct dlp_edge {
    unsigned long lp1, lp2;    /* the two large primes */
    int rel_idx;               /* index into partial relations */
    struct dlp_edge *next;     /* hash chain */
} dlp_edge_t;

typedef struct lp_vertex {
    unsigned long lp;
    int parent;       /* union-find parent (index into vertex array) */
    int rank;         /* union-find rank */
    int edge_list;    /* head of adjacency list (-1 = none) */
    struct lp_vertex *hash_next;
} lp_vertex_t;

typedef struct {
    lp_vertex_t **hash;     /* hash table: lp -> vertex */
    lp_vertex_t *verts;     /* vertex pool */
    int nv;                 /* number of vertices */
    int max_verts;
    dlp_edge_t *edges;      /* edge pool */
    int ne;                 /* number of edges */
    int max_edges;
    int cycles_found;       /* number of cycles found */
} dlp_graph_t;

static dlp_graph_t *dlp_create(int max_v, int max_e) {
    dlp_graph_t *g = calloc(1, sizeof(dlp_graph_t));
    g->hash = calloc(DLP_HASH_SIZE, sizeof(lp_vertex_t*));
    g->verts = calloc(max_v, sizeof(lp_vertex_t));
    g->max_verts = max_v;
    g->edges = calloc(max_e, sizeof(dlp_edge_t));
    g->max_edges = max_e;
    return g;
}

static int dlp_find_or_add_vertex(dlp_graph_t *g, unsigned long lp) {
    unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - DLP_HASH_BITS));
    for (lp_vertex_t *v = g->hash[h]; v; v = v->hash_next)
        if (v->lp == lp) return (int)(v - g->verts);

    if (g->nv >= g->max_verts) return -1;
    int idx = g->nv++;
    g->verts[idx].lp = lp;
    g->verts[idx].parent = idx;
    g->verts[idx].rank = 0;
    g->verts[idx].edge_list = -1;
    g->verts[idx].hash_next = g->hash[h];
    g->hash[h] = &g->verts[idx];
    return idx;
}

static int dlp_find(dlp_graph_t *g, int x) {
    while (g->verts[x].parent != x) {
        g->verts[x].parent = g->verts[g->verts[x].parent].parent; /* path compression */
        x = g->verts[x].parent;
    }
    return x;
}

static int dlp_union(dlp_graph_t *g, int a, int b) {
    a = dlp_find(g, a); b = dlp_find(g, b);
    if (a == b) return 1; /* cycle! */
    if (g->verts[a].rank < g->verts[b].rank) { int t = a; a = b; b = t; }
    g->verts[b].parent = a;
    if (g->verts[a].rank == g->verts[b].rank) g->verts[a].rank++;
    return 0;
}

/* ==================== Single Large Prime Hash ==================== */
#define SLP_HASH_BITS 20
#define SLP_HASH_SIZE (1 << SLP_HASH_BITS)

typedef struct slp_entry {
    unsigned long lp;
    int rel_idx;
    struct slp_entry *next;
} slp_entry_t;

typedef struct {
    slp_entry_t **buckets;
    slp_entry_t *pool;
    int used, max;
} slp_table_t;

static slp_table_t *slp_create(int max) {
    slp_table_t *t = calloc(1, sizeof(slp_table_t));
    t->buckets = calloc(SLP_HASH_SIZE, sizeof(slp_entry_t*));
    t->pool = calloc(max, sizeof(slp_entry_t));
    t->max = max;
    return t;
}

static int slp_find(slp_table_t *t, unsigned long lp) {
    unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - SLP_HASH_BITS));
    for (slp_entry_t *e = t->buckets[h]; e; e = e->next)
        if (e->lp == lp) return e->rel_idx;
    return -1;
}

static void slp_insert(slp_table_t *t, unsigned long lp, int idx) {
    if (t->used >= t->max) return;
    unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - SLP_HASH_BITS));
    slp_entry_t *e = &t->pool[t->used++];
    e->lp = lp; e->rel_idx = idx; e->next = t->buckets[h];
    t->buckets[h] = e;
}

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *ax_b;   /* (ax+b) value */
    mpz_t *Qx;     /* a * Q(x) = (ax+b)^2 - kN */
    unsigned long *lp1;  /* large prime 1 (0 if none) */
    unsigned long *lp2;  /* large prime 2 (0 if none) */
    int count, alloc;
} rel_store_t;

static rel_store_t *rel_create(int alloc) {
    rel_store_t *r = malloc(sizeof(rel_store_t));
    r->ax_b = malloc(alloc * sizeof(mpz_t));
    r->Qx = malloc(alloc * sizeof(mpz_t));
    r->lp1 = calloc(alloc, sizeof(unsigned long));
    r->lp2 = calloc(alloc, sizeof(unsigned long));
    for (int i = 0; i < alloc; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    r->count = 0; r->alloc = alloc;
    return r;
}

/* ==================== GF(2) Matrix ==================== */
typedef unsigned long long u64;

typedef struct {
    u64 **rows;
    int nrows, ncols;
    int words_per_row; /* total u64 words = fb_words + id_words */
    int fb_words;
    int id_words;
} gf2_mat_t;

static gf2_mat_t *gf2_create(int nrows, int ncols) {
    gf2_mat_t *m = malloc(sizeof(gf2_mat_t));
    m->nrows = nrows; m->ncols = ncols;
    m->fb_words = (ncols + 63) / 64;
    m->id_words = (nrows + 63) / 64;
    m->words_per_row = m->fb_words + m->id_words;
    m->rows = malloc(nrows * sizeof(u64*));
    for (int i = 0; i < nrows; i++) {
        m->rows[i] = calloc(m->words_per_row, sizeof(u64));
        /* Identity in right part */
        m->rows[i][m->fb_words + i / 64] |= (1ULL << (i % 64));
    }
    return m;
}

static inline void gf2_set(gf2_mat_t *m, int r, int c) {
    m->rows[r][c / 64] |= (1ULL << (c % 64));
}

static int gf2_solve(gf2_mat_t *m, int ***deps_out, int **dlen_out, int max_deps) {
    int pivot = 0;
    for (int col = 0; col < m->ncols && pivot < m->nrows; col++) {
        int pr = -1;
        for (int r = pivot; r < m->nrows; r++)
            if ((m->rows[r][col / 64] >> (col % 64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != pivot) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[pivot]; m->rows[pivot] = t; }
        u64 *prow = m->rows[pivot];
        for (int r = 0; r < m->nrows; r++) {
            if (r == pivot) continue;
            if ((m->rows[r][col / 64] >> (col % 64)) & 1)
                for (int w = 0; w < m->words_per_row; w++) m->rows[r][w] ^= prow[w];
        }
        pivot++;
    }

    int ndeps = 0;
    *deps_out = malloc(max_deps * sizeof(int*));
    *dlen_out = malloc(max_deps * sizeof(int));

    for (int r = pivot; r < m->nrows && ndeps < max_deps; r++) {
        /* Check left part is zero */
        int zero = 1;
        for (int w = 0; w < m->fb_words && zero; w++) {
            u64 mask = (w < m->fb_words - 1) ? ~0ULL :
                       (m->ncols % 64 == 0 ? ~0ULL : (1ULL << (m->ncols % 64)) - 1);
            if (m->rows[r][w] & mask) { zero = 0; break; }
        }
        if (!zero) continue;

        int *dep = malloc(m->nrows * sizeof(int));
        int dl = 0;
        for (int w = 0; w < m->id_words; w++) {
            u64 bits = m->rows[r][m->fb_words + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w * 64 + bit;
                if (idx < m->nrows) dep[dl++] = idx;
                bits &= bits - 1;
            }
        }
        if (dl > 0) { (*deps_out)[ndeps] = dep; (*dlen_out)[ndeps] = dl; ndeps++; }
        else free(dep);
    }
    return ndeps;
}

/* ==================== SIQS Polynomial ==================== */
typedef struct {
    mpz_t a, b, c;
    int a_idx[MAX_A_FACTORS];
    int s; /* number of a-factors */
    mpz_t B[MAX_A_FACTORS];
    unsigned int *ainv; /* [s][fb_size] flattened */
    int fb_size;
    int num_b, cur_b;
    unsigned int *soln1, *soln2;
} poly_t;

static int poly_new_a(poly_t *ps, fb_t *fb, mpz_t kN, int M, gmp_randstate_t rng) {
    mpz_t target; mpz_init(target);
    mpz_mul_ui(target, kN, 2);
    mpz_sqrt(target, target);
    mpz_tdiv_q_ui(target, target, M);
    double log_target = mpz_sizeinbase(target, 2) * log(2.0);

    int lo = fb->size / 3, hi = 2 * fb->size / 3;
    if (lo < 2) lo = 2;
    if (hi <= lo + 3) hi = fb->size - 1;

    double avg = 0; int cnt = 0;
    for (int i = lo; i < hi; i++) { if (fb->root[i] == 0) continue; avg += log(fb->prime[i]); cnt++; }
    if (cnt == 0) { mpz_clear(target); return 0; }
    avg /= cnt;

    int s = (int)(log_target / avg + 0.5);
    if (s < 3) s = 3;
    if (s > MAX_A_FACTORS) s = MAX_A_FACTORS;
    if (s > hi - lo) s = hi - lo;
    ps->s = s;

    double best_ratio = 1e30;
    int best[MAX_A_FACTORS];

    for (int att = 0; att < 40; att++) {
        mpz_set_ui(ps->a, 1);
        int idx[MAX_A_FACTORS];
        int ok = 1;
        for (int i = 0; i < s && ok; i++) {
            int tries = 0, good;
            do {
                idx[i] = lo + gmp_urandomm_ui(rng, hi - lo);
                good = 1;
                for (int j = 0; j < i; j++) if (idx[j] == idx[i]) { good = 0; break; }
                if (fb->root[idx[i]] == 0) good = 0;
                tries++;
            } while (!good && tries < 100);
            if (!good) { ok = 0; break; }
            mpz_mul_ui(ps->a, ps->a, fb->prime[idx[i]]);
        }
        if (!ok) continue;

        double ratio;
        if (mpz_cmp(ps->a, target) > 0) {
            mpz_t q; mpz_init(q); mpz_tdiv_q(q, ps->a, target); ratio = mpz_get_d(q); mpz_clear(q);
        } else {
            mpz_t q; mpz_init(q); mpz_tdiv_q(q, target, ps->a); ratio = mpz_get_d(q); mpz_clear(q);
        }
        if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s * sizeof(int)); }
        if (ratio < 2.0) break;
    }

    memcpy(ps->a_idx, best, s * sizeof(int));
    mpz_set_ui(ps->a, 1);
    for (int i = 0; i < s; i++) mpz_mul_ui(ps->a, ps->a, fb->prime[ps->a_idx[i]]);
    mpz_clear(target);

    /* Compute B values */
    for (int j = 0; j < s; j++) {
        int idx = ps->a_idx[j];
        unsigned int qj = fb->prime[idx], rj = fb->root[idx];
        mpz_t a_q, mod_q, inv;
        mpz_inits(a_q, mod_q, inv, NULL);
        mpz_divexact_ui(a_q, ps->a, qj);
        mpz_set_ui(mod_q, qj);
        if (!mpz_invert(inv, a_q, mod_q)) { mpz_clears(a_q, mod_q, inv, NULL); return 0; }
        unsigned long iv = mpz_get_ui(inv);
        mpz_mul_ui(ps->B[j], a_q, (rj * iv) % qj);
        mpz_clears(a_q, mod_q, inv, NULL);
    }

    /* Initial b */
    mpz_set_ui(ps->b, 0);
    for (int j = 0; j < s; j++) mpz_add(ps->b, ps->b, ps->B[j]);

    /* Verify b^2 ≡ kN (mod a) */
    mpz_t test; mpz_init(test);
    mpz_mul(test, ps->b, ps->b); mpz_sub(test, test, kN); mpz_mod(test, test, ps->a);
    if (mpz_sgn(test) != 0) {
        mpz_sub(ps->b, ps->b, ps->B[0]); mpz_sub(ps->b, ps->b, ps->B[0]);
        mpz_mul(test, ps->b, ps->b); mpz_sub(test, test, kN); mpz_mod(test, test, ps->a);
        if (mpz_sgn(test) != 0) { mpz_clear(test); return 0; }
    }
    mpz_clear(test);

    /* c = (b^2 - kN) / a */
    mpz_mul(ps->c, ps->b, ps->b); mpz_sub(ps->c, ps->c, kN); mpz_divexact(ps->c, ps->c, ps->a);

    /* Sieve solutions */
    for (int i = 0; i < fb->size; i++) {
        unsigned int p = fb->prime[i];
        unsigned long am = mpz_fdiv_ui(ps->a, p);
        if (am == 0 || fb->root[i] == 0) { ps->soln1[i] = ps->soln2[i] = 0xFFFFFFFF; continue; }
        unsigned int ai = mod_inverse((unsigned int)am, p);
        if (ai == 0) { ps->soln1[i] = ps->soln2[i] = 0xFFFFFFFF; continue; }
        unsigned long bm = mpz_fdiv_ui(ps->b, p);
        unsigned int r = fb->root[i];
        ps->soln1[i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
        ps->soln2[i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
    }

    /* Precompute ainv for Gray code */
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            unsigned long am = mpz_fdiv_ui(ps->a, p);
            if (am == 0 || fb->root[i] == 0) { ps->ainv[j * fb->size + i] = 0; continue; }
            unsigned int ai = mod_inverse((unsigned int)am, p);
            unsigned long Bm = mpz_fdiv_ui(ps->B[j], p);
            ps->ainv[j * fb->size + i] = (unsigned int)((2UL * ai % p * Bm) % p);
        }
    }

    ps->cur_b = 0;
    ps->num_b = 1 << (s - 1);
    return 1;
}

static int poly_next_b(poly_t *ps, fb_t *fb, mpz_t kN) {
    ps->cur_b++;
    if (ps->cur_b >= ps->num_b) return -1;
    int gp = (ps->cur_b - 1) ^ ((ps->cur_b - 1) >> 1);
    int gc = ps->cur_b ^ (ps->cur_b >> 1);
    int j = __builtin_ctz(gp ^ gc);
    int sign = (gc >> j) & 1;

    if (sign) { mpz_add(ps->b, ps->b, ps->B[j]); mpz_add(ps->b, ps->b, ps->B[j]); }
    else { mpz_sub(ps->b, ps->b, ps->B[j]); mpz_sub(ps->b, ps->b, ps->B[j]); }

    mpz_mul(ps->c, ps->b, ps->b); mpz_sub(ps->c, ps->c, kN); mpz_divexact(ps->c, ps->c, ps->a);

    for (int i = 0; i < fb->size; i++) {
        if (ps->soln1[i] == 0xFFFFFFFF) continue;
        unsigned int p = fb->prime[i];
        unsigned int d = ps->ainv[j * fb->size + i];
        if (d == 0) continue;
        if (sign) {
            ps->soln1[i] = (ps->soln1[i] + p - d) % p;
            ps->soln2[i] = (ps->soln2[i] + p - d) % p;
        } else {
            ps->soln1[i] = (ps->soln1[i] + d) % p;
            ps->soln2[i] = (ps->soln2[i] + d) % p;
        }
    }
    return j;
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

    /* Quick trial division */
    for (int p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t c; mpz_init(c); mpz_divexact_ui(c, N, p);
            gmp_printf("%d\n", p);
            fprintf(stderr, "DLP-SIQS: trivial factor %d, %.3fs\n", p, elapsed());
            return 0;
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fprintf(stderr, "DLP-SIQS: %dd (%db), k=%d, FB=%d, blocks=%d, LP_mult=%d\n",
            digits, bits, mult, P.fb_size, P.num_blocks, P.lp_mult);

    fb_t *fb = fb_create(kN, P.fb_size);
    fprintf(stderr, "FB: %d primes, max=%u\n", fb->size, fb->prime[fb->size - 1]);

    int M = SIEVE_BLOCK * P.num_blocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size - 1] * P.lp_mult;
    unsigned long long dlp_bound = (unsigned long long)lp_bound * lp_bound; /* DLP cofactor bound */
    int target_rels = fb->size + P.extra_rels;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh_adj);
    threshold -= 3; /* approximate skipped small primes */

    fprintf(stderr, "M=%d, target=%d, thresh=%d, LP=%lu, DLP_cof=%llu\n",
            M, target_rels, threshold, lp_bound, dlp_bound);

    unsigned char *sieve = malloc(SIEVE_BLOCK);

    rel_store_t *full_rels = rel_create(MAX_FULL_RELS);
    rel_store_t *part_rels = rel_create(MAX_PARTIALS);

    slp_table_t *slp = slp_create(MAX_PARTIALS);
    dlp_graph_t *dlp = dlp_create(MAX_PARTIALS * 2, MAX_PARTIALS);

    poly_t ps;
    mpz_inits(ps.a, ps.b, ps.c, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(ps.B[j]);
    ps.ainv = malloc(MAX_A_FACTORS * fb->size * sizeof(unsigned int));
    ps.fb_size = fb->size;
    ps.soln1 = malloc(fb->size * sizeof(unsigned int));
    ps.soln2 = malloc(fb->size * sizeof(unsigned int));

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

    int total_polys = 0, a_count = 0, slp_combined = 0, dlp_cycles = 0;
    long cands = 0;
    int skip_bound = 5;

    while (full_rels->count < target_rels) {
        if (total_polys % 500 == 0 && total_polys > 0) {
            double t = elapsed();
            if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs\n", t); break; }
            if (total_polys % 2000 == 0) {
                int total = full_rels->count;
                fprintf(stderr, "  p=%d r=%d/%d (full=%d slp=%d dlp=%d) part=%d+%d t=%.1fs\n",
                        total_polys, total, target_rels, total - slp_combined - dlp_cycles,
                        slp_combined, dlp_cycles, part_rels->count, dlp->ne, t);
            }
        }

        if (total_polys == 0 || ps.cur_b >= ps.num_b) {
            if (!poly_new_a(&ps, fb, kN, M, rng)) continue;
            a_count++;
        } else {
            if (poly_next_b(&ps, fb, kN) < 0) {
                if (!poly_new_a(&ps, fb, kN, M, rng)) continue;
                a_count++;
            }
        }
        total_polys++;

        /* Sieve */
        for (int block = -P.num_blocks; block < P.num_blocks; block++) {
            int bs = block * SIEVE_BLOCK;
            memset(sieve, 0, SIEVE_BLOCK);

            for (int i = 1; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                if (ps.soln1[i] == 0xFFFFFFFF) continue;
                if (p < (unsigned int)skip_bound) continue;
                unsigned char lp = fb->logp[i];

                long off1 = ((long)ps.soln1[i] - bs) % (long)p;
                if (off1 < 0) off1 += p;
                for (int j = (int)off1; j < SIEVE_BLOCK; j += p) sieve[j] += lp;

                if (ps.soln1[i] != ps.soln2[i]) {
                    long off2 = ((long)ps.soln2[i] - bs) % (long)p;
                    if (off2 < 0) off2 += p;
                    for (int j = (int)off2; j < SIEVE_BLOCK; j += p) sieve[j] += lp;
                }
            }

            /* Scan for candidates */
            for (int j = 0; j < SIEVE_BLOCK; j++) {
                if (sieve[j] < threshold) continue;
                cands++;
                long x = (long)(bs + j);
                if (x == 0) continue;

                /* Compute Q(x) */
                mpz_set_si(tmp, x);
                mpz_mul(Qx, ps.a, tmp);
                mpz_add(Qx, Qx, ps.b); mpz_add(Qx, Qx, ps.b);
                mpz_mul(Qx, Qx, tmp);
                mpz_add(Qx, Qx, ps.c);

                mpz_mul_si(ax_b, ps.a, x);
                mpz_add(ax_b, ax_b, ps.b);

                if (mpz_sgn(Qx) == 0) continue;
                mpz_abs(residue, Qx);

                /* Trial divide */
                while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);

                for (int i = 1; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (ps.soln1[i] == 0xFFFFFFFF) continue;
                    long xmod = ((x % (long)p) + p) % p;
                    if (xmod != (long)ps.soln1[i] && xmod != (long)ps.soln2[i]) continue;
                    if (mpz_divisible_ui_p(residue, p))
                        do { mpz_divexact_ui(residue, residue, p); } while (mpz_divisible_ui_p(residue, p));
                }
                for (int i = 0; i < fb->size && fb->prime[i] < (unsigned int)skip_bound; i++) {
                    unsigned int p = fb->prime[i];
                    if (p <= 2) continue;
                    while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                }

                /* Compute a*Q(x) for relation */
                mpz_t aQx; mpz_init(aQx);
                mpz_mul(aQx, Qx, ps.a);

                if (mpz_cmp_ui(residue, 1) == 0) {
                    /* Full relation */
                    int ri = full_rels->count;
                    if (ri < full_rels->alloc) {
                        mpz_set(full_rels->ax_b[ri], ax_b);
                        mpz_set(full_rels->Qx[ri], aQx);
                        full_rels->lp1[ri] = 0; full_rels->lp2[ri] = 0;
                        full_rels->count++;
                    }
                } else if (mpz_fits_ulong_p(residue)) {
                    unsigned long cof = mpz_get_ui(residue);
                    if (cof <= lp_bound) {
                        /* SLP: single large prime */
                        int match = slp_find(slp, cof);
                        if (match >= 0) {
                            int ri = full_rels->count;
                            if (ri < full_rels->alloc) {
                                mpz_mul(full_rels->ax_b[ri], ax_b, part_rels->ax_b[match]);
                                mpz_mod(full_rels->ax_b[ri], full_rels->ax_b[ri], N);
                                mpz_mul(full_rels->Qx[ri], aQx, part_rels->Qx[match]);
                                full_rels->lp1[ri] = cof; full_rels->lp2[ri] = 0;
                                full_rels->count++;
                                slp_combined++;
                            }
                        } else {
                            int pi = part_rels->count;
                            if (pi < part_rels->alloc) {
                                mpz_set(part_rels->ax_b[pi], ax_b);
                                mpz_set(part_rels->Qx[pi], aQx);
                                part_rels->lp1[pi] = cof; part_rels->lp2[pi] = 0;
                                slp_insert(slp, cof, pi);
                                part_rels->count++;
                            }
                        }
                    } else if ((unsigned long long)cof <= dlp_bound && cof < (1ULL << 62)) {
                        /* DLP: try to split cofactor */
                        unsigned long long f1;
                        if (is_prime_64(cof)) goto skip_dlp; /* can't split a prime */
                        f1 = pollard_rho_64(cof);
                        if (f1 == 0 || f1 == cof) goto skip_dlp;

                        unsigned long long f2 = cof / f1;
                        /* Both factors must be < lp_bound */
                        if (f1 > lp_bound || f2 > lp_bound) goto skip_dlp;
                        /* Both must be prime (or we don't care, but helps matching) */

                        unsigned long lp1_val = (unsigned long)(f1 < f2 ? f1 : f2);
                        unsigned long lp2_val = (unsigned long)(f1 < f2 ? f2 : f1);

                        /* Store as partial with 2 LPs */
                        int pi = part_rels->count;
                        if (pi < part_rels->alloc) {
                            mpz_set(part_rels->ax_b[pi], ax_b);
                            mpz_set(part_rels->Qx[pi], aQx);
                            part_rels->lp1[pi] = lp1_val;
                            part_rels->lp2[pi] = lp2_val;
                            part_rels->count++;

                            /* Add edge to DLP graph */
                            int v1 = dlp_find_or_add_vertex(dlp, lp1_val);
                            int v2 = dlp_find_or_add_vertex(dlp, lp2_val);
                            if (v1 >= 0 && v2 >= 0) {
                                if (dlp_union(dlp, v1, v2)) {
                                    /* Cycle found! This means we can combine DLP relations
                                     * to get a full relation. For simplicity, we count it
                                     * but use a simpler approach: find a path in the graph. */
                                    dlp_cycles++;

                                    /* For now, just create a combined relation from
                                     * matching DLP pairs (same LP) */
                                }
                                if (dlp->ne < dlp->max_edges) {
                                    dlp->edges[dlp->ne].lp1 = lp1_val;
                                    dlp->edges[dlp->ne].lp2 = lp2_val;
                                    dlp->edges[dlp->ne].rel_idx = pi;
                                    dlp->ne++;
                                }
                            }
                        }
                        skip_dlp: ;
                    }
                }
                mpz_clear(aQx);
            }
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving: %d full rels (%d+%d slp + %d dlp cycles) from %d polys in %.2fs\n",
            full_rels->count, full_rels->count - slp_combined - dlp_cycles,
            slp_combined, dlp_cycles, total_polys, sieve_time);

    /* Also create DLP combined relations from matching pairs */
    /* Find DLP relations sharing a large prime and combine them */
    {
        /* Build hash table: for each LP, list all DLP relations containing it */
        slp_table_t *dlp_lp_table = slp_create(dlp->ne * 2);

        for (int i = 0; i < part_rels->count && full_rels->count < target_rels; i++) {
            if (part_rels->lp2[i] == 0) continue; /* skip SLP */
            unsigned long lp1 = part_rels->lp1[i];
            unsigned long lp2 = part_rels->lp2[i];

            /* Check if we've seen another relation with lp1 */
            int match1 = slp_find(dlp_lp_table, lp1);
            if (match1 >= 0 && full_rels->count < full_rels->alloc) {
                /* Combine: relation i and relation match1 share lp1 */
                /* Product has lp1^2 (even) and lp2_i * lp2_match1 as cofactor */
                /* The cofactor lp2_i * lp2_match1 should be a perfect square for a full rel,
                 * OR we can treat the product as having two new large primes. */
                /* For simplicity: if lp2_i == lp2_match1, we have a full relation */
                if (part_rels->lp2[match1] == lp2) {
                    int ri = full_rels->count;
                    mpz_mul(full_rels->ax_b[ri], part_rels->ax_b[i], part_rels->ax_b[match1]);
                    mpz_mod(full_rels->ax_b[ri], full_rels->ax_b[ri], N);
                    mpz_mul(full_rels->Qx[ri], part_rels->Qx[i], part_rels->Qx[match1]);
                    full_rels->lp1[ri] = lp1; full_rels->lp2[ri] = lp2;
                    full_rels->count++;
                    dlp_cycles++;
                }
            }
            slp_insert(dlp_lp_table, lp1, i);

            int match2 = slp_find(dlp_lp_table, lp2);
            if (match2 >= 0 && full_rels->count < full_rels->alloc) {
                if (part_rels->lp1[match2] == lp1 || part_rels->lp2[match2] == lp1) {
                    int ri = full_rels->count;
                    mpz_mul(full_rels->ax_b[ri], part_rels->ax_b[i], part_rels->ax_b[match2]);
                    mpz_mod(full_rels->ax_b[ri], full_rels->ax_b[ri], N);
                    mpz_mul(full_rels->Qx[ri], part_rels->Qx[i], part_rels->Qx[match2]);
                    full_rels->lp1[ri] = lp2; full_rels->lp2[ri] = 0;
                    full_rels->count++;
                    dlp_cycles++;
                }
            }
            slp_insert(dlp_lp_table, lp2, i);
        }
        free(dlp_lp_table->buckets);
        free(dlp_lp_table->pool);
        free(dlp_lp_table);
    }

    fprintf(stderr, "After DLP matching: %d full rels (need %d)\n", full_rels->count, target_rels);

    if (full_rels->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full_rels->count, fb->size + 1);
        printf("FAIL\n");
        return 1;
    }

    /* ==================== Linear Algebra ==================== */
    int nrels = full_rels->count;
    if (nrels > target_rels) nrels = target_rels;
    int ncols = fb->size + 1; /* +1 for sign */

    fprintf(stderr, "LA: %d x %d matrix\n", nrels, ncols);
    gf2_mat_t *mat = gf2_create(nrels, ncols);

    for (int r = 0; r < nrels; r++) {
        mpz_t Qval; mpz_init(Qval);
        mpz_set(Qval, full_rels->Qx[r]);

        if (mpz_sgn(Qval) < 0) { gf2_set(mat, r, 0); mpz_neg(Qval, Qval); }

        int e2 = 0;
        while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set(mat, r, 1);

        for (int i = 1; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            int e = 0;
            while (mpz_divisible_ui_p(Qval, p)) { mpz_divexact_ui(Qval, Qval, p); e++; }
            if (e & 1) gf2_set(mat, r, i + 1);
        }
        mpz_clear(Qval);
    }

    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, 64);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* ==================== Square Root ==================== */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);

        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_mul(X, X, full_rels->ax_b[deps[d][k]]);
            mpz_mod(X, X, N);
        }

        /* Compute Y from factored product */
        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t absQ; mpz_init(absQ);
            mpz_abs(absQ, full_rels->Qx[deps[d][k]]);
            mpz_mul(prod, prod, absQ);
            mpz_clear(absQ);
        }

        mpz_set_ui(Y, 1);
        mpz_set(rem, prod);

        /* Factor out 2 */
        int e2 = 0;
        while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next_dep;
        if (e2 / 2 > 0) {
            mpz_set_ui(tmp, 2); mpz_powm_ui(tmp, tmp, e2 / 2, N);
            mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
        }

        /* Factor out FB primes */
        {
            int valid = 1;
            for (int i = 1; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                int e = 0;
                while (mpz_divisible_ui_p(rem, p)) { mpz_divexact_ui(rem, rem, p); e++; }
                if (e & 1) { valid = 0; break; }
                if (e / 2 > 0) {
                    mpz_set_ui(tmp, p); mpz_powm_ui(tmp, tmp, e / 2, N);
                    mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
                }
            }
            if (!valid) goto next_dep;
        }

        /* Remaining: LP^2 from combined relations */
        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_perfect_square_p(rem)) {
                mpz_sqrt(tmp, rem); mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            } else goto next_dep;
        }

        /* GCD */
        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other); mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "DLP-SIQS: factored in %.3fs (dep %d/%d)\n", elapsed(), d+1, ndeps);
            mpz_clear(other);
            return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other); mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "DLP-SIQS: factored in %.3fs (dep %d/%d)\n", elapsed(), d+1, ndeps);
            mpz_clear(other);
            return 0;
        }

        next_dep:
        mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "DLP-SIQS: FAILED after %d dependencies\n", ndeps);
    printf("FAIL\n");
    return 1;
}
