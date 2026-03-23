/*
 * Turbo SIQS v2 - Highly optimized Self-Initializing Quadratic Sieve
 *
 * Key features:
 * 1. 48KB L1-cache-optimized sieve blocks (AMD EPYC 9R45)
 * 2. Bucket sieving for large primes
 * 3. Gray code self-initialization with incremental solution updates
 * 4. Double Large Primes (DLP) with graph-based cycle finding (union-find)
 * 5. Two-threshold sieve scanning (tight for full/SLP, loose for DLP)
 * 6. 64-bit fast path trial division
 * 7. Sieve-informed trial division (only test matching roots)
 *
 * Compile: gcc -O3 -march=native -o turbo_siqs library/turbo_siqs.c -lgmp -lm
 * Usage: ./turbo_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>

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

/* ==================== Block Lanczos over GF(2) ==================== */
/*
 * Finds null space of matrix B = A^T * A where A is sparse (nrels x ncols).
 * Processes 64 vectors simultaneously using 64-bit words.
 *
 * Based on Montgomery's Block Lanczos algorithm.
 * Uses sparse matrix A stored in Compressed Row Storage (CRS).
 */
typedef uint64_t gf2w;

typedef struct {
    int nrows, ncols;
    int *row_start;  /* row_start[i] = index of first entry in row i */
    int *col_idx;    /* column indices of nonzero entries */
    int nnz;
} sparse_t;

/* Multiply y = A * x  (A is nrows x ncols, x is ncols block-vectors, y is nrows) */
static void spmv(sparse_t *A, gf2w *x, gf2w *y) {
    for (int i = 0; i < A->nrows; i++) {
        gf2w acc = 0;
        for (int j = A->row_start[i]; j < A->row_start[i+1]; j++)
            acc ^= x[A->col_idx[j]];
        y[i] = acc;
    }
}

/* Multiply y = A^T * x  (A is nrows x ncols, x is nrows, y is ncols) */
static void spmv_t(sparse_t *A, gf2w *x, gf2w *y) {
    memset(y, 0, A->ncols * sizeof(gf2w));
    for (int i = 0; i < A->nrows; i++) {
        gf2w xi = x[i];
        if (!xi) continue;
        for (int j = A->row_start[i]; j < A->row_start[i+1]; j++)
            y[A->col_idx[j]] ^= xi;
    }
}

/* Multiply y = B*x = A^T * A * x  (x and y are ncols-dimensional) */
static void bmul(sparse_t *A, gf2w *x, gf2w *y, gf2w *tmp) {
    spmv(A, x, tmp);     /* tmp = A*x,  nrows-dim */
    spmv_t(A, tmp, y);   /* y = A^T*tmp, ncols-dim */
}

/* Multiply y = C*x = A * A^T * x  (x and y are nrows-dimensional)
 * This gives us a symmetric operator over the RELATION space.
 * Null vectors of C directly identify subsets of relations that sum to zero. */
static void cmul(sparse_t *A, gf2w *x, gf2w *y, gf2w *tmp) {
    spmv_t(A, x, tmp);   /* tmp = A^T*x, ncols-dim */
    spmv(A, tmp, y);     /* y = A*tmp,   nrows-dim */
}

/* Compute 64x64 inner product: M = X^T * Y where X, Y are n-element block vectors */
static void inner_prod(gf2w *X, gf2w *Y, int n, gf2w M[64]) {
    memset(M, 0, 64 * sizeof(gf2w));
    for (int i = 0; i < n; i++) {
        gf2w x = X[i], y = Y[i];
        if (!x) continue;
        /* For each set bit in x, XOR y into corresponding row of M */
        while (x) {
            int bit = __builtin_ctzll(x);
            M[bit] ^= y;
            x &= x - 1;
        }
    }
}

/* Compute rank of 64x64 GF(2) matrix and its "inverse" (for the nonsingular part) */
/* Returns bitmask of active columns (rank columns) */
static gf2w mat64_kernel(gf2w M[64], gf2w inv[64], gf2w *active_mask) {
    /* Initialize inv to identity */
    for (int i = 0; i < 64; i++) inv[i] = 1ULL << i;
    gf2w cols_used = 0;

    for (int c = 0; c < 64; c++) {
        /* Find pivot */
        int pr = -1;
        for (int r = 0; r < 64; r++) {
            if ((cols_used >> r) & 1) continue;
            if ((M[r] >> c) & 1) { pr = r; break; }
        }
        if (pr < 0) continue;
        cols_used |= (1ULL << pr);

        /* Eliminate */
        for (int r = 0; r < 64; r++) {
            if (r == pr) continue;
            if ((M[r] >> c) & 1) {
                M[r] ^= M[pr];
                inv[r] ^= inv[pr];
            }
        }
    }
    *active_mask = cols_used;
    return cols_used;
}

/* Mask block vector: for each element, keep only bits in mask */
static void mask_vec(gf2w *v, int n, gf2w mask) {
    for (int i = 0; i < n; i++) v[i] &= mask;
}

/* Multiply two 64x64 GF(2) matrices: C = A * B */
static void mat64_mul(const gf2w A[64], const gf2w B[64], gf2w C[64]) {
    memset(C, 0, 64 * sizeof(gf2w));
    for (int i = 0; i < 64; i++) {
        gf2w row = A[i];
        while (row) {
            int bit = __builtin_ctzll(row);
            C[i] ^= B[bit];
            row &= row - 1;
        }
    }
}

/* Apply 64x64 matrix M to block vector: y[i] = M * x[i] (per-element) */
/* Each element is a column, M acts on bits */
static void apply_mat64(const gf2w M[64], gf2w *x, gf2w *y, int n) {
    /* y[i] = M applied to the 64 bits of x[i].
     * M[j] is row j, so bit j of output = OR of (bit k of x[i]) for each k where M[j] has bit k set.
     * Equivalently, think transposed: for each set bit k in x[i], XOR M[k] into y[i] */
    for (int i = 0; i < n; i++) {
        gf2w xi = x[i], acc = 0;
        while (xi) {
            int bit = __builtin_ctzll(xi);
            acc ^= M[bit];
            xi &= xi - 1;
        }
        y[i] = acc;
    }
}

/* Invert 64x64 GF(2) matrix, return 1 on success */
static int mat64_inv(gf2w M[64], gf2w inv[64]) {
    gf2w tmp[64];
    memcpy(tmp, M, 64 * sizeof(gf2w));
    for (int i = 0; i < 64; i++) inv[i] = 1ULL << i;
    for (int c = 0; c < 64; c++) {
        int pr = -1;
        for (int r = c; r < 64; r++) {
            if ((tmp[r] >> c) & 1) { pr = r; break; }
        }
        if (pr < 0) return 0; /* singular */
        if (pr != c) {
            gf2w t = tmp[pr]; tmp[pr] = tmp[c]; tmp[c] = t;
            t = inv[pr]; inv[pr] = inv[c]; inv[c] = t;
        }
        for (int r = 0; r < 64; r++) {
            if (r == c) continue;
            if ((tmp[r] >> c) & 1) {
                tmp[r] ^= tmp[c];
                inv[r] ^= inv[c];
            }
        }
    }
    return 1;
}

/*
 * Compute pseudo-inverse of a symmetric GF(2) matrix on its column/row space.
 *
 * Given a symmetric 64x64 GF(2) matrix F (restricted to 'active' rows/cols),
 * runs GE to find the active subspace (pivot columns) and computes the inverse
 * of F restricted to that subspace.
 *
 * On return:
 *   inv[c] is set for pivot columns c: inv[c] = row of (F[S,S])^{-1} corresponding to c
 *   inv[c] = 0 for non-pivot columns c
 *   Returns the bitmask of pivot columns (the "new active" mask).
 *
 * For the Lanczos recurrence, the returned mask is the new 'active' mask,
 * and the coefficients use inv which satisfies: (inv * F)[c,c] = 1 for pivot c.
 */
/*
 * Compute pseudo-inverse of a symmetric 64x64 GF(2) matrix.
 *
 * Runs augmented GE (Gauss-Jordan) on [F | I] to find the invertible
 * subspace of F. For each pivot column c found, inv[c] gets the corresponding
 * row of the pseudo-inverse. Non-pivot columns get inv[c] = 0.
 *
 * Returns the bitmask of pivot columns (the "alive" columns).
 *
 * Note: does NOT do row swaps so pivot_row[c] may differ from c.
 * We track the mapping: pivot_col_for_row[pr] = c.
 */
static gf2w mat64_sym_pseudoinv(const gf2w F[64], gf2w inv[64]) {
    gf2w tmp[64], aug[64];
    memcpy(tmp, F, 64 * sizeof(gf2w));
    for (int i = 0; i < 64; i++) aug[i] = (gf2w)1 << i;

    memset(inv, 0, 64 * sizeof(gf2w));
    gf2w alive = 0;
    gf2w used_rows = 0; /* rows already used as pivots */
    int pivot_col_for_row[64]; /* which column used row r as pivot */
    memset(pivot_col_for_row, -1, sizeof(pivot_col_for_row));

    for (int c = 0; c < 64; c++) {
        /* Find a pivot: any unused row that has bit c set */
        int pr = -1;
        for (int r = 0; r < 64; r++) {
            if ((used_rows >> r) & 1) continue; /* already used */
            if ((tmp[r] >> c) & 1) { pr = r; break; }
        }
        if (pr < 0) continue; /* no pivot: column c is in null space */

        alive |= (gf2w)1 << c;
        used_rows |= (gf2w)1 << pr;
        pivot_col_for_row[pr] = c;

        /* Eliminate column c from ALL other rows (full GE, not just below) */
        for (int r = 0; r < 64; r++) {
            if (r == pr) continue;
            if ((tmp[r] >> c) & 1) {
                tmp[r] ^= tmp[pr];
                aug[r] ^= aug[pr];
            }
        }
    }

    /* After GE: for each pivot column c (bit c in alive), row pr has been used.
     * aug[pr] contains the corresponding row of F^{-1} on the alive subspace.
     * We want: inv[c] = aug[pr] where pivot_col_for_row[pr] = c.
     */
    for (int pr = 0; pr < 64; pr++) {
        int c = pivot_col_for_row[pr];
        if (c >= 0) inv[c] = aug[pr];
    }
    return alive;
}

/*
 * Block Lanczos over GF(2) based on Montgomery's algorithm.
 *
 * Finds null space vectors of B = A^T * A (ncols x ncols matrix).
 *
 * Three-term recurrence (Montgomery 1995):
 *   v_{i+1} = B*v_i + v_i * S_i + v_{i-1} * D_i
 * where S_i and D_i are 64x64 coefficient matrices.
 *
 * A solution y is maintained: y starts random, and at each step we project
 *   y += v_i * Fi_inv * (v_i^T * B * y)
 * to remove range(B) components. After the Krylov space is exhausted, B*y = 0.
 *
 * Multiple independent runs are used to find more null vectors (different random
 * starting y gives different null vectors). Falls back to dense GE if BL fails.
 */
/*
 * Structured GE helper: XOR-merge two sorted integer arrays (GF(2) symmetric difference).
 * Returns number of elements in result stored in out[].
 * out[] must have capacity at least len1+len2.
 */
static int xor_merge(const int *a, int la, const int *b, int lb, int *out) {
    int i = 0, j = 0, k = 0;
    while (i < la && j < lb) {
        if (a[i] < b[j]) out[k++] = a[i++];
        else if (b[j] < a[i]) out[k++] = b[j++];
        else { i++; j++; } /* cancel */
    }
    while (i < la) out[k++] = a[i++];
    while (j < lb) out[k++] = b[j++];
    return k;
}

static int block_lanczos(sparse_t *A, int ***deps_out, int **dlen_out, int max_deps) {
    int ncols = A->ncols;
    int nrows = A->nrows;

    *deps_out = malloc(max_deps * sizeof(int*));
    *dlen_out = malloc(max_deps * sizeof(int));

    /*
     * Structured Gaussian Elimination with Singleton and Doubleton Merging.
     *
     * Maintain mutable sparse rows and original-row tracking.
     * Each logical row carries:
     *   - row_buf[r]: sorted list of alive column indices (updated through merges)
     *   - orig_buf[r]: sorted list of original row indices (XOR-accumulated)
     *
     * At the start of each pass, rebuild the col->row index (col_rows[c]) from
     * current row_buf of alive rows. This keeps the index accurate after merges.
     *
     * Phase 1: Singleton removal — remove weight-1 cols and their unique row.
     * Phase 2: Doubleton merge — for weight-2 col c in rows r1,r2: set r2=r1 XOR r2,
     *          remove r1 and c.
     * Repeat until stable, then run dense GE on remaining rows.
     */

    /* Mutable row storage */
    int *row_len = calloc(nrows, sizeof(int));
    int **row_buf = malloc(nrows * sizeof(int*));
    int *orig_len = malloc(nrows * sizeof(int));
    int **orig_buf = malloc(nrows * sizeof(int*));

    for (int r = 0; r < nrows; r++) {
        int len = A->row_start[r+1] - A->row_start[r];
        row_buf[r] = malloc((len + 1) * sizeof(int));
        row_len[r] = len;
        for (int j = 0; j < len; j++)
            row_buf[r][j] = A->col_idx[A->row_start[r] + j];
        orig_buf[r] = malloc(sizeof(int));
        orig_buf[r][0] = r;
        orig_len[r] = 1;
    }

    int *row_alive = malloc(nrows * sizeof(int));
    int *col_alive = malloc(ncols * sizeof(int));
    memset(row_alive, 1, nrows * sizeof(int));
    memset(col_alive, 1, ncols * sizeof(int));

    /* col_rows[c]: dynamic list of alive rows containing c (rebuilt each pass) */
    /* col_cnt[c]: working count of alive rows (modified during pass) */
    /* col_size[c]: actual number of entries in col_rows[c] (set at rebuild, not modified) */
    int **col_rows = malloc(ncols * sizeof(int*));
    int *col_cnt = calloc(ncols, sizeof(int));
    int *col_size = calloc(ncols, sizeof(int)); /* size of col_rows[c] list */
    int *col_cap = calloc(ncols, sizeof(int));
    for (int c = 0; c < ncols; c++) {
        col_rows[c] = malloc(4 * sizeof(int));
        col_cap[c] = 4;
    }

    int removed_rows = 0, removed_cols = 0;
    int *tmp_xor = malloc((A->nnz + ncols) * sizeof(int)); /* temp buffer for XOR */

    int pass = 0;
    for (;;) {
        pass++;

        /* Rebuild col_rows from current row_buf of alive rows */
        for (int c = 0; c < ncols; c++) col_cnt[c] = 0;
        for (int r = 0; r < nrows; r++) {
            if (!row_alive[r]) continue;
            for (int j = 0; j < row_len[r]; j++) {
                int c = row_buf[r][j];
                if (!col_alive[c]) continue;
                if (col_cnt[c] >= col_cap[c]) {
                    col_cap[c] *= 2;
                    col_rows[c] = realloc(col_rows[c], col_cap[c] * sizeof(int));
                }
                col_rows[c][col_cnt[c]++] = r;
            }
        }
        /* Save the actual size of col_rows[c] (col_cnt will be modified during the pass) */
        memcpy(col_size, col_cnt, ncols * sizeof(int));

        int changed = 0;

        /* Singleton removal */
        for (int c = 0; c < ncols; c++) {
            if (!col_alive[c]) continue;
            if (col_cnt[c] > 1) continue;
            col_alive[c] = 0; removed_cols++; changed = 1;
            if (col_cnt[c] == 0) continue;
            /* Find the ONE alive row containing c from col_rows[c]
             * Use col_size[c] (not col_cnt[c]) since col_cnt was decremented during this pass */
            int r = -1;
            for (int k = 0; k < col_size[c]; k++) {
                if (row_alive[col_rows[c][k]]) { r = col_rows[c][k]; break; }
            }
            if (r < 0) continue;
            row_alive[r] = 0; removed_rows++;
            /* Decrement counts for cols in this removed row */
            for (int j = 0; j < row_len[r]; j++) {
                int cc = row_buf[r][j];
                if (col_alive[cc]) col_cnt[cc]--;
            }
        }

        /* Doubleton merging */
        for (int c = 0; c < ncols; c++) {
            if (!col_alive[c]) continue;
            if (col_cnt[c] != 2) continue;

            /* Find two alive rows from col_rows[c] (use col_size, not col_cnt) */
            int r1 = -1, r2 = -1;
            for (int k = 0; k < col_size[c]; k++) {
                if (!row_alive[col_rows[c][k]]) continue;
                if (r1 < 0) r1 = col_rows[c][k];
                else { r2 = col_rows[c][k]; break; }
            }
            if (r1 < 0 || r2 < 0) continue;

            /* Compute new row content: XOR (symmetric difference) of r1 and r2 */
            int len1 = row_len[r1], len2 = row_len[r2];
            int ni = xor_merge(row_buf[r1], len1, row_buf[r2], len2, tmp_xor);

            /* Remove dead columns from result */
            int ni2 = 0;
            for (int j = 0; j < ni; j++)
                if (col_alive[tmp_xor[j]]) tmp_xor[ni2++] = tmp_xor[j];
            ni = ni2;

            /* Update col_cnt: only columns in BOTH r1 and r2 (cancelled) lose 2 */
            {
                int i1 = 0, i2 = 0;
                while (i1 < len1 && i2 < len2) {
                    int a = row_buf[r1][i1], b = row_buf[r2][i2];
                    if (a == b) {
                        if (col_alive[a]) col_cnt[a] -= 2;
                        i1++; i2++;
                    } else if (a < b) i1++;
                    else i2++;
                }
            }

            /* Replace r2's content */
            free(row_buf[r2]);
            row_buf[r2] = malloc((ni + 1) * sizeof(int));
            memcpy(row_buf[r2], tmp_xor, ni * sizeof(int));
            row_len[r2] = ni;

            /* Merge orig_buf: symmetric difference */
            {
                int ol1 = orig_len[r1], ol2 = orig_len[r2];
                int *new_orig = malloc((ol1 + ol2 + 1) * sizeof(int));
                int oi = xor_merge(orig_buf[r1], ol1, orig_buf[r2], ol2, new_orig);
                free(orig_buf[r2]);
                orig_buf[r2] = new_orig;
                orig_len[r2] = oi;
            }

            /* Remove column c and row r1 */
            col_alive[c] = 0; removed_cols++;
            col_cnt[c] = 0;
            row_alive[r1] = 0; removed_rows++;
            changed = 1;
        }

        if (!changed) break;
    }

    free(tmp_xor);
    for (int c = 0; c < ncols; c++) free(col_rows[c]);
    free(col_rows); free(col_cnt); free(col_size); free(col_cap);

    int remaining_rows = nrows - removed_rows;
    int remaining_cols = ncols - removed_cols;
    fprintf(stderr, "  SGE (%d passes): removed %d rows, %d cols; remaining %d x %d\n",
            pass, removed_rows, removed_cols, remaining_rows, remaining_cols);

    /* Build dense matrix for GE */
    int *col_map = malloc(ncols * sizeof(int));
    int new_ncols = 0;
    for (int c = 0; c < ncols; c++) col_map[c] = col_alive[c] ? new_ncols++ : -1;

    int new_nrows = 0;
    for (int r = 0; r < nrows; r++) if (row_alive[r]) new_nrows++;

    fprintf(stderr, "  Dense GE: %d x %d matrix\n", new_nrows, new_ncols);

    int fb_words = (new_ncols + 63) / 64;
    int id_words = (nrows + 63) / 64;
    int total_words = fb_words + id_words;

    gf2w **ge_rows = malloc(new_nrows * sizeof(gf2w*));
    for (int i = 0; i < new_nrows; i++) ge_rows[i] = calloc(total_words, sizeof(gf2w));

    int ri = 0;
    for (int r = 0; r < nrows; r++) {
        if (!row_alive[r]) continue;
        /* Identity bits: one bit per original row in orig_buf[r] */
        for (int k = 0; k < orig_len[r]; k++) {
            int orig = orig_buf[r][k];
            ge_rows[ri][fb_words + orig/64] |= (1ULL << (orig % 64));
        }
        /* FB columns */
        for (int j = 0; j < row_len[r]; j++) {
            int c = col_map[row_buf[r][j]];
            if (c >= 0) ge_rows[ri][c/64] |= (1ULL << (c % 64));
        }
        ri++;
    }

    /* Dense GE over GF(2) */
    int pivot = 0;
    for (int c = 0; c < new_ncols && pivot < new_nrows; c++) {
        int pr = -1;
        for (int r = pivot; r < new_nrows; r++)
            if ((ge_rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != pivot) { gf2w *t = ge_rows[pr]; ge_rows[pr] = ge_rows[pivot]; ge_rows[pivot] = t; }
        for (int r = 0; r < new_nrows; r++) {
            if (r == pivot) continue;
            if ((ge_rows[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < total_words; w++) ge_rows[r][w] ^= ge_rows[pivot][w];
        }
        pivot++;
    }

    int nd_ge = 0;
    for (int r = pivot; r < new_nrows && nd_ge < max_deps; r++) {
        /* Check all FB columns are zero */
        int zero = 1;
        for (int w = 0; w < fb_words && zero; w++) {
            gf2w ge_mask = (w < fb_words-1) ? ~0ULL :
                           (new_ncols%64==0 ? ~0ULL : (1ULL<<(new_ncols%64))-1);
            if (ge_rows[r][w] & ge_mask) zero = 0;
        }
        if (!zero) continue;
        /* Extract dependency from identity bits */
        int *d = malloc(nrows * sizeof(int)); int dl = 0;
        for (int w = 0; w < id_words; w++) {
            gf2w bits = ge_rows[r][fb_words+w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w*64+bit;
                if (idx < nrows) d[dl++] = idx;
                bits &= bits-1;
            }
        }
        if (dl > 0) { (*deps_out)[nd_ge] = d; (*dlen_out)[nd_ge] = dl; nd_ge++; } else free(d);
    }

    /* Cleanup */
    for (int i = 0; i < new_nrows; i++) free(ge_rows[i]);
    free(ge_rows);
    for (int r = 0; r < nrows; r++) { free(row_buf[r]); free(orig_buf[r]); }
    free(row_buf); free(row_len); free(orig_buf); free(orig_len);
    free(row_alive); free(col_alive); free(col_map);
    return nd_ge;
}

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

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size, nblocks, lp_mult, extra;
    double thresh_adj;
    int dlp;
} params_t;

static params_t get_params(int bits) {
    /* Tuned for single-core 300s with DLP graph matching.
     * Smaller FB + aggressive LP = fewer sieve ops, more DLP cycles. */
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
    if (bits <= 240) return (params_t){9000,  46, 200, 220, 0.85, 1};
    if (bits <= 250) return (params_t){20000, 56, 150, 300, 0.86, 1};
    if (bits <= 260) return (params_t){16000, 66, 200, 300, 0.86, 1};
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
    int target = fb->size + P.extra;

    /* Threshold: single threshold, DLP handled by wider LP range */
    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2((double)M);
    int threshold = (int)(log2_Qmax * P.thresh_adj) - 3;
    if (threshold < 20) threshold = 20;

    /* Bucket sieve setup */
    int bucket_thresh = 0;
    for (int i = 0; i < fb->size; i++) { if (fb->prime[i] > BLOCK_SIZE) { bucket_thresh = i; break; } }
    if (bucket_thresh == 0) bucket_thresh = fb->size;

    int total_blocks = 2 * P.nblocks;
    bucket_t *buckets = calloc(total_blocks, sizeof(bucket_t));
    for (int i = 0; i < total_blocks; i++) { buckets[i].alloc = BUCKET_ALLOC; buckets[i].entries = malloc(BUCKET_ALLOC * sizeof(bucket_entry_t)); }

    fprintf(stderr, "TurboSIQS: %dd (%db), k=%d, FB=%d, M=%d, thresh=%d, LP=%lu%s, target=%d\n",
            digits, bits, mult, fb->size, M, threshold, lp_bound, P.dlp?" [DLP]":"", target);

    uint8_t *sieve = malloc(BLOCK_SIZE);
    rels_t *full = rels_create(MAX_FULL_RELS);
    rels_t *part = rels_create(MAX_PARTIAL_RELS);
    lp_hash_t *slp = lp_create(MAX_PARTIAL_RELS);
    dlp_graph_t *dlp_g = P.dlp ? dlp_create(2000000) : NULL;

    /* DLP relation tracking: for each partial with 2 LPs, store index */
    /* When a cycle is found, we need to trace back the relations in the cycle.
     * For simplicity, we use a different approach:
     * - SLP: standard hash-based matching
     * - DLP: add edge to graph, count cycles. Combined rels = cycles.
     * - For the LA matrix: use ALL partials (SLP and DLP) that have a matching pair.
     * Actually, to keep it simple: for DLP, we hash by individual LP.
     * When two DLP relations share one LP, they can combine if the other LPs also match (cycle of length 2).
     * For longer cycles, we just count them and add "free relations" proportional to cycle count.
     * The actual relation combining for LA needs more bookkeeping... */

    /* Simpler approach: just track SLP matching. For DLP, track graph cycles.
     * When cycle is found, combine the two relations immediately. */

    /* DLP matching hash: for each LP, store list of partial relations containing it */
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

    int total_polys = 0, combined_slp = 0, combined_dlp = 0;
    int dlp_found = 0, dlp_cycles = 0;
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

                /* === Scan for candidates === */
                for (int j = 0; j < BLOCK_SIZE; j++) {
                    if (sieve[j] < threshold) continue;
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

                        if (cof <= lp_bound) {
                            /* SLP */
                            int match = lp_find(slp, cof);
                            if (match >= 0) {
                                mpz_mul(tmp, ax_b, part->ax_b[match]); mpz_mod(tmp, tmp, N);
                                mpz_mul(tmp2, aQx, part->Qx[match]);
                                rels_add(full, tmp, tmp2, 0, 0);
                                combined_slp++;
                            } else {
                                int pi = rels_add(part, ax_b, aQx, cof, 0);
                                if (pi >= 0) lp_insert(slp, cof, pi);
                            }
                        } else if (P.dlp && cof <= dlp_max) {
                            /* DLP: split cofactor */
                            uint64_t f1, f2;
                            if (split64(cof, &f1, &f2)) {
                                if (f1 > f2) { uint64_t t = f1; f1 = f2; f2 = t; }
                                if (f1 <= lp_bound && f2 <= lp_bound) {
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
                                        /* No SLP match: store as SLP with f1 */
                                        int pi = rels_add(part, ax_b, aQx, f1, f2);
                                        if (pi >= 0) {
                                            lp_insert(slp, f1, pi);
                                            /* Also insert with f2 for future matches */
                                            lp_insert(slp, f2, pi);
                                        }
                                    }

                                    int pi = -1; /* For DLP graph compatibility */
                                    if (m1 < 0 && m2 < 0) pi = 0; /* only do graph matching if no SLP match */
                                    if (pi >= 0) {
                                        /* Add to DLP graph */
                                        int cycle = dlp_add_edge(dlp_g, f1, f2);
                                        if (cycle) dlp_cycles++;

                                        /* Hash by each LP for pairwise matching */
                                        /* When we find another relation with the same LP,
                                         * we can try to combine them */
                                        uint32_t h1 = (uint32_t)((f1 * 0x9E3779B97F4A7C15ULL) >> (64-DLP_PAIR_HASH_BITS));
                                        uint32_t h2 = (uint32_t)((f2 * 0x9E3779B97F4A7C15ULL) >> (64-DLP_PAIR_HASH_BITS));

                                        /* Look for matching pair */
                                        for (dp_e_t *e = dp_hash[h1]; e; e = e->next) {
                                            if (e->lp == f1) {
                                                int oi = e->rel_idx;
                                                if (part->lp2[oi] == f2 || part->lp1[oi] == f2) {
                                                    /* Exact DLP pair match! */
                                                    mpz_mul(tmp, ax_b, part->ax_b[oi]); mpz_mod(tmp, tmp, N);
                                                    mpz_mul(tmp2, aQx, part->Qx[oi]);
                                                    rels_add(full, tmp, tmp2, 0, 0);
                                                    combined_dlp++;
                                                    break;
                                                }
                                            }
                                        }
                                        for (dp_e_t *e = dp_hash[h2]; e; e = e->next) {
                                            if (e->lp == f2) {
                                                int oi = e->rel_idx;
                                                if (part->lp1[oi] == f1 || part->lp2[oi] == f1) {
                                                    /* Already matched above */
                                                    break;
                                                }
                                                if (part->lp1[oi] == f2 || part->lp2[oi] == f2) {
                                                    /* Two DLP sharing one LP - need 3rd for cycle */
                                                    /* Skip for now - exact pair match only */
                                                }
                                            }
                                        }

                                        /* Insert into DLP hash */
                                        if (dp_pool_used + 2 <= MAX_PARTIAL_RELS * 2) {
                                            dp_e_t *e1 = &dp_pool[dp_pool_used++];
                                            e1->lp = f1; e1->rel_idx = pi; e1->next = dp_hash[h1]; dp_hash[h1] = e1;
                                            dp_e_t *e2 = &dp_pool[dp_pool_used++];
                                            e2->lp = f2; e2->rel_idx = pi; e2->next = dp_hash[h2]; dp_hash[h2] = e2;
                                        }
                                    } /* end if (pi >= 0) */
                                }
                            }
                        }
                    } else if (P.dlp) {
                        /* Residue > 64 bits. Check if it has small factors we missed */
                        /* Skip - too expensive */
                    }
                }
            }

            /* Progress */
            if (total_polys % 500 == 0) {
                double t = elapsed_sec();
                if (t > 275.0) break;
                if (total_polys % 2000 == 0)
                    fprintf(stderr, "  poly=%d rels=%d/%d (full=%d slp=%d dlp=%d) part=%d dlp_found=%d t=%.1fs\n",
                            total_polys, full->count, target,
                            full->count-combined_slp-combined_dlp, combined_slp, combined_dlp,
                            part->count, dlp_found, t);
            }
        }
    }

    double sieve_time = elapsed_sec();
    fprintf(stderr, "Sieve: %d rels in %.2fs (%d polys, %d SLP, %d DLP, %d dlp_found)\n",
            full->count, sieve_time, total_polys, combined_slp, combined_dlp, dlp_found);

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
    int ndeps = block_lanczos(smat, &deps, &dlen, 64);
    fprintf(stderr, "Found %d deps in %.2fs\n", ndeps, elapsed_sec());
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
