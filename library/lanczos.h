/*
 * lanczos.h - Proper Block Lanczos for GF(2) null space computation
 *
 * Based on Montgomery, "A Block Lanczos Algorithm for Finding
 * Dependencies over GF(2)", EUROCRYPT 1995.
 *
 * For a sparse matrix B over GF(2) (nrows x ncols), finds vectors
 * in the null space of B^T * B (and hence of B itself).
 *
 * Key advantage over Gaussian elimination:
 * - O(weight * ncols * ncols/64) time (exploits sparsity)
 * - vs O(nrows * ncols^2 / 64) for Gauss on dense representation
 * - For QS matrices (weight ~30, ncols ~10000): ~100x faster
 */

#ifndef LANCZOS_H
#define LANCZOS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef LANCZOS_U64_DEFINED
typedef unsigned long long u64;
#define LANCZOS_U64_DEFINED
#endif

/* Sparse matrix in CSR format */
typedef struct {
    int nrows, ncols;
    int *row_off;  /* row_off[r] = start of row r in col */
    int *col;      /* column indices */
    int nnz;
} sp_mat_t;

/* Build sparse matrix from dense GF(2) matrix (row-major, each row is array of u64) */
static sp_mat_t *sp_from_dense(u64 **dense, int nrows, int ncols) {
    sp_mat_t *A = malloc(sizeof(sp_mat_t));
    A->nrows = nrows; A->ncols = ncols;
    A->row_off = malloc((nrows + 1) * sizeof(int));
    int nnz = 0;
    for (int r = 0; r < nrows; r++) {
        A->row_off[r] = nnz;
        int words = (ncols + 63) / 64;
        for (int w = 0; w < words; w++) {
            u64 bits = dense[r][w];
            nnz += __builtin_popcountll(bits);
        }
    }
    A->row_off[nrows] = nnz;
    A->nnz = nnz;
    A->col = malloc(nnz * sizeof(int));
    int idx = 0;
    for (int r = 0; r < nrows; r++) {
        int words = (ncols + 63) / 64;
        for (int w = 0; w < words; w++) {
            u64 bits = dense[r][w];
            while (bits) {
                int b = __builtin_ctzll(bits);
                A->col[idx++] = w * 64 + b;
                bits &= bits - 1;
            }
        }
    }
    return A;
}

/* y = B * x  (sparse multiply) */
static void sp_mul(sp_mat_t *B, u64 *x, u64 *y) {
    for (int r = 0; r < B->nrows; r++) {
        u64 acc = 0;
        for (int j = B->row_off[r]; j < B->row_off[r + 1]; j++)
            acc ^= x[B->col[j]];
        y[r] = acc;
    }
}

/* y = B^T * x  (transpose multiply) */
static void sp_mul_t(sp_mat_t *B, u64 *x, u64 *y) {
    memset(y, 0, B->ncols * sizeof(u64));
    for (int r = 0; r < B->nrows; r++) {
        u64 xr = x[r];
        if (!xr) continue;
        for (int j = B->row_off[r]; j < B->row_off[r + 1]; j++)
            y[B->col[j]] ^= xr;
    }
}

/* y = B^T * B * x  (symmetric multiply, uses temp buffer) */
static void sp_mul_btb(sp_mat_t *B, u64 *x, u64 *y, u64 *tmp) {
    sp_mul(B, x, tmp);    /* tmp = B * x */
    sp_mul_t(B, tmp, y);  /* y = B^T * tmp = B^T * B * x */
}

/* Compute 64x64 inner product: M[i] = sum_k (x[k] >> i & 1) ? y[k] : 0
 * Equivalent to: M = X^T * Y where X, Y are n x 64 bit matrices */
static void inner64(u64 *x, u64 *y, int n, u64 M[64]) {
    memset(M, 0, 64 * sizeof(u64));
    for (int k = 0; k < n; k++) {
        u64 xk = x[k], yk = y[k];
        if (!xk) continue;
        while (xk) {
            int i = __builtin_ctzll(xk);
            M[i] ^= yk;
            xk &= xk - 1;
        }
    }
}

/* 64x64 GF(2) matrix rank and (pseudo-)inverse.
 * Returns rank. Writes row-echelon form to M, inverse to inv.
 * If rank < 64, the inverse is partial (only defined on the row space). */
static int inv64(u64 M[64], u64 inv[64]) {
    u64 m[64], id[64];
    for (int i = 0; i < 64; i++) { m[i] = M[i]; id[i] = 1ULL << i; }
    int rank = 0;
    int piv_col[64];
    for (int c = 0; c < 64; c++) {
        int pr = -1;
        for (int r = rank; r < 64; r++)
            if ((m[r] >> c) & 1) { pr = r; break; }
        if (pr < 0) continue;
        u64 tmp;
        tmp = m[pr]; m[pr] = m[rank]; m[rank] = tmp;
        tmp = id[pr]; id[pr] = id[rank]; id[rank] = tmp;
        for (int r = 0; r < 64; r++) {
            if (r != rank && ((m[r] >> c) & 1)) {
                m[r] ^= m[rank];
                id[r] ^= id[rank];
            }
        }
        piv_col[rank] = c;
        rank++;
    }
    if (inv) memcpy(inv, id, 64 * sizeof(u64));
    return rank;
}

/* Apply 64x64 matrix M to each element of block vector v:
 * For each v[i], replace bit j of result with XOR of (bit k of v[i]) * M[k][j] for all k */
static void apply64(u64 M[64], u64 *v, u64 *out, int n) {
    for (int i = 0; i < n; i++) {
        u64 vi = v[i], res = 0;
        while (vi) {
            int b = __builtin_ctzll(vi);
            res ^= M[b];
            vi &= vi - 1;
        }
        out[i] = res;
    }
}

/*
 * Block Lanczos solver.
 *
 * Finds null vectors of B (nrows x ncols) over GF(2).
 * Returns number of dependencies found.
 *
 * Algorithm:
 * 1. Generate random x
 * 2. Compute A = B^T * B
 * 3. Iterate Lanczos with 64-wide block vectors
 * 4. Accumulate solution in w
 * 5. Extract null vectors from w
 */
static int lanczos_solve(sp_mat_t *B, int ***deps_out, int **dlen_out, int max_deps) {
    int n = B->ncols;

    /* Allocate block vectors */
    u64 *v0 = calloc(n, sizeof(u64));   /* v_{i-1} */
    u64 *v1 = calloc(n, sizeof(u64));   /* v_i */
    u64 *v2 = calloc(n, sizeof(u64));   /* v_{i+1} */
    u64 *Av  = calloc(n, sizeof(u64));  /* A * v_i */
    u64 *AAv = calloc(n, sizeof(u64));  /* A * A * v_i (for computing e) */
    u64 *tmp = malloc(B->nrows * sizeof(u64)); /* temp for B*x */
    u64 *w  = calloc(n, sizeof(u64));   /* solution accumulator */
    u64 *x  = calloc(n, sizeof(u64));   /* initial random vector */

    /* Random initialization */
    srand(42);
    for (int i = 0; i < n; i++)
        x[i] = ((u64)(unsigned)rand() << 32) | (unsigned)rand();

    /* v1 = A * x = B^T * B * x */
    sp_mul_btb(B, x, v1, tmp);

    /* Check if already zero */
    {
        int nz = 0;
        for (int i = 0; i < n; i++) if (v1[i]) { nz = 1; break; }
        if (!nz) { /* x is already in null space! */
            free(v0); free(v1); free(v2); free(Av); free(AAv); free(tmp); free(w); free(x);
            *deps_out = NULL; *dlen_out = NULL;
            return 0;
        }
    }

    u64 d0[64], d1[64]; /* D_{i-1}, D_i */
    memset(d0, 0, sizeof(d0));

    int max_iter = (n + 63) / 64 + 100;
    if (max_iter > 20000) max_iter = 20000;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Av = A * v1 = B^T * B * v1 */
        sp_mul_btb(B, v1, Av, tmp);

        /* D_i = <v1, Av> (64x64 matrix) */
        inner64(v1, Av, n, d1);

        /* Check convergence: D_i == 0 */
        int all_zero = 1;
        for (int i = 0; i < 64; i++) if (d1[i]) { all_zero = 0; break; }
        if (all_zero) {
            if (iter > 0) break;
            /* Try different random start */
            for (int i = 0; i < n; i++) x[i] ^= ((u64)(unsigned)rand() << 32) | (unsigned)rand();
            sp_mul_btb(B, x, v1, tmp);
            continue;
        }

        /* Compute D_i^{-1} */
        u64 d1_inv[64];
        int rank = inv64(d1, d1_inv);
        if (rank == 0) break;

        /* Update w: w += v1 * D_i^{-1} * <v1, x> */
        u64 vx[64];
        inner64(v1, x, n, vx);

        /* coeff = D_i^{-1} * vx */
        u64 coeff[64];
        for (int i = 0; i < 64; i++) {
            u64 acc = 0;
            for (int j = 0; j < 64; j++)
                if ((d1_inv[i] >> j) & 1) acc ^= vx[j];
            coeff[i] = acc;
        }

        /* w += v1 * coeff (apply 64x64 matrix to each block vector element) */
        for (int i = 0; i < n; i++) {
            u64 vi = v1[i], contrib = 0;
            while (vi) {
                int b = __builtin_ctzll(vi);
                contrib ^= coeff[b];
                vi &= vi - 1;
            }
            w[i] ^= contrib;
        }

        /* Compute e_i = D_i^{-1} * <Av, Av> for the recurrence */
        u64 AvAv[64];
        inner64(Av, Av, n, AvAv);
        u64 e[64];
        for (int i = 0; i < 64; i++) {
            u64 acc = 0;
            for (int j = 0; j < 64; j++)
                if ((d1_inv[i] >> j) & 1) acc ^= AvAv[j];
            e[i] = acc;
        }

        /* Compute f_i = D_i^{-1} * D_{i-1} (only if iter > 0) */
        u64 f[64];
        if (iter > 0) {
            for (int i = 0; i < 64; i++) {
                u64 acc = 0;
                for (int j = 0; j < 64; j++)
                    if ((d1_inv[i] >> j) & 1) acc ^= d0[j];
                f[i] = acc;
            }
        }

        /* v2 = Av - v1 * e - v0 * f */
        for (int i = 0; i < n; i++) {
            u64 r = Av[i];
            /* Subtract v1 * e */
            u64 vi = v1[i];
            while (vi) {
                int b = __builtin_ctzll(vi);
                r ^= e[b];
                vi &= vi - 1;
            }
            /* Subtract v0 * f (if iter > 0) */
            if (iter > 0) {
                u64 ui = v0[i];
                while (ui) {
                    int b = __builtin_ctzll(ui);
                    r ^= f[b];
                    ui &= ui - 1;
                }
            }
            v2[i] = r;
        }

        /* Shift: v0 = v1, v1 = v2 */
        u64 *t = v0; v0 = v1; v1 = v2; v2 = t;
        memcpy(d0, d1, 64 * sizeof(u64));

        if (iter > 0 && iter % 50 == 0)
            fprintf(stderr, "  BL iter %d/%d\n", iter, max_iter);
    }

    /* Extract null vectors from w */
    /* w is in the null space of B^T*B, which means B*w should be ~0 */
    int ndeps = 0;
    *deps_out = malloc(max_deps * sizeof(int*));
    *dlen_out = malloc(max_deps * sizeof(int));

    for (int bit = 0; bit < 64 && ndeps < max_deps; bit++) {
        int *dep = malloc(n * sizeof(int));
        int dlen = 0;
        for (int i = 0; i < n; i++)
            if ((w[i] >> bit) & 1) dep[dlen++] = i;
        if (dlen < 2) { free(dep); continue; }

        /* Verify: B * dep_vector should be 0 */
        /* Quick check: compute B * (indicator vector for dep) */
        u64 *test_x = calloc(n, sizeof(u64));
        for (int k = 0; k < dlen; k++) test_x[dep[k]] = 1;
        u64 *test_y = calloc(B->nrows, sizeof(u64));
        sp_mul(B, test_x, test_y);
        int is_null = 1;
        for (int r = 0; r < B->nrows; r++)
            if (test_y[r]) { is_null = 0; break; }
        free(test_x); free(test_y);

        if (is_null) {
            (*deps_out)[ndeps] = dep;
            (*dlen_out)[ndeps] = dlen;
            ndeps++;
        } else {
            free(dep);
        }
    }

    free(v0); free(v1); free(v2); free(Av); free(AAv);
    free(tmp); free(w); free(x);

    return ndeps;
}

#endif /* LANCZOS_H */
