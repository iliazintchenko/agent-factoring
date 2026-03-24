/*
 * Quadratic Sieve (QS) factoring implementation.
 *
 * Handles balanced semiprimes from 30-70 digits.
 *
 * Algorithm:
 *   1. Knuth-Schroeppel multiplier selection
 *   2. Factor base: primes p <= B where (kN|p) = 1
 *   3. Block sieving with log-byte accumulation
 *   4. Trial division on sieve survivors
 *   5. Single large prime (LP) variation
 *   6. GF(2) Gaussian elimination (bit-packed rows + history)
 *   7. X^2 ≡ Y^2 (mod N) and GCD extraction
 *
 * Compile: gcc -O3 -o qs qs.c -lgmp -lm
 * Usage:   ./qs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* =========================================================
 * Limits
 * ========================================================= */
#define MAX_FB        30000
#define BLOCK_SIZE    65536
/* Upper bounds on collected relations/partials */
#define MAX_RELS      60000
#define MAX_PARTS     500000

/* =========================================================
 * Small prime sieve (Eratosthenes)
 * ========================================================= */
static unsigned int *sp;   /* primes array */
static int           nsp;  /* count */

static void gen_primes(int lim) {
    char *sv = calloc(lim + 1, 1);
    int cap = lim < 2 ? 8 : lim / (int)(log((double)lim) - 1) + 64;
    sp  = malloc((size_t)cap * sizeof(unsigned int));
    nsp = 0;
    for (int i = 2; i <= lim; i++) {
        if (!sv[i]) {
            sp[nsp++] = i;
            for (long j = (long)i*i; j <= lim; j += i) sv[j] = 1;
        }
    }
    free(sv);
}

/* =========================================================
 * Modular arithmetic
 * ========================================================= */
static inline unsigned long mulmod(unsigned long a, unsigned long b, unsigned long m) {
    return (unsigned long)((__uint128_t)a * b % m);
}

static unsigned long powmod(unsigned long b, unsigned long e, unsigned long m) {
    unsigned long r = 1; b %= m;
    while (e) { if (e & 1) r = mulmod(r, b, m); b = mulmod(b, b, m); e >>= 1; }
    return r;
}

/* Legendre symbol (n|p), p odd prime */
static int legendre(unsigned long n, unsigned long p) {
    n %= p;
    if (!n) return 0;
    unsigned long v = powmod(n, (p-1)/2, p);
    return (v == 1) ? 1 : -1;
}

static int legendre_mpz(mpz_t n, unsigned long p) {
    return legendre(mpz_fdiv_ui(n, p), p);
}

/* Tonelli-Shanks: sqrt(n) mod p, (n|p) must be 1 */
static unsigned long tonelli(unsigned long n, unsigned long p) {
    n %= p;
    if (!n) return 0;
    if (p == 2) return n & 1;
    if (p % 4 == 3) return powmod(n, (p+1)/4, p);
    unsigned long S = 0, Q = p-1;
    while (!(Q & 1)) { S++; Q >>= 1; }
    unsigned long z = 2;
    while (powmod(z, (p-1)/2, p) != p-1) z++;
    unsigned long M = S;
    unsigned long c = powmod(z, Q, p);
    unsigned long t = powmod(n, Q, p);
    unsigned long R = powmod(n, (Q+1)/2, p);
    for (;;) {
        if (t == 1) return R;
        unsigned long i = 0, tmp = t;
        while (tmp != 1) { tmp = mulmod(tmp, tmp, p); i++; }
        unsigned long b = c;
        for (unsigned long j = 0; j < M-i-1; j++) b = mulmod(b, b, p);
        M = i; c = mulmod(b, b, p); t = mulmod(t, c, p); R = mulmod(R, b, p);
    }
}

/* =========================================================
 * Factor base
 * ========================================================= */
#define LOGSCALE 16  /* log2-units per byte: logval = round(log2(p) * LOGSCALE) */

typedef struct {
    unsigned int   p;
    unsigned int   r1;    /* sqrt(kN) mod p */
    unsigned int   r2;    /* p - r1 */
    unsigned short lv;    /* round(log2(p) * LOGSCALE) */
} FBE;

static FBE fb[MAX_FB];
static int  fbsz;

static void build_fb(mpz_t kN, int B, mpz_t sqkN) {
    mpz_sqrt(sqkN, kN);
    fbsz = 0;
    /* p = 2 */
    fb[0].p  = 2; fb[0].r1 = 0; fb[0].r2 = 0;
    fb[0].lv = (unsigned short)LOGSCALE; /* log2(2)*LOGSCALE = 16 */
    fbsz = 1;
    for (int i = 1; i < nsp && fbsz < MAX_FB; i++) {
        unsigned long p = sp[i];
        if ((int)p > B) break;
        if (legendre_mpz(kN, p) != 1) continue;
        unsigned long r = tonelli(mpz_fdiv_ui(kN, p), p);
        if (!r) continue;
        fb[fbsz].p  = (unsigned int)p;
        fb[fbsz].r1 = (unsigned int)r;
        fb[fbsz].r2 = (unsigned int)(p - r);
        fb[fbsz].lv = (unsigned short)(log2((double)p) * LOGSCALE + 0.5);
        fbsz++;
    }
}

/* =========================================================
 * Multiplier selection (Knuth-Schroeppel)
 * ========================================================= */
static double ks_score(mpz_t N, int k) {
    mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, (unsigned long)k);
    double sc = -0.5 * log2((double)k);
    /* p=2 */
    unsigned long kn8 = mpz_fdiv_ui(kN, 8);
    if      (kn8 == 1) sc += 2.0;
    else if (kn8 == 5) sc += 1.0;
    /* odd primes up to 3000 */
    for (int i = 1; i < nsp && sp[i] < 3000; i++) {
        unsigned long p = sp[i];
        int ls = legendre_mpz(kN, p);
        if      (ls ==  1) sc += 2.0 * log2((double)p) / (double)(p-1);
        else if (ls ==  0) sc +=       log2((double)p) / (double)(p-1);
    }
    mpz_clear(kN);
    return sc;
}

static int pick_mult(mpz_t N) {
    static const int cands[] = {
        1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,
        26,29,30,31,33,35,37,38,39,41,42,43,46,47,
        51,53,55,57,58,59,61,62,65,66,67,69,70,71,
        73,74,77,78,79,82,83,85,86,87,89,91,93,94,
        95,97,101,102,103,105,106,107,109,110,111,
        113,114,115,118,119,122,123,0
    };
    int bk = 1; double bs = -1e18;
    for (int i = 0; cands[i]; i++) {
        double s = ks_score(N, cands[i]);
        if (s > bs) { bs = s; bk = cands[i]; }
    }
    return bk;
}

/* =========================================================
 * Relation storage
 *
 * We store x and full exponent vector.  All mpz_t are
 * heap-initialised once; we never copy the Relation struct
 * (to avoid copying mpz_t internals).
 * ========================================================= */
typedef struct {
    mpz_t         x;      /* (sqkN + delta) [or product for LP pair] */
    int          *fexp;   /* full exponent vector, length fbsz */
    int           sign;   /* parity of negatives */
    unsigned long lpcof;  /* LP cofactor c (for LP-pair relations, else 0) */
    int           poly_col; /* GF2 column for MPQS polynomial A_j (-1 if not MPQS) */
    unsigned long poly_A;   /* A_j value for MPQS (0 if not MPQS) */
    /* mod-2 exponent vector: sign in bit 0, fb exps in bits 1..fbsz */
    /* stored as bit array for matrix */
} Rel;

typedef struct {
    mpz_t         x;
    int          *fexp;
    int           sign;
    unsigned long cof;    /* large prime cofactor */
    int           used;   /* 1 if already matched */
} Part;

static Rel  *rels;
static int   nrels, relcap;

static Part *parts;
static int   nparts, partcap;

static void storage_init(int rc, int pc) {
    relcap  = rc;
    rels    = calloc(rc, sizeof(Rel));
    for (int i = 0; i < rc; i++) { mpz_init(rels[i].x); rels[i].fexp = NULL; rels[i].lpcof = 0; rels[i].poly_col = -1; rels[i].poly_A = 0; }
    partcap = pc;
    parts   = calloc(pc, sizeof(Part));
    for (int i = 0; i < pc; i++) { mpz_init(parts[i].x); parts[i].fexp = NULL; }
    nrels = nparts = 0;
}

/* =========================================================
 * Sieve block
 *
 * Accumulate log2(p)*LOGSCALE in sieve_buf[i] for each
 * position i in [blk_start, blk_start+blk_len) where
 * f(i) = (sqkN + i)^2 - kN ≡ 0 (mod p).
 * ========================================================= */
static unsigned short sieve_buf[BLOCK_SIZE];

static void sieve_block(mpz_t kN, mpz_t sqkN, long bs, int bl) {
    memset(sieve_buf, 0, bl * sizeof(unsigned short));

    /* Pre-compute sqkN mod various values */
    unsigned long sq256 = mpz_fdiv_ui(sqkN, 256);
    unsigned long kn256 = mpz_fdiv_ui(kN,   256);

    for (int fi = 0; fi < fbsz; fi++) {
        unsigned long p  = fb[fi].p;
        unsigned char lv = fb[fi].lv;

        if (p == 2) {
            /* v_2(f(i)) = v_2((sqkN+i)^2 - kN)
             * For threshold estimation, compute low 8 bits of f(i). */
            for (int i = 0; i < bl; i++) {
                /* (sqkN + bs + i) mod 256 */
                long xi = (long)sq256 + bs + i;
                unsigned long xm = ((unsigned long)(xi % 256) + 256) & 255;
                unsigned long fm = (xm*xm - kn256 + 512) & 255;
                int v2 = 0;
                if (fm == 0) v2 = 8;
                else { unsigned long t = fm; while (!(t&1)) { t >>= 1; v2++; } }
                sieve_buf[i] += (unsigned short)(v2 * lv);
            }
            continue;
        }

        /* sqkN mod p */
        unsigned long mp = mpz_fdiv_ui(sqkN, p);

        /* Root 1: x ≡ r1 - mp (mod p), starting offset in block */
        {
            unsigned long r = fb[fi].r1;
            /* off = (r - mp - (bs % p) + 3p) % p */
            long bsp  = ((bs % (long)p) + (long)p) % (long)p;
            long rmmp = ((long)r - (long)mp - bsp + 3L*(long)p) % (long)p;
            for (long i = rmmp; i < bl; i += p) sieve_buf[i] += lv;
        }

        /* Root 2: x ≡ r2 - mp (mod p) */
        if (fb[fi].r2 != fb[fi].r1) {
            unsigned long r = fb[fi].r2;
            long bsp  = ((bs % (long)p) + (long)p) % (long)p;
            long rmmp = ((long)r - (long)mp - bsp + 3L*(long)p) % (long)p;
            for (long i = rmmp; i < bl; i += p) sieve_buf[i] += lv;
        }

        /* p^2 contribution: Hensel-lifted root (for small p, improves accuracy) */
        if (p <= 300) {
            unsigned long p2 = p * p;
            /* Lift each root: r' = r + t*p, where
             *   t = -(f(r)/p) * (2*(sqkN+r))^{-1}  mod p
             *   f(r) = (sqkN+r)^2 - kN
             * We compute (sqkN+r)^2 - kN mod p^2, divide by p to get c,
             * then t = -c * inv(2*(sqkN mod p + r), p) mod p.
             */
            unsigned long mp2  = mpz_fdiv_ui(sqkN, p2);
            unsigned long kn2  = mpz_fdiv_ui(kN,   p2);
            unsigned long bsp2 = ((bs % (long)p2) + p2) % p2;

            for (int ri = 0; ri < 2; ri++) {
                unsigned long r = (ri == 0) ? fb[fi].r1 : fb[fi].r2;
                if (ri == 1 && fb[fi].r2 == fb[fi].r1) break;
                /* (sqkN mod p2 + r)^2 mod p2 */
                unsigned long sr = (mp2 + r) % p2;
                unsigned long fr2 = (sr * sr % p2 + p2 - kn2 % p2) % p2;
                /* fr2 should be divisible by p (since r is a root mod p) */
                if (fr2 % p != 0) continue; /* shouldn't happen */
                long c = (long)(fr2 / p);
                /* derivative: 2*(sqkN mod p + r) mod p */
                unsigned long deriv = (2 * ((mp2 + r) % p)) % p;
                if (deriv == 0) continue;
                unsigned long inv_d = powmod(deriv, p-2, p);
                unsigned long t = (p - (unsigned long)(c % (long)p)) % p * inv_d % p;
                unsigned long r2 = (r + t * p) % p2;
                /* Start offset in block */
                long off2 = ((long)r2 - (long)bsp2 + 3L*(long)p2) % (long)p2;
                for (long i = off2; i < bl; i += p2) sieve_buf[i] += lv;
            }
        }
    }
}

/* =========================================================
 * MPQS: multiple polynomial quadratic sieve support
 *
 * Polynomial: f_j(x) = (A_j * x + B_j)^2 - kN
 *             = A_j * (A_j*x^2 + 2*B_j*x + C_j)  where C_j = (B_j^2-kN)/A_j
 *
 * Roots mod factor-base prime p:
 *   x ≡ (±r_p - B_j mod p) * inv(A_j mod p)  (mod p)
 *   where r_p = fb[fi].r1 = sqrt(kN mod p).
 *
 * The sieve accumulates log(p) for each p dividing h_j(x) = f_j(x)/A_j.
 * Trial division is on h_j(x) (after dividing f_j(x) by A_j).
 * ========================================================= */

/* Per-polynomial sieve roots (computed once per polynomial) */
typedef struct {
    int  r1;   /* first root mod p, offset in sieve */
    int  r2;   /* second root mod p (-1 if r1==r2) */
} MPRoots;

static MPRoots *mpqs_roots = NULL;  /* array of size MAX_FB */

/*
 * Compute MPQS sieve roots for polynomial A_j, B_j.
 * Bj_mod is B_j mod for first block offset (bs=0, so just roots mod p).
 * We write into mpqs_roots[].r1/r2 the starting offsets for the sieve range.
 */
static void compute_mpqs_roots(unsigned long A, mpz_t Bj_mpz) {
    if (!mpqs_roots) mpqs_roots = malloc(MAX_FB * sizeof(MPRoots));
    for (int fi = 0; fi < fbsz; fi++) {
        unsigned long p = fb[fi].p;
        if (p == 2) {
            mpqs_roots[fi].r1 = 0;
            mpqs_roots[fi].r2 = -1;
            continue;
        }
        unsigned long inv_A = powmod(A % p, p-2, p);
        unsigned long Bm = mpz_fdiv_ui(Bj_mpz, p);
        /* Root 1: (r1 - Bm) * inv_A mod p */
        unsigned long r1 = fb[fi].r1;
        unsigned long x1 = ((r1 + p - Bm) % p * inv_A) % p;
        /* Root 2: (r2 - Bm) * inv_A mod p */
        unsigned long r2 = fb[fi].r2;
        unsigned long x2 = ((r2 + p - Bm) % p * inv_A) % p;
        mpqs_roots[fi].r1 = (int)x1;
        mpqs_roots[fi].r2 = (x1 == x2) ? -1 : (int)x2;
    }
}

/*
 * Sieve one block for MPQS polynomial.
 * bs = block start offset (in sieve coordinates, e.g. -M to M).
 * bl = block length.
 */
static void sieve_block_mpqs(long bs, int bl) {
    memset(sieve_buf, 0, bl * sizeof(unsigned short));
    /* p=2: add log2(2)*LOGSCALE = 16 to every even position (where h_j is even).
     * h_j(x) has a fixed parity pattern; just add to every position as average
     * (adds 8 on average; actual positions needing p=2 get more from trial div). */
    unsigned short lv2 = fb[0].lv;  /* fb[0].lv = LOGSCALE = 16 */
    for (int i = 0; i < bl; i++) sieve_buf[i] += lv2;  /* add 16 to all */

    for (int fi = 1; fi < fbsz; fi++) {
        unsigned long p = fb[fi].p;
        unsigned short lv = fb[fi].lv;
        /* Root 1 */
        int x1 = mpqs_roots[fi].r1;
        long bsp = ((bs % (long)p) + (long)p) % (long)p;
        long off1 = ((long)x1 - bsp + (long)p) % (long)p;
        for (long i = off1; i < bl; i += p) sieve_buf[i] += lv;
        /* Root 2 */
        if (mpqs_roots[fi].r2 >= 0) {
            int x2 = mpqs_roots[fi].r2;
            long off2 = ((long)x2 - bsp + (long)p) % (long)p;
            for (long i = off2; i < bl; i += p) sieve_buf[i] += lv;
        }
    }
}

/* =========================================================
 * Trial division
 * Returns 1 if fully smooth (cofactor = 1), 0 otherwise.
 * fexp[] is filled with exponents; on return val = cofactor.
 * ========================================================= */
static int trial_div(mpz_t val, int *fexp) {
    memset(fexp, 0, fbsz * sizeof(int));
    for (int fi = 0; fi < fbsz; fi++) {
        unsigned long p = fb[fi].p;
        while (mpz_divisible_ui_p(val, p)) {
            mpz_divexact_ui(val, val, p);
            fexp[fi]++;
        }
    }
    return mpz_cmp_ui(val, 1) == 0;
}

/* =========================================================
 * GF(2) matrix with bit-packed rows and history
 * ========================================================= */
typedef unsigned long long u64;
#define WB 64

typedef struct {
    u64 *data;
    int  nr, nc, mw, hw, tw;
} GF2M;

static GF2M *gf2_new(int nr, int nc) {
    GF2M *m = malloc(sizeof(GF2M));
    m->nr = nr; m->nc = nc;
    m->mw = (nc + WB-1) / WB;
    m->hw = (nr + WB-1) / WB;
    m->tw = m->mw + m->hw;
    m->data = calloc((size_t)nr * m->tw, sizeof(u64));
    /* history = identity: row i has bit i set in history words */
    for (int i = 0; i < nr; i++)
        m->data[(size_t)i * m->tw + m->mw + i/WB] |= (1ULL << (i%WB));
    return m;
}

static inline void gf2_set(GF2M *m, int r, int c) {
    m->data[(size_t)r * m->tw + c/WB] |= (1ULL << (c%WB));
}

static inline int gf2_get(GF2M *m, int r, int c) {
    return (int)((m->data[(size_t)r * m->tw + c/WB] >> (c%WB)) & 1);
}

static void gf2_xor(GF2M *m, int d, int s) {
    u64 *dp = m->data + (size_t)d * m->tw;
    u64 *sp2 = m->data + (size_t)s * m->tw;
    for (int w = 0; w < m->tw; w++) dp[w] ^= sp2[w];
}

/* Gaussian elimination over GF(2).
 * Puts matrix in reduced row echelon form.
 * deps[] receives row indices of null-space rows.
 * Returns count of null-space rows. */
static int gf2_elim(GF2M *m, int *deps) {
    int nd = 0;
    int *piv = calloc(m->nr, sizeof(int)); /* piv[r] = 1 if row r is a pivot */
    int *pcol = calloc(m->nr, sizeof(int)); /* pcol[r] = pivot column of row r (if piv[r]) */

    for (int c = 0; c < m->nc; c++) {
        /* Find first non-pivot row with bit c set */
        int pr = -1;
        for (int r = 0; r < m->nr; r++) {
            if (!piv[r] && gf2_get(m, r, c)) { pr = r; break; }
        }
        if (pr < 0) continue;
        piv[pr] = 1; pcol[pr] = c;
        /* Eliminate column c from all other rows */
        for (int r = 0; r < m->nr; r++) {
            if (r != pr && gf2_get(m, r, c)) gf2_xor(m, r, pr);
        }
    }

    /* Collect zero rows (non-pivot rows with all-zero matrix part) */
    for (int r = 0; r < m->nr; r++) {
        if (piv[r]) continue;
        int z = 1;
        for (int w = 0; w < m->mw; w++) if (m->data[(size_t)r * m->tw + w]) { z = 0; break; }
        if (z) deps[nd++] = r;
    }
    free(piv); free(pcol);
    return nd;
}

static void gf2_free(GF2M *m) { free(m->data); free(m); }

/* =========================================================
 * LP matching: sort by cofactor, combine matching pairs
 *
 * IMPORTANT: Part structs contain mpz_t, which cannot be
 * safely copied by qsort.  We sort an index array instead.
 * ========================================================= */
static int cmp_part_idx(const void *a, const void *b) {
    int ia = *(const int *)a, ib = *(const int *)b;
    if (parts[ia].cof < parts[ib].cof) return -1;
    if (parts[ia].cof > parts[ib].cof) return  1;
    return 0;
}

static void match_lp(mpz_t N) {
    if (nparts < 2) return;
    int *idx = malloc(nparts * sizeof(int));
    for (int i = 0; i < nparts; i++) idx[i] = i;
    qsort(idx, nparts, sizeof(int), cmp_part_idx);

    int matched = 0;
    for (int i = 0; i+1 < nparts && nrels < relcap; i++) {
        int ia = idx[i], ib = idx[i+1];
        if (parts[ia].cof != parts[ib].cof || parts[ia].cof <= 1) continue;

        /* Skip if either was already matched in a previous pass.
         * Don't skip the outer loop (i++) so ib can still pair with idx[i+2]. */
        if (parts[ia].used || parts[ib].used) { continue; }

        int *ce = malloc(fbsz * sizeof(int));
        for (int j = 0; j < fbsz; j++)
            ce[j] = parts[ia].fexp[j] + parts[ib].fexp[j];
        int csign = parts[ia].sign ^ parts[ib].sign;

        /* Skip trivially-zero combinations (all exponents even AND sign even).
         * These produce X ≡ Y (mod N) and yield only trivial GCDs. */
        if (!csign) {
            int any_odd = 0;
            for (int j = 0; j < fbsz; j++) if (ce[j] & 1) { any_odd = 1; break; }
            if (!any_odd) { free(ce); i++; continue; }
        }

        int ri = nrels++;
        mpz_mul(rels[ri].x, parts[ia].x, parts[ib].x);
        mpz_mod(rels[ri].x, rels[ri].x, N);
        rels[ri].fexp  = ce;
        rels[ri].sign  = csign;
        rels[ri].lpcof = parts[ia].cof;
        parts[ia].used = 1;
        parts[ib].used = 1;
        matched++;
        i++; /* skip partner */
    }
    free(idx);
    fprintf(stderr, "LP matched: %d pairs from %d partials, total rels=%d\n",
            matched, nparts, nrels);
}

/* =========================================================
 * Factor extraction from a dependency row
 * ========================================================= */
static int try_dep(GF2M *mat, int dep, mpz_t N) {
    mpz_t X, Y, g;
    mpz_init_set_ui(X, 1);
    mpz_init_set_ui(Y, 1);
    mpz_init(g);

    int *se = calloc(fbsz, sizeof(int));
    int ss = 0;

    /* Track LP cofactors: each appears twice (once per matched partial).
     * Their product is a perfect square.  We need to include each LP in Y.
     * Since the same LP appears twice in the combined relation, the net
     * contribution to the exponent is 2*e_LP.  For the sqrt, we need LP^e_LP.
     * In a combined LP pair: lpcof appears ONCE in the combined relation
     * (it was squared in the combining, but we store lpcof once).
     * So for each LP relation in this dependency, multiply Y by lpcof.
     */
    mpz_t lp_product; mpz_init_set_ui(lp_product, 1);

    /* Track MPQS polynomial A_j occurrence counts for Y computation */
    /* We use a simple hash map: poly_col -> (count, A_j value) */
    /* Max poly_col index = fbsz + 1 + nrels (bounded; use a small array) */
    typedef struct { int col; int count; unsigned long A; } PolyCount;
    PolyCount *pcounts = NULL;
    int npoly_counts = 0;
    int poly_count_cap = 64;
    pcounts = malloc(poly_count_cap * sizeof(PolyCount));

    u64 *hist = mat->data + (size_t)dep * mat->tw + mat->mw;
    for (int i = 0; i < nrels; i++) {
        if (!((hist[i/WB] >> (i%WB)) & 1)) continue;
        mpz_mul(X, X, rels[i].x);
        mpz_mod(X, X, N);
        for (int j = 0; j < fbsz; j++) se[j] += rels[i].fexp[j];
        ss += rels[i].sign;
        /* Accumulate LP cofactors */
        if (rels[i].lpcof > 1) {
            mpz_mul_ui(lp_product, lp_product, rels[i].lpcof);
        }
        /* Track MPQS poly A_j occurrence */
        if (rels[i].poly_col >= 0 && rels[i].poly_A > 0) {
            int found_pc = 0;
            for (int k = 0; k < npoly_counts; k++) {
                if (pcounts[k].col == rels[i].poly_col) {
                    pcounts[k].count++;
                    found_pc = 1;
                    break;
                }
            }
            if (!found_pc) {
                if (npoly_counts >= poly_count_cap) {
                    poly_count_cap *= 2;
                    pcounts = realloc(pcounts, poly_count_cap * sizeof(PolyCount));
                }
                pcounts[npoly_counts].col   = rels[i].poly_col;
                pcounts[npoly_counts].count = 1;
                pcounts[npoly_counts].A     = rels[i].poly_A;
                npoly_counts++;
            }
        }
    }

    /* All exponents must be even and sign parity must be even
     * (guaranteed by the GF(2) elimination zero-row condition) */
    if (ss & 1) { free(se); mpz_clear(X); mpz_clear(Y); mpz_clear(g); mpz_clear(lp_product); return 0; }
    for (int j = 0; j < fbsz; j++) {
        if (se[j] & 1) { free(se); mpz_clear(X); mpz_clear(Y); mpz_clear(g); mpz_clear(lp_product); return 0; }
    }

    /* Y = (product p^(e/2)) * lp_product
     *
     * For each LP-combined relation (lpcof = c) in this dependency:
     *   X contribution: x1 * x2 where x1^2 = f1*c, x2^2 = f2*c
     *   So (x1*x2)^2 = f1*f2*c^2
     *   Y contribution from smooth part: sqrt(f1*f2) [tracked in se[]]
     *   Y contribution from LP: c (NOT sqrt(c^2) — we need to multiply by c)
     * Thus: Y = sqrt(smooth product) * product(lpcof_i for LP rels in dep)
     */
    for (int j = 0; j < fbsz; j++) {
        if (se[j] > 0) {
            mpz_t pw; mpz_init(pw);
            mpz_ui_pow_ui(pw, (unsigned long)fb[j].p, (unsigned long)(se[j]/2));
            mpz_mul(Y, Y, pw);
            mpz_mod(Y, Y, N);
            mpz_clear(pw);
        }
    }
    /* Multiply Y by all LP cofactors */
    if (mpz_cmp_ui(lp_product, 1) > 0) {
        mpz_mod(lp_product, lp_product, N);
        mpz_mul(Y, Y, lp_product);
        mpz_mod(Y, Y, N);
    }
    mpz_clear(lp_product);

    /* Multiply Y by A_j^(count/2) for each MPQS polynomial in the dependency.
     * Since all poly_col columns are zero in the dependency, each A_j appears
     * an even number of times.  Y needs A_j^(count/2). */
    for (int k = 0; k < npoly_counts; k++) {
        int cnt = pcounts[k].count;
        if (cnt & 1) {
            /* Odd count: dependency is invalid (shouldn't happen with correct GF2) */
            free(se); free(pcounts); mpz_clear(X); mpz_clear(Y); mpz_clear(g); return 0;
        }
        unsigned long Ak = pcounts[k].A;
        /* Multiply Y by Ak^(cnt/2) */
        mpz_t pw; mpz_init(pw);
        mpz_ui_pow_ui(pw, Ak, (unsigned long)(cnt/2));
        mpz_mod(pw, pw, N);
        mpz_mul(Y, Y, pw);
        mpz_mod(Y, Y, N);
        mpz_clear(pw);
    }
    free(pcounts);
    free(se);

    int found = 0;
    /* GCD(X-Y, N) */
    mpz_sub(g, X, Y); mpz_gcd(g, g, N);
    if (mpz_cmp_ui(g,1)>0 && mpz_cmp(g,N)<0) {
        gmp_printf("FACTOR: %Zd\n", g); fflush(stdout); found = 1;
    }
    if (!found) {
        /* GCD(X+Y, N) */
        mpz_add(g, X, Y); mpz_gcd(g, g, N);
        if (mpz_cmp_ui(g,1)>0 && mpz_cmp(g,N)<0) {
            gmp_printf("FACTOR: %Zd\n", g); fflush(stdout); found = 1;
        }
    }

    mpz_clear(X); mpz_clear(Y); mpz_clear(g);
    return found;
}

/* =========================================================
 * Main
 * ========================================================= */
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr,"Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N;
    mpz_init(N);
    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr,"Invalid N\n"); return 1;
    }

    int dig  = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);
    fprintf(stderr, "QS: %d digits (%d bits)\n", dig, bits);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* --- Trial division for small factors --- */
    gen_primes(1000000);
    for (int i = 0; i < nsp; i++) {
        if (mpz_divisible_ui_p(N, sp[i])) {
            printf("FACTOR: %u\n", sp[i]); fflush(stdout); return 0;
        }
    }

    /* Perfect square check */
    if (mpz_perfect_square_p(N)) {
        mpz_t r; mpz_init(r); mpz_sqrt(r, N);
        gmp_printf("FACTOR: %Zd\n", r); fflush(stdout); mpz_clear(r); return 0;
    }

    /* --- Parameter selection --- */
    double lnN  = (double)bits * log(2.0);
    double L    = exp(sqrt(lnN * log(lnN)));

    double cB;
    if      (dig <= 30) cB = 0.55;
    else if (dig <= 40) cB = 0.52;
    else if (dig <= 50) cB = 0.48;
    else if (dig <= 60) cB = 0.46;
    else                cB = 0.44;

    int B = (int)pow(L, cB);
    if (B < 500)     B = 500;
    if (B > 2000000) B = 2000000;

    /*
     * Initial sieve half-width M.
     *
     * Larger M means more sieving per pass before LP matching.  LP matching
     * efficiency scales as P^2, so it's beneficial to accumulate many partials
     * before matching.  We set M large enough that after one sieve pass we have
     * enough partials for good LP matching (target P ≈ sqrt(target) after pass 0).
     *
     * Rule of thumb: M = B * mult where mult is chosen so 2M * (smooth_rate_per_pos)
     * ≈ sqrt(target).  For 50-digit: smooth_rate ≈ 6/M, target=3563, so M ≈ sqrt(3563)/6
     * per M which gives 2M ≈ sqrt(3563)/6 * 1M → M ≈ 10M.  Use B*150.
     */
    long M;  /* half-sieve width (total: 2*M positions per pass) */
    if      (dig <= 30) M = (long)B * 30;
    else if (dig <= 40) M = (long)B * 60;
    else if (dig <= 50) M = (long)B * 150;
    else if (dig <= 60) M = (long)B * 100;
    else                M = (long)B * 150;
    if (M < 300000)    M = 300000;
    if (M > 500000000) M = 500000000;

    /*
     * LP bound: we use lp_bound = B * lp_mult, NOT B^2.
     *
     * With B^2, the matching probability is ~P^2/(2*π(B^2)) which is tiny
     * (need P ≈ sqrt(2*π(B^2)) ≈ sqrt(2*B^2/ln(B^2)) ≈ B/sqrt(2*ln(B)) partials
     * before even 1 match — for B=70K that's ~17K partials, but we only collect
     * a few thousand per pass).
     *
     * With lp_bound = lp_mult * B, π(lp_bound) ≈ lp_mult*B/ln(lp_mult*B), so
     * the expected pairs from P partials is P^2*ln(lp_mult*B)/(2*lp_mult*B).
     * For dig=50, B=70K, lp_mult=100, P=50K:
     *   pairs ≈ 50K^2 * ln(7M) / (2 * 7M) ≈ 2.5B * 15.8 / 14M ≈ 2821 pairs!
     *
     * Tradeoff: smaller lp_bound means fewer LP numbers exist, but we match them.
     * Empirically lp_mult in [50, 200] works best.
     */
    /*
     * LP multiplier: lp_bound = B * lp_mult.
     *
     * Larger lp_mult: more LP candidates captured but rarer collisions.
     * Smaller lp_mult: fewer LP candidates but more frequent collisions.
     *
     * Expected LP pairs from P partials uniformly over primes < lp_bound:
     *   E[pairs] = P^2 / (2 * π(lp_bound)) ≈ P^2 * ln(lp_bound) / (2 * lp_bound)
     *
     * We tune lp_mult so that at P ≈ sqrt(fbsz) partials we start getting matches:
     * target E[pairs] ≥ 1 when P = sqrt(target) / 2:
     *   P^2 * ln(B*lp_mult) / (2 * B * lp_mult) ≥ 1
     *   target * ln(B*lp_mult) / (8 * B * lp_mult) ≥ 1
     *   lp_mult ≈ target * ln(B*lp_mult) / (8 * B)
     * For B=70K, target=3500: lp_mult ≈ 3500*ln(7M)/8/70K ≈ 3500*15.7/560K ≈ 98
     * For B=120K, target=5700: lp_mult ≈ 5700*ln(12M)/8/120K ≈ 5700*16.3/960K ≈ 97
     *
     * So lp_mult ≈ 100 works across a range of sizes.
     */
    /*
     * LP multiplier: lp_bound = B * lp_mult.
     *
     * Optimal lp_mult maximizes LP pairs per unit sieve length.
     * LP pairs ∝ P^2/lp_bound where P ∝ log(lp_mult).
     * Optimal: lp_mult = e^2 ≈ 7.
     *
     * In practice, using lp_mult in [5, 20] works well.  For smaller N (30-40
     * digits), we use larger multipliers since the threshold is less critical.
     */
    int lp_mult;
    if      (dig <= 30) lp_mult = 30;
    else if (dig <= 40) lp_mult = 20;
    else if (dig <= 50) lp_mult = 10;
    else if (dig <= 60) lp_mult = 100;
    else                lp_mult = 100;

    unsigned long lp_bound = (unsigned long)B * (unsigned long)lp_mult;
    if (lp_bound > 2000000000UL) lp_bound = 2000000000UL;

    /* --- Multiplier and factor base --- */
    int k = pick_mult(N);
    fprintf(stderr, "k=%d B=%d M=%ld lp_bound=%lu\n", k, B, M, lp_bound);

    /* Ensure small_primes is large enough */
    if ((unsigned int)B > sp[nsp-1]) { free(sp); gen_primes(B + 100000); }

    mpz_t kN, sqkN;
    mpz_init(kN); mpz_init(sqkN);
    mpz_mul_ui(kN, N, (unsigned long)k);
    build_fb(kN, B, sqkN);
    fprintf(stderr, "Factor base: %d primes\n", fbsz);

    int target = fbsz + 80;

    /* Allocate relation storage */
    int rc = target * 5 + 200;
    if (rc > MAX_RELS) rc = MAX_RELS;
    storage_init(rc, MAX_PARTS);

    /*
     * Sieve threshold:
     *   f(delta) at delta=M has size ~ 2*M*sqrt(kN).
     *   In log2-bits: ~= 1 + log2(M) + 0.5*(log2(k) + bits)
     *   A B-smooth number of this size has log2 "used up" = log2(size).
     *   We accept positions that accumulated ~70% of log2(size) in sieve.
     *
     *   threshold (in log2*LOGSCALE bytes) = 0.70 * log2(2*M*sqrt(kN)) * LOGSCALE
     */
    /* Sieve threshold:
     * We want to accept positions where f(delta) has a B^2-smooth part >= T.
     * For a position with a single large prime < B^2, the smooth part is >= f / B^2.
     * The accumulated sieve value ≈ log2(smooth_part) * LOGSCALE.
     * threshold ≈ (log2(f(M)) - 2*log2(B) - slack) * LOGSCALE
     * where slack accounts for p^k>1 under-counting and small primes.
     *
     * Use: thr = max(0, log2(f(M)) - 2*log2(B) - 2) * LOGSCALE * 0.95
     */
    double log2_fM   = 1.0 + log2((double)M) + 0.5*(log2((double)k) + (double)bits);
    double log2_lpb  = log2((double)lp_bound);
    /*
     * Sieve threshold calculation:
     *
     * A "useful" LP position has f(delta) = smooth_part * lp_cofactor
     * with lp_cofactor < lp_bound.  The sieve accumulates roughly
     * log2(smooth_part)*LOGSCALE.  The minimum smooth accumulation for a
     * useful position is (log2_fM - log2_lpb)*LOGSCALE.
     *
     * We set threshold = (log2_fM - log2_lpb) * LOGSCALE * factor,
     * where factor < 1 accounts for sieve under-counting (prime powers,
     * approximate p=2 handling, etc.).  factor ≈ 0.80-0.85 works well.
     */
    double lp_smooth_bits = log2_fM - log2_lpb;
    if (lp_smooth_bits < 8.0) lp_smooth_bits = 8.0;
    /*
     * Threshold factor: accounts for sieve under-counting.
     * The sieve accumulates log2(p)*LOGSCALE per occurrence of prime p,
     * but misses prime power contributions (p^k for k>1 with p>300) and
     * uses an approximation for p=2.  For larger N, f(delta) has more
     * prime power contributions on average, so we use a smaller factor.
     *
     * Empirically tuned: 0.85 for 30-digit, 0.82 for 40-digit, 0.76 for
     * 50-digit, 0.68 for 60-digit, 0.62 for 70-digit.
     */
    double thr_factor;
    if      (dig <= 30) thr_factor = 0.88;
    else if (dig <= 40) thr_factor = 0.86;
    else if (dig <= 50) thr_factor = 0.82;
    else if (dig <= 60) thr_factor = 0.78;
    else                thr_factor = 0.73;
    double thr_f    = lp_smooth_bits * LOGSCALE * thr_factor;
    unsigned short thr = (unsigned short)(thr_f + 0.5);
    if (thr < 20) thr = 20;
    fprintf(stderr, "Threshold: lp_smooth_bits=%.1f thr=%u (log2_fM=%.1f log2_lpb=%.1f factor=%.2f)\n",
            lp_smooth_bits, thr, log2_fM, log2_lpb, thr_factor);

    /* Sieve loop */
    mpz_t xv, fv;
    mpz_init(xv); mpz_init(fv);

    /* MPQS state (initialized null for non-MPQS cases) */
    int npolys = 0;
    unsigned long *poly_A = NULL;
    mpz_t *poly_B = NULL;

    long sieved  = 0;
    long cands   = 0;
    int  found   = 0;
    int  pass    = 0;
    long curM    = M;
    long prevM   = 0;   /* how far we've already sieved */

    while (nrels < target && !found) {
        /* Positive and negative deltas: sieve [prevM, curM) in each direction */
        for (int dir = 0; dir <= 1 && nrels < target; dir++) {
            for (long bs = prevM; bs < curM && nrels < target; bs += BLOCK_SIZE) {
                long blk_start, blk_len;
                if (dir == 0) {
                    blk_start = bs;
                    blk_len   = BLOCK_SIZE;
                    if (blk_start + blk_len > curM) blk_len = curM - blk_start;
                } else {
                    /* negative direction: sieve [-curM, -prevM) */
                    blk_start = -(bs + BLOCK_SIZE);
                    blk_len   = BLOCK_SIZE;
                    if (bs + blk_len > curM) {
                        blk_len = (int)(curM - bs);
                        blk_start = -curM;
                    }
                    if (blk_len <= 0) break;
                }
                if (blk_len <= 0) break;

                sieve_block(kN, sqkN, blk_start, (int)blk_len);
                sieved += blk_len;

                /* Scan survivors */
                for (int i = 0; i < (int)blk_len; i++) {
                    if (sieve_buf[i] < thr) continue;
                    long delta = blk_start + i;
                    if (delta == 0) continue;
                    cands++;

                    /* x = sqkN + delta */
                    if (delta >= 0) mpz_add_ui(xv, sqkN, (unsigned long)delta);
                    else {
                        unsigned long ud = (unsigned long)(-delta);
                        if (mpz_cmp_ui(sqkN, ud) < 0) continue;
                        mpz_sub_ui(xv, sqkN, ud);
                    }

                    /* fv = x^2 - kN */
                    mpz_mul(fv, xv, xv);
                    mpz_sub(fv, fv, kN);
                    int sign = 0;
                    if (mpz_sgn(fv) < 0) { mpz_neg(fv, fv); sign = 1; }
                    if (mpz_sgn(fv) == 0) continue;

                    /* Trial division (modifies fv to cofactor) */
                    int *fe = malloc(fbsz * sizeof(int));
                    int smooth = trial_div(fv, fe);

                    if (smooth) {
                        if (nrels < relcap) {
                            int ri = nrels++;
                            mpz_set(rels[ri].x, xv);
                            rels[ri].fexp  = fe;
                            rels[ri].sign  = sign;
                            rels[ri].lpcof = 0;
                        } else free(fe);
                    } else {
                        /* Check LP */
                        unsigned long cof = 0;
                        if (mpz_fits_ulong_p(fv)) cof = mpz_get_ui(fv);
                        if (cof > 1 && cof < lp_bound) {
                            if (nparts < partcap) {
                                int pi = nparts++;
                                mpz_set(parts[pi].x, xv);
                                parts[pi].fexp = fe;
                                parts[pi].sign = sign;
                                parts[pi].cof  = cof;
                            } else free(fe);
                        } else free(fe);
                    }
                }

                /* Progress */
                if ((sieved % (BLOCK_SIZE * 32)) < BLOCK_SIZE) {
                    clock_gettime(CLOCK_MONOTONIC, &t1);
                    double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
                    fprintf(stderr, "\r  sieved=%.2fM smooth=%d part=%d cands=%ld t=%.1fs  ",
                            sieved/1e6, nrels, nparts, cands, el);
                }
            }
        }

        fprintf(stderr, "\nPass %d: smooth=%d part=%d sieved=%.2fM\n",
                pass, nrels, nparts, sieved/1e6);

        /* LP matching */
        match_lp(N);
        if (nrels >= target) break;

        /* Timeout check */
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
        /* For large N (60+digit), switch to MPQS quickly if QS is not converging */
        double phase1_limit = (dig >= 55) ? 5.0 : 270.0;
        if (el > phase1_limit && dig >= 55) { fprintf(stderr,"Switching to MPQS\n"); break; }
        if (el > 270.0) { fprintf(stderr,"Timeout\n"); break; }

        /* Update prevM: we've now sieved up to curM */
        prevM = curM;

        /* Extend range if not progressing */
        if (nrels < target) {
            long newM = (long)(curM * 1.5);
            if (newM > 1500000000L) { fprintf(stderr,"Range too large\n"); break; }
            fprintf(stderr, "Extending M to %ld\n", newM);
            curM = newM;
        }
        pass++;
        if (pass > 25) break;
    }

    fprintf(stderr,"\nAfter single-poly QS: %d (need %d), partials: %d\n", nrels, target, nparts);

    /* =========================================================
     * MPQS phase: for large N where single-poly QS is too slow.
     *
     * Generate multiple MPQS polynomials f_j(x) = (A_j*x+B_j)^2 - kN.
     * Each polynomial contributes relations with the "large prime" A_j.
     * Two relations from the same polynomial combine (LP-matching style)
     * to give a valid GF2 relation without needing A_j in the factor base.
     *
     * A_j is stored as poly_col in each Rel; the GF2 matrix includes
     * one extra column per distinct A_j used.
     * ========================================================= */
    if (nrels < target && dig >= 55) {
        /* MPQS parameters */
        long mpqs_M = M;   /* half-range per polynomial */
        if (mpqs_M > 5000000) mpqs_M = 5000000;  /* keep each poly sieve small */
        if (mpqs_M < 100000) mpqs_M = 100000;

        /* Target A size: sqrt(sqrt(2kN) / mpqs_M) */
        /* log2(target_A) = 0.25 * (bits + log2(k)) - 0.5 * log2(mpqs_M) */
        double log2_target_A = 0.25 * ((double)bits + log2((double)k)) - 0.5 * log2((double)mpqs_M);
        unsigned long target_A = (unsigned long)pow(2.0, log2_target_A);
        if (target_A < 1000000) target_A = 1000000;

        /* MPQS threshold:
         * h_j(x) = C_j + 2*B_j*x + A_j*x^2
         * At |x| = mpqs_M: h_j(M) ≈ 2*B_j*M + A_j*M^2.
         * With B_j ≈ sqrt(kN) and A_j = sqrt(sqrt(kN)/M), we have:
         *   A_j*M^2 = sqrt(kN/M) * M^2 = M^{3/2} * kN^{1/4}
         *   2*B_j*M ≈ 2*M*sqrt(kN)
         * For typical parameters, 2*B_j*M dominates, so log2_hM ≈ log2(2*M*sqrt(kN))
         * which equals the single-poly log2_fM (with mpqs_M replacing M).
         */
        double log2_hM = 1.0 + log2((double)mpqs_M) + 0.5*(log2((double)k) + (double)bits);
        if (log2_hM < 8.0) log2_hM = 8.0;
        double mpqs_lp_smooth = log2_hM - log2_lpb;
        if (mpqs_lp_smooth < 8.0) mpqs_lp_smooth = 8.0;
        double thr_mpqs_f = mpqs_lp_smooth * LOGSCALE * thr_factor;
        unsigned short thr_mpqs = (unsigned short)(thr_mpqs_f + 0.5);
        if (thr_mpqs < 20) thr_mpqs = 20;
        fprintf(stderr, "MPQS: M=%ld A≈%lu log2_hM=%.1f thr_mpqs=%u\n",
                mpqs_M, target_A, log2_hM, thr_mpqs);

        /* Allocate polynomial storage */
        int max_polys = 2000;
        poly_A = malloc(max_polys * sizeof(unsigned long));
        poly_B = malloc(max_polys * sizeof(mpz_t));
        for (int j = 0; j < max_polys; j++) mpz_init(poly_B[j]);

        /* Generate MPQS polynomials by searching for primes A near target_A
         * with (kN|A) == 1. */
        mpz_t Bj, Bj2, Cj, fj, hj, yj;
        mpz_init(Bj); mpz_init(Bj2); mpz_init(Cj);
        mpz_init(fj); mpz_init(hj); mpz_init(yj);

        unsigned long A_search = target_A;
        /* Find primes with (kN|A) == 1, and not in factor base.
         *
         * MPQS requires B_j to be LIFTED so that B_j ≈ sqrt(kN).
         * Starting from Bval_small = tonelli(kN mod A, A) ∈ [0, A),
         * we lift: B_j = Bval_small + k * A where k = (sqkN - Bval_small) / A
         * so that B_j ≡ sqrt(kN) (mod A) and B_j ≈ sqkN.
         *
         * With this lift:
         *   C_j = (B_j^2 - kN) / A_j ≈ (sqkN + delta)^2 - kN) / A_j ≈ small
         *   h_j(x) = A_j*x^2 + 2*B_j*x + C_j has values bounded by ~2*M*sqkN/A_j^???
         *
         * Actually: h_j(0) = C_j ≈ 0, h_j(±M) ≈ A_j*M^2 ± 2*sqkN*M
         * For A_j = sqrt(sqkN / M): h_j(M) ≈ sqkN/M * M^2 + 2*sqkN*M ≈ 3*M*sqkN
         * So |h_j| ≤ 3*M*sqkN ≈ same order as single-poly f(M) ≈ 2*M*sqkN.
         * The advantage is diversity across polynomials (independent values).
         */
        mpz_t Bj_lifted, kN_div_A, tmp_A;
        mpz_init(Bj_lifted); mpz_init(kN_div_A); mpz_init(tmp_A);

        while (A_search < 2*target_A + 100000) {
            /* Check if A_search is prime */
            int is_prime = 1;
            for (int ii = 0; ii < nsp && (long long)sp[ii]*sp[ii] <= A_search; ii++) {
                if (A_search % sp[ii] == 0) { is_prime = 0; break; }
            }
            if (!is_prime || (int)A_search <= B) { A_search++; continue; }
            if (legendre_mpz(kN, A_search) != 1) { A_search++; continue; }

            /* Valid A found: compute B_j = sqrt(kN mod A) */
            unsigned long Bval = tonelli(mpz_fdiv_ui(kN, A_search), A_search);
            if (Bval == 0) { A_search++; continue; }

            /* Lift B_j to be approximately sqkN.
             * B_j_lifted = Bval + k*A where k = floor((sqkN - Bval) / A) */
            mpz_set_ui(Bj_lifted, Bval);
            if (mpz_cmp(sqkN, Bj_lifted) > 0) {
                mpz_sub(tmp_A, sqkN, Bj_lifted);
                mpz_set_ui(kN_div_A, A_search);
                mpz_tdiv_q(tmp_A, tmp_A, kN_div_A);  /* k = floor((sqkN - Bval) / A) */
                mpz_mul_ui(tmp_A, tmp_A, A_search);   /* k * A */
                mpz_add(Bj_lifted, Bj_lifted, tmp_A);
            }
            /* Check: B_j^2 ≡ kN (mod A)? */
            /* (sqkN ≈ B_j_lifted, so yes by construction) */

            /* Store polynomial with lifted B */
            poly_A[npolys] = A_search;
            mpz_set(poly_B[npolys], Bj_lifted);
            npolys++;
            if (npolys >= max_polys) break;
            A_search++;
        }
        mpz_clear(Bj_lifted); mpz_clear(kN_div_A); mpz_clear(tmp_A);
        fprintf(stderr, "MPQS: generated %d polynomials, A in [%lu, %lu)\n",
                npolys, target_A, A_search);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double el_start = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;

        /* For each polynomial, sieve and collect relations */
        for (int j = 0; j < npolys && nrels < target && !found; j++) {
            unsigned long Aj = poly_A[j];

            /* Precompute MPQS roots for this polynomial (uses full mpz_t B_j) */
            compute_mpqs_roots(Aj, poly_B[j]);

            /* Sieve [-mpqs_M, mpqs_M] */
            for (int dir = 0; dir <= 1 && nrels < target; dir++) {
                for (long bs = 0; bs < mpqs_M && nrels < target; bs += BLOCK_SIZE) {
                    long blk_start, blk_len;
                    if (dir == 0) {
                        blk_start = bs;
                        blk_len = BLOCK_SIZE;
                        if (bs + blk_len > mpqs_M) blk_len = mpqs_M - bs;
                    } else {
                        blk_start = -(bs + BLOCK_SIZE);
                        blk_len = BLOCK_SIZE;
                        if (bs + blk_len > mpqs_M) {
                            blk_len = mpqs_M - bs;
                            blk_start = -mpqs_M;
                        }
                    }
                    if (blk_len <= 0) break;

                    sieve_block_mpqs(blk_start, (int)blk_len);
                    sieved += blk_len;

                    /* Scan survivors */
                    for (int i = 0; i < (int)blk_len; i++) {
                        if (sieve_buf[i] < thr_mpqs) continue;
                        long x = blk_start + i;
                        cands++;

                        /* y = A_j * x + B_j */
                        mpz_set_si(yj, x);
                        mpz_mul_ui(yj, yj, Aj);
                        mpz_add(yj, yj, poly_B[j]);  /* B_j is lifted, ≈ sqrt(kN) */

                        /* fj = yj^2 - kN. Use |yj| since (-y)^2 = y^2. */
                        if (mpz_sgn(yj) < 0) mpz_neg(yj, yj);
                        mpz_mul(fj, yj, yj);
                        mpz_sub(fj, fj, kN);
                        int sign = 0;
                        if (mpz_sgn(fj) < 0) { mpz_neg(fj, fj); sign = 1; }
                        if (mpz_sgn(fj) == 0) continue;

                        /* Divide out A_j: hj = fj / Aj */
                        if (!mpz_divisible_ui_p(fj, Aj)) continue;
                        mpz_divexact_ui(hj, fj, Aj);

                        /* Trial divide hj over factor base */
                        int *fe = malloc(fbsz * sizeof(int));
                        mpz_set(fv, hj);
                        int smooth = trial_div(fv, fe);

                        if (smooth) {
                            if (nrels < relcap) {
                                int ri = nrels++;
                                mpz_set(rels[ri].x, yj);
                                rels[ri].fexp     = fe;
                                rels[ri].sign     = sign;
                                rels[ri].lpcof    = 0;
                                rels[ri].poly_col = fbsz + 1 + j; /* column for A_j */
                                rels[ri].poly_A   = Aj;
                            } else free(fe);
                        } else {
                            /* LP check on hj residual */
                            unsigned long cof = 0;
                            if (mpz_fits_ulong_p(fv)) cof = mpz_get_ui(fv);
                            if (cof > 1 && cof < lp_bound) {
                                if (nparts < partcap) {
                                    int pi = nparts++;
                                    mpz_set(parts[pi].x, yj);
                                    parts[pi].fexp = fe;
                                    parts[pi].sign = sign;
                                    parts[pi].cof  = cof;
                                } else free(fe);
                            } else free(fe);
                        }
                    }
                }
            }

            /* Progress every 10 polynomials */
            if ((j+1) % 10 == 0) {
                clock_gettime(CLOCK_MONOTONIC, &t1);
                double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
                fprintf(stderr, "  MPQS poly %d/%d: rels=%d sieved=%.1fM t=%.1fs\n",
                        j+1, npolys, nrels, sieved/1e6, el);
            }
        }

        /* LP matching for any LP partials collected during MPQS */
        match_lp(N);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double el_end = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
        fprintf(stderr, "MPQS done: %d rels in %.1fs (%.1fs for MPQS phase)\n",
                nrels, el_end, el_end - el_start);

        mpz_clear(Bj); mpz_clear(Bj2); mpz_clear(Cj);
        mpz_clear(fj); mpz_clear(hj); mpz_clear(yj);
    }

    fprintf(stderr,"\nRelations: %d (need %d), partials: %d\n", nrels, target, nparts);

    if (nrels < fbsz + 2) {
        fprintf(stderr,"Insufficient relations\n"); return 1;
    }

    /* --- GF(2) matrix --- */
    int ur   = nrels;
    /* Columns: 0 = sign, 1..fbsz = prime exps, fbsz+1..fbsz+npolys = MPQS poly A_j */
    int nc   = fbsz + 1 + npolys;
    fprintf(stderr,"Matrix: %d x %d\n", ur, nc);

    GF2M *mat = gf2_new(ur, nc);
    for (int i = 0; i < ur; i++) {
        if (rels[i].sign & 1) gf2_set(mat, i, 0);
        for (int j = 0; j < fbsz; j++)
            if (rels[i].fexp[j] & 1) gf2_set(mat, i, j+1);
        /* Set MPQS polynomial column if applicable */
        if (rels[i].poly_col >= 0)
            gf2_set(mat, i, rels[i].poly_col);
    }

    /* --- Elimination --- */
    int *deps = malloc(ur * sizeof(int));
    int ndeps = gf2_elim(mat, deps);
    fprintf(stderr,"Dependencies: %d\n", ndeps);

    /* --- Factor extraction --- */
    for (int d = 0; d < ndeps && !found; d++) {
        found = try_dep(mat, deps[d], N);
        if (found) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
            fprintf(stderr,"Factor found at dep %d/%d, t=%.1fs\n", d, ndeps, el);
        }
    }
    if (!found) fprintf(stderr,"No factor found\n");

    free(deps);
    gf2_free(mat);
    /* Cleanup MPQS polynomial storage */
    if (poly_B) {
        for (int j = 0; j < 2000; j++) mpz_clear(poly_B[j]);
        free(poly_B);
    }
    if (poly_A) free(poly_A);
    if (mpqs_roots) { free(mpqs_roots); mpqs_roots = NULL; }
    mpz_clear(N); mpz_clear(kN); mpz_clear(sqkN);
    mpz_clear(xv); mpz_clear(fv);
    return found ? 0 : 1;
}
