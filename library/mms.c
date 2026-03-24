/*
 * Multi-Multiplier Sieve (MMS) — novel factoring approach.
 *
 * For each sieve point x, evaluates f_k(x) = (x + floor(sqrt(kN)))^2 - kN
 * for multiple multipliers k=1..K.  Since (x+m_k)^2 ≡ f_k(x) (mod N),
 * smooth values give congruences of squares.  Cross-multiplier large-prime
 * collisions amplify relation collection.
 *
 * Usage: ./mms <N>
 * Output: <p> <q>
 * Build: gcc -O3 -o mms mms.c -lgmp -lm
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned int u32;
typedef unsigned long long u64;

/* ── Limits ── */
#define MAX_FB    16384
#define MAX_K     12
#define MAX_REL   (MAX_FB + 1024)
#define MAX_PART  524288

/* ── Factor base ── */
typedef struct {
    u32 p;
    int logp;
    int nroots[MAX_K];
    u32 roots[MAX_K][2];
} fb_t;

/* ── Relation: stores exponent vector mod 2 and the LHS value ── */
typedef struct {
    u64  vec[(MAX_FB + 64) / 64];   /* bit vector: bit 0 = sign, bit i+1 = fb[i] */
    mpz_t x_val;                     /* the (x + m_k) value */
    mpz_t rhs_abs;                   /* |f_k(x)| for square-root reconstruction */
    int   active;
    int   usable;                    /* 1 = full or combined (safe for linalg) */
} rel_t;

/* ── Globals ── */
static mpz_t g_N;
static int   g_digits;
static int   g_K;
static u32   g_B;
static mpz_t g_m[MAX_K], g_kN[MAX_K];
static fb_t  g_fb[MAX_FB];
static int   g_fb_size;
static int   g_vec_words;
static rel_t g_rels[MAX_REL];
static int   g_nrels;

/* ── Partial relations (stored separately from usable ones) ── */
typedef struct {
    u64  vec[(MAX_FB + 64) / 64];
    mpz_t x_val;
    mpz_t rhs_abs;
} partial_t;

static partial_t g_parts[MAX_PART];
static int       g_nparts;

/* ── Large-prime matching ── */
#define LP_HASH_SIZE 131072
#define LP_HASH_MASK (LP_HASH_SIZE - 1)
typedef struct lp_node { u64 lp; int idx; struct lp_node *next; } lp_node_t;
static lp_node_t *g_lp_ht[LP_HASH_SIZE];
static lp_node_t  g_lp_pool[MAX_PART * 2];
static int        g_lp_next;

/* ── Tonelli-Shanks ── */
static u32 mod_sqrt(u32 a, u32 p) {
    if (a == 0) return 0;
    if (p == 2) return a & 1;
    u64 aa = a % p, pp = p;
    /* Euler criterion */
    u64 pw = 1, b = aa;
    for (u64 e = (p - 1) / 2; e; e >>= 1) {
        if (e & 1) pw = pw * b % pp;
        b = b * b % pp;
    }
    if (pw != 1) return 0;
    u32 Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }
    if (S == 1) {
        b = aa; pw = 1;
        for (u64 e = (u64)(p + 1) / 4; e; e >>= 1) {
            if (e & 1) pw = pw * b % pp;
            b = b * b % pp;
        }
        return (u32)pw;
    }
    u64 z = 2;
    for (;;) {
        pw = 1; b = z;
        for (u64 e = (p - 1) / 2; e; e >>= 1) {
            if (e & 1) pw = pw * b % pp;
            b = b * b % pp;
        }
        if (pw == pp - 1) break;
        z++;
    }
    u64 M = S;
    u64 c = 1; b = z;
    for (u64 e = (u64)Q; e; e >>= 1) { if (e & 1) c = c * b % pp; b = b * b % pp; }
    u64 t = 1; b = aa;
    for (u64 e = (u64)Q; e; e >>= 1) { if (e & 1) t = t * b % pp; b = b * b % pp; }
    u64 R = 1; b = aa;
    for (u64 e = (u64)(Q + 1) / 2; e; e >>= 1) { if (e & 1) R = R * b % pp; b = b * b % pp; }
    for (;;) {
        if (t == 1) return (u32)R;
        u64 tmp = t; u64 i = 0;
        while (tmp != 1) { tmp = tmp * tmp % pp; i++; }
        u64 bv = c;
        for (u64 j = 0; j < M - i - 1; j++) bv = bv * bv % pp;
        M = i; c = bv * bv % pp; t = t * c % pp; R = R * bv % pp;
    }
}

/* ── Factor base construction ── */
static int build_factor_base(void) {
    char *sieve = calloc(g_B + 1, 1);
    for (u32 i = 2; i <= g_B; i++) sieve[i] = 1;
    for (u32 i = 2; (u64)i * i <= g_B; i++)
        if (sieve[i])
            for (u32 j = i * i; j <= g_B; j += i) sieve[j] = 0;

    int n = 0;
    /* p = 2 */
    g_fb[n].p = 2; g_fb[n].logp = 1;
    for (int k = 0; k < g_K; k++) {
        u32 r = mpz_fdiv_ui(g_kN[k], 2);
        g_fb[n].nroots[k] = 1;
        g_fb[n].roots[k][0] = r;
    }
    n++;

    for (u32 p = 3; p <= g_B && n < MAX_FB; p += 2) {
        if (!sieve[p]) continue;
        int any = 0;
        g_fb[n].p = p;
        g_fb[n].logp = (int)(log2((double)p) + 0.5);
        for (int k = 0; k < g_K; k++) {
            u32 r = mpz_fdiv_ui(g_kN[k], p);
            u32 s = mod_sqrt(r, p);
            if (s || r == 0) {
                if (r == 0) {
                    g_fb[n].nroots[k] = 1;
                    g_fb[n].roots[k][0] = 0;
                } else {
                    g_fb[n].nroots[k] = 2;
                    g_fb[n].roots[k][0] = s;
                    g_fb[n].roots[k][1] = p - s;
                }
                any = 1;
            } else {
                g_fb[n].nroots[k] = 0;
            }
        }
        if (any) n++;
    }
    free(sieve);
    return n;
}

/* ── Trial-divide val by factor base, fill exponent vector (mod 2).
      Returns: 0 = not smooth enough
               1 = fully smooth (lp=0)
               1 = one large prime ≤ B^2 (lp = the prime)
               2 = two large primes (lp1, lp2, each ≤ B^2) ── */
static int trial_div(mpz_t val, u64 *vec, u64 *lp1, u64 *lp2, int k_idx) {
    mpz_t rem;
    mpz_init_set(rem, val);
    *lp1 = *lp2 = 0;
    memset(vec, 0, g_vec_words * sizeof(u64));

    int neg = mpz_sgn(rem) < 0;
    if (neg) { mpz_neg(rem, rem); vec[0] |= 1ULL; }

    for (int i = 0; i < g_fb_size; i++) {
        u32 p = g_fb[i].p;
        int exp = 0;
        while (mpz_fdiv_ui(rem, p) == 0) {
            mpz_fdiv_q_ui(rem, rem, p);
            exp++;
        }
        if (exp & 1) {
            int bit = i + 1;
            vec[bit / 64] ^= 1ULL << (bit % 64);
        }
    }

    if (mpz_cmp_ui(rem, 1) == 0) {
        mpz_clear(rem);
        return 1;  /* fully smooth */
    }

    u64 LP_BOUND = (u64)g_B * (u64)g_B;

    if (mpz_sizeinbase(rem, 2) <= 52) {
        u64 r = mpz_get_ui(rem);
        if (r <= LP_BOUND) {
            *lp1 = r;
            mpz_clear(rem);
            return 1;  /* single large prime */
        }
    }

    /* 2-LP disabled for now: combining 2-LP+1-LP leaves unsquared LP */

    mpz_clear(rem);
    return 0;
}

/* ── Store a relation ── */
static void store_rel(u64 *vec, mpz_t x_val, mpz_t rhs_abs, int usable) {
    if (g_nrels >= MAX_REL) return;
    rel_t *r = &g_rels[g_nrels];
    memcpy(r->vec, vec, g_vec_words * sizeof(u64));
    mpz_init_set(r->x_val, x_val);
    mpz_init_set(r->rhs_abs, rhs_abs);
    r->active = 1;
    r->usable = usable;
    g_nrels++;
}

/* ── Combine two partials sharing a large prime into a usable relation ── */
static void combine_parts(int p1, int p2) {
    if (g_nrels >= MAX_REL) return;
    partial_t *r1 = &g_parts[p1], *r2 = &g_parts[p2];

    u64 vec[MAX_FB / 64 + 2];
    memset(vec, 0, sizeof(vec));
    int bits1 = 0, bits2 = 0, bitsR = 0;
    for (int w = 0; w < g_vec_words; w++) {
        bits1 += __builtin_popcountll(r1->vec[w]);
        bits2 += __builtin_popcountll(r2->vec[w]);
        vec[w] = r1->vec[w] ^ r2->vec[w];
        bitsR += __builtin_popcountll(vec[w]);
    }
    mpz_t xv, rv;
    mpz_init(xv); mpz_init(rv);
    mpz_mul(xv, r1->x_val, r2->x_val);
    mpz_mod(xv, xv, g_N);
    mpz_mul(rv, r1->rhs_abs, r2->rhs_abs);

    store_rel(vec, xv, rv, 1);

    mpz_clear(xv); mpz_clear(rv);
}

/* ── Try to extract factor from a GF(2) null-space vector ── */
static int try_factor(int *rows, int nrows_used) {
    mpz_t x, y, rhs_prod, g, tmp;
    mpz_init_set_ui(x, 1);
    mpz_init(y);
    mpz_init_set_ui(rhs_prod, 1);
    mpz_init(g);
    mpz_init(tmp);

    for (int i = 0; i < nrows_used; i++) {
        rel_t *r = &g_rels[rows[i]];
        mpz_mul(x, x, r->x_val);
        mpz_mod(x, x, g_N);
        mpz_mul(rhs_prod, rhs_prod, r->rhs_abs);
    }

    /* y = sqrt(rhs_prod) — it should be a perfect square */
    if (!mpz_perfect_square_p(rhs_prod)) {
        mpz_clear(x); mpz_clear(y); mpz_clear(rhs_prod);
        mpz_clear(g); mpz_clear(tmp);
        return 0;
    }
    mpz_sqrt(y, rhs_prod);
    mpz_mod(y, y, g_N);

    /* Try gcd(x - y, N) and gcd(x + y, N) */
    mpz_sub(tmp, x, y);
    mpz_gcd(g, tmp, g_N);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, g_N) < 0) {
        mpz_fdiv_q(tmp, g_N, g);
        gmp_printf("%Zd %Zd\n", g, tmp);
        mpz_clear(x); mpz_clear(y); mpz_clear(rhs_prod);
        mpz_clear(g); mpz_clear(tmp);
        return 1;
    }
    mpz_add(tmp, x, y);
    mpz_gcd(g, tmp, g_N);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, g_N) < 0) {
        mpz_fdiv_q(tmp, g_N, g);
        gmp_printf("%Zd %Zd\n", g, tmp);
        mpz_clear(x); mpz_clear(y); mpz_clear(rhs_prod);
        mpz_clear(g); mpz_clear(tmp);
        return 1;
    }

    mpz_clear(x); mpz_clear(y); mpz_clear(rhs_prod);
    mpz_clear(g); mpz_clear(tmp);
    return 0;
}

/* ── GF(2) Gaussian elimination and factor extraction ── */
static int solve_and_factor(void) {
    /* Collect only usable relations (full + combined, not raw partials) */
    int *use_idx = malloc(g_nrels * sizeof(int));
    int n = 0;
    for (int i = 0; i < g_nrels; i++)
        if (g_rels[i].active && g_rels[i].usable)
            use_idx[n++] = i;

    /* Count trivially-zero vectors (sanity check) */
    int nzero = 0;
    for (int i = 0; i < n; i++) {
        int z = 1;
        for (int w = 0; w < g_vec_words; w++)
            if (g_rels[use_idx[i]].vec[w]) { z = 0; break; }
        if (z) nzero++;
    }
    fprintf(stderr, "LinAlg: %d usable rels, %d cols, %d trivially-zero vecs\n",
            n, g_fb_size + 1, nzero);
    if (n < g_fb_size + 1) { free(use_idx); return 0; }

    int ncols = g_fb_size + 1;
    int nwords = (ncols + 63) / 64;

    u64 **mat = malloc(n * sizeof(u64*));
    for (int i = 0; i < n; i++) {
        mat[i] = malloc(nwords * sizeof(u64));
        memcpy(mat[i], g_rels[use_idx[i]].vec, nwords * sizeof(u64));
    }

    int hwords = (n + 63) / 64;
    u64 **hist = malloc(n * sizeof(u64*));
    for (int i = 0; i < n; i++) {
        hist[i] = calloc(hwords, sizeof(u64));
        hist[i][i / 64] |= 1ULL << (i % 64);
    }

    /* Proper GF(2) Gaussian elimination: each row used as pivot at most once */
    char *is_pivot = calloc(n, 1);  /* marks rows already used as pivots */
    int *pivots = malloc(ncols * sizeof(int));
    memset(pivots, -1, ncols * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        int w = col / 64, b = col % 64;
        /* Find a non-pivot row with a 1 in this column */
        int piv = -1;
        for (int row = 0; row < n; row++) {
            if (!is_pivot[row] && (mat[row][w] & (1ULL << b))) {
                piv = row;
                break;
            }
        }
        if (piv < 0) continue;  /* empty column */

        pivots[col] = piv;
        is_pivot[piv] = 1;

        /* Eliminate this column from ALL other rows (including earlier pivot rows) */
        for (int row = 0; row < n; row++) {
            if (row == piv) continue;
            if (mat[row][w] & (1ULL << b)) {
                for (int ww = 0; ww < nwords; ww++)
                    mat[row][ww] ^= mat[piv][ww];
                for (int ww = 0; ww < hwords; ww++)
                    hist[row][ww] ^= hist[piv][ww];
            }
        }
    }

    int npivots = 0;
    for (int c = 0; c < ncols; c++) if (pivots[c] >= 0) npivots++;
    fprintf(stderr, "Rank: %d, expect %d dependencies\n", npivots, n - npivots);
    free(is_pivot);

    /* Find zero rows — these are dependencies */
    int found = 0, ndeps = 0;
    for (int row = 0; row < n; row++) {
        int zero = 1;
        for (int w = 0; w < nwords; w++)
            if (mat[row][w]) { zero = 0; break; }
        if (!zero) continue;

        /* Collect the ORIGINAL relation indices used in this dependency */
        int *used = malloc(n * sizeof(int));
        int nu = 0;
        for (int r = 0; r < n; r++) {
            if (hist[row][r / 64] & (1ULL << (r % 64)))
                used[nu++] = use_idx[r];  /* map back to original index */
        }

        ndeps++;
        if (nu > 0 && !found) {
            found = try_factor(used, nu);
        }
        free(used);
    }
    fprintf(stderr, "Found %d dependencies, factor extracted: %d\n", ndeps, found);

    for (int i = 0; i < n; i++) { free(mat[i]); free(hist[i]); }
    free(mat); free(hist); free(pivots); free(use_idx);
    return found;
}

/* ── Sieve one block ── */
#define BLOCK 65536
static signed char g_sieve[BLOCK];

static void sieve_block(int ki, long base_x, int len, int negative) {
    memset(g_sieve, 0, len);
    for (int i = 0; i < g_fb_size; i++) {
        fb_t *f = &g_fb[i];
        u32 p = f->p;
        if (f->nroots[ki] == 0) continue;
        long m_mod = (long)mpz_fdiv_ui(g_m[ki], p);

        for (int ri = 0; ri < f->nroots[ki]; ri++) {
            long r = f->roots[ki][ri];
            long start;
            if (!negative) {
                /* Need (base_x + j + m_k) ≡ r (mod p) */
                start = ((r - m_mod - base_x) % (long)p + (long)p) % (long)p;
            } else {
                /* f_k(-(base_x+j+1)) = (m_k - base_x - j - 1)^2 - kN
                   Need (m_k - base_x - j - 1) ≡ r (mod p)
                   j ≡ m_k - 1 - r - base_x (mod p) */
                start = ((m_mod - 1 - r - base_x) % (long)p + (long)p) % (long)p;
            }
            for (long j = start; j < len; j += p)
                g_sieve[j] += f->logp;

            /* Also sieve for p^2 (adds another logp) */
            if ((u64)p * p <= g_B) {
                u64 pp = (u64)p * p;
                /* Roots of kN mod p^2 via Hensel lift */
                /* Skip for simplicity in v1 */
            }
        }
    }
}

/* ── Main factoring ── */
static int factor(void) {
    /* Parameters */
    double logN = g_digits * log(10.0);
    double loglogN = log(logN);
    double Lhalf = exp(0.5 * sqrt(logN * loglogN));

    g_B = (u32)(Lhalf);
    if (g_B < 100) g_B = 100;
    if (g_B > 400000) g_B = 400000;

    g_K = 3 + g_digits / 20;
    if (g_K > MAX_K) g_K = MAX_K;
    if (g_K < 2) g_K = 2;

    fprintf(stderr, "MMS: %d digits, B=%u, K=%d\n", g_digits, g_B, g_K);

    /* Init multipliers — use square-free values only to avoid redundant parity vectors */
    static const int sqfree[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23};
    for (int k = 0; k < g_K; k++) {
        mpz_init(g_m[k]);
        mpz_init(g_kN[k]);
        mpz_mul_ui(g_kN[k], g_N, sqfree[k]);
        mpz_sqrt(g_m[k], g_kN[k]);
    }

    /* Build factor base */
    g_fb_size = build_factor_base();
    g_vec_words = (g_fb_size + 1 + 63) / 64;
    fprintf(stderr, "FB: %d primes up to %u\n", g_fb_size, g_fb[g_fb_size - 1].p);

    int target = g_fb_size + 32;
    if (target > MAX_REL) target = MAX_REL;

    /* Sieve threshold: we want log2(|f_k(x)|) minus sieve_sum to be < log2(B^2)
       f_k(x) ≈ 2x * sqrt(kN), log2 ≈ log2(2x) + digits*3.32/2
       At x = M, we want the sieve sum to account for most of the log.
       threshold ≈ log2(2*M*sqrt(N)) - 2*log2(B) */
    int half_bits = (int)(g_digits * 3.321928 / 2.0);
    int log2B = (int)(log2((double)g_B));
    /* We'll use an adaptive threshold per block based on x position */

    g_nrels = 0;
    g_nparts = 0;
    g_lp_next = 0;
    memset(g_lp_ht, 0, sizeof(g_lp_ht));

    mpz_t val, xv, absval;
    mpz_init(val); mpz_init(xv); mpz_init(absval);

    int full_rels = 0, partials = 0, combined = 0;

    for (long bx = 0; bx < (1L << 32) && (full_rels + combined) < target; bx += BLOCK) {
        int blen = BLOCK;

        /* Compute threshold for this block position */
        double log2_val = log2(2.0 * (bx + BLOCK / 2)) + half_bits;
        int thresh = (int)(log2_val - 2.0 * log2B);
        if (thresh < 10) thresh = 10;
        if (thresh > 120) thresh = 120;

        for (int ki = 0; ki < g_K; ki++) {
            for (int neg = 0; neg <= 1; neg++) {
                sieve_block(ki, bx, blen, neg);

                for (int j = 0; j < blen; j++) {
                    if (g_sieve[j] < thresh) continue;

                    long x = bx + j;
                    if (!neg) {
                        mpz_set_ui(xv, (unsigned long)x);
                        mpz_add(xv, xv, g_m[ki]);
                        mpz_mul(val, xv, xv);
                        mpz_sub(val, val, g_kN[ki]);
                    } else {
                        mpz_set(xv, g_m[ki]);
                        mpz_sub_ui(xv, xv, (unsigned long)(x + 1));
                        mpz_mul(val, xv, xv);
                        mpz_sub(val, val, g_kN[ki]);
                    }

                    mpz_abs(absval, val);
                    u64 vec_buf[(MAX_FB + 64) / 64];
                    u64 lp1 = 0, lp2 = 0;

                    int td = trial_div(val, vec_buf, &lp1, &lp2, ki);
                    if (td > 0) {
                        mpz_mod(xv, xv, g_N);

                        if (lp1 == 0) {
                            /* Fully smooth */
                            store_rel(vec_buf, xv, absval, 1);
                            full_rels++;
                        } else if (g_nparts < MAX_PART) {
                            /* Store partial (1-LP or 2-LP) */
                            partial_t *pp = &g_parts[g_nparts];
                            memcpy(pp->vec, vec_buf, g_vec_words * sizeof(u64));
                            mpz_init_set(pp->x_val, xv);
                            mpz_init_set(pp->rhs_abs, absval);
                            int my_idx = g_nparts;
                            g_nparts++;
                            partials++;

                            /* Match single large prime */
                            u32 h = (u32)(lp1 * 2654435761ULL) & LP_HASH_MASK;
                            int match = -1;
                            for (lp_node_t *nd = g_lp_ht[h]; nd; nd = nd->next) {
                                if (nd->lp == lp1) { match = nd->idx; break; }
                            }
                            if (match >= 0) {
                                combine_parts(match, my_idx);
                                combined++;
                            }
                            if (g_lp_next < MAX_PART * 2) {
                                lp_node_t *nd = &g_lp_pool[g_lp_next++];
                                nd->lp = lp1;
                                nd->idx = my_idx;
                                nd->next = g_lp_ht[h];
                                g_lp_ht[h] = nd;
                            }
                        }
                    }
                }
            }
        }

        if ((bx / BLOCK) % 20 == 0 && bx > 0) {
            fprintf(stderr, "  x=%ld rels=%d (full=%d part=%d comb=%d) / %d\n",
                    bx, g_nrels, full_rels, partials, combined, target);
        }
    }

    fprintf(stderr, "Collected %d rels (full=%d part=%d comb=%d)\n",
            g_nrels, full_rels, partials, combined);

    mpz_clear(val); mpz_clear(xv); mpz_clear(absval);

    if (g_nrels < g_fb_size + 1) {
        fprintf(stderr, "Not enough relations (%d < %d)\n", g_nrels, g_fb_size + 1);
        return 0;
    }

    /* Solve */
    int ok = solve_and_factor();

    for (int k = 0; k < g_K; k++) {
        mpz_clear(g_m[k]);
        mpz_clear(g_kN[k]);
    }
    for (int i = 0; i < g_nrels; i++) {
        if (g_rels[i].active) {
            mpz_clear(g_rels[i].x_val);
            mpz_clear(g_rels[i].rhs_abs);
        }
    }
    return ok;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }
    mpz_init_set_str(g_N, argv[1], 10);
    g_digits = (int)mpz_sizeinbase(g_N, 10);

    /* Quick trial division up to 10^6 */
    {
        u32 primes_small[] = {2, 3, 5, 7, 11, 13};
        for (int i = 0; i < 6; i++) {
            if (mpz_fdiv_ui(g_N, primes_small[i]) == 0) {
                mpz_t q; mpz_init(q);
                mpz_fdiv_q_ui(q, g_N, primes_small[i]);
                gmp_printf("%u %Zd\n", primes_small[i], q);
                mpz_clear(q); mpz_clear(g_N);
                return 0;
            }
        }
        for (u32 p = 17; p < 1000000; p += 2) {
            if (mpz_fdiv_ui(g_N, p) == 0) {
                mpz_t q; mpz_init(q);
                mpz_fdiv_q_ui(q, g_N, p);
                gmp_printf("%u %Zd\n", p, q);
                mpz_clear(q); mpz_clear(g_N);
                return 0;
            }
        }
    }

    if (!factor()) {
        fprintf(stderr, "FAILED\n");
        mpz_clear(g_N);
        return 1;
    }
    mpz_clear(g_N);
    return 0;
}
