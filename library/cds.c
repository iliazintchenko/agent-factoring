/*
 * Cofactor Descent Sieve (CDS) — novel factoring approach.
 *
 * Extends the MMS approach with ECM-assisted cofactor splitting:
 * when sieving finds Q(x) = smooth_part * cofactor (cofactor > B),
 * attempt to split the cofactor using a few ECM curves. If the
 * cofactor splits into medium-sized primes, the relation becomes
 * usable with multi-LP matching.
 *
 * Novel elements:
 * 1. Multi-multiplier sieving (from MMS)
 * 2. ECM cofactor descent: cofactor C → C1 * C2 via ECM
 * 3. Extended GF(2) matrix: each unique medium prime gets a column
 *    (no need for LP matching — primes become part of the factor base)
 *
 * Build: gcc -O3 -o cds cds.c -lgmp -L/usr/local/lib -lecm -lm
 * Usage: ./cds <N>
 */

#include <gmp.h>
#include <ecm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned int u32;
typedef unsigned long long u64;

/* ── Limits ── */
#define MAX_FB    16384
#define MAX_K     8
#define BLOCK     65536
#define MAX_REL   (MAX_FB + 4096)  /* extra space for extended LP columns */
#define MAX_PART  524288

/* ── Factor base ── */
typedef struct {
    u32 p;
    int logp;
    int nroots[MAX_K];
    u32 roots[MAX_K][2];
} fb_t;

/* ── Extended relation: tracks both FB primes and medium primes ── */
/* Medium primes get columns beyond g_fb_size in the exponent vector */
#define MAX_MED_PRIMES 8192
typedef struct {
    u64 *vec;           /* packed bit vector */
    mpz_t x_val;       /* LHS: product of (x + m_k) */
    mpz_t rhs_abs;     /* |product of f_k(x)| */
    int active;
} rel_t;

/* ── Globals ── */
static mpz_t g_N;
static int g_digits, g_K;
static u32 g_B;
static mpz_t g_m[MAX_K], g_kN[MAX_K];
static fb_t g_fb[MAX_FB];
static int g_fb_size;

static rel_t g_rels[MAX_REL];
static int g_nrels;

/* Extended column tracking: medium primes that appear as cofactors */
static u64 g_med_primes[MAX_MED_PRIMES];
static int g_nmed;
static int g_total_cols;   /* fb_size + 1 + nmed */
static int g_vec_words;

/* Hash table for medium prime → column index */
#define MED_HT_SIZE 32768
#define MED_HT_MASK (MED_HT_SIZE - 1)
typedef struct med_entry { u64 p; int col; struct med_entry *next; } med_entry_t;
static med_entry_t *g_med_ht[MED_HT_SIZE];
static med_entry_t g_med_pool[MAX_MED_PRIMES];

/* ── Tonelli-Shanks ── */
static u32 mod_sqrt(u32 a, u32 p) {
    if (a == 0) return 0;
    if (p == 2) return a & 1;
    u64 aa = a % p, pp = p;
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
    u64 M = S, c = 1;
    b = z;
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

/* ── Factor base ── */
static int build_factor_base(void) {
    char *sieve = calloc(g_B + 1, 1);
    for (u32 i = 2; i <= g_B; i++) sieve[i] = 1;
    for (u32 i = 2; (u64)i * i <= g_B; i++)
        if (sieve[i])
            for (u32 j = i * i; j <= g_B; j += i) sieve[j] = 0;

    int n = 0;
    g_fb[n].p = 2; g_fb[n].logp = 1;
    for (int k = 0; k < g_K; k++) {
        g_fb[n].nroots[k] = 1;
        g_fb[n].roots[k][0] = mpz_fdiv_ui(g_kN[k], 2);
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

/* ── Medium prime column management ── */
static int get_or_add_med_col(u64 p) {
    u32 h = (u32)(p * 2654435761ULL) & MED_HT_MASK;
    for (med_entry_t *e = g_med_ht[h]; e; e = e->next)
        if (e->p == p) return e->col;
    if (g_nmed >= MAX_MED_PRIMES) return -1;
    int col = g_fb_size + 1 + g_nmed;  /* +1 for sign */
    g_med_primes[g_nmed] = p;
    med_entry_t *e = &g_med_pool[g_nmed];
    e->p = p;
    e->col = col;
    e->next = g_med_ht[h];
    g_med_ht[h] = e;
    g_nmed++;
    /* Update total columns and vec words */
    g_total_cols = g_fb_size + 1 + g_nmed;
    g_vec_words = (g_total_cols + 63) / 64;
    return col;
}

/* ── Trial division + ECM cofactor descent ── */
static int trial_div_with_descent(mpz_t val, u64 **vec_out,
                                   mpz_t x_val, mpz_t rhs_abs) {
    /* Returns 1 if relation is usable, 0 otherwise.
       vec_out is allocated and filled on success. */
    mpz_t rem;
    mpz_init_set(rem, val);

    /* Temporary vec — we'll reallocate to proper size later */
    int max_words = (MAX_FB + MAX_MED_PRIMES + 64 + 63) / 64;
    u64 *vec = calloc(max_words, sizeof(u64));

    int neg = mpz_sgn(rem) < 0;
    if (neg) { mpz_neg(rem, rem); vec[0] |= 1ULL; }

    /* Trial division by factor base */
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
        /* Fully smooth */
        *vec_out = vec;
        mpz_clear(rem);
        return 1;
    }

    /* Cofactor descent: try to split rem using ECM */
    u64 LP_BOUND = (u64)g_B * (u64)g_B;
    int nbits = mpz_sizeinbase(rem, 2);

    if (nbits <= 64) {
        u64 r = mpz_get_ui(rem);
        if (r <= LP_BOUND) {
            /* Single medium prime — add as extended column */
            int col = get_or_add_med_col(r);
            if (col >= 0) {
                vec[col / 64] ^= 1ULL << (col % 64);
                *vec_out = vec;
                mpz_clear(rem);
                return 1;
            }
        }
    }

    /* Try ECM with small B1 to split the cofactor */
    if (nbits >= 20 && nbits <= 128) {
        mpz_t factor_found;
        mpz_init(factor_found);
        ecm_params params;
        ecm_init(params);
        params->B1done = 1.0;

        /* Try a few curves with small B1 */
        int ecm_B1_values[] = {100, 300, 1000, 3000};
        int n_ecm = (nbits < 40) ? 2 : 4;

        for (int ei = 0; ei < n_ecm; ei++) {
            mpz_set_ui(params->sigma, 42 + ei);  /* seed 42 */
            int ret = ecm_factor(factor_found, rem, (double)ecm_B1_values[ei], params);
            if (ret > 0 && mpz_cmp_ui(factor_found, 1) > 0 &&
                mpz_cmp(factor_found, rem) < 0) {
                /* Found a factor! Split rem = factor_found * cofactor */
                mpz_t cof;
                mpz_init(cof);
                mpz_fdiv_q(cof, rem, factor_found);

                /* Check both parts are small enough for medium primes */
                int ok = 1;
                mpz_t parts[2];
                mpz_init_set(parts[0], factor_found);
                mpz_init_set(parts[1], cof);

                for (int pi = 0; pi < 2 && ok; pi++) {
                    /* Each part might still be composite — trial divide more */
                    while (mpz_sizeinbase(parts[pi], 2) > 64) {
                        /* Try to split further with tiny ECM */
                        mpz_t ff2;
                        mpz_init(ff2);
                        mpz_set_ui(params->sigma, 42 + 10 + pi);
                        int r2 = ecm_factor(ff2, parts[pi], 100.0, params);
                        if (r2 > 0 && mpz_cmp_ui(ff2, 1) > 0 &&
                            mpz_cmp(ff2, parts[pi]) < 0) {
                            /* Record this sub-factor */
                            u64 sf = mpz_get_ui(ff2);
                            int col = get_or_add_med_col(sf);
                            if (col >= 0) vec[col / 64] ^= 1ULL << (col % 64);
                            else { ok = 0; }
                            mpz_fdiv_q(parts[pi], parts[pi], ff2);
                        } else {
                            ok = 0;  /* can't split further */
                        }
                        mpz_clear(ff2);
                    }
                    if (ok && mpz_cmp_ui(parts[pi], 1) > 0) {
                        u64 pv = mpz_get_ui(parts[pi]);
                        int col = get_or_add_med_col(pv);
                        if (col >= 0) vec[col / 64] ^= 1ULL << (col % 64);
                        else ok = 0;
                    }
                }

                mpz_clear(parts[0]); mpz_clear(parts[1]); mpz_clear(cof);

                if (ok) {
                    *vec_out = vec;
                    mpz_clear(factor_found);
                    ecm_clear(params);
                    mpz_clear(rem);
                    return 1;
                }
                break;  /* ECM found something but couldn't fully resolve */
            }
        }
        mpz_clear(factor_found);
        ecm_clear(params);
    }

    free(vec);
    mpz_clear(rem);
    return 0;
}

/* ── Sieve ── */
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
            if (!negative)
                start = ((r - m_mod - base_x) % (long)p + (long)p) % (long)p;
            else
                start = ((m_mod - 1 - r - base_x) % (long)p + (long)p) % (long)p;
            for (long j = start; j < len; j += p)
                g_sieve[j] += f->logp;
        }
    }
}

/* ── Factor extraction ── */
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

    if (!mpz_perfect_square_p(rhs_prod)) {
        mpz_clear(x); mpz_clear(y); mpz_clear(rhs_prod);
        mpz_clear(g); mpz_clear(tmp);
        return 0;
    }
    mpz_sqrt(y, rhs_prod);
    mpz_mod(y, y, g_N);

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

/* ── GF(2) elimination ── */
static int solve_and_factor(void) {
    int n = g_nrels;
    int ncols = g_total_cols;
    int nwords = (ncols + 63) / 64;

    fprintf(stderr, "LinAlg: %d rels, %d cols (%d FB + %d med)\n",
            n, ncols, g_fb_size + 1, g_nmed);
    if (n < ncols + 1) return 0;

    u64 **mat = malloc(n * sizeof(u64*));
    for (int i = 0; i < n; i++) {
        mat[i] = calloc(nwords, sizeof(u64));
        /* Copy only the words that exist in the relation's vec */
        int rel_words = g_vec_words;
        if (rel_words > nwords) rel_words = nwords;
        memcpy(mat[i], g_rels[i].vec, rel_words * sizeof(u64));
    }

    int hwords = (n + 63) / 64;
    u64 **hist = malloc(n * sizeof(u64*));
    for (int i = 0; i < n; i++) {
        hist[i] = calloc(hwords, sizeof(u64));
        hist[i][i / 64] |= 1ULL << (i % 64);
    }

    char *is_pivot = calloc(n, 1);
    for (int col = 0; col < ncols; col++) {
        int w = col / 64, b = col % 64;
        int piv = -1;
        for (int row = 0; row < n; row++) {
            if (!is_pivot[row] && (mat[row][w] & (1ULL << b))) {
                piv = row;
                break;
            }
        }
        if (piv < 0) continue;
        is_pivot[piv] = 1;
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

    int found = 0, ndeps = 0;
    for (int row = 0; row < n; row++) {
        int zero = 1;
        for (int w = 0; w < nwords; w++)
            if (mat[row][w]) { zero = 0; break; }
        if (!zero) continue;

        int *used = malloc(n * sizeof(int));
        int nu = 0;
        for (int r = 0; r < n; r++)
            if (hist[row][r / 64] & (1ULL << (r % 64)))
                used[nu++] = r;
        ndeps++;
        if (nu > 0 && !found)
            found = try_factor(used, nu);
        free(used);
    }

    fprintf(stderr, "Found %d dependencies, factor extracted: %d\n", ndeps, found);

    for (int i = 0; i < n; i++) { free(mat[i]); free(hist[i]); }
    free(mat); free(hist); free(is_pivot);
    return found;
}

/* ── Main factoring ── */
static int factor(void) {
    double logN = g_digits * log(10.0);
    double loglogN = log(logN);

    /* Smaller B than MMS — rely on ECM descent for cofactors */
    g_B = (u32)(exp(0.45 * sqrt(logN * loglogN)));
    if (g_B < 100) g_B = 100;
    if (g_B > 200000) g_B = 200000;

    g_K = 3 + g_digits / 20;
    if (g_K > MAX_K) g_K = MAX_K;
    if (g_K < 2) g_K = 2;

    static const int sqfree[] = {1,2,3,5,6,7,10,11,13,14,15,17};
    for (int k = 0; k < g_K; k++) {
        mpz_init(g_m[k]);
        mpz_init(g_kN[k]);
        mpz_mul_ui(g_kN[k], g_N, sqfree[k]);
        mpz_sqrt(g_m[k], g_kN[k]);
    }

    g_fb_size = build_factor_base();
    g_nmed = 0;
    g_total_cols = g_fb_size + 1;
    g_vec_words = (g_total_cols + MAX_MED_PRIMES + 63) / 64;
    memset(g_med_ht, 0, sizeof(g_med_ht));

    fprintf(stderr, "CDS: %d digits, B=%u, K=%d, FB=%d\n",
            g_digits, g_B, g_K, g_fb_size);

    int target = g_fb_size + 32;
    int half_bits = (int)(g_digits * 3.321928 / 2.0);
    int log2B = (int)(log2((double)g_B));

    g_nrels = 0;
    int full_rels = 0, ecm_rels = 0, failed = 0;

    mpz_t val, xv, absval;
    mpz_init(val); mpz_init(xv); mpz_init(absval);

    for (long bx = 0; bx < (1L << 32) && g_nrels < target; bx += BLOCK) {
        int blen = BLOCK;
        double log2_val = log2(2.0 * (bx + BLOCK / 2)) + half_bits;
        /* Lower threshold to get more candidates for ECM descent */
        int thresh = (int)(log2_val - 3.0 * log2B);
        if (thresh < 8) thresh = 8;
        if (thresh > 120) thresh = 120;

        for (int ki = 0; ki < g_K; ki++) {
            for (int neg = 0; neg <= 1; neg++) {
                sieve_block(ki, bx, blen, neg);

                for (int j = 0; j < blen; j++) {
                    if (g_sieve[j] < thresh) continue;
                    if (g_nrels >= MAX_REL) break;

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
                    mpz_mod(xv, xv, g_N);

                    u64 *vec = NULL;
                    if (trial_div_with_descent(val, &vec, xv, absval)) {
                        rel_t *r = &g_rels[g_nrels];
                        r->vec = vec;
                        mpz_init_set(r->x_val, xv);
                        mpz_init_set(r->rhs_abs, absval);
                        r->active = 1;
                        g_nrels++;
                        full_rels++;
                    }
                }
            }
        }

        if ((bx / BLOCK) % 20 == 0 && bx > 0) {
            g_total_cols = g_fb_size + 1 + g_nmed;
            fprintf(stderr, "  x=%ld rels=%d med=%d / target=%d\n",
                    bx, g_nrels, g_nmed, g_fb_size + 1 + g_nmed);
            target = g_total_cols + 32;
        }
    }

    g_total_cols = g_fb_size + 1 + g_nmed;
    g_vec_words = (g_total_cols + 63) / 64;
    fprintf(stderr, "Collected %d rels, %d med primes, %d total cols\n",
            g_nrels, g_nmed, g_total_cols);

    mpz_clear(val); mpz_clear(xv); mpz_clear(absval);

    if (g_nrels < g_total_cols + 1) {
        fprintf(stderr, "Not enough relations\n");
        return 0;
    }

    int ok = solve_and_factor();

    for (int k = 0; k < g_K; k++) {
        mpz_clear(g_m[k]);
        mpz_clear(g_kN[k]);
    }
    for (int i = 0; i < g_nrels; i++) {
        if (g_rels[i].active) {
            free(g_rels[i].vec);
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

    /* Quick trial division */
    u32 small_primes[] = {2,3,5,7,11,13};
    for (int i = 0; i < 6; i++) {
        if (mpz_fdiv_ui(g_N, small_primes[i]) == 0) {
            mpz_t q; mpz_init(q);
            mpz_fdiv_q_ui(q, g_N, small_primes[i]);
            gmp_printf("%u %Zd\n", small_primes[i], q);
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

    srand(42);

    if (!factor()) {
        fprintf(stderr, "FAILED\n");
        mpz_clear(g_N);
        return 1;
    }
    mpz_clear(g_N);
    return 0;
}
