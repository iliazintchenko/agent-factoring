/*
 * mpqs.c - Multiple Polynomial Quadratic Sieve (SIQS variant)
 * Optimized for balanced semiprimes, 30-80 digits
 *
 * Compile: gcc -O3 -march=native library/mpqs.c -o mpqs -lgmp -lm
 * Usage: ./mpqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>

#define SIEVE_SIZE 65536
#define MAX_FB 50000
#define MAX_RELS 120000
#define LP_HASH_SIZE (1 << 18)
#define LP_HASH_MASK (LP_HASH_SIZE - 1)

typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned char u8;

/* ==================== Parameters ==================== */
typedef struct {
    int fb_target;
    int sieve_half;   /* sieve from -M to M */
    int nfactors;     /* # primes in a */
    int lp_mult;
    int thresh_adj;   /* subtract from theoretical threshold */
} params_t;

static params_t get_params(int bits) {
    /* s chosen so that product of s mid-range FB primes ≈ sqrt(2N)/M */
    if (bits <= 100) return (params_t){100,  16384, 4, 30, 5};
    if (bits <= 120) return (params_t){200,  32768, 5, 30, 7};
    if (bits <= 140) return (params_t){400,  32768, 6, 40, 9};
    if (bits <= 160) return (params_t){900,  65536, 7, 50, 11};
    if (bits <= 180) return (params_t){2000, 65536, 8, 50, 13};
    if (bits <= 200) return (params_t){4000, 98304, 8, 60, 15};
    if (bits <= 220) return (params_t){7000, 131072, 9, 70, 17};
    if (bits <= 240) return (params_t){12000, 196608, 10, 80, 19};
    if (bits <= 260) return (params_t){22000, 262144, 11, 90, 21};
    if (bits <= 280) return (params_t){40000, 327680, 11, 100, 23};
    return (params_t){70000, 524288, 12, 110, 25};
}

/* ==================== Factor Base ==================== */
static u32 fb_p[MAX_FB];     /* primes */
static u32 fb_r[MAX_FB];     /* sqrt(N) mod p */
static u8  fb_log[MAX_FB];
static int fb_size;

static u32 tonelli(mpz_t N, u32 p) {
    u32 n = mpz_fdiv_ui(N, p);
    if (n == 0) return 0;
    if (p == 2) return n & 1;

    /* Euler criterion */
    mpz_t b, e, m, r;
    mpz_inits(b, e, m, r, NULL);
    mpz_set_ui(b, n); mpz_set_ui(m, p);
    mpz_set_ui(e, (p-1)/2);
    mpz_powm(r, b, e, m);
    if (mpz_cmp_ui(r, 1) != 0) { mpz_clears(b,e,m,r,NULL); return 0; }

    if (p % 4 == 3) {
        mpz_set_ui(e, (p+1)/4);
        mpz_powm(r, b, e, m);
        u32 ret = mpz_get_ui(r);
        mpz_clears(b,e,m,r,NULL);
        return ret;
    }

    /* Full Tonelli-Shanks */
    u32 Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    u32 z = 2;
    while (1) {
        mpz_set_ui(b, z); mpz_set_ui(e, (p-1)/2);
        mpz_powm(r, b, e, m);
        if (mpz_cmp_ui(r, p-1) == 0) break;
        z++;
    }
    mpz_t M, c, t, R, bb, tmp;
    mpz_inits(M, c, t, R, bb, tmp, NULL);
    mpz_set_ui(M, S);
    mpz_set_ui(c, z); mpz_set_ui(e, Q); mpz_powm(c, c, e, m);
    mpz_set_ui(t, n); mpz_powm(t, t, e, m);
    mpz_set_ui(R, n); mpz_set_ui(e, (Q+1)/2); mpz_powm(R, R, e, m);
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            u32 ret = mpz_get_ui(R);
            mpz_clears(b,e,m,r,M,c,t,R,bb,tmp,NULL);
            return ret;
        }
        int i = 0; mpz_set(tmp, t);
        while (mpz_cmp_ui(tmp, 1) != 0) { mpz_mul(tmp,tmp,tmp); mpz_mod(tmp,tmp,m); i++; }
        mpz_set(bb, c);
        for (int j = 0; j < (int)mpz_get_ui(M) - i - 1; j++) { mpz_mul(bb,bb,bb); mpz_mod(bb,bb,m); }
        mpz_set_ui(M, i);
        mpz_mul(c, bb, bb); mpz_mod(c, c, m);
        mpz_mul(t, t, c); mpz_mod(t, t, m);
        mpz_mul(R, R, bb); mpz_mod(R, R, m);
    }
}

static void build_fb(mpz_t N, int target) {
    fb_p[0] = 2; fb_r[0] = 1; fb_log[0] = 1; fb_size = 1;

    int bound = target * 20 + 50000;
    char *sieve = calloc(bound, 1);
    for (int i = 2; (long)i*i < bound; i++)
        if (!sieve[i]) for (int j = i*i; j < bound; j += i) sieve[j] = 1;

    for (int p = 3; p < bound && fb_size < target; p += 2) {
        if (sieve[p]) continue;
        u32 r = tonelli(N, p);
        if (r == 0) continue;
        fb_p[fb_size] = p;
        fb_r[fb_size] = r;
        fb_log[fb_size] = (u8)(log2((double)p) + 0.5);
        fb_size++;
    }
    free(sieve);
}

/* ==================== Relations ==================== */
typedef struct {
    mpz_t val;     /* ax+b */
    u8 *exps;      /* exponents mod 2 (bit vector stored as bytes) */
    int *full_exps; /* full exponents (for sqrt step) */
    u64 lp;        /* large prime (0 if full) */
} rel_t;

static rel_t rels[MAX_RELS];
static int nrels;
static int nfull;

/* LP hash */
typedef struct lp_node { u64 lp; int idx; struct lp_node *next; } lp_node_t;
static lp_node_t *lp_hash[LP_HASH_SIZE];
static lp_node_t lp_pool[MAX_RELS];
static int lp_used;

static int lp_match(u64 lp) {
    u32 h = (u32)(lp * 2654435761ULL) & LP_HASH_MASK;
    for (lp_node_t *e = lp_hash[h]; e; e = e->next)
        if (e->lp == lp) return e->idx;
    return -1;
}

static void lp_insert(u64 lp, int idx) {
    u32 h = (u32)(lp * 2654435761ULL) & LP_HASH_MASK;
    lp_node_t *e = &lp_pool[lp_used++];
    e->lp = lp; e->idx = idx; e->next = lp_hash[h]; lp_hash[h] = e;
}

/* ==================== Modular inverse ==================== */
static u32 modinv(u32 a, u32 m) {
    long long old_r = a, r = m, old_s = 1, s = 0;
    while (r) {
        long long q = old_r / r;
        long long tmp = r; r = old_r - q*r; old_r = tmp;
        tmp = s; s = old_s - q*s; old_s = tmp;
    }
    return (u32)((old_s % (long long)m + m) % m);
}

/* ==================== SIQS Core ==================== */

static mpz_t gN;
static u8 *sieve_arr;
static int *soln1, *soln2;

/* Polynomial: Q(x) = (a*x + b)^2 - N = a^2*x^2 + 2abx + (b^2-N)
 * We sieve g(x) = a*x^2 + 2*b*x + c where c = (b^2-N)/a
 * A smooth g(x) means (ax+b)^2 ≡ g(x)*a (mod N)
 */

static mpz_t poly_a, poly_b, poly_c;
static mpz_t Bvals[20];
static int a_idx[20];
static int n_afactors;
static int poly_count;

static void choose_a(const params_t *par) {
    mpz_t target, two_n;
    mpz_inits(target, two_n, NULL);
    mpz_mul_ui(two_n, gN, 2);
    mpz_sqrt(target, two_n);
    mpz_tdiv_q_ui(target, target, par->sieve_half);

    n_afactors = par->nfactors;

    /* target_a^(1/s) = ideal per-factor size */
    double log_target = mpz_sizeinbase(target, 2);
    double log_per = log_target / n_afactors;
    u32 ideal_prime = (u32)(1UL << (int)log_per);

    /* Find FB primes near ideal size */
    int center = 1;
    for (int i = 1; i < fb_size; i++) {
        if (fb_p[i] >= ideal_prime) { center = i; break; }
    }
    if (center < n_afactors + 2) center = n_afactors + 2;
    if (center + n_afactors >= fb_size) center = fb_size - n_afactors - 1;

    /* Deterministic pseudo-random selection */
    srand(42 + poly_count);

    int range = fb_size / 4;
    if (range < 2 * n_afactors) range = 2 * n_afactors;
    int lo = center - range/2;
    if (lo < 1) lo = 1;
    int hi = lo + range;
    if (hi >= fb_size) { hi = fb_size - 1; lo = hi - range; if (lo < 1) lo = 1; }

    int *used = calloc(fb_size, sizeof(int));
    mpz_set_ui(poly_a, 1);
    int sel = 0;

    for (int k = 0; k < n_afactors && sel < 20; k++) {
        int idx;
        int tries = 0;
        do {
            idx = lo + (rand() % (hi - lo));
            tries++;
        } while (used[idx] && tries < 1000);
        if (tries >= 1000) break;
        used[idx] = 1;
        a_idx[sel++] = idx;
        mpz_mul_ui(poly_a, poly_a, fb_p[idx]);
    }
    n_afactors = sel;
    free(used);

    /* Compute B values */
    mpz_t adq;
    mpz_init(adq);

    for (int j = 0; j < n_afactors; j++) {
        mpz_init(Bvals[j]);
        u32 qj = fb_p[a_idx[j]];
        mpz_tdiv_q_ui(adq, poly_a, qj);
        u32 adq_mod = mpz_fdiv_ui(adq, qj);
        u32 inv = modinv(adq_mod, qj);
        u32 gamma = (u32)(((u64)inv * fb_r[a_idx[j]]) % qj);
        if (gamma > qj/2) gamma = qj - gamma;
        mpz_mul_ui(Bvals[j], adq, gamma);
    }

    /* b = sum of Bvals */
    mpz_set_ui(poly_b, 0);
    for (int j = 0; j < n_afactors; j++)
        mpz_add(poly_b, poly_b, Bvals[j]);

    /* c = (b^2 - N) / a */
    mpz_mul(poly_c, poly_b, poly_b);
    mpz_sub(poly_c, poly_c, gN);
    if (!mpz_divisible_p(poly_c, poly_a)) {
        /* b^2 - N not divisible by a: b is wrong */
        fprintf(stderr, "WARNING: b^2 - N not divisible by a\n");
    }
    mpz_tdiv_q(poly_c, poly_c, poly_a);

    mpz_clear(adq);
    mpz_clears(target, two_n, NULL);
}

/* Compute sieve start positions */
static void compute_roots(int M) {
    for (int i = 1; i < fb_size; i++) {
        u32 p = fb_p[i];

        /* Check if p divides a */
        int divides_a = 0;
        for (int j = 0; j < n_afactors; j++) {
            if (a_idx[j] == i) { divides_a = 1; break; }
        }

        if (divides_a) {
            /* Single root: x = -c * (2b)^(-1) mod p */
            u32 b_mod = mpz_fdiv_ui(poly_b, p);
            u32 c_mod = mpz_fdiv_ui(poly_c, p);
            if (b_mod == 0) { soln1[i] = soln2[i] = -1; continue; }
            u32 inv2b = modinv((2 * b_mod) % p, p);
            u32 x0 = (u32)(((u64)(p - c_mod) * inv2b) % p);
            /* Map to [0, 2M) array index */
            soln1[i] = (int)((x0 + (u64)M) % p);
            soln2[i] = -1;
            continue;
        }

        u32 a_mod = mpz_fdiv_ui(poly_a, p);
        u32 b_mod = mpz_fdiv_ui(poly_b, p);
        u32 inv_a = modinv(a_mod, p);

        /* roots of a*x^2 + 2*b*x + c = 0 mod p are:
         * x = (-b +/- sqrt(N)) / a mod p  (since a*c + b^2 = N, sqrt(disc) = sqrt(N)) */
        u32 r1 = fb_r[i];
        u32 r2 = p - r1;

        /* x1 = (r1 - b) * a^(-1) mod p */
        u32 x1 = (u32)(((u64)((r1 + p - b_mod) % p) * inv_a) % p);
        u32 x2 = (u32)(((u64)((r2 + p - b_mod) % p) * inv_a) % p);

        /* Map to array coords [0, 2M): sieve[i] corresponds to x = i - M */
        soln1[i] = (int)(((u64)x1 + M) % p);
        soln2[i] = (int)(((u64)x2 + M) % p);
    }
}

/* Sieve a block of the array */
static void sieve_block(int start, int len) {
    memset(sieve_arr, 0, len);

    for (int i = 1; i < fb_size; i++) {
        u32 p = fb_p[i];
        u8 lp = fb_log[i];

        if (soln1[i] >= 0) {
            int pos = soln1[i] - start;
            if (pos < 0) pos += p * ((-pos + p - 1) / p);
            for (; pos < len; pos += p)
                sieve_arr[pos] += lp;
        }

        if (soln2[i] >= 0) {
            int pos = soln2[i] - start;
            if (pos < 0) pos += p * ((-pos + p - 1) / p);
            for (; pos < len; pos += p)
                sieve_arr[pos] += lp;
        }
    }
}

/* Try to create a relation from sieve hit at position x (relative to sieve center) */
static int process_candidate(int x, const params_t *par) {
    /* g(x) = a*x^2 + 2*b*x + c */
    mpz_t gx, ax_b, rem;
    mpz_inits(gx, ax_b, rem, NULL);

    mpz_set_si(gx, x);
    mpz_mul(gx, gx, poly_a);
    mpz_mul_ui(gx, gx, x);  /* a*x^2 ... wait, need 2b*x */

    /* Actually: g(x) = a*x^2 + 2*b*x + c */
    mpz_set_si(ax_b, x);
    mpz_mul(ax_b, ax_b, poly_a);  /* a*x */
    mpz_add(ax_b, ax_b, poly_b);  /* a*x + b = the value whose square we compute */

    /* g(x) = (ax+b)^2 - N, divided by a = a*x^2 + 2bx + c */
    mpz_mul(gx, ax_b, ax_b);
    mpz_sub(gx, gx, gN);
    /* We actually want g(x)/a for factoring, but the relation involves ax+b */
    /* Actually Q(x) = (ax+b)^2 - N. We need Q(x) to be smooth (possibly times a) */
    /* Q(x) = a * g(x) where g(x) = a*x^2 + 2bx + c. So Q(x) is smooth if g(x) is smooth (up to factors of a) */

    /* Let's factor g(x) = a*x^2 + 2bx + c */
    mpz_set_si(gx, x);
    mpz_mul_si(gx, gx, x);
    mpz_mul(gx, gx, poly_a);    /* a*x^2 */
    mpz_t bx;
    mpz_init(bx);
    mpz_mul_si(bx, poly_b, 2*x); /* 2*b*x */
    mpz_add(gx, gx, bx);
    mpz_add(gx, gx, poly_c);    /* + c */
    mpz_clear(bx);

    int sign = 0;
    if (mpz_sgn(gx) < 0) { sign = 1; mpz_neg(gx, gx); }
    if (mpz_sgn(gx) == 0) {
        mpz_clears(gx, ax_b, rem, NULL);
        return 0;
    }

    if (nrels >= MAX_RELS) { mpz_clears(gx, ax_b, rem, NULL); return 0; }

    rel_t *rel = &rels[nrels];
    rel->exps = calloc(fb_size, 1);
    rel->full_exps = calloc(fb_size, sizeof(int));
    rel->exps[0] = sign; /* -1 factor */
    rel->full_exps[0] = sign;
    mpz_init_set(rel->val, ax_b);
    rel->lp = 0;

    /* Include the 'a' coefficient factors in the exponent vector.
     * Since (ax+b)^2 = N + a*g(x), the relation is:
     * (ax+b)^2 ≡ a * g(x) (mod N)
     * So we need the exponent vector for a * g(x), not just g(x). */
    for (int j = 0; j < n_afactors; j++) {
        int idx = a_idx[j];
        rel->exps[idx] ^= 1;   /* a contributes 1 to each factor */
        rel->full_exps[idx]++;
    }

    /* Trial divide g(x) by FB */
    for (int i = 1; i < fb_size; i++) {
        u32 p = fb_p[i];
        while (mpz_fdiv_ui(gx, p) == 0) {
            mpz_tdiv_q_ui(gx, gx, p);
            rel->exps[i] ^= 1; /* toggle parity */
            rel->full_exps[i]++;
        }
    }

    /* Dedup: check if this val already exists */
    for (int k = 0; k < nrels; k++) {
        if (rels[k].exps && mpz_cmp(rels[k].val, ax_b) == 0) {
            /* Duplicate - skip */
            mpz_clear(rel->val);
            free(rel->exps);
            free(rel->full_exps);
            mpz_clears(gx, ax_b, rem, NULL);
            return 0;
        }
    }

    /* Check residue */
    int ok = 0;
    if (mpz_cmp_ui(gx, 1) == 0) {
        /* Full relation - verify correctness */
        if (nrels < 3) {
            /* Verify: val^2 mod N == product of p^full_exps mod N */
            mpz_t v2, prod_check, ptmp;
            mpz_inits(v2, prod_check, ptmp, NULL);
            mpz_mul(v2, rel->val, rel->val);
            mpz_mod(v2, v2, gN);
            mpz_set_ui(prod_check, 1);
            for (int j = 1; j < fb_size; j++) {
                if (rel->full_exps[j] > 0) {
                    mpz_ui_pow_ui(ptmp, fb_p[j], rel->full_exps[j]);
                    mpz_mul(prod_check, prod_check, ptmp);
                }
            }
            mpz_mod(prod_check, prod_check, gN);
            if (mpz_cmp(v2, prod_check) != 0) {
                gmp_fprintf(stderr, "REL VERIFY FAIL: val=%Zd, val^2 mod N=%Zd, prod mod N=%Zd\n",
                           rel->val, v2, prod_check);
            } else {
                fprintf(stderr, "REL OK: x=%d, %d factors in a*g(x)\n", x, rel->full_exps[0] ? -1 : 1);
            }
            mpz_clears(v2, prod_check, ptmp, NULL);
        }
        nrels++;
        nfull++;
        ok = 1;
    } else if (mpz_sizeinbase(gx, 2) <= 52) {
        u64 lp = mpz_get_ui(gx);
        u64 bound = (u64)fb_p[fb_size-1] * par->lp_mult;
        if (lp <= bound) {
            rel->lp = lp;
            nrels++;
            if (lp_match(lp) >= 0) nfull++;
            else lp_insert(lp, nrels - 1);
            ok = 1;
        }
    }

    if (!ok) {
        mpz_clear(rel->val);
        free(rel->exps);
        rel->exps = NULL;
    }

    mpz_clears(gx, ax_b, rem, NULL);
    return ok;
}

/* ==================== Linear Algebra (Gaussian Elimination over GF(2)) ==================== */

static int g_dep_lens[64];
static int *g_deps;

static int gauss_elim(int max_deps) {
    int nr = nrels, nc = fb_size;
    int fw = (nc + 63) / 64, iw = (nr + 63) / 64;

    u64 **mat = malloc(nr * sizeof(u64*));
    for (int r = 0; r < nr; r++) {
        mat[r] = calloc(fw + iw, sizeof(u64));
        if (!rels[r].exps) continue;
        for (int c = 0; c < nc; c++)
            if (rels[r].exps[c]) mat[r][c/64] |= (1ULL << (c%64));
        mat[r][fw + r/64] |= (1ULL << (r%64));
    }

    int *piv = calloc(nc, sizeof(int));
    for (int c = 0; c < nc; c++) piv[c] = -1;

    for (int c = 0; c < nc; c++) {
        int pr = -1;
        for (int r = 0; r < nr; r++) {
            if (!((mat[r][c/64] >> (c%64)) & 1)) continue;
            int used = 0;
            for (int k = 0; k < c; k++) if (piv[k] == r) { used = 1; break; }
            if (!used) { pr = r; break; }
        }
        if (pr < 0) continue;
        piv[c] = pr;
        for (int r = 0; r < nr; r++) {
            if (r == pr) continue;
            if ((mat[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < fw + iw; w++) mat[r][w] ^= mat[pr][w];
        }
    }

    int ndeps = 0;
    g_deps = malloc(max_deps * nr * sizeof(int));

    for (int r = 0; r < nr && ndeps < max_deps; r++) {
        int zero = 1;
        for (int w = 0; w < fw && zero; w++) {
            u64 mask = (w == fw-1 && nc%64) ? ((1ULL << (nc%64))-1) : ~0ULL;
            if (mat[r][w] & mask) zero = 0;
        }
        if (!zero) continue;

        int *d = &g_deps[ndeps * nr];
        int cnt = 0;
        for (int w = 0; w < iw; w++) {
            u64 bits = mat[r][fw + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w*64 + bit;
                if (idx < nr) d[cnt++] = idx;
                bits &= bits - 1;
            }
        }
        g_dep_lens[ndeps] = cnt;
        ndeps++;
    }

    for (int r = 0; r < nr; r++) free(mat[r]);
    free(mat); free(piv);
    return ndeps;
}

/* ==================== Square Root ==================== */

static int try_sqrt(int *dep, int dlen) {
    mpz_t X, Y, tmp, g;
    mpz_inits(X, Y, tmp, g, NULL);
    mpz_set_ui(X, 1);
    mpz_set_ui(Y, 1);

    /*
     * We have: for each relation i, (val_i)^2 ≡ a_i * g(x_i) (mod N)
     * where g(x_i) = product of fb_primes^full_exps * lp
     * The dependency means product of g(x_i) is a perfect square.
     * So: (product val_i)^2 ≡ (product a_i) * (product g(x_i)) (mod N)
     *
     * For SIQS: val_i = a_i*x + b_i, and (a_i*x + b_i)^2 - N = a_i * g(x_i)
     * So (val_i)^2 ≡ a_i * g(x_i) (mod N)
     *
     * Product: (prod val_i)^2 ≡ (prod a_i) * (prod g(x_i)) (mod N)
     * If prod g(x_i) = S^2, then (prod val_i)^2 ≡ (prod a_i) * S^2 (mod N)
     *
     * Actually, since we're doing SIQS with varying a's, this is tricky.
     * Simpler approach: compute val_i^2 - N for each relation, factor that,
     * and build the congruence X^2 ≡ Y^2 (mod N).
     *
     * val^2 - N = a * g(x). So the relation is val^2 ≡ N (mod nothing useful)...
     *
     * Standard QS approach: (ax+b)^2 ≡ g(x)*a (mod N) where g(x) is smooth.
     * But we need (ax+b)^2 ≡ a*g(x) (mod N).
     * The matrix should track exponents of ALL prime factors including a.
     *
     * Actually the simpler fix: just compute val^2 mod N for each relation in the dependency,
     * take the product, and also compute the product of (val^2 - N) = a*g(x), factor it, take sqrt.
     */

    /* Use stored full exponents for a*g(x) */
    int *total = calloc(fb_size, sizeof(int));

    for (int i = 0; i < dlen; i++) {
        int ri = dep[i];
        if (!rels[ri].full_exps) continue;

        /* X = product of val_i mod N */
        mpz_mul(X, X, rels[ri].val);
        mpz_mod(X, X, gN);

        /* Sum exponents from stored full_exps (which includes a*g(x) factors) */
        for (int j = 0; j < fb_size; j++)
            total[j] += rels[ri].full_exps[j];
    }

    /* Debug: verify product matches for small deps */
    static int sqrt_dbg = 0;
    if (sqrt_dbg < 1 && dlen <= 4) {
        sqrt_dbg++;
        mpz_t prod_ag;
        mpz_init_set_ui(prod_ag, 1);
        for (int i = 0; i < dlen; i++) {
            int ri = dep[i];
            /* Compute a*g(x) = val^2 - N */
            mpz_t ag;
            mpz_init(ag);
            mpz_mul(ag, rels[ri].val, rels[ri].val);
            mpz_sub(ag, ag, gN);
            gmp_fprintf(stderr, "  rel[%d]: val=%Zd, a*g(x)=%Zd\n", ri, rels[ri].val, ag);
            mpz_mul(prod_ag, prod_ag, ag);
            mpz_clear(ag);
        }
        gmp_fprintf(stderr, "  product(a*g)=%Zd\n", prod_ag);

        /* Now compute product from exponents */
        mpz_t prod_exp;
        mpz_init_set_ui(prod_exp, 1);
        if (total[0] % 2) mpz_neg(prod_exp, prod_exp); /* sign */
        for (int j = 1; j < fb_size; j++) {
            if (total[j] > 0) {
                mpz_t pp;
                mpz_init(pp);
                mpz_ui_pow_ui(pp, fb_p[j], total[j]);
                mpz_mul(prod_exp, prod_exp, pp);
                mpz_clear(pp);
            }
        }
        gmp_fprintf(stderr, "  product(p^e)=%Zd\n", prod_exp);
        if (mpz_cmp(prod_ag, prod_exp) != 0) {
            fprintf(stderr, "  MISMATCH: product of a*g(x) != product of p^exp\n");
        } else {
            fprintf(stderr, "  MATCH: products agree\n");
        }
        mpz_clear(prod_ag);
        mpz_clear(prod_exp);
    }

    /* Check evenness */
    int all_even = 1;
    for (int j = 0; j < fb_size; j++) {
        if (total[j] % 2 != 0) {
            fprintf(stderr, "  ODD exp at fb[%d]=%u: %d\n", j, j < fb_size ? fb_p[j] : 0, total[j]);
            all_even = 0;
        }
    }
    if (!all_even) {
        fprintf(stderr, "  Dependency has odd exponents, skipping\n");
        free(total);
        mpz_clears(X, Y, tmp, g, NULL);
        return 0;
    }

    /* Y = product of p^(e/2) */
    mpz_set_ui(Y, 1);
    for (int j = 1; j < fb_size; j++) {
        int half = total[j] / 2;
        if (half > 0) {
            mpz_ui_pow_ui(tmp, fb_p[j], half);
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, gN);
        }
    }

    /* Verify: X^2 ≡ Y^2 (mod N)? */
    mpz_t X2, Y2;
    mpz_inits(X2, Y2, NULL);
    mpz_mul(X2, X, X); mpz_mod(X2, X2, gN);
    mpz_mul(Y2, Y, Y); mpz_mod(Y2, Y2, gN);
    if (mpz_cmp(X2, Y2) != 0) {
        fprintf(stderr, "  ERROR: X^2 != Y^2 mod N!\n");
        gmp_fprintf(stderr, "  X=%Zd\n  Y=%Zd\n  X2=%Zd\n  Y2=%Zd\n", X, Y, X2, Y2);
    }
    mpz_clears(X2, Y2, NULL);

    /* gcd(X - Y, N) and gcd(X + Y, N) */
    mpz_sub(tmp, X, Y);
    mpz_gcd(g, tmp, gN);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
        mpz_t cof;
        mpz_init(cof);
        mpz_divexact(cof, gN, g);
        gmp_printf("%Zd %Zd\n", g, cof);
        mpz_clear(cof);
        free(total);
        mpz_clears(X, Y, tmp, g, NULL);
        return 1;
    }

    mpz_add(tmp, X, Y);
    mpz_gcd(g, tmp, gN);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
        mpz_t cof;
        mpz_init(cof);
        mpz_divexact(cof, gN, g);
        gmp_printf("%Zd %Zd\n", g, cof);
        mpz_clear(cof);
        free(total);
        mpz_clears(X, Y, tmp, g, NULL);
        return 1;
    }

    free(total);
    mpz_clears(X, Y, tmp, g, NULL);
    return 0;
}

/* ==================== Main ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec tstart;
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    mpz_init_set_str(gN, argv[1], 10);
    int bits = mpz_sizeinbase(gN, 2);
    int digits = mpz_sizeinbase(gN, 10);
    params_t par = get_params(bits);

    fprintf(stderr, "MPQS: %d digits (%d bits), FB=%d, M=%d, s=%d\n",
            digits, bits, par.fb_target, par.sieve_half, par.nfactors);

    build_fb(gN, par.fb_target);
    fprintf(stderr, "FB: %d primes, largest=%u\n", fb_size, fb_p[fb_size-1]);

    int target = fb_size + 30;

    sieve_arr = aligned_alloc(64, 2 * par.sieve_half + 64);
    soln1 = malloc(fb_size * sizeof(int));
    soln2 = malloc(fb_size * sizeof(int));
    mpz_inits(poly_a, poly_b, poly_c, NULL);
    memset(lp_hash, 0, sizeof(lp_hash));

    /* Compute sieve threshold */
    /* Theoretical: log2(sqrt(N) * M * sqrt(2/e)) */
    double log2N = bits;
    int threshold = (int)(log2N / 2.0 + log2(par.sieve_half)) - par.thresh_adj;
    if (threshold < 20) threshold = 20;
    fprintf(stderr, "Threshold: %d bits\n", threshold);

    poly_count = 0;
    int total_candidates = 0;
    int sieve_len = 2 * par.sieve_half;

    while (nfull < target) {
        /* Timeout check */
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - tstart.tv_sec) + (now.tv_nsec - tstart.tv_nsec) / 1e9;
        if (elapsed > 290.0) {
            fprintf(stderr, "TIMEOUT at %.1fs: %d full, %d total rels\n", elapsed, nfull, nrels);
            return 1;
        }

        choose_a(&par);
        poly_count++;
        compute_roots(par.sieve_half);

        /* Sieve the interval [-M, M) in blocks */
        int block_size = 32768;
        for (int block_start = 0; block_start < sieve_len; block_start += block_size) {
            int blen = block_size;
            if (block_start + blen > sieve_len) blen = sieve_len - block_start;

            sieve_block(block_start, blen);

            /* Scan for candidates */
            for (int j = 0; j < blen; j++) {
                if (sieve_arr[j] >= threshold) {
                    int x = block_start + j - par.sieve_half;
                    total_candidates++;
                    process_candidate(x, &par);
                }
            }
        }

        if (poly_count % 50 == 0) {
            struct timespec now2;
            clock_gettime(CLOCK_MONOTONIC, &now2);
            double el = (now2.tv_sec - tstart.tv_sec) + (now2.tv_nsec - tstart.tv_nsec) / 1e9;
            fprintf(stderr, "%.1fs: poly=%d cand=%d rels=%d/%d (full=%d)\n",
                    el, poly_count, total_candidates, nrels, target, nfull);
        }
    }

    struct timespec sieve_end;
    clock_gettime(CLOCK_MONOTONIC, &sieve_end);
    double sieve_time = (sieve_end.tv_sec - tstart.tv_sec) + (sieve_end.tv_nsec - tstart.tv_nsec) / 1e9;
    fprintf(stderr, "Sieve done: %.1fs, %d rels (%d full) from %d candidates, %d polys\n",
            sieve_time, nrels, nfull, total_candidates, poly_count);

    /* Linear algebra */
    int ndeps = gauss_elim(64);

    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* Verify first dependency: check that sum of exps mod 2 = 0 for each column */
    if (ndeps > 0) {
        int *check = calloc(fb_size, sizeof(int));
        for (int i = 0; i < g_dep_lens[0]; i++) {
            int ri = g_deps[i];
            if (!rels[ri].exps) { fprintf(stderr, "DEP[0] includes null rel %d!\n", ri); continue; }
            for (int j = 0; j < fb_size; j++)
                check[j] ^= rels[ri].exps[j];
        }
        int ok = 1;
        for (int j = 0; j < fb_size; j++) {
            if (check[j]) { fprintf(stderr, "DEP[0] NOT null at col %d\n", j); ok = 0; }
        }
        if (ok) fprintf(stderr, "DEP[0] verified: %d rels, null vector OK\n", g_dep_lens[0]);
        free(check);
    }

    for (int d = 0; d < ndeps; d++) {
        if (try_sqrt(&g_deps[d * nrels], g_dep_lens[d])) {
            struct timespec tend;
            clock_gettime(CLOCK_MONOTONIC, &tend);
            double total = (tend.tv_sec - tstart.tv_sec) + (tend.tv_nsec - tstart.tv_nsec) / 1e9;
            fprintf(stderr, "Total: %.3fs\n", total);
            return 0;
        }
    }

    fprintf(stderr, "FAIL: no factor from %d deps\n", ndeps);
    return 1;
}
