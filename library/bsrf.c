/*
 * bsrf.c - Quadratic Sieve with a=q² polynomials
 *
 * For polynomial Q_a(x) = ((q²·x+b)² - N) / q²
 * where q is prime, b² ≡ N (mod q²) (Hensel lifted), and q² ≈ sqrt(2N)/M.
 *
 * Since q² is a perfect square, the relation
 *   (q²·x+b)² ≡ q² · Q_a(x)  (mod N)
 * has q² contributing even exponents on the RHS. Each relation is independently
 * usable without pairing from the same polynomial.
 *
 * Sieve roots for prime p in factor base:
 *   x ≡ (-b ± r_p) · (q²)^{-1} (mod p)
 *   where r_p = sqrt(N mod p).
 *
 * Usage: ./bsrf <N>
 * Output: "p q" where N = p*q
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <ecm.h>

static int SIEVE_SIZE = 131072;
static int HALF_M     = 65536;

static mpz_t N_g, temp_g;

/* ===== Factor base ===== */

typedef struct {
    unsigned long p;
    unsigned long r;   /* sqrt(N mod p) */
    double logp;
} fb_entry_t;

static fb_entry_t *fb   = NULL;
static int         fb_size = 0;
static int         fb_cap  = 0;
static unsigned long smooth_bound;
static unsigned long max_lp_val;

/* ===== Modular arithmetic ===== */

static unsigned long mod_inv(unsigned long a, unsigned long m) {
    long or_ = (long)a, r = (long)m, os = 1, s = 0;
    while (r) {
        long q = or_ / r, t = r;
        r = or_ - q * r; or_ = t;
        t = s; s = os - q * s; os = t;
    }
    return os < 0 ? (unsigned long)(os + (long)m) : (unsigned long)os;
}

static int legendre_sym(unsigned long a, unsigned long p) {
    unsigned long val = a % p;
    if (val == 0) return 0;
    mpz_t base, ev, mod, res;
    mpz_init_set_ui(base, val); mpz_init_set_ui(mod, p);
    mpz_init(ev); mpz_init(res);
    mpz_sub_ui(ev, mod, 1); mpz_tdiv_q_2exp(ev, ev, 1);
    mpz_powm(res, base, ev, mod);
    int ret = (mpz_cmp_ui(res, 1) == 0) ? 1 : (mpz_cmp_ui(res, 0) == 0 ? 0 : -1);
    mpz_clear(base); mpz_clear(ev); mpz_clear(mod); mpz_clear(res);
    return ret;
}

static unsigned long sqrt_mod_p(unsigned long nv, unsigned long p) {
    if (p == 2) return nv & 1;
    nv %= p;
    if (nv == 0) return 0;
    if (p % 4 == 3) {
        mpz_t nn, pp, r, e;
        mpz_init_set_ui(nn, nv); mpz_init_set_ui(pp, p);
        mpz_init(r); mpz_init_set_ui(e, (p + 1) / 4);
        mpz_powm(r, nn, e, pp);
        unsigned long res = mpz_get_ui(r);
        mpz_clear(nn); mpz_clear(pp); mpz_clear(r); mpz_clear(e);
        return res;
    }
    unsigned long q = p - 1, s = 0;
    while (q % 2 == 0) { q /= 2; s++; }
    unsigned long z = 2;
    while (legendre_sym(z, p) != -1) z++;
    mpz_t pp, nn, Mv, cv, tv, Rv, bv, tmp;
    mpz_init_set_ui(pp, p); mpz_init_set_ui(nn, nv); mpz_init_set_ui(Mv, s);
    mpz_init(cv); mpz_init(tv); mpz_init(Rv); mpz_init(bv); mpz_init(tmp);
    mpz_set_ui(tmp, q); mpz_set_ui(cv, z); mpz_powm(cv, cv, tmp, pp);
    mpz_set_ui(tmp, (q + 1) / 2); mpz_powm(Rv, nn, tmp, pp);
    mpz_set_ui(tmp, q); mpz_powm(tv, nn, tmp, pp);
    while (mpz_cmp_ui(tv, 1) != 0) {
        mpz_set(tmp, tv);
        unsigned long i = 0;
        while (mpz_cmp_ui(tmp, 1) != 0) { mpz_mul(tmp, tmp, tmp); mpz_mod(tmp, tmp, pp); i++; }
        unsigned long mv = mpz_get_ui(Mv);
        mpz_set(bv, cv);
        for (unsigned long j = 0; j < mv - i - 1; j++) { mpz_mul(bv, bv, bv); mpz_mod(bv, bv, pp); }
        mpz_set_ui(Mv, i);
        mpz_mul(cv, bv, bv); mpz_mod(cv, cv, pp);
        mpz_mul(tv, tv, cv); mpz_mod(tv, tv, pp);
        mpz_mul(Rv, Rv, bv); mpz_mod(Rv, Rv, pp);
    }
    unsigned long res = mpz_get_ui(Rv);
    mpz_clear(pp); mpz_clear(nn); mpz_clear(Mv); mpz_clear(cv);
    mpz_clear(tv); mpz_clear(Rv); mpz_clear(bv); mpz_clear(tmp);
    return res;
}

static int build_factor_base(unsigned long bound) {
    char *sieve = calloc(bound + 1, 1);
    for (unsigned long i = 2; i <= bound; i++) sieve[i] = 1;
    for (unsigned long i = 2; i * i <= bound; i++)
        if (sieve[i])
            for (unsigned long j = i * i; j <= bound; j += i) sieve[j] = 0;
    fb_cap  = 16384;
    fb      = malloc(fb_cap * sizeof(fb_entry_t));
    fb_size = 0;
    fb[fb_size].p = 2; fb[fb_size].r = mpz_get_ui(N_g) & 1; fb[fb_size].logp = 1.0; fb_size++;
    for (unsigned long p = 3; p <= bound; p++) {
        if (!sieve[p]) continue;
        unsigned long nm = mpz_fdiv_ui(N_g, p);
        if (legendre_sym(nm, p) < 0) continue;
        if (fb_size >= fb_cap) { fb_cap *= 2; fb = realloc(fb, fb_cap * sizeof(fb_entry_t)); }
        fb[fb_size].p    = p;
        fb[fb_size].logp = log2((double)p);
        fb[fb_size].r    = (nm == 0) ? 0 : sqrt_mod_p(nm, p);
        fb_size++;
    }
    free(sieve);
    return fb_size;
}

/* ===== Bit matrix ===== */

typedef struct { unsigned long *rows; int nrows, ncols, wpr; } bm_t;

static bm_t *bm_alloc(int nr, int nc) {
    bm_t *m = malloc(sizeof(bm_t)); m->nrows = nr; m->ncols = nc;
    m->wpr = (nc + 63) / 64; m->rows = calloc((size_t)nr * m->wpr, 8); return m;
}
static void bm_free(bm_t *m) { free(m->rows); free(m); }
static inline void bm_flip(bm_t *m, int r, int c) {
    m->rows[(size_t)r * m->wpr + c / 64] ^= 1UL << (c % 64);
}
static inline int bm_get(bm_t *m, int r, int c) {
    return (m->rows[(size_t)r * m->wpr + c / 64] >> (c % 64)) & 1;
}
static inline void bm_xor_row(bm_t *m, int d, int s) {
    unsigned long *dp = m->rows + (size_t)d * m->wpr;
    unsigned long *sp = m->rows + (size_t)s * m->wpr;
    for (int i = 0; i < m->wpr; i++) dp[i] ^= sp[i];
}

static int gauss_elim(bm_t *M, int **deps_out, int *ndeps_out) {
    int nr = M->nrows, nc = M->ncols;
    bm_t *H = bm_alloc(nr, nr);
    for (int i = 0; i < nr; i++) bm_flip(H, i, i);
    int rank = 0;
    for (int c = 0; c < nc && rank < nr; c++) {
        int piv = -1;
        for (int r = rank; r < nr; r++) if (bm_get(M, r, c)) { piv = r; break; }
        if (piv < 0) continue;
        if (piv != rank) {
            size_t wsz = (size_t)M->wpr * 8;
            unsigned long *t = malloc(wsz);
            memcpy(t, M->rows + (size_t)piv * M->wpr, wsz);
            memcpy(M->rows + (size_t)piv * M->wpr, M->rows + (size_t)rank * M->wpr, wsz);
            memcpy(M->rows + (size_t)rank * M->wpr, t, wsz); free(t);
            wsz = (size_t)H->wpr * 8; t = malloc(wsz);
            memcpy(t, H->rows + (size_t)piv * H->wpr, wsz);
            memcpy(H->rows + (size_t)piv * H->wpr, H->rows + (size_t)rank * H->wpr, wsz);
            memcpy(H->rows + (size_t)rank * H->wpr, t, wsz); free(t);
        }
        for (int r = 0; r < nr; r++) if (r != rank && bm_get(M, r, c)) { bm_xor_row(M, r, rank); bm_xor_row(H, r, rank); }
        rank++;
    }
    *ndeps_out = 0; *deps_out = NULL;
    for (int r = 0; r < nr; r++) {
        int zero = 1;
        for (int w = 0; w < M->wpr && zero; w++) if (M->rows[(size_t)r * M->wpr + w]) zero = 0;
        if (!zero) continue;
        int cnt = 0; for (int j = 0; j < nr; j++) if (bm_get(H, r, j)) cnt++;
        if (cnt < 2) continue;
        *deps_out = realloc(*deps_out, (*ndeps_out + 1) * (nr + 1) * sizeof(int));
        int *dep  = *deps_out + (*ndeps_out) * (nr + 1); dep[0] = cnt; int idx = 1;
        for (int j = 0; j < nr; j++) if (bm_get(H, r, j)) dep[idx++] = j;
        (*ndeps_out)++;
        if (*ndeps_out >= 128) break;
    }
    bm_free(H); return rank;
}

/* ===== Relations ===== */

/*
 * For polynomial with parameter q (prime), a = q²:
 *   (q²·x + b)² ≡ q² · Q_a(x)  (mod N)
 *
 * Relation: x_lhs = (q²·x + b) mod N
 *           Q_a(x) is B-smooth (possibly with 1 LP)
 *
 * For GF(2) matrix: column 0 = sign, columns 1..fb_size = mod-2 exponents of Q_a(x) over fb.
 * For Y reconstruction: use full_exp (actual exponent of each fb prime in Q_a(x)).
 *
 * In a valid dependency, for each fb prime p:
 *   XOR of mod-2 exponents = 0  =>  sum of actual exponents is even  =>  half-exponent is integer.
 * Y = (product of q) * (product of p^(total_exp/2)) * (product of LP for LP-merged pairs).
 * X = product of x_lhs mod N.
 * Then gcd(X ± Y, N) gives a factor.
 */
typedef struct {
    unsigned long  q;
    mpz_t          x_lhs;
    int            sign;
    unsigned char *exponents;    /* mod-2 exponents over fb */
    unsigned short *full_exp;   /* actual exponents over fb (for Y reconstruction) */
    int            nlp;
    unsigned long  lp;
} rel_t;

static rel_t *all_rels     = NULL;
static int    n_all_rels   = 0;
static int    all_rels_cap = 0;

/* "Full" matrix relation: direct (nlp==0) or LP-merged pair */
typedef struct { int ri1, ri2; } full_rel_t;

static full_rel_t *full_rels     = NULL;
static int         n_full_rels   = 0;
static int         full_rels_cap = 0;

static void add_full_rel(int ri1, int ri2) {
    if (n_full_rels >= full_rels_cap) {
        full_rels_cap = full_rels_cap ? full_rels_cap * 2 : 4096;
        full_rels = realloc(full_rels, full_rels_cap * sizeof(full_rel_t));
    }
    full_rels[n_full_rels].ri1 = ri1;
    full_rels[n_full_rels].ri2 = ri2;
    n_full_rels++;
}

/* LP hash table */
#define LP_HASH_SIZE (1 << 20)
typedef struct lp_node { unsigned long lp; int ri; struct lp_node *next; } lp_node_t;
static lp_node_t **lp_hash = NULL;

static void lp_hash_init(void) { lp_hash = calloc(LP_HASH_SIZE, sizeof(lp_node_t *)); }
static void lp_hash_free(void) {
    for (int i = 0; i < LP_HASH_SIZE; i++) {
        lp_node_t *n = lp_hash[i];
        while (n) { lp_node_t *nx = n->next; free(n); n = nx; }
    }
    free(lp_hash); lp_hash = NULL;
}
static int lp_lookup_insert(unsigned long lp, int ri) {
    unsigned int h = (unsigned int)((lp * 2654435761UL) >> 12) & (LP_HASH_SIZE - 1);
    for (lp_node_t *n = lp_hash[h]; n; n = n->next) if (n->lp == lp) return n->ri;
    lp_node_t *nn = malloc(sizeof(lp_node_t)); nn->lp = lp; nn->ri = ri; nn->next = lp_hash[h]; lp_hash[h] = nn;
    return -1;
}

static void collect_candidate(unsigned long q_val, mpz_t b_val, int x_off, mpz_t Q_mpz) {
    if (n_all_rels >= all_rels_cap) {
        all_rels_cap = all_rels_cap ? all_rels_cap * 2 : 8192;
        all_rels = realloc(all_rels, all_rels_cap * sizeof(rel_t));
    }
    rel_t *r = &all_rels[n_all_rels];
    r->q         = q_val;
    r->nlp       = 0; r->lp = 0;
    r->exponents = calloc(fb_size, 1);
    r->full_exp  = calloc(fb_size, sizeof(unsigned short));
    mpz_init(r->x_lhs);

    /* x_lhs = (q²·x + b) mod N */
    mpz_t a_mpz; mpz_init_set_ui(a_mpz, q_val); mpz_mul(a_mpz, a_mpz, a_mpz);
    mpz_set_si(r->x_lhs, x_off);
    mpz_mul(r->x_lhs, r->x_lhs, a_mpz);
    mpz_add(r->x_lhs, r->x_lhs, b_val);
    mpz_mod(r->x_lhs, r->x_lhs, N_g);
    mpz_clear(a_mpz);

    /* Factor |Q_mpz| over fb */
    mpz_t rem; mpz_init(rem); mpz_abs(rem, Q_mpz);
    r->sign = (mpz_sgn(Q_mpz) < 0) ? 1 : 0;
    int done_early = 0;
    for (int i = 0; i < fb_size; i++) {
        unsigned long p = fb[i].p; int exp = 0;
        while (mpz_divisible_ui_p(rem, p)) { mpz_divexact_ui(rem, rem, p); exp++; }
        r->exponents[i] = exp & 1;
        r->full_exp[i]  = (unsigned short)(exp < 65535 ? exp : 65535);
        /* Check early exit every 32 primes to amortize mpz overhead */
        if (exp > 0 && (i & 31) == 31) {
            if (mpz_cmp_ui(rem, 1) == 0) { done_early = 1; break; }
            if (mpz_sizeinbase(rem, 2) <= 32) { done_early = 1; break; }
        }
    }

    if (mpz_cmp_ui(rem, 1) == 0) {
        /* Fully smooth */
        n_all_rels++;
        add_full_rel(n_all_rels - 1, -1);
        mpz_clear(rem); return;
    }
    if (mpz_fits_ulong_p(rem)) {
        unsigned long cof = mpz_get_ui(rem);
        if (cof <= max_lp_val && mpz_probab_prime_p(rem, 10)) {
            r->nlp = 1; r->lp = cof;
            n_all_rels++;
            int other = lp_lookup_insert(cof, n_all_rels - 1);
            if (other >= 0) add_full_rel(other, n_all_rels - 1);
            mpz_clear(rem); return;
        }
    }
    /* Not usable */
    free(r->exponents); r->exponents = NULL;
    free(r->full_exp);  r->full_exp  = NULL;
    mpz_clear(r->x_lhs);
    mpz_clear(rem);
}

/* Try a dependency: fr_indices[0..cnt-1] are indices into full_rels[]. */
static int try_full_rels(int *fr_indices, int cnt) {
    mpz_t X, Y, g;
    mpz_init_set_ui(X, 1); mpz_init_set_ui(Y, 1); mpz_init(g);
    int total_sign = 0;

    /* Accumulate exponent sums using int to avoid overflow */
    int *exp_sum = calloc(fb_size, sizeof(int));

    for (int k = 0; k < cnt; k++) {
        full_rel_t *fr = &full_rels[fr_indices[k]];
        int rlist[2] = {fr->ri1, fr->ri2};
        int rcount   = (fr->ri2 >= 0) ? 2 : 1;
        for (int j = 0; j < rcount; j++) {
            rel_t *r = &all_rels[rlist[j]];
            mpz_mul(X, X, r->x_lhs); mpz_mod(X, X, N_g);
            total_sign ^= r->sign;
            mpz_mul_ui(Y, Y, r->q); mpz_mod(Y, Y, N_g);
            for (int i = 0; i < fb_size; i++) exp_sum[i] += r->full_exp[i];
        }
    }

    int ret = 0;
    if (total_sign != 0) goto done;

    /* Y *= product p^(exp_sum[i]/2); exponents must all be even */
    for (int i = 0; i < fb_size; i++) {
        if (exp_sum[i] == 0) continue;
        if (exp_sum[i] & 1) goto done;  /* should not happen for valid dep */
        int half = exp_sum[i] / 2;
        for (int e = 0; e < half; e++) { mpz_mul_ui(Y, Y, fb[i].p); mpz_mod(Y, Y, N_g); }
    }

    /* For LP-merged full_rels: multiply Y by lp (lp appears twice total, sqrt = lp^1) */
    for (int k = 0; k < cnt; k++) {
        full_rel_t *fr = &full_rels[fr_indices[k]];
        if (fr->ri2 < 0) continue;
        unsigned long lp = all_rels[fr->ri1].lp;
        mpz_mul_ui(Y, Y, lp); mpz_mod(Y, Y, N_g);
    }

    /* gcd(X ± Y, N) */
    mpz_sub(temp_g, X, Y); mpz_gcd(g, temp_g, N_g);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_g) < 0) {
        mpz_t cf; mpz_init(cf); mpz_divexact(cf, N_g, g);
        if (mpz_cmp(g, cf) > 0) mpz_swap(g, cf);
        gmp_printf("%Zd %Zd\n", g, cf); mpz_clear(cf); ret = 1; goto done;
    }
    mpz_add(temp_g, X, Y); mpz_gcd(g, temp_g, N_g);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_g) < 0) {
        mpz_t cf; mpz_init(cf); mpz_divexact(cf, N_g, g);
        if (mpz_cmp(g, cf) > 0) mpz_swap(g, cf);
        gmp_printf("%Zd %Zd\n", g, cf); mpz_clear(cf); ret = 1; goto done;
    }
done:
    free(exp_sum); mpz_clear(X); mpz_clear(Y); mpz_clear(g);
    return ret;
}

/* ===== Polynomial generation: compute q and b for a=q² ===== */
/*
 * Choose q such that q² ≈ sqrt(2N) / M
 * b = Hensel lift: b0 = sqrt(N mod q), b = b0 + q·((N-b0²)/q · inv(2b0) mod q)
 * Center: if b > q²/2, use b = q²-b.
 */
static int gen_poly(unsigned long q_target, mpz_t b_out, unsigned long *q_out) {
    mpz_t cand; mpz_init_set_ui(cand, q_target);
    if (mpz_get_ui(cand) % 2 == 0) mpz_add_ui(cand, cand, 1);
    for (int iters = 0; iters < 10000; iters++) {
        if (mpz_probab_prime_p(cand, 10) > 0) {
            unsigned long q = mpz_get_ui(cand);
            unsigned long nm = mpz_fdiv_ui(N_g, q);
            if (nm > 0 && legendre_sym(nm, q) == 1) {
                unsigned long b0 = sqrt_mod_p(nm, q);
                /* Hensel lift: b = b0 + q * ((N-b0²)/q * inv(2b0) mod q) */
                mpz_t nb0sq; mpz_init(nb0sq);
                mpz_set_ui(nb0sq, b0); mpz_mul(nb0sq, nb0sq, nb0sq);
                mpz_sub(nb0sq, N_g, nb0sq);
                mpz_divexact_ui(nb0sq, nb0sq, q);
                unsigned long two_b0 = (2 * b0) % q;
                if (two_b0 == 0) { mpz_clear(nb0sq); mpz_add_ui(cand, cand, 2); iters++; continue; }
                unsigned long inv2b0 = mod_inv(two_b0, q);
                unsigned long lift = (unsigned long)(mpz_fdiv_ui(nb0sq, q) * (unsigned long long)inv2b0 % q);
                mpz_clear(nb0sq);
                /* b = b0 + q*lift */
                mpz_set_ui(b_out, q); mpz_mul_ui(b_out, b_out, lift); mpz_add_ui(b_out, b_out, b0);
                /* Verify b² ≡ N (mod q²) */
                mpz_t q2, bsq, chk, Nq2;
                mpz_init_set_ui(q2, q); mpz_mul(q2, q2, q2);
                mpz_init(bsq); mpz_mul(bsq, b_out, b_out);
                mpz_init(chk); mpz_mod(chk, bsq, q2);
                mpz_init(Nq2); mpz_mod(Nq2, N_g, q2);
                int ok = (mpz_cmp(chk, Nq2) == 0);
                if (!ok) {
                    /* Try complementary root: b = q²-b */
                    mpz_sub(b_out, q2, b_out);
                    mpz_mul(bsq, b_out, b_out); mpz_mod(chk, bsq, q2);
                    ok = (mpz_cmp(chk, Nq2) == 0);
                }
                mpz_clear(q2); mpz_clear(bsq); mpz_clear(chk); mpz_clear(Nq2);
                if (ok) {
                    /* Center: if b > q²/2, use q²-b */
                    mpz_t q2c, hq2;
                    mpz_init_set_ui(q2c, q); mpz_mul(q2c, q2c, q2c);
                    mpz_init(hq2); mpz_tdiv_q_2exp(hq2, q2c, 1);
                    if (mpz_cmp(b_out, hq2) > 0) mpz_sub(b_out, q2c, b_out);
                    mpz_clear(q2c); mpz_clear(hq2);
                    *q_out = q; mpz_clear(cand); return 1;
                }
            }
        }
        mpz_add_ui(cand, cand, 2);
    }
    mpz_clear(cand); return 0;
}

/* ===== Main factoring ===== */

static int factor_bsrf(void) {
    size_t n_bits   = mpz_sizeinbase(N_g, 2);
    size_t n_digits = mpz_sizeinbase(N_g, 10);

    double ln_n    = (double)n_bits * log(2.0);
    double ln_ln_n = log(ln_n);
    smooth_bound   = (unsigned long)exp(0.5 * sqrt(ln_n * ln_ln_n));
    if (smooth_bound < 500)     smooth_bound = 500;
    if (smooth_bound > 8000000) smooth_bound = 8000000;

    max_lp_val = smooth_bound * 100UL;
    if (max_lp_val > (1UL << 32)) max_lp_val = (1UL << 32);

    fprintf(stderr, "BSRF a=q²: %zu-digit (%zu-bit), B=%lu, LP_max=%lu\n",
            n_digits, n_bits, smooth_bound, max_lp_val);

    int fb_count = build_factor_base(smooth_bound);
    fprintf(stderr, "Factor base: %d primes\n", fb_count);

    /* Adaptive sieve size: larger for bigger numbers to find more smooths per poly */
    if (n_digits >= 65) {
        SIEVE_SIZE = 4 * 1024 * 1024;   /* 4M for 65+ digits */
    } else if (n_digits >= 60) {
        SIEVE_SIZE = 2 * 1024 * 1024;   /* 2M for 60-64 digits */
    } else if (n_digits >= 55) {
        SIEVE_SIZE = 1024 * 1024;        /* 1M for 55-59 digits */
    } else if (n_digits >= 45) {
        SIEVE_SIZE = 512 * 1024;         /* 512K for 45-54 digits */
    } else {
        SIEVE_SIZE = 131072;             /* 128K for smaller */
    }
    HALF_M = SIEVE_SIZE / 2;

    /* q² ≈ sqrt(2N)/M  =>  q ≈ N^{1/4} / sqrt(M/2) */
    double N_d     = mpz_get_d(N_g);
    double q_ideal = pow(N_d, 0.25) / sqrt((double)SIEVE_SIZE / 2.0);
    if (q_ideal < 10.0) q_ideal = 10.0;

    fprintf(stderr, "Target q ≈ %.3e (a=q²≈%.3e, M=%d)\n",
            q_ideal, q_ideal * q_ideal, SIEVE_SIZE);

    lp_hash_init();
    float *sieve_log = malloc(SIEVE_SIZE * sizeof(float));
    if (!sieve_log) return 0;

    /* LP allowance in sieve threshold */
    float lp_allow = (float)log2((double)max_lp_val);

    int target_full = fb_count + 200;
    if (target_full > 500000) target_full = 500000;
    fprintf(stderr, "Target full relations: %d\n", target_full);

    int polys_used = 0;
    unsigned long q_up   = (unsigned long)q_ideal;
    unsigned long q_down = (unsigned long)q_ideal;
    int go_up = 1;

    while (n_full_rels < target_full && polys_used < 1000000) {
        unsigned long q_try;
        unsigned long step = (unsigned long)(q_ideal * 0.0005 + 3);
        if (go_up) {
            q_try = q_up; q_up += step;
        } else {
            if (q_down <= step + 3) { go_up = !go_up; continue; }
            q_try = q_down; q_down -= step;
        }
        go_up = !go_up;

        mpz_t b_val; mpz_init(b_val);
        unsigned long q_val = 0;
        if (!gen_poly(q_try, b_val, &q_val)) { mpz_clear(b_val); continue; }
        polys_used++;

        /* a = q², c = (b²-N)/q² */
        mpz_t q2_mpz; mpz_init_set_ui(q2_mpz, q_val); mpz_mul(q2_mpz, q2_mpz, q2_mpz);
        mpz_t c_coeff; mpz_init(c_coeff);
        mpz_mul(c_coeff, b_val, b_val); mpz_sub(c_coeff, c_coeff, N_g);
        mpz_divexact(c_coeff, c_coeff, q2_mpz);

        /* Sieve threshold */
        double ad   = (double)q_val * (double)q_val;
        double bd   = mpz_get_d(b_val);
        double hm   = (double)HALF_M;
        double edge = ad * hm * hm + 2.0 * fabs(bd) * hm + fabs(mpz_get_d(c_coeff));
        float cur_thr = (float)log2(edge) - lp_allow;
        if (cur_thr < 2.0f) cur_thr = 2.0f;

        memset(sieve_log, 0, SIEVE_SIZE * sizeof(float));

        for (int fi = 0; fi < fb_size; fi++) {
            unsigned long p   = fb[fi].p;
            float logp        = (float)fb[fi].logp;
            unsigned long r_p = fb[fi].r;

            if (p == 2) {
                /* Q_a(x) = q²x² + 2bx + c; mod 2: (q²≡1 if q odd) x² + c mod 2.
                 * x² ≡ x mod 2, so Q_a ≡ x + c mod 2. Root: x ≡ c mod 2. */
                unsigned long cmod2 = (unsigned long)(mpz_get_ui(c_coeff) & 1);
                long start = ((long)cmod2 + HALF_M) % 2;
                if (start < 0) start += 2;
                for (long pos = start; pos < SIEVE_SIZE; pos += 2) sieve_log[pos] += logp;
                continue;
            }

            unsigned long qmodp  = q_val % p;
            unsigned long q2modp = (unsigned long)((unsigned long long)qmodp * qmodp % p);
            unsigned long bmodp  = mpz_fdiv_ui(b_val, p);

            if (q2modp == 0) {
                /* p | q: linear: 2b·x + c ≡ 0 (mod p) */
                unsigned long two_b = (2 * bmodp) % p;
                if (two_b == 0) continue;
                unsigned long cmodp = mpz_fdiv_ui(c_coeff, p);
                unsigned long x0    = (unsigned long)((p - cmodp % p) % p * (unsigned long long)mod_inv(two_b, p) % p);
                long st = ((long)x0 + HALF_M) % (long)p; if (st < 0) st += (long)p;
                for (long pos = st; pos < SIEVE_SIZE; pos += (long)p) sieve_log[pos] += logp;
                continue;
            }

            /* Two roots: x ≡ (-b ± r_p) · (q²)^{-1} (mod p) */
            unsigned long inv_q2 = mod_inv(q2modp, p);
            /* root 1: (-b + r_p) mod p */
            unsigned long diff1 = (r_p + p - bmodp) % p;
            unsigned long x1    = (unsigned long)((unsigned long long)diff1 * inv_q2 % p);
            long st1 = ((long)x1 + HALF_M) % (long)p; if (st1 < 0) st1 += (long)p;
            for (long pos = st1; pos < SIEVE_SIZE; pos += (long)p) sieve_log[pos] += logp;

            if (r_p == 0) continue;  /* double root */

            /* root 2: (-b - r_p) mod p */
            unsigned long diff2 = (p - r_p + p - bmodp) % p;
            unsigned long x2    = (unsigned long)((unsigned long long)diff2 * inv_q2 % p);
            if (x2 == x1) continue;
            long st2 = ((long)x2 + HALF_M) % (long)p; if (st2 < 0) st2 += (long)p;
            for (long pos = st2; pos < SIEVE_SIZE; pos += (long)p) sieve_log[pos] += logp;
        }

        /* Scan sieve */
        for (int i = 0; i < SIEVE_SIZE; i++) {
            if (sieve_log[i] < cur_thr) continue;
            int x = i - HALF_M;
            /* Q_a(x) = q²x² + 2bx + c */
            mpz_t qv, bv2; mpz_init(qv); mpz_init_set(bv2, b_val); mpz_mul_ui(bv2, bv2, 2);
            mpz_set_si(qv, x); mpz_mul(qv, qv, q2_mpz); mpz_add(qv, qv, bv2);
            mpz_mul_si(qv, qv, x); mpz_add(qv, qv, c_coeff);
            mpz_clear(bv2);
            collect_candidate(q_val, b_val, x, qv);
            mpz_clear(qv);
        }

        mpz_clear(q2_mpz); mpz_clear(c_coeff); mpz_clear(b_val);

        if (polys_used % 500 == 0)
            fprintf(stderr, "  polys=%d full_rels=%d all_rels=%d\n", polys_used, n_full_rels, n_all_rels);
    }

    free(sieve_log);
    fprintf(stderr, "Collected: %d full_rels, %d all_rels, polys=%d\n", n_full_rels, n_all_rels, polys_used);

    if (n_full_rels < fb_count + 10) {
        fprintf(stderr, "Too few full relations (%d)\n", n_full_rels);
        lp_hash_free(); return 0;
    }

    /* Build GF(2) matrix: rows = full_rels, cols = sign + fb primes.
     * Use bm_flip (XOR toggle) so LP-merged pairs combine correctly. */
    int nrows = n_full_rels, ncols = 1 + fb_size;
    fprintf(stderr, "Building matrix %d x %d...\n", nrows, ncols);
    bm_t *M = bm_alloc(nrows, ncols);
    for (int row = 0; row < nrows; row++) {
        full_rel_t *fr = &full_rels[row];
        int rlist[2] = {fr->ri1, fr->ri2};
        int rcount   = (fr->ri2 >= 0) ? 2 : 1;
        for (int j = 0; j < rcount; j++) {
            rel_t *r = &all_rels[rlist[j]];
            if (r->sign) bm_flip(M, row, 0);                          /* XOR sign */
            for (int i = 0; i < fb_size; i++)
                if (r->exponents[i]) bm_flip(M, row, 1 + i);          /* XOR exponents */
        }
    }

    int *deps = NULL, ndeps = 0;
    fprintf(stderr, "Gaussian elimination %d x %d...\n", nrows, ncols);
    gauss_elim(M, &deps, &ndeps);
    fprintf(stderr, "Dependencies found: %d\n", ndeps);

    int factored = 0;
    for (int d = 0; d < ndeps && !factored; d++) {
        int *dep = deps + d * (nrows + 1);
        factored = try_full_rels(dep + 1, dep[0]);
    }

    free(deps); bm_free(M); lp_hash_free();
    for (int i = 0; i < n_all_rels; i++) {
        mpz_clear(all_rels[i].x_lhs);
        free(all_rels[i].exponents);
        free(all_rels[i].full_exp);
    }
    free(all_rels); all_rels = NULL; n_all_rels = 0; all_rels_cap = 0;
    free(full_rels); full_rels = NULL; n_full_rels = 0; full_rels_cap = 0;
    free(fb); fb = NULL; fb_size = 0; fb_cap = 0;
    return factored;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    mpz_init(N_g); mpz_init(temp_g);
    if (mpz_set_str(N_g, argv[1], 10) != 0) { fprintf(stderr, "Invalid: %s\n", argv[1]); return 1; }
    /* Trial division up to 1M */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N_g, p)) {
            mpz_t f, cof; mpz_init_set_ui(f, p); mpz_init(cof); mpz_divexact(cof, N_g, f);
            if (mpz_cmp_ui(cof, 1) > 0) { gmp_printf("%Zd %Zd\n", f, cof); mpz_clear(f); mpz_clear(cof); goto cleanup; }
            mpz_clear(f); mpz_clear(cof);
        }
    }
    if (factor_bsrf()) goto cleanup;
    fprintf(stderr, "BSRF failed\n");
    mpz_clear(N_g); mpz_clear(temp_g); return 1;
cleanup:
    mpz_clear(N_g); mpz_clear(temp_g); return 0;
}
