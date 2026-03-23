/*
 * MCFRAC: Multi-Continued-Fraction Factoring
 *
 * Novel approach: Instead of using the single CF expansion of sqrt(N),
 * use MULTIPLE CF expansions simultaneously:
 *   sqrt(k*N) for various small multipliers k
 *
 * Each expansion generates small residues |p_i^2 - k*N*q_i^2|.
 * These residues from different k can be combined in a single GF(2) matrix
 * because they all produce congruences mod N.
 *
 * Key insight: different multipliers k expose different factor base primes
 * (only primes where k*N is a QR). By using multiple k simultaneously,
 * we effectively use a larger factor base and find more smooth residues.
 *
 * Hypothesis: the combined smooth-finding rate from M multipliers
 * should exceed M times the single-multiplier rate because:
 * 1. Different k give different "easy" primes
 * 2. The combined FB is larger, increasing smoothness probability
 * 3. CF residues are naturally small (bounded by 2*sqrt(kN))
 *
 * Compile: gcc -O3 -march=native -o mcfrac library/mcfrac.c -lgmp -lm
 * Usage: ./mcfrac <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

#define SEED 42
#define MAX_FB 50000
#define MAX_RELS 200000
#define MAX_PARTIALS 500000
#define MAX_MULTIPLIERS 32

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static uint32_t sqrt_mod_p(uint32_t n, uint32_t p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    uint64_t b = n % p, e = (p - 1) / 2, r = 1, m = p;
    { uint64_t bb = b, ee = e; while (ee) { if (ee & 1) r = (r * bb) % m; bb = (bb * bb) % m; ee >>= 1; } }
    if (r != 1) return 0;
    if (p % 4 == 3) { b = n % p; e = (p + 1) / 4; r = 1; while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; } return (uint32_t)r; }
    uint32_t Q = p - 1, S = 0; while (Q % 2 == 0) { Q /= 2; S++; }
    uint32_t z = 2; while (1) { b = z; e = (p-1)/2; r = 1; while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; } if (r == (uint64_t)(p-1)) break; z++; }
    uint64_t M_val = S; b = z; e = Q; uint64_t c = 1; while (e) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = Q; uint64_t t = 1; while (e) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; uint64_t R = 1; while (e) { if (e & 1) R = (R * b) % m; b = (b * b) % m; e >>= 1; }
    while (1) { if (t == 1) return (uint32_t)R; int i = 0; uint64_t tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; }
        uint64_t bb = c; for (int j = 0; j < (int)M_val - i - 1; j++) bb = (bb * bb) % p; M_val = i; c = (bb * bb) % p; t = (t * c) % p; R = (R * bb) % p; }
}

/* ==================== Factor Base ==================== */
typedef struct {
    uint32_t *prime;
    uint8_t *logp;
    int size;
} fb_t;

/* Build unified FB for multiple multipliers */
static fb_t *fb_create_multi(mpz_t N, int *multipliers, int num_mult, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(uint32_t));
    fb->logp = malloc(alloc * sizeof(uint8_t));
    fb->prime[0] = 2; fb->logp[0] = 1; fb->size = 1;

    int bound = target * 20 + 100000;
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i * i <= bound; i++)
        if (!sv[i]) for (int j = i * i; j <= bound; j += i) sv[j] = 1;

    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        /* Accept prime if kN is QR mod p for ANY multiplier k */
        int accept = 0;
        for (int ki = 0; ki < num_mult && !accept; ki++) {
            mpz_t kN; mpz_init(kN);
            mpz_mul_ui(kN, N, multipliers[ki]);
            unsigned long nm = mpz_fdiv_ui(kN, i);
            if (nm == 0 || sqrt_mod_p((uint32_t)nm, i)) accept = 1;
            mpz_clear(kN);
        }
        if (!accept) continue;
        fb->prime[fb->size] = i;
        fb->logp[fb->size] = (uint8_t)(log2(i) + 0.5);
        fb->size++;
    }
    free(sv);
    return fb;
}

/* ==================== LP Hash ==================== */
#define LP_HASH_BITS 20
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { uint64_t lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t = calloc(1, sizeof(lp_t)); t->b = calloc(LP_HASH_SIZE, sizeof(lp_e_t*)); t->pool = calloc(m, sizeof(lp_e_t)); t->max = m; return t; }
static int lp_find(lp_t *t, uint64_t lp) { uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); for (lp_e_t *e = t->b[h]; e; e = e->next) if (e->lp == lp) return e->idx; return -1; }
static void lp_insert(lp_t *t, uint64_t lp, int idx) { if (t->used >= t->max) return; uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); lp_e_t *e = &t->pool[t->used++]; e->lp = lp; e->idx = idx; e->next = t->b[h]; t->b[h] = e; }

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *value;  /* the x value (product of p_i or ax+b values) */
    mpz_t *cofact; /* the Q value (product) */
    uint64_t *lp;
    int count, alloc;
} rels_t;
static rels_t *rels_create(int n) {
    rels_t *r = malloc(sizeof(rels_t));
    r->value = malloc(n * sizeof(mpz_t));
    r->cofact = malloc(n * sizeof(mpz_t));
    r->lp = calloc(n, sizeof(uint64_t));
    for (int i = 0; i < n; i++) { mpz_init(r->value[i]); mpz_init(r->cofact[i]); }
    r->count = 0; r->alloc = n;
    return r;
}

/* ==================== GF(2) Matrix ==================== */
typedef uint64_t u64;
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
        int pr = -1; for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) { if (r == piv) continue; if ((m->rows[r][c/64] >> (c%64)) & 1) for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w]; }
        piv++;
    }
    int nd = 0; *deps = malloc(max * sizeof(int*)); *dlen = malloc(max * sizeof(int));
    for (int r = piv; r < m->nr && nd < max; r++) {
        int z = 1; for (int w = 0; w < m->fbw && z; w++) { u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL<<(m->nc%64))-1); if (m->rows[r][w] & mask) z = 0; }
        if (!z) continue;
        int *d = malloc(m->nr * sizeof(int)); int dl = 0;
        for (int w = 0; w < m->idw; w++) { u64 bits = m->rows[r][m->fbw+w]; while (bits) { int bit = __builtin_ctzll(bits); int idx = w*64+bit; if (idx < m->nr) d[dl++] = idx; bits &= bits-1; } }
        if (dl > 0) { (*deps)[nd] = d; (*dlen)[nd] = dl; nd++; } else free(d);
    }
    return nd;
}

/* ==================== CF expansion state ==================== */
typedef struct {
    int k;           /* multiplier */
    mpz_t kN;        /* k * N */
    mpz_t sqrtkN;    /* floor(sqrt(kN)) */
    mpz_t P, Q0, Q1; /* CF state */
    mpz_t p_prev, p_curr;  /* convergent numerators */
    int step;        /* iteration count */
} cf_state_t;

static void cf_init(cf_state_t *s, mpz_t N, int k) {
    s->k = k;
    mpz_init(s->kN); mpz_mul_ui(s->kN, N, k);
    mpz_init(s->sqrtkN); mpz_sqrt(s->sqrtkN, s->kN);
    mpz_init_set(s->P, s->sqrtkN);
    mpz_init_set_ui(s->Q0, 1);
    mpz_init(s->Q1);
    mpz_mul(s->Q1, s->sqrtkN, s->sqrtkN);
    mpz_sub(s->Q1, s->kN, s->Q1);
    if (mpz_sgn(s->Q1) == 0) mpz_set_ui(s->Q1, 1); /* perfect square */
    mpz_init_set_ui(s->p_prev, 1);
    mpz_init_set(s->p_curr, s->sqrtkN);
    s->step = 0;
}

/* Advance CF one step. Returns the Q value (residue) and the convergent p */
static void cf_step(cf_state_t *s, mpz_t Q_out, mpz_t p_out) {
    /* b = floor((sqrtkN + P) / Q1) */
    mpz_t b, Pnew, Qnew, p_new;
    mpz_inits(b, Pnew, Qnew, p_new, NULL);

    mpz_add(b, s->sqrtkN, s->P);
    mpz_tdiv_q(b, b, s->Q1);

    /* P_new = b * Q1 - P */
    mpz_mul(Pnew, b, s->Q1);
    mpz_sub(Pnew, Pnew, s->P);

    /* Q_new = Q0 + b * (P - P_new) */
    mpz_sub(Qnew, s->P, Pnew);
    mpz_mul(Qnew, Qnew, b);
    mpz_add(Qnew, Qnew, s->Q0);

    /* p_new = b * p_curr + p_prev */
    mpz_mul(p_new, b, s->p_curr);
    mpz_add(p_new, p_new, s->p_prev);

    /* Update state */
    mpz_set(s->Q0, s->Q1);
    mpz_set(s->Q1, Qnew);
    mpz_set(s->P, Pnew);
    mpz_set(s->p_prev, s->p_curr);
    mpz_set(s->p_curr, p_new);
    s->step++;

    /* Output: Q = Q1 (the denominator), p = p_curr (convergent numerator) */
    /* Actually, the residue is (-1)^step * Q1 such that p_curr^2 ≡ (-1)^step * Q1 (mod kN) */
    mpz_set(Q_out, s->Q1);
    mpz_set(p_out, s->p_curr);

    mpz_clears(b, Pnew, Qnew, p_new, NULL);
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size;
    int lp_mult;
    int extra;
    int num_multipliers;
} params_t;

static params_t get_params(int bits) {
    if (bits <= 100) return (params_t){100, 30, 40, 8};
    if (bits <= 120) return (params_t){200, 40, 50, 12};
    if (bits <= 140) return (params_t){400, 50, 60, 16};
    if (bits <= 160) return (params_t){800, 60, 80, 20};
    if (bits <= 180) return (params_t){1500, 70, 100, 24};
    if (bits <= 200) return (params_t){3000, 80, 120, 28};
    if (bits <= 220) return (params_t){6000, 90, 150, 32};
    if (bits <= 240) return (params_t){12000, 100, 200, 32};
    if (bits <= 260) return (params_t){25000, 110, 250, 32};
    if (bits <= 280) return (params_t){45000, 120, 300, 32};
    return (params_t){70000, 140, 400, 32};
}

/* ==================== Main ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Trial division */
    for (unsigned long p = 2; p < 100000; p++) {
        if (p > 2 && p % 2 == 0) continue;
        if (mpz_divisible_ui_p(N, p)) {
            printf("%lu\n", p);
            return 0;
        }
    }

    params_t P = get_params(bits);

    /* Select multipliers: squarefree numbers 1,2,3,5,6,7,10,... */
    static const int candidate_mults[] = {
        1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,
        33,34,35,37,38,39,41,42,43,46,47,51,53,55,57,58,59,0
    };
    int multipliers[MAX_MULTIPLIERS];
    int num_mult = 0;
    for (int i = 0; candidate_mults[i] && num_mult < P.num_multipliers; i++) {
        multipliers[num_mult++] = candidate_mults[i];
    }

    fb_t *fb = fb_create_multi(N, multipliers, num_mult, P.fb_size);
    uint64_t lp_bound = (uint64_t)fb->prime[fb->size - 1] * P.lp_mult;
    int target = fb->size + P.extra;

    fprintf(stderr, "MCFRAC: %dd (%db), FB=%d, mults=%d, LP=%lu, target=%d\n",
            digits, bits, fb->size, num_mult, (unsigned long)lp_bound, target);

    /* Initialize CF expansions */
    cf_state_t *cfs = malloc(num_mult * sizeof(cf_state_t));
    for (int i = 0; i < num_mult; i++)
        cf_init(&cfs[i], N, multipliers[i]);

    /* Relation storage */
    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_t *lpt = lp_create(MAX_PARTIALS);

    mpz_t Q_val, p_val, residue, tmp;
    mpz_inits(Q_val, p_val, residue, tmp, NULL);

    int total_steps = 0, combined = 0;

    /* Main loop: round-robin through CF expansions */
    while (full->count < target) {
        if (elapsed() > 280) {
            fprintf(stderr, "TIMEOUT at %.1fs\n", elapsed());
            break;
        }

        for (int mi = 0; mi < num_mult && full->count < target; mi++) {
            cf_step(&cfs[mi], Q_val, p_val);
            total_steps++;

            /* The relation is: p_val^2 ≡ (-1)^step * Q_val (mod k*N)
             * Which gives: p_val^2 ≡ (-1)^step * Q_val (mod N)
             * We need Q_val (possibly * k) to be smooth over the FB */

            mpz_abs(residue, Q_val);

            /* Also need to handle the factor k:
             * p^2 = kN * q^2 + (-1)^step * Q_val
             * => p^2 ≡ (-1)^step * Q_val (mod N)  [since kN ≡ 0 mod N]
             * But we also need: p^2 - (-1)^step * Q_val = kN * q^2
             * So the smooth value is Q_val itself
             * And the congruence is p^2 ≡ (-1)^step * Q_val (mod N)
             * For the GF(2) matrix, we factor Q_val over FB
             * and track the (-1)^step sign bit */

            /* Trial divide Q_val by FB primes */
            while (mpz_even_p(residue))
                mpz_tdiv_q_2exp(residue, residue, 1);
            for (int i = 1; i < fb->size; i++) {
                uint32_t p = fb->prime[i];
                while (mpz_divisible_ui_p(residue, p))
                    mpz_divexact_ui(residue, residue, p);
            }

            /* Divide out multiplier k */
            int k = cfs[mi].k;
            while (k > 1) {
                for (int i = 0; i < fb->size; i++) {
                    uint32_t p = fb->prime[i];
                    if (k % p == 0) {
                        k /= p;
                        while (mpz_divisible_ui_p(residue, p))
                            mpz_divexact_ui(residue, residue, p);
                        break;
                    }
                }
                if (k == cfs[mi].k) break; /* couldn't reduce further */
            }

            if (mpz_cmp_ui(residue, 1) == 0) {
                /* Full relation */
                int ri = full->count;
                if (ri < full->alloc) {
                    mpz_set(full->value[ri], p_val);
                    /* Store (-1)^step * Q_val for the matrix */
                    if (cfs[mi].step & 1)
                        mpz_neg(full->cofact[ri], Q_val);
                    else
                        mpz_set(full->cofact[ri], Q_val);
                    full->lp[ri] = 0;
                    full->count++;
                }
            } else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) <= lp_bound) {
                uint64_t lp_val = mpz_get_ui(residue);
                int match = lp_find(lpt, lp_val);
                if (match >= 0) {
                    int ri = full->count;
                    if (ri < full->alloc) {
                        mpz_mul(full->value[ri], p_val, part->value[match]);
                        mpz_mod(full->value[ri], full->value[ri], N);
                        mpz_t Q_signed; mpz_init(Q_signed);
                        if (cfs[mi].step & 1)
                            mpz_neg(Q_signed, Q_val);
                        else
                            mpz_set(Q_signed, Q_val);
                        mpz_mul(full->cofact[ri], Q_signed, part->cofact[match]);
                        full->lp[ri] = lp_val;
                        full->count++;
                        combined++;
                        mpz_clear(Q_signed);
                    }
                } else {
                    int pi = part->count;
                    if (pi < part->alloc) {
                        mpz_set(part->value[pi], p_val);
                        if (cfs[mi].step & 1)
                            mpz_neg(part->cofact[pi], Q_val);
                        else
                            mpz_set(part->cofact[pi], Q_val);
                        part->lp[pi] = lp_val;
                        lp_insert(lpt, lp_val, pi);
                        part->count++;
                    }
                }
            }
        }

        if (total_steps % 100000 == 0) {
            fprintf(stderr, "  steps=%d rels=%d/%d (full=%d+%d) part=%d t=%.1fs\n",
                    total_steps, full->count, target,
                    full->count - combined, combined, part->count, elapsed());
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "CF done: %d rels (%d full + %d combined) from %d steps in %.2fs\n",
            full->count, full->count - combined, combined, total_steps, sieve_time);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "MCFRAC: FAIL (not enough relations)\n");
        printf("FAIL\n");
        return 1;
    }

    /* === Linear Algebra === */
    int nrels = full->count;
    if (nrels > target) nrels = target;
    int ncols = fb->size + 1;
    gf2_t *mat = gf2_create(nrels, ncols);

    for (int r = 0; r < nrels; r++) {
        mpz_t Qval;
        mpz_init(Qval);
        mpz_set(Qval, full->cofact[r]);
        if (mpz_sgn(Qval) < 0) { gf2_set(mat, r, 0); mpz_neg(Qval, Qval); }
        int e2 = 0;
        while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set(mat, r, 1);
        for (int i = 1; i < fb->size; i++) {
            uint32_t p = fb->prime[i]; int e = 0;
            while (mpz_divisible_ui_p(Qval, p)) { mpz_divexact_ui(Qval, Qval, p); e++; }
            if (e & 1) gf2_set(mat, r, i + 1);
        }
        mpz_clear(Qval);
    }

    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, 64);
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n",
            ndeps, nrels, ncols, elapsed() - sieve_time);

    /* === Square Root === */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);

        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_mul(X, X, full->value[deps[d][k]]);
            mpz_mod(X, X, N);
        }

        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq; mpz_init(aq);
            mpz_abs(aq, full->cofact[deps[d][k]]);
            mpz_mul(prod, prod, aq);
            mpz_clear(aq);
        }

        mpz_set(rem, prod);
        int e2 = 0;
        while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next_dep;

        mpz_set_ui(Y, 1);
        if (e2 / 2 > 0) {
            mpz_set_ui(tmp, 2);
            mpz_powm_ui(tmp, tmp, e2 / 2, N);
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, N);
        }

        {
            int valid = 1;
            for (int i = 1; i < fb->size && valid; i++) {
                uint32_t p = fb->prime[i]; int e = 0;
                while (mpz_divisible_ui_p(rem, p)) { mpz_divexact_ui(rem, rem, p); e++; }
                if (e & 1) { valid = 0; break; }
                if (e / 2 > 0) {
                    mpz_set_ui(tmp, p);
                    mpz_powm_ui(tmp, tmp, e / 2, N);
                    mpz_mul(Y, Y, tmp);
                    mpz_mod(Y, Y, N);
                }
            }
            if (!valid) goto next_dep;
        }

        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_perfect_square_p(rem)) {
                mpz_sqrt(tmp, rem); mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            } else goto next_dep;
        }

        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "MCFRAC: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "MCFRAC: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }

        next_dep:
        mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "MCFRAC: FAILED\n");
    printf("FAIL\n");
    return 1;
}
