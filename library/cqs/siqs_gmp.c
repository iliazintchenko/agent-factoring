/*
 * Custom SIQS for balanced semiprimes, using GMP.
 * Build: gcc -O2 -march=native -o siqs_gmp siqs_gmp.c -lgmp -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define SIEVE_BLOCK  32768   /* 32KB sieve block */
#define MAX_FB       100000
#define MAX_RELS     300000
#define MAX_FACTORS  300
#define LP_HASH_SIZE (1 << 20)
#define SEED         42

/* ============== Parameters per bit size ============== */
typedef struct { int bits; uint32_t fb_target; uint32_t fb_bound; int nb; } param_row_t;
static const param_row_t PARAM_TABLE[] = {
    {130, 400,    5000,    2},
    {140, 600,    8000,    2},
    {150, 900,    13000,   3},
    {160, 1200,   18000,   3},
    {170, 1800,   30000,   4},
    {180, 2500,   45000,   5},
    {190, 3500,   65000,   6},
    {200, 5000,   90000,   7},
    {210, 6500,   120000,  8},
    {220, 8500,   170000,  10},
    {230, 12000,  260000,  12},
    {240, 16000,  360000,  14},
    {250, 22000,  500000,  16},
    {260, 30000,  700000,  18},
    {270, 40000,  900000,  20},
    {280, 55000,  1300000, 22},
    {290, 70000,  1600000, 24},
    {300, 90000,  2200000, 26},
};
#define NPARAM (sizeof(PARAM_TABLE)/sizeof(PARAM_TABLE[0]))

static void get_params(int bits, uint32_t *fb_target, uint32_t *fb_bound, int *nb) {
    if (bits <= PARAM_TABLE[0].bits) {
        *fb_target = PARAM_TABLE[0].fb_target;
        *fb_bound = PARAM_TABLE[0].fb_bound;
        *nb = PARAM_TABLE[0].nb;
        return;
    }
    for (int i = 0; i < (int)NPARAM - 1; i++) {
        if (bits >= PARAM_TABLE[i].bits && bits < PARAM_TABLE[i+1].bits) {
            double t = (double)(bits - PARAM_TABLE[i].bits) /
                       (PARAM_TABLE[i+1].bits - PARAM_TABLE[i].bits);
            *fb_target = PARAM_TABLE[i].fb_target + (uint32_t)(t * (PARAM_TABLE[i+1].fb_target - PARAM_TABLE[i].fb_target));
            *fb_bound = PARAM_TABLE[i].fb_bound + (uint32_t)(t * (PARAM_TABLE[i+1].fb_bound - PARAM_TABLE[i].fb_bound));
            *nb = PARAM_TABLE[i].nb + (int)(t * (PARAM_TABLE[i+1].nb - PARAM_TABLE[i].nb));
            if (*nb < 2) *nb = 2;
            return;
        }
    }
    *fb_target = PARAM_TABLE[NPARAM-1].fb_target;
    *fb_bound = PARAM_TABLE[NPARAM-1].fb_bound;
    *nb = PARAM_TABLE[NPARAM-1].nb;
}

/* ============== Utilities ============== */
static double elapsed_sec(struct timespec *s) {
    struct timespec n; clock_gettime(CLOCK_MONOTONIC, &n);
    return (n.tv_sec - s->tv_sec) + (n.tv_nsec - s->tv_nsec) / 1e9;
}

static uint32_t *sieve_primes(uint32_t limit, uint32_t *count) {
    uint8_t *siv = calloc(limit + 1, 1);
    uint32_t *p = malloc((limit / 2 + 100) * sizeof(uint32_t));
    uint32_t n = 0;
    for (uint32_t i = 2; i <= limit; i++) {
        if (!siv[i]) {
            p[n++] = i;
            if ((uint64_t)i * i <= limit)
                for (uint64_t j = (uint64_t)i * i; j <= limit; j += i) siv[j] = 1;
        }
    }
    free(siv); *count = n; return p;
}

static uint32_t tonelli_shanks(uint64_t n_mod_p, uint32_t p) {
    if (p == 2) return n_mod_p & 1;
    uint64_t a = n_mod_p % p;
    if (a == 0) return 0;
    if (p % 4 == 3) {
        /* r = a^((p+1)/4) mod p */
        uint64_t r = 1, b = a, e = (p + 1) / 4;
        for (; e; e >>= 1) { if (e & 1) r = r * b % p; b = b * b % p; }
        return (uint32_t)r;
    }
    uint32_t s = 0, q = p - 1;
    while (!(q & 1)) { s++; q >>= 1; }
    uint64_t z = 2;
    for (;;z++) {
        uint64_t r = 1, b = z, e = (p - 1) / 2;
        for (; e; e >>= 1) { if (e & 1) r = r * b % p; b = b * b % p; }
        if (r == p - 1) break;
    }
    uint64_t M = s, c, t, R;
    { uint64_t r=1,b=z,e=q; for(;e;e>>=1){if(e&1)r=r*b%p;b=b*b%p;} c=r; }
    { uint64_t r=1,b=a,e=q; for(;e;e>>=1){if(e&1)r=r*b%p;b=b*b%p;} t=r; }
    { uint64_t r=1,b=a,e=(q+1)/2; for(;e;e>>=1){if(e&1)r=r*b%p;b=b*b%p;} R=r; }
    for (;;) {
        if (t == 1) return (uint32_t)R;
        uint64_t i = 1, tt = t * t % p;
        while (tt != 1) { tt = tt * tt % p; i++; }
        uint64_t b = c;
        for (uint64_t j = 0; j < M - i - 1; j++) b = b * b % p;
        M = i; c = b * b % p; t = t * c % p; R = R * b % p;
    }
}

static uint32_t mod_inv(uint32_t a, uint32_t m) {
    int64_t g = m, x = 0, y = 1, a1 = a;
    while (a1) { int64_t q = g / a1, t; t = g - q*a1; g=a1; a1=t; t = x - q*y; x=y; y=t; }
    return (uint32_t)((x % (int64_t)m + m) % m);
}

/* ============== State ============== */
typedef struct {
    uint32_t p; uint32_t r; uint8_t logp;
} fb_t;

typedef struct {
    mpz_t y;              /* a*x + B */
    uint32_t *factors;    /* FB indices in exponent vector (including a's primes) */
    uint32_t nfactors;
    uint32_t large_prime; /* 0 = full, else LP value */
} rel_t;

typedef struct {
    mpz_t N, kN, sqrtKN;
    uint32_t k;
    int nbits;
    fb_t *fb; uint32_t fb_size;
    uint32_t fb_bound, lp_bound;
    uint8_t *sieve;
    int num_blocks, sieve_half; /* sieve covers [-sieve_half, sieve_half) */
    uint8_t sieve_init;

    /* Current polynomial */
    mpz_t poly_a, poly_b, poly_c;
    uint32_t *a_fb_idx; int a_nfact; /* FB indices of primes in 'a' */
    uint32_t *root1, *root2;

    /* Relations */
    rel_t *rels; uint32_t nrels;
    uint32_t n_full, n_partial, n_combined;

    /* LP hash */
    uint32_t *lp_first;  /* hash -> first chain index (1-based, 0=empty) */
    uint32_t *lp_next;   /* chain: next[i] -> next in chain */
    uint32_t *lp_val;    /* LP value */
    uint32_t *lp_rel;    /* relation index */
    uint32_t lp_count;

    gmp_randstate_t rng;
    struct timespec t0;
    uint64_t poly_count;
} state_t;

/* ============== Multiplier ============== */
static uint32_t choose_k(mpz_t N) {
    static const uint32_t ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,37,38,39,41,42,43};
    uint32_t np; uint32_t *prm = sieve_primes(300, &np);
    double best = -1e30; uint32_t bestk = 1;
    mpz_t kN; mpz_init(kN);
    for (int i = 0; i < (int)(sizeof(ks)/sizeof(ks[0])); i++) {
        mpz_mul_ui(kN, N, ks[i]);
        double sc = -0.5 * log((double)ks[i]);
        for (uint32_t j = 0; j < np; j++) {
            uint32_t p = prm[j];
            if (p == 2) {
                uint32_t m8 = mpz_fdiv_ui(kN, 8);
                if (m8==1) sc += 2*log(2.); else if (m8==5) sc += 1.5*log(2.); else if (m8&1) sc += log(2.);
            } else {
                int kr = mpz_kronecker_ui(kN, p);
                if (kr==1) { double c = 2*log((double)p)/(p-1); if(ks[i]%p==0) c*=0.5; sc+=c; }
                else if (kr==0) sc += log((double)p)/p;
            }
        }
        if (sc > best) { best = sc; bestk = ks[i]; }
    }
    mpz_clear(kN); free(prm);
    return bestk;
}

/* ============== Factor Base ============== */
static void build_fb(state_t *S) {
    uint32_t target, bound; int nb;
    get_params(S->nbits, &target, &bound, &nb);
    S->num_blocks = nb;
    S->sieve_half = nb * SIEVE_BLOCK;

    uint32_t np; uint32_t *prm = sieve_primes(bound + 1000, &np);
    S->fb = malloc(MAX_FB * sizeof(fb_t));
    S->fb_size = 0;
    /* fb[0] = sentinel for sign (-1) */
    S->fb[S->fb_size++] = (fb_t){1, 0, 0};
    for (uint32_t i = 0; i < np && S->fb_size < target; i++) {
        uint32_t p = prm[i];
        int kr = mpz_kronecker_ui(S->kN, p);
        if (kr != 1 && p != 2) continue;
        uint32_t r = tonelli_shanks(mpz_fdiv_ui(S->kN, p), p);
        uint8_t lp = (uint8_t)(log2((double)p) + 0.5);
        if (!lp) lp = 1;
        S->fb[S->fb_size++] = (fb_t){p, r, lp};
    }
    free(prm);
    S->fb_bound = S->fb[S->fb_size - 1].p;
    S->lp_bound = (uint64_t)S->fb_bound * 40;
    if (S->lp_bound > 0x7FFFFFFFU) S->lp_bound = 0x7FFFFFFFU;
    fprintf(stderr, "FB: %u primes, bound=%u, LP=%u, blocks=%d, half=%d\n",
            S->fb_size, S->fb_bound, S->lp_bound, S->num_blocks, S->sieve_half);
}

/* ============== Polynomial ============== */
static void new_poly(state_t *S) {
    mpz_t target_a, cur;
    mpz_init(target_a); mpz_init(cur);
    /* a ≈ sqrt(2*kN) / M */
    mpz_mul_ui(target_a, S->kN, 2);
    mpz_sqrt(target_a, target_a);
    mpz_fdiv_q_ui(target_a, target_a, S->sieve_half);

    double log_target = mpz_sizeinbase(target_a, 2) * log(2.0);
    /* Pick primes from middle range of FB */
    uint32_t lo = S->fb_size / 4, hi = S->fb_size * 3 / 4;
    if (lo < 10) lo = 10; if (hi >= S->fb_size) hi = S->fb_size - 1;

    double avg_lp = 0;
    for (uint32_t i = lo; i <= hi; i++) avg_lp += log((double)S->fb[i].p);
    avg_lp /= (hi - lo + 1);
    int nf = (int)(log_target / avg_lp + 0.5);
    if (nf < 2) nf = 2; if (nf > 15) nf = 15;

    S->a_nfact = nf;
    S->a_fb_idx = realloc(S->a_fb_idx, nf * sizeof(uint32_t));
    int *used = calloc(S->fb_size, sizeof(int));
    mpz_set_ui(cur, 1);
    double rem = log_target;

    for (int i = 0; i < nf; i++) {
        double want = rem / (nf - i);
        /* Collect candidates within factor of 2 of target */
        uint32_t cands[64]; int ncands = 0;
        for (uint32_t idx = lo; idx <= hi && ncands < 64; idx++) {
            if (used[idx] || !S->fb[idx].p) continue;
            double lp = log((double)S->fb[idx].p);
            if (fabs(lp - want) < want * 0.5)
                cands[ncands++] = idx;
        }
        if (ncands == 0) {
            /* Fallback: pick any unused prime */
            for (uint32_t idx = lo; idx <= hi; idx++) {
                if (!used[idx] && S->fb[idx].p) { cands[ncands++] = idx; break; }
            }
        }
        /* Random selection from candidates */
        uint32_t pick = cands[gmp_urandomm_ui(S->rng, ncands)];
        used[pick] = 1;
        S->a_fb_idx[i] = pick;
        mpz_mul_ui(cur, cur, S->fb[pick].p);
        rem -= log((double)S->fb[pick].p);
    }
    mpz_set(S->poly_a, cur);
    free(used);

    /* Compute B via CRT: B ≡ sqrt(kN) (mod a), B^2 ≡ kN (mod a) */
    mpz_set_ui(S->poly_b, 0);
    mpz_t adq, tmp;
    mpz_init(adq); mpz_init(tmp);
    for (int i = 0; i < nf; i++) {
        uint32_t qi = S->fb[S->a_fb_idx[i]].p;
        uint32_t ri = S->fb[S->a_fb_idx[i]].r;
        mpz_fdiv_q_ui(adq, S->poly_a, qi);
        uint32_t adq_mod = mpz_fdiv_ui(adq, qi);
        uint32_t inv = mod_inv(adq_mod, qi);
        uint64_t coeff = (uint64_t)ri * inv % qi;
        mpz_set_ui(tmp, coeff);
        mpz_mul(tmp, tmp, adq);
        mpz_add(S->poly_b, S->poly_b, tmp);
    }
    mpz_mod(S->poly_b, S->poly_b, S->poly_a);
    /* Ensure B is odd for uniqueness */
    if (mpz_even_p(S->poly_b)) mpz_sub(S->poly_b, S->poly_a, S->poly_b);

    /* C = (B^2 - kN) / a */
    mpz_mul(S->poly_c, S->poly_b, S->poly_b);
    mpz_sub(S->poly_c, S->poly_c, S->kN);
    mpz_divexact(S->poly_c, S->poly_c, S->poly_a);

    /* Compute sieve roots for each FB prime */
    for (uint32_t i = 1; i < S->fb_size; i++) {
        uint32_t p = S->fb[i].p;
        if (!p) { S->root1[i] = S->root2[i] = UINT32_MAX; continue; }
        /* Check if p divides a */
        int divs_a = 0;
        for (int j = 0; j < nf; j++) if (S->a_fb_idx[j] == i) { divs_a = 1; break; }
        if (divs_a) {
            /* For p | a: Q(x) = 2*B*x + C (mod p), single root x = -C/(2B) mod p */
            uint32_t bm = mpz_fdiv_ui(S->poly_b, p);
            uint32_t cm = mpz_fdiv_ui(S->poly_c, p);
            uint32_t inv_2b = mod_inv((2ULL * bm) % p, p);
            uint32_t r1 = (uint64_t)(p - cm) * inv_2b % p;
            r1 = (r1 + (uint32_t)(S->sieve_half % p)) % p;
            S->root1[i] = r1;
            S->root2[i] = UINT32_MAX; /* only one root */
            continue;
        }

        uint32_t r = S->fb[i].r;
        uint32_t bm = mpz_fdiv_ui(S->poly_b, p);
        uint32_t am = mpz_fdiv_ui(S->poly_a, p);
        uint32_t ainv = mod_inv(am, p);
        /* x = (±r - B) * a^{-1} mod p */
        uint32_t r1 = (uint64_t)((r + p - bm) % p) * ainv % p;
        uint32_t r2 = (uint64_t)((p - r + p - bm) % p) * ainv % p;
        /* Shift from x-space to sieve position space: pos = x + sieve_half */
        S->root1[i] = (r1 + (uint32_t)(S->sieve_half % p)) % p;
        S->root2[i] = (r2 + (uint32_t)(S->sieve_half % p)) % p;
    }

    /* Compute sieve init value = log2(Q_max) where Q_max ≈ sqrt(2*kN)*M */
    double l2q = mpz_sizeinbase(S->kN, 2) / 2.0 + log2((double)S->sieve_half) + 0.5;
    S->sieve_init = (uint8_t)(l2q < 240 ? l2q : 240);

    mpz_clear(adq); mpz_clear(tmp); mpz_clear(target_a); mpz_clear(cur);
    S->poly_count++;
}

/* ============== Sieve one block ============== */
static void do_sieve_block(state_t *S, uint32_t block_start) {
    memset(S->sieve, S->sieve_init, SIEVE_BLOCK);
    for (uint32_t i = 2; i < S->fb_size; i++) {
        if (S->root1[i] == UINT32_MAX && S->root2[i] == UINT32_MAX) continue;
        uint32_t p = S->fb[i].p;
        uint8_t lp = S->fb[i].logp;
        /* Find starting positions in this block */
        uint32_t s1 = UINT32_MAX, s2 = UINT32_MAX;
        if (S->root1[i] != UINT32_MAX) {
            if (S->root1[i] >= block_start) s1 = S->root1[i] - block_start;
            else { uint32_t k = (block_start - S->root1[i] + p - 1) / p; s1 = S->root1[i] + k * p - block_start; }
        }
        if (S->root2[i] != UINT32_MAX) {
            if (S->root2[i] >= block_start) s2 = S->root2[i] - block_start;
            else { uint32_t k = (block_start - S->root2[i] + p - 1) / p; s2 = S->root2[i] + k * p - block_start; }
        }
        if (p < SIEVE_BLOCK) {
            if (s1 != UINT32_MAX && s2 != UINT32_MAX) {
                while (s1 < SIEVE_BLOCK && s2 < SIEVE_BLOCK) {
                    S->sieve[s1] -= lp; S->sieve[s2] -= lp; s1 += p; s2 += p;
                }
            }
            if (s1 != UINT32_MAX) while (s1 < SIEVE_BLOCK) { S->sieve[s1] -= lp; s1 += p; }
            if (s2 != UINT32_MAX) while (s2 < SIEVE_BLOCK) { S->sieve[s2] -= lp; s2 += p; }
        } else {
            if (s1 != UINT32_MAX && s1 < SIEVE_BLOCK) S->sieve[s1] -= lp;
            if (s2 != UINT32_MAX && s2 < SIEVE_BLOCK) S->sieve[s2] -= lp;
        }
    }
}

/* ============== Trial division ============== */
/* Returns: 0=not smooth, 1=full, 2=partial(LP). Fills factors[] with FB indices. */
static int trial_div(state_t *S, int32_t x, uint32_t *factors, uint32_t *nf, uint32_t *lp_out) {
    mpz_t Q, tmp;
    mpz_init(Q); mpz_init(tmp);
    /* Q(x) = a*x^2 + 2*B*x + C */
    mpz_set_si(tmp, x);
    mpz_mul(Q, S->poly_a, tmp);   /* a*x */
    mpz_mul_ui(tmp, S->poly_b, 2);
    mpz_add(Q, Q, tmp);           /* a*x + 2B */
    mpz_set_si(tmp, x);
    mpz_mul(Q, Q, tmp);           /* (a*x + 2B)*x = a*x^2 + 2Bx */
    mpz_add(Q, Q, S->poly_c);    /* + C */

    *nf = 0;

    /* Include factors of 'a' in the exponent vector */
    for (int i = 0; i < S->a_nfact; i++)
        factors[(*nf)++] = S->a_fb_idx[i];

    /* Handle sign: if Q < 0, factor out -1 (index 0) */
    if (mpz_sgn(Q) < 0) { factors[(*nf)++] = 0; mpz_neg(Q, Q); }
    if (mpz_sgn(Q) == 0) { mpz_clear(Q); mpz_clear(tmp); return 0; }

    /* Trial divide by FB primes */
    for (uint32_t i = 1; i < S->fb_size && *nf < MAX_FACTORS - 2; i++) {
        uint32_t p = S->fb[i].p;
        if (!p) continue;
        while (mpz_divisible_ui_p(Q, p)) {
            mpz_divexact_ui(Q, Q, p);
            factors[(*nf)++] = i;
        }
        if (mpz_cmp_ui(Q, 1) == 0) { mpz_clear(Q); mpz_clear(tmp); return 1; }
    }

    /* Check large prime */
    if (mpz_fits_ulong_p(Q)) {
        uint64_t cof = mpz_get_ui(Q);
        if (cof > 1 && cof <= S->lp_bound) {
            *lp_out = (uint32_t)cof;
            mpz_clear(Q); mpz_clear(tmp);
            return 2;
        }
    }
    mpz_clear(Q); mpz_clear(tmp);
    return 0;
}

/* ============== Store relation ============== */
static void store_rel(state_t *S, int32_t x, uint32_t *factors, uint32_t nf,
                      int type, uint32_t lp) {
    if (S->nrels >= MAX_RELS) return;
    rel_t *r = &S->rels[S->nrels];
    mpz_init(r->y);
    /* y = a*x + B mod N */
    mpz_mul_si(r->y, S->poly_a, x);
    mpz_add(r->y, r->y, S->poly_b);
    mpz_mod(r->y, r->y, S->N);

    r->factors = malloc(nf * sizeof(uint32_t));
    memcpy(r->factors, factors, nf * sizeof(uint32_t));
    r->nfactors = nf;
    r->large_prime = (type == 2) ? lp : 0;

    if (type == 1) {
        S->n_full++;
    } else {
        S->n_partial++;
        /* Check for LP match */
        uint32_t h = lp % LP_HASH_SIZE;
        uint32_t idx = S->lp_first[h];
        while (idx) {
            if (S->lp_val[idx] == lp) { S->n_combined++; break; }
            idx = S->lp_next[idx];
        }
        if (!idx && S->lp_count + 1 < LP_HASH_SIZE) {
            S->lp_count++;
            S->lp_val[S->lp_count] = lp;
            S->lp_rel[S->lp_count] = S->nrels;
            S->lp_next[S->lp_count] = S->lp_first[h];
            S->lp_first[h] = S->lp_count;
        }
    }
    S->nrels++;
}

/* ============== Gaussian Elimination ============== */
static int do_linalg(state_t *S, mpz_t factor) {
    /* Build matrix using only FULL relations for simplicity */
    uint32_t ncols = S->fb_size;
    uint32_t max_rows = S->n_full + S->n_combined;
    if (max_rows < ncols + 1) {
        fprintf(stderr, "  LA: not enough rows (%u < %u)\n", max_rows, ncols+1);
        return 0;
    }

    /* Collect full relations and combined partials */
    uint32_t nrows = 0;
    uint32_t *rmap = malloc(max_rows * sizeof(uint32_t));
    uint32_t *rmap2 = malloc(max_rows * sizeof(uint32_t)); /* second rel for combined */

    for (uint32_t i = 0; i < S->nrels && nrows < max_rows; i++) {
        rel_t *r = &S->rels[i];
        if (r->large_prime == 0) {
            rmap[nrows] = i; rmap2[nrows] = UINT32_MAX; nrows++;
        }
    }

    /* Also include combined partials */
    /* For each LP that appears twice, combine them */
    for (uint32_t h = 0; h < LP_HASH_SIZE && nrows < max_rows; h++) {
        uint32_t idx = S->lp_first[h];
        while (idx) {
            uint32_t lp = S->lp_val[idx];
            uint32_t r1_idx = S->lp_rel[idx];
            /* Find second relation with same LP */
            uint32_t idx2 = S->lp_next[idx];
            while (idx2) {
                if (S->lp_val[idx2] == lp) {
                    rmap[nrows] = r1_idx;
                    rmap2[nrows] = S->lp_rel[idx2];
                    nrows++;
                    break;
                }
                idx2 = S->lp_next[idx2];
            }
            idx = S->lp_next[idx];
        }
    }

    if (nrows < ncols + 1) {
        free(rmap); free(rmap2);
        return 0;
    }
    if (nrows > ncols + 128) nrows = ncols + 128; /* limit matrix size */

    fprintf(stderr, "  LA: %u rows x %u cols\n", nrows, ncols);

    uint32_t nwords = (ncols + 63) / 64;
    uint64_t **mat = malloc(nrows * sizeof(uint64_t *));
    for (uint32_t i = 0; i < nrows; i++) mat[i] = calloc(nwords, sizeof(uint64_t));

    uint32_t hist_words = (nrows + 63) / 64;
    uint64_t **hist = malloc(nrows * sizeof(uint64_t *));
    for (uint32_t i = 0; i < nrows; i++) {
        hist[i] = calloc(hist_words, sizeof(uint64_t));
        hist[i][i/64] |= (1ULL << (i%64));
    }

    /* Fill matrix: exponent vector mod 2 */
    for (uint32_t row = 0; row < nrows; row++) {
        uint32_t *exps = calloc(ncols, sizeof(uint32_t));
        /* First relation */
        rel_t *r = &S->rels[rmap[row]];
        for (uint32_t j = 0; j < r->nfactors; j++)
            if (r->factors[j] < ncols) exps[r->factors[j]]++;
        /* Second relation (for combined partials) */
        if (rmap2[row] != UINT32_MAX) {
            rel_t *r2 = &S->rels[rmap2[row]];
            for (uint32_t j = 0; j < r2->nfactors; j++)
                if (r2->factors[j] < ncols) exps[r2->factors[j]]++;
            /* LP appears in both, adds 2 (even) - no contribution mod 2 */
        }
        for (uint32_t j = 0; j < ncols; j++)
            if (exps[j] & 1) mat[row][j/64] |= (1ULL << (j%64));
        free(exps);
    }

    /* Save original rmap before elimination (history references pre-swap indices) */
    uint32_t *orig_rmap = malloc(nrows * sizeof(uint32_t));
    uint32_t *orig_rmap2 = malloc(nrows * sizeof(uint32_t));
    memcpy(orig_rmap, rmap, nrows * sizeof(uint32_t));
    memcpy(orig_rmap2, rmap2, nrows * sizeof(uint32_t));

    /* Gaussian elimination */
    for (uint32_t col = 0; col < ncols; col++) {
        /* Find pivot row */
        uint32_t piv = UINT32_MAX;
        for (uint32_t r = col; r < nrows; r++) {  /* start from col for efficiency */
            if (mat[r][col/64] & (1ULL << (col%64))) { piv = r; break; }
        }
        if (piv == UINT32_MAX) continue;
        /* Swap to position col (if needed) */
        if (piv != col) {
            uint64_t *t; t = mat[col]; mat[col] = mat[piv]; mat[piv] = t;
            t = hist[col]; hist[col] = hist[piv]; hist[piv] = t;
            uint32_t u; u = rmap[col]; rmap[col] = rmap[piv]; rmap[piv] = u;
            u = rmap2[col]; rmap2[col] = rmap2[piv]; rmap2[piv] = u;
        }
        /* Eliminate */
        for (uint32_t r = 0; r < nrows; r++) {
            if (r == col) continue;
            if (mat[r][col/64] & (1ULL << (col%64))) {
                for (uint32_t w = 0; w < nwords; w++) mat[r][w] ^= mat[col][w];
                for (uint32_t w = 0; w < hist_words; w++) hist[r][w] ^= hist[col][w];
            }
        }
    }

    /* Mask off extra bits beyond ncols in the last word */
    uint64_t last_mask = (ncols % 64 == 0) ? ~0ULL : ((1ULL << (ncols % 64)) - 1);

    /* Find null space vectors */
    int found = 0;
    mpz_t X, Y, g;
    mpz_init(X); mpz_init(Y); mpz_init(g);

    int deps_found = 0, deps_tried = 0, deps_odd = 0;
    for (uint32_t row = 0; row < nrows && !found; row++) {
        int zero = 1;
        for (uint32_t w = 0; w < nwords - 1; w++) if (mat[row][w]) { zero=0; break; }
        if (zero && (mat[row][nwords-1] & last_mask)) zero = 0;
        if (!zero) continue;
        deps_found++;

        /* Combine relations in this dependency */
        mpz_set_ui(X, 1);
        uint32_t *total_exp = calloc(ncols, sizeof(uint32_t));

        for (uint32_t i = 0; i < nrows; i++) {
            if (!(hist[row][i/64] & (1ULL << (i%64)))) continue;
            /* Multiply y values */
            rel_t *r = &S->rels[orig_rmap[i]];
            mpz_mul(X, X, r->y);
            mpz_mod(X, X, S->N);
            for (uint32_t j = 0; j < r->nfactors; j++)
                if (r->factors[j] < ncols) total_exp[r->factors[j]]++;
            if (orig_rmap2[i] != UINT32_MAX) {
                rel_t *r2 = &S->rels[orig_rmap2[i]];
                mpz_mul(X, X, r2->y);
                mpz_mod(X, X, S->N);
                for (uint32_t j = 0; j < r2->nfactors; j++)
                    if (r2->factors[j] < ncols) total_exp[r2->factors[j]]++;
                /* LP^2 contribution: LP * LP = LP^2 */
                /* The LP itself is not in the FB. We need to include it in Y. */
                uint32_t lp = S->rels[rmap[i]].large_prime;
                if (lp > 1) {
                    /* lp appears with total exponent 2 (once in each relation) */
                    mpz_t lp_mpz; mpz_init_set_ui(lp_mpz, lp);
                    /* We'll multiply Y by lp after computing the FB part */
                    /* For now, just note it */
                    mpz_clear(lp_mpz);
                }
            }
        }

        deps_tried++;
        /* Y = product of p^(e/2) mod N */
        mpz_set_ui(Y, 1);
        for (uint32_t j = 0; j < ncols; j++) {
            if (total_exp[j] == 0) continue;
            /* Verify exponent is even */
            if (total_exp[j] & 1) { deps_odd++; goto next_dep; }  /* shouldn't happen */
            uint32_t e2 = total_exp[j] / 2;
            if (j == 0) continue; /* sign: -1^(2k) = 1 */
            if (j < S->fb_size) {
                mpz_t pw; mpz_init(pw);
                mpz_ui_pow_ui(pw, S->fb[j].p, e2);
                mpz_mul(Y, Y, pw);
                mpz_mod(Y, Y, S->N);
                mpz_clear(pw);
            }
        }

        /* Include LP contributions for combined partials */
        for (uint32_t i = 0; i < nrows; i++) {
            if (!(hist[row][i/64] & (1ULL << (i%64)))) continue;
            if (orig_rmap2[i] != UINT32_MAX) {
                uint32_t lp = S->rels[orig_rmap[i]].large_prime;
                if (lp > 1) {
                    /* LP appears twice -> LP^1 in Y */
                    mpz_mul_ui(Y, Y, lp);
                    mpz_mod(Y, Y, S->N);
                }
            }
        }

        /* Check gcd(X-Y, N) and gcd(X+Y, N) */
        mpz_sub(g, X, Y); mpz_gcd(g, g, S->N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, S->N) < 0) {
            mpz_set(factor, g); found = 1;
        }
        if (!found) {
            mpz_add(g, X, Y); mpz_gcd(g, g, S->N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, S->N) < 0) {
                mpz_set(factor, g); found = 1;
            }
        }
next_dep:
        free(total_exp);
    }

    fprintf(stderr, "  LA: deps_found=%d, deps_tried=%d, deps_odd=%d\n",
            deps_found, deps_tried, deps_odd);
    mpz_clear(X); mpz_clear(Y); mpz_clear(g);
    for (uint32_t i = 0; i < nrows; i++) { free(mat[i]); free(hist[i]); }
    free(mat); free(hist); free(rmap); free(rmap2);
    free(orig_rmap); free(orig_rmap2);
    return found;
}

/* ============== Main ============== */
static int factor_siqs(mpz_t N, mpz_t fac) {
    state_t S = {0};
    clock_gettime(CLOCK_MONOTONIC, &S.t0);
    mpz_init_set(S.N, N);
    mpz_init(S.kN); mpz_init(S.sqrtKN);
    mpz_init(S.poly_a); mpz_init(S.poly_b); mpz_init(S.poly_c);
    gmp_randinit_mt(S.rng); gmp_randseed_ui(S.rng, SEED);
    S.nbits = mpz_sizeinbase(N, 2);
    S.k = choose_k(N);
    mpz_mul_ui(S.kN, N, S.k);
    mpz_sqrt(S.sqrtKN, S.kN);
    fprintf(stderr, "SIQS: %d bits, k=%u\n", S.nbits, S.k);

    build_fb(&S);
    S.sieve = aligned_alloc(64, SIEVE_BLOCK + 64);
    S.root1 = calloc(S.fb_size, sizeof(uint32_t));
    S.root2 = calloc(S.fb_size, sizeof(uint32_t));
    S.rels = calloc(MAX_RELS, sizeof(rel_t));
    S.lp_first = calloc(LP_HASH_SIZE, sizeof(uint32_t));
    S.lp_next = calloc(LP_HASH_SIZE, sizeof(uint32_t));
    S.lp_val = calloc(LP_HASH_SIZE, sizeof(uint32_t));
    S.lp_rel = calloc(LP_HASH_SIZE, sizeof(uint32_t));
    S.a_fb_idx = NULL;

    uint32_t target = S.fb_size + 64;
    int found = 0;
    uint32_t factors[MAX_FACTORS]; uint32_t nf;

    while (!found) {
        new_poly(&S);
        uint8_t cutoff = (uint8_t)(log2((double)S.lp_bound) + 8);

        for (int blk = 0; blk < S.num_blocks * 2; blk++) {
            uint32_t bstart = blk * SIEVE_BLOCK;
            do_sieve_block(&S, bstart);
            for (uint32_t j = 0; j < SIEVE_BLOCK; j++) {
                if (S.sieve[j] < cutoff) {
                    int32_t x = (int32_t)(bstart + j) - S.sieve_half;
                    uint32_t lp_val = 0;
                    int res = trial_div(&S, x, factors, &nf, &lp_val);
                    if (res == 1) {
                        store_rel(&S, x, factors, nf, 1, 0);
                    } else if (res == 2) {
                        store_rel(&S, x, factors, nf, 2, lp_val);
                    }
                }
            }
        }

        uint32_t eff = S.n_full + S.n_combined;
        if (S.poly_count % 50 == 0 || eff >= target) {
            fprintf(stderr, "\r[%.1fs] poly=%lu full=%u partial=%u combined=%u eff=%u/%u (%.0f%%)   ",
                    elapsed_sec(&S.t0), S.poly_count, S.n_full, S.n_partial,
                    S.n_combined, eff, target, 100.0*eff/target);
        }
        if (eff >= target) {
            fprintf(stderr, "\nLinear algebra...\n");
            found = do_linalg(&S, fac);
            if (!found) { target += 32; fprintf(stderr, "LA failed, need more.\n"); }
        }
        if (elapsed_sec(&S.t0) > 290.0) break;
    }
    fprintf(stderr, "\nTime: %.2fs\n", elapsed_sec(&S.t0));

    /* Cleanup */
    for (uint32_t i = 0; i < S.nrels; i++) { mpz_clear(S.rels[i].y); free(S.rels[i].factors); }
    free(S.rels); free(S.sieve); free(S.root1); free(S.root2);
    free(S.fb); free(S.a_fb_idx);
    free(S.lp_first); free(S.lp_next); free(S.lp_val); free(S.lp_rel);
    mpz_clear(S.N); mpz_clear(S.kN); mpz_clear(S.sqrtKN);
    mpz_clear(S.poly_a); mpz_clear(S.poly_b); mpz_clear(S.poly_c);
    gmp_randclear(S.rng);
    return found;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    mpz_t N, fac, cof;
    mpz_init(N); mpz_init(fac); mpz_init(cof);
    mpz_set_str(N, argv[1], 10);

    /* Quick trial division */
    uint32_t np; uint32_t *sp = sieve_primes(1000000, &np);
    for (uint32_t i = 0; i < np; i++) {
        if (mpz_divisible_ui_p(N, sp[i])) {
            mpz_divexact_ui(cof, N, sp[i]);
            gmp_printf("%u\n%Zd\n", sp[i], cof);
            free(sp); return 0;
        }
    }
    free(sp);

    if (factor_siqs(N, fac)) {
        mpz_divexact(cof, N, fac);
        if (mpz_cmp(fac, cof) > 0) mpz_swap(fac, cof);
        gmp_printf("%Zd\n%Zd\n", fac, cof);
    } else {
        fprintf(stderr, "FAIL\n"); return 1;
    }
    mpz_clear(N); mpz_clear(fac); mpz_clear(cof);
    return 0;
}
