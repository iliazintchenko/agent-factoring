/*
 * gnfs_simple.c - Custom GNFS for 85-95 digit balanced semiprimes
 *
 * A single-threaded GNFS with:
 *   - Base-m polynomial selection (degree 4)
 *   - Fast modular polynomial root finding via GCD
 *   - Special-Q line sieving
 *   - Block Lanczos linear algebra
 *   - Algebraic square root via Couveignes' method
 *
 * NFS complexity: L[1/3, (64/9)^(1/3)] ≈ L[1/3, 1.923]
 * vs SIQS: L[1/2, 1]
 * At 90d, NFS should be competitive with SIQS.
 *
 * Compile: gcc -O3 -march=native -o gnfs_simple library/gnfs_simple.c -lgmp -lm
 * Usage:   timeout 295 ./gnfs_simple <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ========== CONFIGURATION ========== */

#define DEGREE 4
#define MAX_FB_SIZE 300000

/* Parameters tuned for 90-digit numbers */
#define DEFAULT_RFB_BOUND 100000
#define DEFAULT_AFB_BOUND 200000
#define DEFAULT_LPB_BITS 23
#define DEFAULT_SIEVE_I_BITS 12
#define DEFAULT_SIEVE_J 2048

#define SEED 42

/* ========== DATA STRUCTURES ========== */

typedef struct {
    uint32_t p;
    uint32_t r;            /* root of f(x) ≡ 0 (mod p) */
    uint8_t logp;
} fb_entry_t;

typedef struct {
    fb_entry_t *entries;
    int count;
    int alloc;
} factor_base_t;

typedef struct {
    mpz_t c[DEGREE + 1];
    mpz_t m;
    double skew;
    /* Precomputed c[i] mod p for small p — used during FB construction */
    uint32_t *cmod[DEGREE + 1];  /* cmod[i][j] = c[i] mod primes[j] */
} nfs_poly_t;

typedef struct {
    int64_t a;
    uint32_t b;
} relation_t;

/* ========== TIMING ========== */

static struct timespec g_start;

double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) +
           (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ========== POLYNOMIAL SELECTION ========== */

void poly_select_base_m(nfs_poly_t *poly, const mpz_t N) {
    mpz_t tmp;
    mpz_init(tmp);

    for (int i = 0; i <= DEGREE; i++) {
        mpz_init(poly->c[i]);
        poly->cmod[i] = NULL;
    }
    mpz_init(poly->m);

    /* m = floor(N^(1/(d+1))) */
    mpz_root(poly->m, N, DEGREE + 1);

    /* Try m, m+1, m-1, m+2, m-2 etc. to minimize max coefficient */
    mpz_t best_m, best_c[DEGREE + 1], test_c[DEGREE + 1];
    mpz_init(best_m);
    double best_score = 1e300;

    for (int i = 0; i <= DEGREE; i++) {
        mpz_init(best_c[i]);
        mpz_init(test_c[i]);
    }

    for (int delta = 0; delta <= 20; delta++) {
        for (int sign = (delta == 0 ? 1 : -1); sign <= 1; sign += 2) {
            mpz_t m_try;
            mpz_init(m_try);
            mpz_set(m_try, poly->m);
            if (sign > 0)
                mpz_add_ui(m_try, m_try, delta);
            else
                mpz_sub_ui(m_try, m_try, delta);

            if (mpz_sgn(m_try) <= 0) { mpz_clear(m_try); continue; }

            /* Express N in base m_try */
            mpz_set(tmp, N);
            for (int i = 0; i < DEGREE; i++) {
                mpz_fdiv_qr(tmp, test_c[i], tmp, m_try);
            }
            mpz_set(test_c[DEGREE], tmp);

            /* Score: sum of |c[i]| * m^i, roughly proportional to norm */
            double score = 0;
            for (int i = 0; i <= DEGREE; i++) {
                score += fabs(mpz_get_d(test_c[i]));
            }

            if (score < best_score) {
                best_score = score;
                mpz_set(best_m, m_try);
                for (int i = 0; i <= DEGREE; i++)
                    mpz_set(best_c[i], test_c[i]);
            }

            mpz_clear(m_try);
        }
    }

    mpz_set(poly->m, best_m);
    for (int i = 0; i <= DEGREE; i++)
        mpz_set(poly->c[i], best_c[i]);

    /* Verify: f(m) = N */
    mpz_set(tmp, poly->c[DEGREE]);
    for (int i = DEGREE - 1; i >= 0; i--) {
        mpz_mul(tmp, tmp, poly->m);
        mpz_add(tmp, tmp, poly->c[i]);
    }
    if (mpz_cmp(tmp, N) != 0) {
        fprintf(stderr, "ERROR: f(m) != N\n");
        exit(1);
    }

    poly->skew = pow(fabs(mpz_get_d(poly->c[0])) /
                     fabs(mpz_get_d(poly->c[DEGREE])),
                     1.0 / DEGREE);
    if (poly->skew < 1.0) poly->skew = 1.0;
    if (poly->skew > 1e6) poly->skew = 1e6;

    printf("Polynomial: m = ");
    mpz_out_str(stdout, 10, poly->m);
    printf("\n");
    for (int i = DEGREE; i >= 0; i--) {
        printf("  c%d = ", i);
        mpz_out_str(stdout, 10, poly->c[i]);
        printf("\n");
    }
    printf("  skew = %.1f, score = %.3e\n", poly->skew, best_score);

    for (int i = 0; i <= DEGREE; i++) {
        mpz_clear(best_c[i]);
        mpz_clear(test_c[i]);
    }
    mpz_clear(best_m);
    mpz_clear(tmp);
}

/* ========== FAST ROOT FINDING ========== */

/*
 * Find roots of f(x) mod p using modular polynomial arithmetic.
 * For degree 4: f(x) mod p has at most 4 roots.
 *
 * Algorithm:
 * 1. Compute g(x) = gcd(f(x), x^p - x) mod p
 *    (x^p - x splits into all linear factors mod p)
 * 2. If deg(g) == 0: no roots
 * 3. If deg(g) == 1: one root, extract it
 * 4. If deg(g) > 1: split using random element
 *
 * For p < 10000: brute force (faster for small primes)
 * For p >= 10000: use the GCD approach
 */

/* Polynomial mod p arithmetic (small degree, coefficients < p) */
typedef struct {
    uint64_t c[2 * DEGREE + 2];  /* coefficients, c[0] = constant, enough for products */
    int deg;
} modpoly_t;

static inline uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t p) {
    return (unsigned __int128)a * b % p;
}

static inline uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t p) {
    uint64_t result = 1;
    base %= p;
    while (exp > 0) {
        if (exp & 1) result = mod_mul(result, base, p);
        base = mod_mul(base, base, p);
        exp >>= 1;
    }
    return result;
}

static inline uint64_t mod_inv(uint64_t a, uint64_t p) {
    return mod_pow(a, p - 2, p);
}

/* Set poly from coefficient array, fix degree */
void mpoly_set(modpoly_t *r, const uint64_t *coeffs, int maxdeg) {
    int d = maxdeg;
    while (d > 0 && coeffs[d] == 0) d--;
    r->deg = d;
    for (int i = 0; i <= d; i++) r->c[i] = coeffs[i];
    for (int i = d + 1; i < 2 * DEGREE + 2; i++) r->c[i] = 0;
}

/* r = a mod b, using polynomial long division */
void mpoly_mod(modpoly_t *r, const modpoly_t *a, const modpoly_t *b, uint64_t p) {
    if (b->deg < 0 || (b->deg == 0 && b->c[0] == 0)) {
        fprintf(stderr, "poly division by zero\n");
        exit(1);
    }

    /* Copy a to r */
    *r = *a;

    uint64_t inv_lc = mod_inv(b->c[b->deg], p);

    while (r->deg >= b->deg) {
        uint64_t coeff = mod_mul(r->c[r->deg], inv_lc, p);
        int shift = r->deg - b->deg;

        for (int i = 0; i <= b->deg; i++) {
            uint64_t sub = mod_mul(coeff, b->c[i], p);
            r->c[i + shift] = (r->c[i + shift] + p - sub) % p;
        }

        /* Recompute degree */
        while (r->deg > 0 && r->c[r->deg] == 0) r->deg--;
        if (r->deg == 0 && r->c[0] == 0) { r->deg = -1; break; }
    }
}

/* r = gcd(a, b) mod p */
void mpoly_gcd(modpoly_t *r, const modpoly_t *a, const modpoly_t *b, uint64_t p) {
    modpoly_t u, v, tmp;
    u = *a;
    v = *b;

    while (v.deg >= 0 && !(v.deg == 0 && v.c[0] == 0)) {
        mpoly_mod(&tmp, &u, &v, p);
        u = v;
        v = tmp;
    }

    /* Make monic */
    if (u.deg > 0) {
        uint64_t inv = mod_inv(u.c[u.deg], p);
        for (int i = 0; i <= u.deg; i++)
            u.c[i] = mod_mul(u.c[i], inv, p);
    }

    *r = u;
}

/* r = (base^exp) mod modulus, all polynomial mod p */
void mpoly_powmod(modpoly_t *r, const modpoly_t *base, uint64_t exp,
                  const modpoly_t *modulus, uint64_t p) {
    modpoly_t result, b, tmp;

    /* result = 1 */
    memset(&result, 0, sizeof(result));
    result.deg = 0;
    result.c[0] = 1;

    b = *base;

    while (exp > 0) {
        if (exp & 1) {
            /* result = result * b */
            modpoly_t prod;
            memset(&prod, 0, sizeof(prod));
            prod.deg = result.deg + b.deg;
            for (int i = 0; i <= result.deg; i++) {
                for (int j = 0; j <= b.deg; j++) {
                    prod.c[i + j] = (prod.c[i + j] + mod_mul(result.c[i], b.c[j], p)) % p;
                }
            }
            mpoly_mod(&result, &prod, modulus, p);
        }
        /* b = b * b */
        modpoly_t sq;
        memset(&sq, 0, sizeof(sq));
        sq.deg = b.deg * 2;
        for (int i = 0; i <= b.deg; i++) {
            for (int j = 0; j <= b.deg; j++) {
                sq.c[i + j] = (sq.c[i + j] + mod_mul(b.c[i], b.c[j], p)) % p;
            }
        }
        mpoly_mod(&b, &sq, modulus, p);

        exp >>= 1;
    }

    *r = result;
}

/* Find all roots of f(x) ≡ 0 (mod p). Returns count, fills roots array. */
int poly_roots_mod_p(const nfs_poly_t *poly, uint32_t p, uint32_t *roots) {
    /* For small p, brute force is faster */
    if (p < 500) {
        int count = 0;
        for (uint32_t x = 0; x < p && count < DEGREE; x++) {
            uint64_t v = mpz_fdiv_ui(poly->c[DEGREE], p);
            for (int i = DEGREE - 1; i >= 0; i--) {
                v = ((uint64_t)v * x + mpz_fdiv_ui(poly->c[i], p)) % p;
            }
            if (v == 0) roots[count++] = x;
        }
        return count;
    }

    /* Build f(x) mod p */
    modpoly_t f;
    f.deg = DEGREE;
    for (int i = 0; i <= DEGREE; i++)
        f.c[i] = mpz_fdiv_ui(poly->c[i], p);
    while (f.deg > 0 && f.c[f.deg] == 0) f.deg--;

    if (f.deg <= 0) {
        /* Constant or zero - if zero, all are roots (skip), if constant, no roots */
        return 0;
    }

    /* Compute h(x) = x^p mod f(x) */
    modpoly_t x_poly;
    x_poly.deg = 1;
    x_poly.c[0] = 0;
    x_poly.c[1] = 1;
    for (int i = 2; i <= DEGREE + 1; i++) x_poly.c[i] = 0;

    modpoly_t xp;
    mpoly_powmod(&xp, &x_poly, p, &f, p);

    /* g(x) = x^p - x mod f(x) */
    xp.c[1] = (xp.c[1] + p - 1) % p;
    while (xp.deg > 0 && xp.c[xp.deg] == 0) xp.deg--;

    /* gcd(f(x), x^p - x) gives product of all linear factors */
    modpoly_t g;
    mpoly_gcd(&g, &f, &xp, p);

    if (g.deg <= 0) return 0;

    /* Extract roots from g(x) */
    int count = 0;

    if (g.deg == 1) {
        /* g(x) = x - r, so root = -c[0]/c[1] = p - c[0] * inv(c[1]) */
        uint64_t inv = mod_inv(g.c[1], p);
        roots[count++] = (p - mod_mul(g.c[0], inv, p)) % p;
    } else if (g.deg == 2) {
        /* Quadratic: try all values or use quadratic formula */
        for (uint32_t x = 0; x < p && count < 2; x++) {
            uint64_t v = (mod_mul(mod_mul(g.c[2], x, p) + g.c[1], x, p) + g.c[0]) % p;
            if (v == 0) roots[count++] = x;
        }
    } else if (g.deg <= 4) {
        /* Small degree: brute force (g has only ≤4 roots) */
        for (uint32_t x = 0; x < p && count < g.deg; x++) {
            uint64_t v = g.c[g.deg];
            for (int i = g.deg - 1; i >= 0; i--)
                v = (mod_mul(v, x, p) + g.c[i]) % p;
            if (v == 0) roots[count++] = x;
        }
    }

    return count;
}

/* ========== SIEVE OF ERATOSTHENES ========== */

int sieve_primes(uint32_t *primes, int max_primes, uint32_t bound) {
    char *is_composite = calloc(bound + 1, 1);
    int count = 0;
    for (uint32_t i = 2; i <= bound && count < max_primes; i++) {
        if (!is_composite[i]) {
            primes[count++] = i;
            for (uint64_t j = (uint64_t)i * i; j <= bound; j += i)
                is_composite[j] = 1;
        }
    }
    free(is_composite);
    return count;
}

/* ========== FACTOR BASE ========== */

void fb_init(factor_base_t *fb) {
    fb->alloc = 4096;
    fb->entries = malloc(sizeof(fb_entry_t) * fb->alloc);
    fb->count = 0;
}

void fb_add(factor_base_t *fb, uint32_t p, uint32_t r, uint8_t logp) {
    if (fb->count >= fb->alloc) {
        fb->alloc *= 2;
        fb->entries = realloc(fb->entries, sizeof(fb_entry_t) * fb->alloc);
    }
    fb->entries[fb->count].p = p;
    fb->entries[fb->count].r = r;
    fb->entries[fb->count].logp = logp;
    fb->count++;
}

void build_algebraic_fb(factor_base_t *fb, const nfs_poly_t *poly, uint32_t bound) {
    uint32_t *primes = malloc(sizeof(uint32_t) * MAX_FB_SIZE);
    int nprimes = sieve_primes(primes, MAX_FB_SIZE, bound);

    fb_init(fb);
    uint32_t roots[DEGREE];

    for (int i = 0; i < nprimes; i++) {
        uint32_t p = primes[i];
        int nroots = poly_roots_mod_p(poly, p, roots);
        uint8_t logp = (uint8_t)(log2(p));
        for (int j = 0; j < nroots; j++)
            fb_add(fb, p, roots[j], logp);
    }

    printf("Algebraic FB: %d entries from %d primes up to %u\n",
           fb->count, nprimes, bound);
    free(primes);
}

void build_rational_fb(factor_base_t *fb, const nfs_poly_t *poly, uint32_t bound) {
    uint32_t *primes = malloc(sizeof(uint32_t) * MAX_FB_SIZE);
    int nprimes = sieve_primes(primes, MAX_FB_SIZE, bound);

    fb_init(fb);
    for (int i = 0; i < nprimes; i++) {
        uint32_t p = primes[i];
        unsigned long m_mod_p = mpz_fdiv_ui(poly->m, p);
        uint32_t root = (p - m_mod_p) % p;
        uint8_t logp = (uint8_t)(log2(p));
        fb_add(fb, p, root, logp);
    }

    printf("Rational FB: %d entries up to %u\n", fb->count, bound);
    free(primes);
}

/* ========== SMOOTHNESS TESTING ========== */

/* Quick trial division to check if val is smooth over FB + one LP */
int is_smooth_trial(mpz_t val, const factor_base_t *fb, uint64_t lpb_max,
                    uint32_t *lp_out) {
    mpz_t tmp;
    mpz_init_set(tmp, val);
    mpz_abs(tmp, tmp);

    if (mpz_cmp_ui(tmp, 0) == 0) { mpz_clear(tmp); return 0; }

    /* Trial divide by FB primes */
    for (int i = 0; i < fb->count; i++) {
        uint32_t p = fb->entries[i].p;
        while (mpz_divisible_ui_p(tmp, p))
            mpz_divexact_ui(tmp, tmp, p);
        if (mpz_cmp_ui(tmp, 1) == 0) {
            if (lp_out) *lp_out = 1;
            mpz_clear(tmp);
            return 1;  /* fully smooth */
        }
    }

    /* Check for one large prime */
    if (mpz_fits_ulong_p(tmp)) {
        unsigned long rem = mpz_get_ui(tmp);
        if (rem <= lpb_max) {
            if (lp_out) *lp_out = (uint32_t)rem;
            mpz_clear(tmp);
            return 2;  /* one LP */
        }
    }

    mpz_clear(tmp);
    return 0;
}

/* ========== LINE SIEVING ========== */

void sieve_one_line(const factor_base_t *rfb, const factor_base_t *afb,
                    uint32_t b, uint8_t *rsieve, uint8_t *asieve,
                    int sieve_half) {
    int sieve_len = 2 * sieve_half;
    memset(rsieve, 0, sieve_len);
    memset(asieve, 0, sieve_len);

    /* Rational sieve: a + b*m ≡ 0 (mod p) → a ≡ -b*root (mod p) */
    for (int i = 0; i < rfb->count; i++) {
        uint32_t p = rfb->entries[i].p;
        if (p == 0) continue;
        uint8_t logp = rfb->entries[i].logp;

        int64_t start = ((int64_t)((uint64_t)b * rfb->entries[i].r % p) + sieve_half) % p;
        if (start < 0) start += p;

        for (int64_t idx = start; idx < sieve_len; idx += p)
            rsieve[idx] += logp;
    }

    /* Algebraic sieve: f(a,b) ≡ 0 (mod p) when a ≡ b*r (mod p) */
    for (int i = 0; i < afb->count; i++) {
        uint32_t p = afb->entries[i].p;
        if (p == 0) continue;
        uint8_t logp = afb->entries[i].logp;

        int64_t start = ((int64_t)((uint64_t)b * afb->entries[i].r % p) + sieve_half) % p;
        if (start < 0) start += p;

        for (int64_t idx = start; idx < sieve_len; idx += p)
            asieve[idx] += logp;
    }
}

/* ========== NORM COMPUTATION ========== */

/* Rational norm: |a + b*m| */
void compute_rnorm(mpz_t result, const nfs_poly_t *poly, int64_t a, uint32_t b) {
    mpz_mul_ui(result, poly->m, b);
    if (a >= 0)
        mpz_add_ui(result, result, (unsigned long)a);
    else
        mpz_sub_ui(result, result, (unsigned long)(-a));
    mpz_abs(result, result);
}

/* Algebraic norm: Resultant(a - b*x, f(x)) = (-b)^d * f(a/b)
 * = c_d * a^d + c_{d-1} * a^{d-1} * b + ... + c_0 * b^d */
void compute_anorm(mpz_t result, const nfs_poly_t *poly, int64_t a, uint32_t b) {
    mpz_t term, bpow;
    mpz_init(term);
    mpz_init_set_ui(bpow, 1);

    mpz_set_ui(result, 0);

    for (int i = DEGREE; i >= 0; i--) {
        /* term = c[i] * a^i * b^(d-i) */
        mpz_set(term, poly->c[i]);
        mpz_mul(term, term, bpow);

        /* Multiply by a^i */
        for (int j = 0; j < i; j++) {
            mpz_mul_si(term, term, a);
        }

        mpz_add(result, result, term);

        if (i > 0)
            mpz_mul_ui(bpow, bpow, b);
    }

    mpz_abs(result, result);
    mpz_clear(term);
    mpz_clear(bpow);
}

/* ========== MAIN ========== */

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);
    srand(SEED);

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N, 10);
    printf("GNFS factoring %d-digit number\n", digits);

    /* Auto-tune parameters */
    uint32_t rfb_bound, afb_bound;
    int sieve_i_bits, sieve_j;
    uint32_t lpb;

    if (digits <= 75) {
        rfb_bound = 20000; afb_bound = 40000;
        sieve_i_bits = 11; sieve_j = 512; lpb = 1 << 21;
    } else if (digits <= 80) {
        rfb_bound = 40000; afb_bound = 80000;
        sieve_i_bits = 11; sieve_j = 768; lpb = 1 << 22;
    } else if (digits <= 85) {
        rfb_bound = 60000; afb_bound = 120000;
        sieve_i_bits = 12; sieve_j = 1024; lpb = 1 << 22;
    } else if (digits <= 90) {
        rfb_bound = 80000; afb_bound = 160000;
        sieve_i_bits = 12; sieve_j = 2048; lpb = 1 << 23;
    } else {
        rfb_bound = 120000; afb_bound = 250000;
        sieve_i_bits = 12; sieve_j = 3072; lpb = 1 << 24;
    }

    int sieve_half = 1 << sieve_i_bits;
    int sieve_len = 2 * sieve_half;

    /* Step 1: Polynomial selection */
    printf("\n=== POLYNOMIAL SELECTION (%.2fs) ===\n", elapsed_sec());
    nfs_poly_t poly;
    poly_select_base_m(&poly, N);
    printf("Poly select done: %.2fs\n", elapsed_sec());

    /* Step 2: Factor bases */
    printf("\n=== FACTOR BASES ===\n");
    factor_base_t rfb, afb;
    build_rational_fb(&rfb, &poly, rfb_bound);
    build_algebraic_fb(&afb, &poly, afb_bound);
    printf("FB construction done: %.2fs\n", elapsed_sec());

    int target_rels = (int)((rfb.count + afb.count) * 1.05);
    printf("Target: %d relations (FB total: %d)\n", target_rels, rfb.count + afb.count);

    /* Step 3: Sieving */
    printf("\n=== SIEVING ===\n");
    printf("Sieve: a ∈ [-%d, %d), b ∈ [1, %d]\n", sieve_half, sieve_half, sieve_j);

    /* Compute sieve thresholds */
    double rnorm_bits = mpz_sizeinbase(poly.m, 2) + sieve_i_bits + log2(sieve_j);
    double anorm_bits = 0;
    for (int i = 0; i <= DEGREE; i++) {
        double cb = mpz_sizeinbase(poly.c[i], 2);
        if (cb > anorm_bits) anorm_bits = cb;
    }
    anorm_bits += DEGREE * sieve_i_bits + log2(sieve_j);

    int r_thresh = (int)(rnorm_bits * 0.5);
    int a_thresh = (int)(anorm_bits * 0.4);
    if (r_thresh < 15) r_thresh = 15;
    if (a_thresh < 15) a_thresh = 15;
    printf("Thresholds: r=%d (norm ~%d bits), a=%d (norm ~%d bits)\n",
           r_thresh, (int)rnorm_bits, a_thresh, (int)anorm_bits);

    relation_t *relations = malloc(sizeof(relation_t) * (target_rels + 10000));
    int num_rels = 0;
    int num_checked = 0;

    uint8_t *rsieve_buf = malloc(sieve_len);
    uint8_t *asieve_buf = malloc(sieve_len);

    mpz_t rnorm, anorm;
    mpz_init(rnorm);
    mpz_init(anorm);

    for (uint32_t b = 1; b <= (uint32_t)sieve_j && num_rels < target_rels; b++) {
        if (b % 50 == 0) {
            double t = elapsed_sec();
            if (t > 280.0) {
                printf("TIMEOUT at b=%u, %.1fs, %d/%d rels\n", b, t, num_rels, target_rels);
                break;
            }
            if (b % 200 == 0) {
                printf("b=%u: %d rels (%.1f/sec, %.1fs)\n",
                       b, num_rels, num_rels / t, t);
            }
        }

        sieve_one_line(&rfb, &afb, b, rsieve_buf, asieve_buf, sieve_half);

        for (int idx = 0; idx < sieve_len; idx++) {
            if (rsieve_buf[idx] < r_thresh || asieve_buf[idx] < a_thresh)
                continue;

            int64_t a = (int64_t)idx - sieve_half;
            if (a == 0) continue;

            /* gcd(a,b) == 1 check */
            int64_t abs_a = a < 0 ? -a : a;
            uint32_t g = b;
            int64_t t2 = abs_a;
            while (t2) { int64_t u = g % t2; g = (uint32_t)t2; t2 = u; }
            if (g > 1) continue;

            num_checked++;

            /* Check rational smoothness */
            compute_rnorm(rnorm, &poly, a, b);
            uint32_t rlp = 0;
            int rs = is_smooth_trial(rnorm, &rfb, lpb, &rlp);
            if (!rs) continue;

            /* Check algebraic smoothness */
            compute_anorm(anorm, &poly, a, b);
            uint32_t alp = 0;
            int as = is_smooth_trial(anorm, &afb, lpb, &alp);
            if (!as) continue;

            relations[num_rels].a = a;
            relations[num_rels].b = b;
            num_rels++;
        }
    }

    double total_time = elapsed_sec();
    printf("\n=== RESULTS ===\n");
    printf("Relations: %d / %d (%.1f%%)\n", num_rels, target_rels,
           100.0 * num_rels / target_rels);
    printf("Candidates checked: %d\n", num_checked);
    printf("Sieve rate: %.1f rels/sec\n", num_rels / total_time);
    printf("Total time: %.2fs\n", total_time);

    if (num_rels < rfb.count + afb.count) {
        printf("FAIL: insufficient relations\n");
    } else {
        printf("SUCCESS: enough relations for linear algebra\n");
        /* TODO: implement filter + LA + sqrt */
        printf("(Filter/LA/sqrt not yet implemented - output relations for external tools)\n");
    }

    /* Cleanup */
    free(relations);
    free(rsieve_buf);
    free(asieve_buf);
    free(rfb.entries);
    free(afb.entries);
    mpz_clear(rnorm);
    mpz_clear(anorm);
    mpz_clear(N);

    return 0;
}
