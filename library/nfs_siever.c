/*
 * nfs_siever.c - Custom NFS lattice siever for 85-95 digit semiprimes
 *
 * Implements the lattice sieving step of the General Number Field Sieve.
 * Designed for degree-4 polynomials, optimized for single-core.
 *
 * Output format: GGNFS-compatible relations
 *   a,b:hex_alg_primes:hex_rat_primes
 *
 * Compile:
 *   gcc -O3 -march=native -mavx512bw -o nfs_siever library/nfs_siever.c -lgmp -lm
 *
 * Usage:
 *   ./nfs_siever -f <startq> -c <qrange> -o <outfile> -a <jobfile>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>

/* Maximum polynomial degree */
#define MAX_DEGREE 6

/* Sieve parameters */
#define SIEVE_I_BITS 11      /* I = 2^11 = 2048 */
#define SIEVE_I (1 << SIEVE_I_BITS)
#define SIEVE_J (SIEVE_I / 2)
#define SIEVE_SIZE (SIEVE_I * SIEVE_J)

/* Factor base limits */
#define MAX_FB_SIZE 200000   /* Max factor base entries per side */
#define MAX_BUCKET_ALLOC (1 << 24)

/* Large prime bound */
#define MAX_LPB_BITS 27

typedef struct {
    uint32_t p;       /* prime */
    uint32_t r;       /* root of polynomial mod p */
    uint8_t logp;     /* floor(log2(p)) */
} fb_entry_t;

typedef struct {
    mpz_t n;          /* number to factor */
    int degree;       /* polynomial degree (4 for 90d) */
    mpz_t c[MAX_DEGREE+1]; /* algebraic poly coeffs c[0]..c[degree] */
    mpz_t Y0, Y1;     /* rational poly: Y1*x - Y0 (so m = Y0/Y1) */
    double skew;       /* polynomial skewness */

    uint32_t rlim;     /* rational factor base bound */
    uint32_t alim;     /* algebraic factor base bound */
    uint32_t lpbr;     /* rational large prime bound (bits) */
    uint32_t lpba;     /* algebraic large prime bound (bits) */
    double rlambda;    /* rational sieve lambda */
    double alambda;    /* algebraic sieve lambda */
} nfs_poly_t;

typedef struct {
    fb_entry_t *entries;
    uint32_t count;
    uint32_t alloc;
} factor_base_t;

static nfs_poly_t poly;
static factor_base_t rat_fb;   /* rational factor base */
static factor_base_t alg_fb;   /* algebraic factor base */

/* Sieve arrays */
static uint8_t *rat_sieve;    /* rational side sieve */
static uint8_t *alg_sieve;    /* algebraic side sieve */

/* Relation output */
static FILE *outfp;
static uint64_t total_rels = 0;

/*-----------------------------------------------------------------
 * Polynomial evaluation: compute f(a,b) = sum c[i] * a^i * b^(d-i)
 * This is the homogeneous form.
 *-----------------------------------------------------------------*/
static void eval_alg_poly(mpz_t result, const nfs_poly_t *p, int64_t a, uint32_t b) {
    mpz_t tmp, apow, bpow;
    mpz_inits(tmp, apow, bpow, NULL);

    mpz_set_ui(result, 0);
    mpz_set_si(apow, 1);          /* a^0 */
    mpz_ui_pow_ui(bpow, b, p->degree); /* b^d */

    for (int i = 0; i <= p->degree; i++) {
        /* term = c[i] * a^i * b^(d-i) */
        mpz_mul(tmp, p->c[i], apow);
        mpz_mul(tmp, tmp, bpow);
        mpz_add(result, result, tmp);

        /* update powers */
        mpz_mul_si(apow, apow, a);
        if (i < p->degree && b > 0) {
            mpz_divexact_ui(bpow, bpow, b);
        }
    }

    mpz_clears(tmp, apow, bpow, NULL);
}

/*-----------------------------------------------------------------
 * Rational norm: |a - b*m| where m = Y0/Y1
 * Actually: |Y1*a - Y0*b| (to stay in integers)
 * But for sieving, we want |a + b*(Y0/Y1)| on the rational side.
 * Actually the rational polynomial is g(x) = Y1*x + Y0 (linear)
 * So g(a,b) = Y1*a + Y0*b (homogeneous form)
 *-----------------------------------------------------------------*/
static void eval_rat_poly(mpz_t result, const nfs_poly_t *p, int64_t a, uint32_t b) {
    mpz_t tmp;
    mpz_init(tmp);

    mpz_mul_si(result, p->Y1, a);
    mpz_mul_ui(tmp, p->Y0, b);
    mpz_add(result, result, tmp);
    mpz_abs(result, result);

    mpz_clear(tmp);
}

/*-----------------------------------------------------------------
 * Parse a GGNFS-format job file
 *-----------------------------------------------------------------*/
static int parse_job_file(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Cannot open job file: %s\n", filename);
        return -1;
    }

    mpz_init(poly.n);
    for (int i = 0; i <= MAX_DEGREE; i++)
        mpz_init_set_ui(poly.c[i], 0);
    mpz_inits(poly.Y0, poly.Y1, NULL);
    poly.degree = 0;
    poly.skew = 1.0;
    poly.rlim = 500000;
    poly.alim = 1000000;
    poly.lpbr = 25;
    poly.lpba = 26;
    poly.rlambda = 2.4;
    poly.alambda = 2.4;

    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;

        char key[64];
        char val[960];
        if (sscanf(line, "%63[^:]: %959s", key, val) != 2) continue;

        if (strcmp(key, "n") == 0) mpz_set_str(poly.n, val, 10);
        else if (strcmp(key, "skew") == 0) poly.skew = atof(val);
        else if (key[0] == 'c' && key[1] >= '0' && key[1] <= '6') {
            int idx = key[1] - '0';
            mpz_set_str(poly.c[idx], val, 10);
            if (idx > poly.degree) poly.degree = idx;
        }
        else if (strcmp(key, "Y0") == 0) mpz_set_str(poly.Y0, val, 10);
        else if (strcmp(key, "Y1") == 0) mpz_set_str(poly.Y1, val, 10);
        else if (strcmp(key, "rlim") == 0) poly.rlim = atoi(val);
        else if (strcmp(key, "alim") == 0) poly.alim = atoi(val);
        else if (strcmp(key, "lpbr") == 0) poly.lpbr = atoi(val);
        else if (strcmp(key, "lpba") == 0) poly.lpba = atoi(val);
        else if (strcmp(key, "rlambda") == 0) poly.rlambda = atof(val);
        else if (strcmp(key, "alambda") == 0) poly.alambda = atof(val);
    }

    fclose(f);

    gmp_printf("n = %Zd\n", poly.n);
    printf("degree = %d, skew = %.2f\n", poly.degree, poly.skew);
    printf("rlim = %u, alim = %u, lpbr = %u, lpba = %u\n",
           poly.rlim, poly.alim, poly.lpbr, poly.lpba);

    return 0;
}

/*-----------------------------------------------------------------
 * Build factor bases: find roots of polynomials mod each prime p
 *-----------------------------------------------------------------*/
static uint32_t mod_poly_eval(const nfs_poly_t *p, uint32_t x, uint32_t prime) {
    uint64_t result = 0;
    uint64_t xpow = 1;
    for (int i = 0; i <= p->degree; i++) {
        /* mpz_fdiv_ui returns non-negative remainder for any sign */
        uint64_t ci = mpz_fdiv_ui(p->c[i], prime);
        result = (result + ci * xpow) % prime;
        xpow = (xpow * x) % prime;
    }
    return (uint32_t)result;
}

static void fb_init(factor_base_t *fb, uint32_t initial_size) {
    fb->alloc = initial_size;
    fb->count = 0;
    fb->entries = malloc(fb->alloc * sizeof(fb_entry_t));
}

static void fb_add(factor_base_t *fb, uint32_t p, uint32_t r) {
    if (fb->count >= fb->alloc) {
        fb->alloc *= 2;
        fb->entries = realloc(fb->entries, fb->alloc * sizeof(fb_entry_t));
    }
    fb->entries[fb->count].p = p;
    fb->entries[fb->count].r = r;
    fb->entries[fb->count].logp = (uint8_t)(log2(p) + 0.5);
    fb->count++;
}

/* Simple primality test for small numbers */
static int is_prime(uint32_t n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (uint32_t i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i+2) == 0) return 0;
    }
    return 1;
}

/* Polynomial arithmetic mod p for root finding */
typedef struct {
    uint32_t *coeff;  /* coefficients */
    int deg;          /* degree (-1 for zero poly) */
} poly_modp_t;

static poly_modp_t poly_alloc(int max_deg) {
    poly_modp_t r;
    r.coeff = calloc(max_deg + 1, sizeof(uint32_t));
    r.deg = -1;
    return r;
}

static void poly_free(poly_modp_t *p) { free(p->coeff); p->coeff = NULL; p->deg = -1; }

static void poly_normalize(poly_modp_t *p) {
    while (p->deg >= 0 && p->coeff[p->deg] == 0) p->deg--;
}

/* Compute (a * b) mod f mod p using schoolbook multiplication */
static poly_modp_t poly_mulmod(const poly_modp_t *a, const poly_modp_t *b,
                                const poly_modp_t *f, uint32_t p) {
    if (a->deg < 0 || b->deg < 0) {
        poly_modp_t z = poly_alloc(0);
        return z;
    }
    int prod_deg = a->deg + b->deg;
    poly_modp_t r = poly_alloc(prod_deg);

    for (int i = 0; i <= a->deg; i++) {
        if (a->coeff[i] == 0) continue;
        for (int j = 0; j <= b->deg; j++) {
            uint64_t t = (uint64_t)a->coeff[i] * b->coeff[j] + r.coeff[i+j];
            r.coeff[i+j] = t % p;
        }
    }
    r.deg = prod_deg;
    poly_normalize(&r);

    /* Reduce mod f */
    if (f->deg >= 0) {
        uint64_t f_inv = 0;
        /* Compute inverse of leading coeff of f */
        {
            int64_t aa = f->coeff[f->deg], bb = p, old_s = 1, s = 0;
            while (bb) { int64_t q = aa/bb, tmp = aa - q*bb; aa = bb; bb = tmp; tmp = old_s - q*s; old_s = s; s = tmp; }
            f_inv = ((old_s % (int64_t)p) + p) % p;
        }
        while (r.deg >= f->deg) {
            uint64_t scale = (uint64_t)r.coeff[r.deg] * f_inv % p;
            int shift = r.deg - f->deg;
            for (int i = 0; i <= f->deg; i++) {
                uint64_t sub = (uint64_t)scale * f->coeff[i] % p;
                r.coeff[i + shift] = (r.coeff[i + shift] + p - sub) % p;
            }
            poly_normalize(&r);
        }
    }
    return r;
}

/* Compute x^n mod f mod p using repeated squaring */
static poly_modp_t poly_powmod(uint64_t n, const poly_modp_t *f, uint32_t p) {
    poly_modp_t result = poly_alloc(f->deg);
    result.coeff[0] = 1; result.deg = 0;  /* result = 1 */

    poly_modp_t base = poly_alloc(f->deg);
    base.coeff[1] = 1; base.deg = 1;  /* base = x */

    while (n > 0) {
        if (n & 1) {
            poly_modp_t tmp = poly_mulmod(&result, &base, f, p);
            poly_free(&result);
            result = tmp;
        }
        n >>= 1;
        if (n > 0) {
            poly_modp_t tmp = poly_mulmod(&base, &base, f, p);
            poly_free(&base);
            base = tmp;
        }
    }
    poly_free(&base);
    return result;
}

/* GCD of two polynomials mod p */
static poly_modp_t poly_gcd(poly_modp_t a, poly_modp_t b, uint32_t p) {
    while (b.deg >= 0) {
        /* a = a mod b */
        uint64_t b_inv;
        {
            int64_t aa = b.coeff[b.deg], bb = p, old_s = 1, s = 0;
            while (bb) { int64_t q = aa/bb, tmp = aa - q*bb; aa = bb; bb = tmp; tmp = old_s - q*s; old_s = s; s = tmp; }
            b_inv = ((old_s % (int64_t)p) + p) % p;
        }
        while (a.deg >= b.deg) {
            uint64_t scale = (uint64_t)a.coeff[a.deg] * b_inv % p;
            int shift = a.deg - b.deg;
            for (int i = 0; i <= b.deg; i++) {
                uint64_t sub = (uint64_t)scale * b.coeff[i] % p;
                a.coeff[i + shift] = (a.coeff[i + shift] + p - sub) % p;
            }
            poly_normalize(&a);
        }
        /* swap */
        poly_modp_t tmp = a; a = b; b = tmp;
    }
    return a;
}

/* Find all roots of algebraic polynomial mod p using GCD method.
 * For large p, compute gcd(x^p - x, f(x)) mod p.
 * For small p (< 1000), use brute force. */
static int find_alg_roots(uint32_t *roots, uint32_t p) {
    if (p < 500) {
        /* Brute force for small primes */
        int count = 0;
        for (uint32_t x = 0; x < p && count < MAX_DEGREE; x++) {
            if (mod_poly_eval(&poly, x, p) == 0) {
                roots[count++] = x;
            }
        }
        return count;
    }

    /* Build f(x) mod p - mpz_fdiv_ui handles negative coefficients correctly */
    poly_modp_t f = poly_alloc(poly.degree);
    for (int i = 0; i <= poly.degree; i++) {
        f.coeff[i] = mpz_fdiv_ui(poly.c[i], p);
    }
    f.deg = poly.degree;
    poly_normalize(&f);

    if (f.deg < 0) { poly_free(&f); return 0; }

    /* Compute x^p mod f(x) mod p */
    poly_modp_t xp = poly_powmod((uint64_t)p, &f, p);

    /* Compute g = gcd(x^p - x, f) mod p */
    /* xp - x */
    if (xp.deg < 1) {
        poly_modp_t tmp = poly_alloc(1);
        for (int i = 0; i <= xp.deg; i++) tmp.coeff[i] = xp.coeff[i];
        tmp.deg = xp.deg;
        poly_free(&xp);
        xp = tmp;
    }
    xp.coeff[1] = (xp.coeff[1] + p - 1) % p;
    poly_normalize(&xp);

    /* Copy f for GCD */
    poly_modp_t f_copy = poly_alloc(f.deg);
    for (int i = 0; i <= f.deg; i++) f_copy.coeff[i] = f.coeff[i];
    f_copy.deg = f.deg;

    poly_modp_t g = poly_gcd(xp, f_copy, p);

    int nroots = 0;
    if (g.deg <= 0) {
        /* No roots */
    } else if (g.deg == 1) {
        /* One root: g = a1*x + a0, root = -a0/a1 mod p */
        uint64_t a1_inv;
        {
            int64_t aa = g.coeff[1], bb = p, old_s = 1, s = 0;
            while (bb) { int64_t q = aa/bb, tmp = aa - q*bb; aa = bb; bb = tmp; tmp = old_s - q*s; old_s = s; s = tmp; }
            a1_inv = ((old_s % (int64_t)p) + p) % p;
        }
        roots[nroots++] = (uint32_t)((uint64_t)(p - g.coeff[0]) * a1_inv % p);
    } else {
        /* Multiple roots: use Cantor-Zassenhaus splitting */
        /* Split g into linear factors by random (x+a)^((p-1)/2) - 1 */
        poly_modp_t factors[MAX_DEGREE];
        int nfactors = 0;
        factors[nfactors] = poly_alloc(g.deg);
        for (int i = 0; i <= g.deg; i++) factors[nfactors].coeff[i] = g.coeff[i];
        factors[nfactors].deg = g.deg;
        nfactors = 1;

        uint32_t rng_state = p;  /* Simple RNG seeded by p */

        for (int attempts = 0; attempts < 30 && nfactors < g.deg; attempts++) {
            /* Pick random a */
            rng_state = rng_state * 1103515245 + 12345;
            uint32_t a = rng_state % p;

            /* Try to split each non-linear factor */
            int new_nfactors = nfactors;
            for (int fi = 0; fi < nfactors; fi++) {
                if (factors[fi].deg <= 1) continue;

                /* Compute (x+a)^((p-1)/2) mod factors[fi] mod p */
                poly_modp_t xa = poly_alloc(factors[fi].deg);
                xa.coeff[0] = a; xa.coeff[1] = 1; xa.deg = 1;

                /* Need to compute xa^((p-1)/2) mod factors[fi] */
                poly_modp_t result = poly_alloc(factors[fi].deg);
                result.coeff[0] = 1; result.deg = 0;
                uint64_t exp = ((uint64_t)p - 1) / 2;

                while (exp > 0) {
                    if (exp & 1) {
                        poly_modp_t tmp = poly_mulmod(&result, &xa, &factors[fi], p);
                        poly_free(&result);
                        result = tmp;
                    }
                    exp >>= 1;
                    if (exp > 0) {
                        poly_modp_t tmp = poly_mulmod(&xa, &xa, &factors[fi], p);
                        poly_free(&xa);
                        xa = tmp;
                    }
                }
                poly_free(&xa);

                /* result - 1 */
                result.coeff[0] = (result.coeff[0] + p - 1) % p;
                poly_normalize(&result);

                /* Compute gcd(result, factors[fi]) */
                poly_modp_t fi_copy = poly_alloc(factors[fi].deg);
                for (int k = 0; k <= factors[fi].deg; k++) fi_copy.coeff[k] = factors[fi].coeff[k];
                fi_copy.deg = factors[fi].deg;

                poly_modp_t h = poly_gcd(result, fi_copy, p);

                if (h.deg > 0 && h.deg < factors[fi].deg) {
                    /* Successful split! h is one factor, factors[fi]/h is the other */
                    /* Compute quotient factors[fi] / h */
                    poly_modp_t quot = poly_alloc(factors[fi].deg);
                    /* Simple polynomial division */
                    poly_modp_t rem = poly_alloc(factors[fi].deg);
                    for (int k = 0; k <= factors[fi].deg; k++) rem.coeff[k] = factors[fi].coeff[k];
                    rem.deg = factors[fi].deg;

                    uint64_t h_inv;
                    {
                        int64_t aa2 = h.coeff[h.deg], bb2 = p, old_s2 = 1, s2 = 0;
                        while (bb2) { int64_t q2 = aa2/bb2, tmp = aa2 - q2*bb2; aa2 = bb2; bb2 = tmp; tmp = old_s2 - q2*s2; old_s2 = s2; s2 = tmp; }
                        h_inv = ((old_s2 % (int64_t)p) + p) % p;
                    }
                    while (rem.deg >= h.deg) {
                        uint64_t scale = (uint64_t)rem.coeff[rem.deg] * h_inv % p;
                        int qdeg = rem.deg - h.deg;
                        quot.coeff[qdeg] = (uint32_t)scale;
                        if (qdeg > quot.deg) quot.deg = qdeg;
                        for (int k = 0; k <= h.deg; k++) {
                            uint64_t sub = (uint64_t)scale * h.coeff[k] % p;
                            rem.coeff[k + qdeg] = (rem.coeff[k + qdeg] + p - sub) % p;
                        }
                        poly_normalize(&rem);
                    }
                    poly_free(&rem);

                    /* Replace factors[fi] with h, add quot */
                    poly_free(&factors[fi]);
                    factors[fi] = h;
                    if (new_nfactors < MAX_DEGREE) {
                        factors[new_nfactors] = quot;
                        new_nfactors++;
                    } else {
                        poly_free(&quot);
                    }
                } else {
                    poly_free(&h);
                }
            }
            nfactors = new_nfactors;
        }

        /* Extract roots from linear factors */
        for (int fi = 0; fi < nfactors && nroots < MAX_DEGREE; fi++) {
            if (factors[fi].deg == 1) {
                uint64_t a1_inv;
                {
                    int64_t aa = factors[fi].coeff[1], bb = p, old_s = 1, s = 0;
                    while (bb) { int64_t q = aa/bb, tmp = aa - q*bb; aa = bb; bb = tmp; tmp = old_s - q*s; old_s = s; s = tmp; }
                    a1_inv = ((old_s % (int64_t)p) + p) % p;
                }
                roots[nroots++] = (uint32_t)((uint64_t)(p - factors[fi].coeff[0]) * a1_inv % p);
            }
            poly_free(&factors[fi]);
        }
    }

    poly_free(&f);
    poly_free(&g);
    return nroots;
}

static void build_factor_bases(void) {
    printf("Building factor bases...\n");

    fb_init(&rat_fb, MAX_FB_SIZE);
    fb_init(&alg_fb, MAX_FB_SIZE);

    /* Rational side: g(x) = Y1*x + Y0, root = -Y0/Y1 mod p */
    for (uint32_t p = 2; p <= poly.rlim; p++) {
        if (!is_prime(p)) continue;

        uint64_t y1_mod = mpz_fdiv_ui(poly.Y1, p);
        if (y1_mod == 0) {
            /* Projective root - handle specially */
            fb_add(&rat_fb, p, p); /* sentinel for projective */
            continue;
        }

        /* r = -Y0 * Y1^(-1) mod p */
        uint64_t y0_mod = mpz_fdiv_ui(poly.Y0, p);
        /* modular inverse of Y1 mod p using extended GCD */
        int64_t g, s, t;
        {
            int64_t a = y1_mod, b = p;
            int64_t old_s = 1, old_t = 0;
            s = 0; t = 1;
            while (b != 0) {
                int64_t q = a / b;
                int64_t tmp = a - q * b; a = b; b = tmp;
                tmp = old_s - q * s; old_s = s; s = tmp;
                tmp = old_t - q * t; old_t = t; t = tmp;
            }
            g = a; s = old_s; /* t = old_t; */
        }
        if (g != 1) continue; /* shouldn't happen for primes */
        int64_t y1_inv = ((s % (int64_t)p) + p) % p;

        uint64_t r = (uint64_t)(p - y0_mod) * y1_inv % p;
        fb_add(&rat_fb, p, (uint32_t)r);
    }

    /* Algebraic side: find roots of f(x) mod p */
    for (uint32_t p = 2; p <= poly.alim; p++) {
        if (!is_prime(p)) continue;

        uint32_t roots[MAX_DEGREE];
        int nroots = find_alg_roots(roots, p);
        for (int i = 0; i < nroots; i++) {
            fb_add(&alg_fb, p, roots[i]);
        }
    }

    printf("Rational FB: %u entries (up to %u)\n", rat_fb.count, poly.rlim);
    printf("Algebraic FB: %u entries (up to %u)\n", alg_fb.count, poly.alim);
}

/*-----------------------------------------------------------------
 * Lattice sieve for a given special-q
 *
 * For special-q prime q with root r (f(r) ≡ 0 mod q):
 * We sieve over lattice L = {(a,b) : a ≡ r*b mod q}
 *
 * Lattice basis: e0 = (q, 0), e1 = (r, 1)
 * After reduction: e0, e1 are short vectors.
 *
 * Map lattice to sieve coordinates: (a,b) = i*e0 + j*e1
 * Sieve over i in [-I/2, I/2), j in [0, J)
 *-----------------------------------------------------------------*/

typedef struct {
    int64_t x, y;
} vec2_t;

/* Simple 2D lattice reduction (partial LLL for 2 vectors) */
static void reduce_lattice(vec2_t *e0, vec2_t *e1) {
    /* Gauss/Lagrange reduction for 2D lattice */
    double n0 = (double)e0->x * e0->x + (double)e0->y * e0->y;
    double n1 = (double)e1->x * e1->x + (double)e1->y * e1->y;

    if (n0 < n1) {
        vec2_t tmp = *e0; *e0 = *e1; *e1 = tmp;
        double t = n0; n0 = n1; n1 = t;
    }

    for (int iter = 0; iter < 100; iter++) {
        double dot = (double)e0->x * e1->x + (double)e0->y * e1->y;
        int64_t q = (int64_t)round(dot / n1);
        if (q == 0) break;

        e0->x -= q * e1->x;
        e0->y -= q * e1->y;
        n0 = (double)e0->x * e0->x + (double)e0->y * e0->y;

        if (n0 >= n1) break;

        vec2_t tmp = *e0; *e0 = *e1; *e1 = tmp;
        double t = n0; n0 = n1; n1 = t;
    }

    /* e1 should be the shorter vector */
    if (n0 < n1) {
        vec2_t tmp = *e0; *e0 = *e1; *e1 = tmp;
    }
}

/* Convert sieve coordinates (i,j) to (a,b) using lattice basis */
static inline void ij_to_ab(int64_t *a, int64_t *b,
                             int32_t i, int32_t j,
                             const vec2_t *e0, const vec2_t *e1) {
    *a = (int64_t)i * e0->x + (int64_t)j * e1->x;
    *b = (int64_t)i * e0->y + (int64_t)j * e1->y;
}

/*-----------------------------------------------------------------
 * Sieve one special-q
 *-----------------------------------------------------------------*/
static void sieve_special_q(uint32_t q, uint32_t qroot) {
    /* Set up lattice */
    vec2_t e0 = {q, 0};
    vec2_t e1 = {qroot, 1};
    reduce_lattice(&e0, &e1);

    /* Initialize sieve arrays */
    memset(alg_sieve, 0, SIEVE_SIZE);
    memset(rat_sieve, 0, SIEVE_SIZE);

    /* Compute target sieve thresholds */
    /* For algebraic side: log2(|F(a,b)|) ≈ half_degree * log2(I*skew) + ... */
    double log2_alg_norm = 0;
    {
        /* Rough estimate of algebraic norm in sieve region */
        double max_a = SIEVE_I * fmax(fabs(e0.x), fabs(e1.x));
        double max_b = SIEVE_J * fmax(fabs(e0.y), fabs(e1.y));
        if (max_b < 1) max_b = 1;
        log2_alg_norm = 0;
        for (int i = 0; i <= poly.degree; i++) {
            double ci = fabs(mpz_get_d(poly.c[i]));
            if (ci > 0) {
                double term = log2(ci) + i * log2(fmax(max_a, 1)) + (poly.degree - i) * log2(max_b);
                if (term > log2_alg_norm) log2_alg_norm = term;
            }
        }
    }
    double log2_rat_norm = 0;
    {
        double max_a = SIEVE_I * fmax(fabs(e0.x), fabs(e1.x));
        double max_b = SIEVE_J * fmax(fabs(e0.y), fabs(e1.y));
        double y1 = fabs(mpz_get_d(poly.Y1));
        double y0 = fabs(mpz_get_d(poly.Y0));
        log2_rat_norm = log2(y1 * max_a + y0 * max_b + 1);
    }

    uint8_t alg_thresh = (uint8_t)(log2_alg_norm * poly.alambda / (poly.alambda + 1));
    uint8_t rat_thresh = (uint8_t)(log2_rat_norm * poly.rlambda / (poly.rlambda + 1));
    /* Subtract special-q contribution */
    uint8_t logq = (uint8_t)(log2(q) + 0.5);
    if (alg_thresh > logq) alg_thresh -= logq;

    /* Sieve algebraic side by factor base primes */
    for (uint32_t fi = 0; fi < alg_fb.count; fi++) {
        uint32_t p = alg_fb.entries[fi].p;
        uint32_t r = alg_fb.entries[fi].r;
        uint8_t logp = alg_fb.entries[fi].logp;

        if (p == q) continue; /* Skip special-q */

        /* Find sieve hits: we need (a,b) such that a ≡ r*b mod p
         * In lattice coords (i,j):
         *   a = i*e0.x + j*e1.x
         *   b = i*e0.y + j*e1.y
         *   (i*e0.x + j*e1.x) ≡ r*(i*e0.y + j*e1.y) mod p
         *   i*(e0.x - r*e0.y) ≡ -j*(e1.x - r*e1.y) mod p
         */
        int64_t a0 = ((int64_t)e0.x - (int64_t)r * e0.y) % (int64_t)p;
        int64_t a1 = ((int64_t)e1.x - (int64_t)r * e1.y) % (int64_t)p;
        if (a0 < 0) a0 += p;
        if (a1 < 0) a1 += p;

        /* For each j, find starting i: i ≡ -j * a1 * a0^(-1) mod p */
        if (a0 == 0) {
            /* All i values at j positions where a1*j ≡ 0 mod p */
            if (a1 == 0) {
                /* All positions hit - degenerate, skip */
                continue;
            }
            /* Only j=0 hits (and j multiples of p) */
            for (int32_t j = 0; j < SIEVE_J; j += p) {
                for (int32_t i = -(SIEVE_I/2); i < SIEVE_I/2; i++) {
                    uint32_t idx = (j * SIEVE_I) + (i + SIEVE_I/2);
                    if (idx < SIEVE_SIZE)
                        alg_sieve[idx] += logp;
                }
            }
            continue;
        }

        /* Compute a0_inv = a0^(-1) mod p */
        uint64_t a0_inv;
        {
            int64_t aa = a0, bb = p;
            int64_t old_s = 1, s = 0;
            while (bb != 0) {
                int64_t qq = aa / bb;
                int64_t tmp = aa - qq * bb; aa = bb; bb = tmp;
                tmp = old_s - qq * s; old_s = s; s = tmp;
            }
            a0_inv = ((old_s % (int64_t)p) + p) % p;
        }

        /* step = a0_inv mod p (stride in i for consecutive j) */
        int64_t di = (int64_t)(a1 * a0_inv % p);

        for (int32_t j = 0; j < SIEVE_J; j++) {
            /* Starting i for this j:
             * i_start ≡ -j * a1 * a0_inv mod p */
            int64_t i_off = (-(int64_t)j * di % (int64_t)p);
            if (i_off < 0) i_off += p;

            /* Convert to sieve array index: sieve[i] where i in [0, SIEVE_I)
             * represents i_actual = i - SIEVE_I/2 in lattice coords.
             * We need i_actual ≡ i_off mod p, so i ≡ (i_off + SIEVE_I/2) mod p */
            uint32_t i_sieve = (uint32_t)((i_off + SIEVE_I/2) % p);

            for (uint32_t i = i_sieve; i < (uint32_t)SIEVE_I; i += p) {
                uint32_t idx = j * SIEVE_I + i;
                alg_sieve[idx] += logp;
            }
        }
    }

    /* Sieve rational side by factor base primes (similar logic) */
    for (uint32_t fi = 0; fi < rat_fb.count; fi++) {
        uint32_t p = rat_fb.entries[fi].p;
        uint32_t r = rat_fb.entries[fi].r;
        uint8_t logp = rat_fb.entries[fi].logp;

        if (p == 1) continue;

        /* Rational polynomial: g(a,b) = Y1*a + Y0*b = 0 mod p
         * So a ≡ -Y0/Y1 * b mod p, i.e., a ≡ r*b mod p (where r is the root)
         * Same lattice sieve logic as algebraic side */
        int64_t a0 = ((int64_t)e0.x - (int64_t)r * e0.y) % (int64_t)p;
        int64_t a1 = ((int64_t)e1.x - (int64_t)r * e1.y) % (int64_t)p;
        if (a0 < 0) a0 += p;
        if (a1 < 0) a1 += p;

        if (a0 == 0) {
            if (a1 == 0) continue;
            for (int32_t j = 0; j < SIEVE_J; j += p) {
                for (int32_t i = 0; i < SIEVE_I; i++) {
                    uint32_t idx = j * SIEVE_I + i;
                    rat_sieve[idx] += logp;
                }
            }
            continue;
        }

        uint64_t a0_inv;
        {
            int64_t aa = a0, bb = p;
            int64_t old_s = 1, s = 0;
            while (bb != 0) {
                int64_t qq = aa / bb;
                int64_t tmp = aa - qq * bb; aa = bb; bb = tmp;
                tmp = old_s - qq * s; old_s = s; s = tmp;
            }
            a0_inv = ((old_s % (int64_t)p) + p) % p;
        }

        int64_t di = (int64_t)(a1 * a0_inv % p);

        for (int32_t j = 0; j < SIEVE_J; j++) {
            int64_t i_off = (-(int64_t)j * di % (int64_t)p);
            if (i_off < 0) i_off += p;

            uint32_t i_sieve = (uint32_t)((i_off + SIEVE_I/2) % p);

            for (uint32_t i = i_sieve; i < (uint32_t)SIEVE_I; i += p) {
                uint32_t idx = j * SIEVE_I + i;
                rat_sieve[idx] += logp;
            }
        }
    }

    /* Scan sieve for candidates and trial divide */
    mpz_t anorm, rnorm, cofactor, tmp;
    mpz_inits(anorm, rnorm, cofactor, tmp, NULL);

    uint64_t candidates = 0;
    uint64_t rels_this_q = 0;

    for (int32_t j = 0; j < SIEVE_J; j++) {
        for (int32_t i = 0; i < SIEVE_I; i++) {
            uint32_t idx = j * SIEVE_I + i;

            /* Check if both sides exceed threshold */
            if (alg_sieve[idx] < alg_thresh || rat_sieve[idx] < rat_thresh)
                continue;

            /* Convert to (a,b) */
            int64_t a, b_val;
            ij_to_ab(&a, &b_val, i - SIEVE_I/2, j, &e0, &e1);

            if (b_val == 0) continue;
            /* Convention: b must be positive. Negate both if b < 0 */
            if (b_val < 0) { a = -a; b_val = -b_val; }
            if (a == 0) continue;

            candidates++;

            /* Compute algebraic norm: f(a,b) */
            eval_alg_poly(anorm, &poly, a, (uint32_t)b_val);
            mpz_abs(anorm, anorm);

            /* Compute rational norm: |Y1*a + Y0*b| */
            eval_rat_poly(rnorm, &poly, a, (uint32_t)b_val);

            /* Trial divide algebraic norm */
            mpz_set(cofactor, anorm);

            /* Remove special-q */
            while (mpz_divisible_ui_p(cofactor, q))
                mpz_divexact_ui(cofactor, cofactor, q);

            /* Buffer for prime factors */
            uint32_t alg_primes[256];
            int alg_nprimes = 0;
            alg_primes[alg_nprimes++] = q;

            /* Trial divide by algebraic FB */
            for (uint32_t fi = 0; fi < alg_fb.count && mpz_cmp_ui(cofactor, 1) > 0; fi++) {
                uint32_t p = alg_fb.entries[fi].p;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    if (alg_nprimes < 256)
                        alg_primes[alg_nprimes++] = p;
                }
            }

            /* Check algebraic cofactor: must be 1 or a large prime */
            int alg_ok = 0;
            if (mpz_cmp_ui(cofactor, 1) == 0) {
                alg_ok = 1;
            } else if (mpz_sizeinbase(cofactor, 2) <= poly.lpba) {
                /* Large prime on algebraic side */
                alg_ok = 1;
                if (alg_nprimes < 256) {
                    alg_primes[alg_nprimes++] = mpz_get_ui(cofactor);
                }
            }

            if (!alg_ok) continue;

            /* Trial divide rational norm */
            mpz_set(cofactor, rnorm);

            uint32_t rat_primes[256];
            int rat_nprimes = 0;

            for (uint32_t fi = 0; fi < rat_fb.count && mpz_cmp_ui(cofactor, 1) > 0; fi++) {
                uint32_t p = rat_fb.entries[fi].p;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    if (rat_nprimes < 256)
                        rat_primes[rat_nprimes++] = p;
                }
            }

            int rat_ok = 0;
            if (mpz_cmp_ui(cofactor, 1) == 0) {
                rat_ok = 1;
            } else if (mpz_sizeinbase(cofactor, 2) <= poly.lpbr) {
                rat_ok = 1;
                if (rat_nprimes < 256) {
                    rat_primes[rat_nprimes++] = mpz_get_ui(cofactor);
                }
            }

            if (!rat_ok) continue;

            /* Output relation in GGNFS format: a,b:alg_hex:rat_hex */
            fprintf(outfp, "%ld,%lu:", a, (unsigned long)b_val);
            for (int k = 0; k < alg_nprimes; k++) {
                if (k > 0) fprintf(outfp, ",");
                fprintf(outfp, "%x", alg_primes[k]);
            }
            fprintf(outfp, ":");
            for (int k = 0; k < rat_nprimes; k++) {
                if (k > 0) fprintf(outfp, ",");
                fprintf(outfp, "%x", rat_primes[k]);
            }
            fprintf(outfp, "\n");

            rels_this_q++;
            total_rels++;
        }
    }

    mpz_clears(anorm, rnorm, cofactor, tmp, NULL);

    if (rels_this_q > 0) {
        printf("q=%u: %lu candidates, %lu relations (total: %lu)\n",
               q, candidates, rels_this_q, total_rels);
    }
}

/*-----------------------------------------------------------------
 * Main
 *-----------------------------------------------------------------*/
int main(int argc, char *argv[]) {
    uint32_t start_q = 0, q_range = 0;
    char *job_file = NULL;
    char *out_file = NULL;

    /* Parse command line */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i+1 < argc) start_q = atoi(argv[++i]);
        else if (strcmp(argv[i], "-c") == 0 && i+1 < argc) q_range = atoi(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i+1 < argc) out_file = argv[++i];
        else if (strcmp(argv[i], "-a") == 0 && i+1 < argc) job_file = argv[++i];
    }

    if (!job_file || !start_q || !q_range) {
        fprintf(stderr, "Usage: %s -f <startq> -c <qrange> -o <outfile> -a <jobfile>\n", argv[0]);
        return 1;
    }

    if (!out_file) out_file = "rels.out";
    outfp = fopen(out_file, "w");
    if (!outfp) {
        fprintf(stderr, "Cannot open output file: %s\n", out_file);
        return 1;
    }

    /* Parse polynomial */
    if (parse_job_file(job_file) < 0) return 1;

    /* Build factor bases */
    build_factor_bases();

    /* Allocate sieve arrays */
    alg_sieve = calloc(SIEVE_SIZE, 1);
    rat_sieve = calloc(SIEVE_SIZE, 1);
    if (!alg_sieve || !rat_sieve) {
        fprintf(stderr, "Failed to allocate sieve arrays\n");
        return 1;
    }

    printf("Sieving q from %u to %u (I=%d, J=%d)\n",
           start_q, start_q + q_range, SIEVE_I, SIEVE_J);

    clock_t start_time = clock();

    /* Sieve over special-q range */
    for (uint32_t q = start_q; q < start_q + q_range; q++) {
        if (!is_prime(q)) continue;

        /* Find roots of algebraic polynomial mod q */
        uint32_t roots[MAX_DEGREE];
        int nroots = find_alg_roots(roots, q);

        for (int ri = 0; ri < nroots; ri++) {
            sieve_special_q(q, roots[ri]);
        }
    }

    double elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("\nTotal: %lu relations in %.1f seconds (%.1f rels/sec)\n",
           total_rels, elapsed, total_rels / (elapsed + 0.001));

    fclose(outfp);
    free(alg_sieve);
    free(rat_sieve);

    return 0;
}
