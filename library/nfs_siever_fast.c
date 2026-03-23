/*
 * nfs_siever_fast.c - High-performance NFS lattice siever for 90-95 digit semiprimes
 *
 * Key optimizations over nfs_siever.c:
 *   1. Bucket sieve for large primes (p > bucket_thresh)
 *   2. Small primes skipped in sieve, handled by trial division
 *   3. 32KB cache-line blocks for sieve array
 *   4. AVX512 for fast sieve-array scanning
 *   5. Sieve of Eratosthenes for factor base generation
 *   6. Precomputed per-j starting positions via recurrence
 *
 * Output format: GGNFS-compatible relations
 *   a,b:hex_alg_primes:hex_rat_primes
 *
 * Compile:
 *   gcc -O3 -march=native -mavx512bw -o nfs_siever_fast library/nfs_siever_fast.c -lgmp -lm
 *
 * Usage:
 *   ./nfs_siever_fast -f <startq> -c <qrange> -o <outfile> -a <jobfile>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <immintrin.h>
#include <gmp.h>

/* ============================================================
 * Configuration constants
 * ============================================================ */

#define MAX_DEGREE 6

/* Sieve geometry — I=4096, J=2048 for 90d numbers */
#define SIEVE_I_BITS 12
#define SIEVE_I (1 << SIEVE_I_BITS)
#define SIEVE_J (SIEVE_I / 2)
#define SIEVE_AREA ((uint32_t)SIEVE_I * SIEVE_J)

/* Cache-line block size for sieve: 32KB = 2^15 */
#define BLOCK_SIZE (1 << 15)
#define BLOCK_MASK (BLOCK_SIZE - 1)
#define NUM_BLOCKS ((SIEVE_AREA + BLOCK_SIZE - 1) / BLOCK_SIZE)

/* Factor base limits */
#define MAX_FB_SIZE 300000

/* Small prime threshold: primes below this are not sieved, handled
 * by trial division during cofactorization.  For I=4096 these primes
 * hit a large fraction of positions and are cheap to trial-divide. */
#define SMALL_PRIME_CUTOFF 64

/* Bucket sieve threshold: primes >= bucket_thresh are sieved via
 * bucket-sort.  A good heuristic is sqrt(BLOCK_SIZE) so that each
 * prime touches at most ~1 position per block.  For 32KB blocks,
 * sqrt(32768) ≈ 181.  We round up a bit. */
#define BUCKET_THRESH 256

/* Maximum prime to sieve: primes above this are not sieved but are
 * still used during trial division.  This dramatically reduces bucket
 * fill time.  Set to 0 to sieve all primes up to alim/rlim.
 * For 90d numbers with I=4096, a good value is ~200K-500K.
 * Each doubling of this value roughly doubles sieve time but only
 * adds ~1 bit to the sieve accuracy. */
/* Default: sieve primes up to 200K. Can be overridden with -S flag. */
static uint32_t max_sieve_prime = 200000;

/* Maximum bucket entries per block.  With ~100K FB primes and
 * NUM_BLOCKS blocks the average is manageable; we use a generous
 * allocation and fall back if full. */
#define MAX_BUCKET_ENTRIES (1 << 20)  /* 1M entries per block */

/* Large prime bound */
#define MAX_LPB_BITS 27

/* Maximum prime factors stored per relation */
#define MAX_REL_PRIMES 256

/* ============================================================
 * Data structures
 * ============================================================ */

typedef struct {
    uint32_t p;       /* prime */
    uint32_t r;       /* root of polynomial mod p */
    uint8_t  logp;    /* floor(log2(p)) */
} fb_entry_t;

typedef struct {
    fb_entry_t *entries;
    uint32_t count;
    uint32_t alloc;
} factor_base_t;

/* Bucket entry: packed as (position_in_block : 15 bits, logp : 8 bits)
 * but for simplicity we store (uint16_t pos, uint8_t logp) or pack into 32 bits. */
typedef struct {
    uint16_t pos;     /* position within 32KB block (0..32767) */
    uint8_t  logp;
} __attribute__((packed)) bucket_entry_t;

typedef struct {
    bucket_entry_t *entries;
    uint32_t count;
    uint32_t alloc;
} bucket_t;

typedef struct {
    mpz_t n;
    int degree;
    mpz_t c[MAX_DEGREE + 1];   /* algebraic coefficients */
    mpz_t Y0, Y1;              /* rational linear poly */
    double skew;

    uint32_t rlim, alim;
    uint32_t lpbr, lpba;
    uint32_t mfbr, mfba;       /* double-large-prime bounds (bits) */
    double rlambda, alambda;
} nfs_poly_t;

typedef struct {
    int64_t x, y;
} vec2_t;

/* ============================================================
 * Globals
 * ============================================================ */

static nfs_poly_t poly;
static factor_base_t rat_fb, alg_fb;

/* Sieve arrays (flat, SIEVE_AREA bytes each) */
static uint8_t *alg_sieve;
static uint8_t *rat_sieve;

/* Bucket arrays — one per block, per side */
static bucket_t *alg_buckets;   /* [NUM_BLOCKS] */
static bucket_t *rat_buckets;

/* Prime sieve for factor-base generation */
static uint8_t *prime_sieve_bits;
static uint32_t prime_sieve_limit;

static FILE *outfp;
static uint64_t total_rels = 0;

/* ============================================================
 * Prime sieve (Eratosthenes)
 * ============================================================ */

static void init_prime_sieve(uint32_t limit) {
    prime_sieve_limit = limit;
    uint32_t bytes = (limit / 2) + 1;
    prime_sieve_bits = calloc(bytes, 1);  /* 0 = prime, 1 = composite */
    for (uint32_t i = 3; (uint64_t)i * i <= limit; i += 2) {
        if (prime_sieve_bits[i / 2] == 0) {
            for (uint32_t j = i * i; j <= limit; j += 2 * i)
                prime_sieve_bits[j / 2] = 1;
        }
    }
}

static inline int is_prime(uint32_t n) {
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    if (n <= prime_sieve_limit)
        return prime_sieve_bits[n / 2] == 0;
    /* Fallback for numbers beyond sieve */
    for (uint32_t d = 3; (uint64_t)d * d <= n; d += 2)
        if (n % d == 0) return 0;
    return 1;
}

/* ============================================================
 * Modular arithmetic helpers
 * ============================================================ */

static inline uint32_t mod_inverse(uint32_t a, uint32_t p) {
    int64_t aa = a, bb = p, old_s = 1, s = 0;
    while (bb) {
        int64_t q = aa / bb;
        int64_t tmp;
        tmp = aa - q * bb; aa = bb; bb = tmp;
        tmp = old_s - q * s; old_s = s; s = tmp;
    }
    return (uint32_t)(((old_s % (int64_t)p) + p) % p);
}

/* ============================================================
 * Polynomial evaluation
 * ============================================================ */

static uint32_t mod_poly_eval(const nfs_poly_t *p, uint32_t x, uint32_t prime) {
    uint64_t result = 0, xpow = 1;
    for (int i = 0; i <= p->degree; i++) {
        uint64_t ci = mpz_fdiv_ui(p->c[i], prime);
        result = (result + ci * xpow) % prime;
        xpow = (xpow * x) % prime;
    }
    return (uint32_t)result;
}

/* Homogeneous algebraic: F(a,b) = sum c[i] * a^i * b^(d-i) */
static void eval_alg_poly(mpz_t result, const nfs_poly_t *p, int64_t a, int64_t b) {
    mpz_t tmp, apow, bpow;
    mpz_inits(tmp, apow, bpow, NULL);

    mpz_set_ui(result, 0);
    mpz_set_si(apow, 1);
    mpz_ui_pow_ui(bpow, (unsigned long)labs(b), p->degree);
    if (b < 0) {
        /* b^d with sign: if degree is odd and b<0, bpow is negative */
        if (p->degree % 2 == 1) mpz_neg(bpow, bpow);
    }

    for (int i = 0; i <= p->degree; i++) {
        mpz_mul(tmp, p->c[i], apow);
        mpz_mul(tmp, tmp, bpow);
        mpz_add(result, result, tmp);

        mpz_mul_si(apow, apow, a);
        if (i < p->degree && b != 0)
            mpz_divexact_ui(bpow, bpow, (unsigned long)labs(b));
    }
    mpz_clears(tmp, apow, bpow, NULL);
}

/* Rational norm: |Y1*a + Y0*b| */
static void eval_rat_poly(mpz_t result, const nfs_poly_t *p, int64_t a, int64_t b) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul_si(result, p->Y1, a);
    mpz_mul_si(tmp, p->Y0, b);
    mpz_add(result, result, tmp);
    mpz_abs(result, result);
    mpz_clear(tmp);
}

/* ============================================================
 * Job file parser
 * ============================================================ */

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
    poly.rlim = 1200000;
    poly.alim = 1200000;
    poly.lpbr = 25;
    poly.lpba = 25;
    poly.mfbr = 50;
    poly.mfba = 50;
    poly.rlambda = 2.5;
    poly.alambda = 2.5;

    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char key[64], val[960];
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
        else if (strcmp(key, "rlim") == 0) poly.rlim = (uint32_t)atol(val);
        else if (strcmp(key, "alim") == 0) poly.alim = (uint32_t)atol(val);
        else if (strcmp(key, "lpbr") == 0) poly.lpbr = atoi(val);
        else if (strcmp(key, "lpba") == 0) poly.lpba = atoi(val);
        else if (strcmp(key, "mfbr") == 0) poly.mfbr = atoi(val);
        else if (strcmp(key, "mfba") == 0) poly.mfba = atoi(val);
        else if (strcmp(key, "rlambda") == 0) poly.rlambda = atof(val);
        else if (strcmp(key, "alambda") == 0) poly.alambda = atof(val);
    }
    fclose(f);

    gmp_printf("n = %Zd\n", poly.n);
    printf("degree = %d, skew = %.2f\n", poly.degree, poly.skew);
    printf("rlim = %u, alim = %u, lpbr = %u, lpba = %u\n",
           poly.rlim, poly.alim, poly.lpbr, poly.lpba);
    printf("mfbr = %u, mfba = %u\n", poly.mfbr, poly.mfba);

    return 0;
}

/* ============================================================
 * Factor base construction
 * ============================================================ */

static void fb_init(factor_base_t *fb, uint32_t initial_size) {
    fb->alloc = initial_size;
    fb->count = 0;
    fb->entries = malloc(fb->alloc * sizeof(fb_entry_t));
}

static inline void fb_add(factor_base_t *fb, uint32_t p, uint32_t r) {
    if (fb->count >= fb->alloc) {
        fb->alloc = fb->alloc * 3 / 2;
        fb->entries = realloc(fb->entries, fb->alloc * sizeof(fb_entry_t));
    }
    fb_entry_t *e = &fb->entries[fb->count++];
    e->p = p;
    e->r = r;
    e->logp = (uint8_t)(log2((double)p) + 0.5);
}

/* ---- Polynomial root-finding mod p (Cantor-Zassenhaus) ---- */

typedef struct {
    uint32_t *coeff;
    int deg;
} poly_modp_t;

static poly_modp_t poly_alloc(int max_deg) {
    poly_modp_t r;
    r.coeff = calloc(max_deg + 2, sizeof(uint32_t));
    r.deg = -1;
    return r;
}

static void poly_free(poly_modp_t *p) {
    free(p->coeff);
    p->coeff = NULL;
    p->deg = -1;
}

static void poly_normalize(poly_modp_t *p) {
    while (p->deg >= 0 && p->coeff[p->deg] == 0) p->deg--;
}

static poly_modp_t poly_mulmod(const poly_modp_t *a, const poly_modp_t *b,
                                const poly_modp_t *f, uint32_t p) {
    if (a->deg < 0 || b->deg < 0) return poly_alloc(0);
    int prod_deg = a->deg + b->deg;
    poly_modp_t r = poly_alloc(prod_deg);
    for (int i = 0; i <= a->deg; i++) {
        if (a->coeff[i] == 0) continue;
        for (int j = 0; j <= b->deg; j++) {
            uint64_t t = (uint64_t)a->coeff[i] * b->coeff[j] + r.coeff[i + j];
            r.coeff[i + j] = t % p;
        }
    }
    r.deg = prod_deg;
    poly_normalize(&r);

    if (f->deg >= 0) {
        uint64_t f_inv = mod_inverse(f->coeff[f->deg], p);
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

static poly_modp_t poly_powmod(uint64_t n, const poly_modp_t *f, uint32_t p) {
    poly_modp_t result = poly_alloc(f->deg);
    result.coeff[0] = 1; result.deg = 0;
    poly_modp_t base = poly_alloc(f->deg);
    base.coeff[1] = 1; base.deg = 1;
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

static poly_modp_t poly_gcd(poly_modp_t a, poly_modp_t b, uint32_t p) {
    while (b.deg >= 0) {
        uint64_t b_inv = mod_inverse(b.coeff[b.deg], p);
        while (a.deg >= b.deg) {
            uint64_t scale = (uint64_t)a.coeff[a.deg] * b_inv % p;
            int shift = a.deg - b.deg;
            for (int i = 0; i <= b.deg; i++) {
                uint64_t sub = (uint64_t)scale * b.coeff[i] % p;
                a.coeff[i + shift] = (a.coeff[i + shift] + p - sub) % p;
            }
            poly_normalize(&a);
        }
        poly_modp_t tmp = a; a = b; b = tmp;
    }
    return a;
}

static int find_alg_roots(uint32_t *roots, uint32_t p) {
    if (p < 500) {
        int count = 0;
        for (uint32_t x = 0; x < p && count < MAX_DEGREE; x++) {
            if (mod_poly_eval(&poly, x, p) == 0)
                roots[count++] = x;
        }
        return count;
    }

    poly_modp_t f = poly_alloc(poly.degree);
    for (int i = 0; i <= poly.degree; i++)
        f.coeff[i] = mpz_fdiv_ui(poly.c[i], p);
    f.deg = poly.degree;
    poly_normalize(&f);
    if (f.deg < 0) { poly_free(&f); return 0; }

    poly_modp_t xp = poly_powmod((uint64_t)p, &f, p);
    /* Make room for degree 1 if xp is constant */
    if (xp.deg < 1) {
        poly_modp_t tmp = poly_alloc(1);
        for (int i = 0; i <= xp.deg && i <= 1; i++) tmp.coeff[i] = xp.coeff[i];
        tmp.deg = xp.deg;
        poly_free(&xp);
        xp = tmp;
    }
    xp.coeff[1] = (xp.coeff[1] + p - 1) % p;
    poly_normalize(&xp);

    poly_modp_t f_copy = poly_alloc(f.deg);
    for (int i = 0; i <= f.deg; i++) f_copy.coeff[i] = f.coeff[i];
    f_copy.deg = f.deg;

    poly_modp_t g = poly_gcd(xp, f_copy, p);
    int nroots = 0;

    if (g.deg <= 0) {
        /* no roots */
    } else if (g.deg == 1) {
        uint64_t a1_inv = mod_inverse(g.coeff[1], p);
        roots[nroots++] = (uint32_t)((uint64_t)(p - g.coeff[0]) * a1_inv % p);
    } else {
        /* Cantor-Zassenhaus splitting */
        poly_modp_t factors[MAX_DEGREE];
        int nfactors = 0;
        factors[nfactors] = poly_alloc(g.deg);
        for (int i = 0; i <= g.deg; i++) factors[nfactors].coeff[i] = g.coeff[i];
        factors[nfactors].deg = g.deg;
        nfactors = 1;

        uint32_t rng = p;
        for (int att = 0; att < 40 && nfactors < g.deg; att++) {
            rng = rng * 1103515245 + 12345;
            uint32_t a = rng % p;

            int new_nf = nfactors;
            for (int fi = 0; fi < nfactors; fi++) {
                if (factors[fi].deg <= 1) continue;

                poly_modp_t xa = poly_alloc(factors[fi].deg);
                xa.coeff[0] = a; xa.coeff[1] = 1; xa.deg = 1;

                poly_modp_t res = poly_alloc(factors[fi].deg);
                res.coeff[0] = 1; res.deg = 0;
                uint64_t exp = ((uint64_t)p - 1) / 2;

                while (exp > 0) {
                    if (exp & 1) {
                        poly_modp_t tmp = poly_mulmod(&res, &xa, &factors[fi], p);
                        poly_free(&res);
                        res = tmp;
                    }
                    exp >>= 1;
                    if (exp > 0) {
                        poly_modp_t tmp = poly_mulmod(&xa, &xa, &factors[fi], p);
                        poly_free(&xa);
                        xa = tmp;
                    }
                }
                poly_free(&xa);

                res.coeff[0] = (res.coeff[0] + p - 1) % p;
                poly_normalize(&res);

                poly_modp_t fi_copy = poly_alloc(factors[fi].deg);
                for (int k = 0; k <= factors[fi].deg; k++) fi_copy.coeff[k] = factors[fi].coeff[k];
                fi_copy.deg = factors[fi].deg;

                poly_modp_t h = poly_gcd(res, fi_copy, p);

                if (h.deg > 0 && h.deg < factors[fi].deg) {
                    /* Polynomial division to get quotient */
                    poly_modp_t quot = poly_alloc(factors[fi].deg);
                    poly_modp_t rem = poly_alloc(factors[fi].deg);
                    for (int k = 0; k <= factors[fi].deg; k++) rem.coeff[k] = factors[fi].coeff[k];
                    rem.deg = factors[fi].deg;

                    uint64_t h_inv = mod_inverse(h.coeff[h.deg], p);
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

                    poly_free(&factors[fi]);
                    factors[fi] = h;
                    if (new_nf < MAX_DEGREE) {
                        factors[new_nf] = quot;
                        new_nf++;
                    } else {
                        poly_free(&quot);
                    }
                } else {
                    poly_free(&h);
                }
            }
            nfactors = new_nf;
        }

        for (int fi = 0; fi < nfactors && nroots < MAX_DEGREE; fi++) {
            if (factors[fi].deg == 1) {
                uint64_t a1_inv = mod_inverse(factors[fi].coeff[1], p);
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
            fb_add(&rat_fb, p, p); /* projective root sentinel */
            continue;
        }
        uint64_t y0_mod = mpz_fdiv_ui(poly.Y0, p);
        uint32_t y1_inv = mod_inverse((uint32_t)y1_mod, p);
        uint64_t r = (uint64_t)(p - y0_mod % p) % p * y1_inv % p;
        fb_add(&rat_fb, p, (uint32_t)r);
    }

    /* Algebraic side */
    for (uint32_t p = 2; p <= poly.alim; p++) {
        if (!is_prime(p)) continue;
        uint32_t roots[MAX_DEGREE];
        int nroots = find_alg_roots(roots, p);
        for (int i = 0; i < nroots; i++)
            fb_add(&alg_fb, p, roots[i]);
    }

    printf("Rational FB: %u entries (up to %u)\n", rat_fb.count, poly.rlim);
    printf("Algebraic FB: %u entries (up to %u)\n", alg_fb.count, poly.alim);
}

/* ============================================================
 * Lattice reduction (Gauss/Lagrange for 2D)
 * ============================================================ */

static void reduce_lattice(vec2_t *e0, vec2_t *e1) {
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

    /* Ensure e0 is longer (e1 shorter) for convention */
    n0 = (double)e0->x * e0->x + (double)e0->y * e0->y;
    n1 = (double)e1->x * e1->x + (double)e1->y * e1->y;
    if (n0 < n1) {
        vec2_t tmp = *e0; *e0 = *e1; *e1 = tmp;
    }
}

static inline void ij_to_ab(int64_t *a, int64_t *b,
                              int32_t i, int32_t j,
                              const vec2_t *e0, const vec2_t *e1) {
    *a = (int64_t)i * e0->x + (int64_t)j * e1->x;
    *b = (int64_t)i * e0->y + (int64_t)j * e1->y;
}

/* ============================================================
 * Bucket operations
 * ============================================================ */

static void bucket_init(bucket_t *bk, uint32_t alloc) {
    bk->entries = malloc(alloc * sizeof(bucket_entry_t));
    bk->count = 0;
    bk->alloc = alloc;
}

static inline void bucket_add(bucket_t *bk, uint16_t pos, uint8_t logp) {
    if (__builtin_expect(bk->count < bk->alloc, 1)) {
        bucket_entry_t *e = &bk->entries[bk->count++];
        e->pos = pos;
        e->logp = logp;
    }
}

static void bucket_reset(bucket_t *bk) {
    bk->count = 0;
}

/* ============================================================
 * Core sieve routines
 * ============================================================ */

/*
 * For a factor-base prime p with root r (i.e. f(r) ≡ 0 mod p),
 * the sieve-hit condition in lattice coordinates (i,j) is:
 *
 *   a ≡ r*b (mod p)
 *   where a = i*e0.x + j*e1.x,  b = i*e0.y + j*e1.y
 *
 * This gives:  i * (e0.x - r*e0.y) ≡ -j * (e1.x - r*e1.y) (mod p)
 *
 * Let α0 = (e0.x - r*e0.y) mod p,  α1 = (e1.x - r*e1.y) mod p.
 *
 * If α0 ≠ 0: for each j, starting i = (-j * α1 * α0^{-1}) mod p,
 *             then step by p.
 *   Use recurrence: i_start[j+1] = (i_start[j] - di) mod p
 *     where di = α1 * α0^{-1} mod p.
 *
 * For the "line sieve" (medium primes, SMALL_PRIME_CUTOFF <= p < BUCKET_THRESH):
 *   We iterate over j, and for each j walk through positions at stride p.
 *
 * For the "bucket sieve" (large primes, p >= BUCKET_THRESH):
 *   We compute linear positions in the flat sieve array and store
 *   (block_index, offset_in_block, logp) in bucket lists.
 */

/* Compute sieve-hit parameters for a given FB entry in the lattice.
 * Returns 0 if degenerate (skip this prime), 1 otherwise. */
static int compute_sieve_params(uint32_t p, uint32_t r,
                                 const vec2_t *e0, const vec2_t *e1,
                                 uint32_t *out_di, uint32_t *out_i0) {
    int64_t a0 = ((int64_t)e0->x - (int64_t)r * e0->y) % (int64_t)p;
    int64_t a1 = ((int64_t)e1->x - (int64_t)r * e1->y) % (int64_t)p;
    if (a0 < 0) a0 += p;
    if (a1 < 0) a1 += p;

    if (a0 == 0) {
        /* Degenerate — either all positions hit or only specific j multiples of p.
         * For medium primes this is extremely rare; for simplicity we skip. */
        return 0;
    }

    uint32_t a0_inv = mod_inverse((uint32_t)a0, p);
    *out_di = (uint32_t)((uint64_t)a1 * a0_inv % p);
    /* Starting i for j=0: we want i ≡ 0 mod p in the range [0, I).
     * The sieve array maps i in [-I/2, I/2) to index i + I/2.
     * Starting index for j=0: (I/2) % p (first hit in row). */
    *out_i0 = (uint32_t)((uint64_t)(SIEVE_I / 2) % p);
    return 1;
}

/* Line sieve for medium primes: SMALL_PRIME_CUTOFF <= p < SIEVE_I.
 * This writes directly into the sieve array, one row at a time. */
static void line_sieve(uint8_t *sieve, const factor_base_t *fb,
                        const vec2_t *e0, const vec2_t *e1,
                        uint32_t skip_q) {
    for (uint32_t fi = 0; fi < fb->count; fi++) {
        uint32_t p = fb->entries[fi].p;
        if (p < SMALL_PRIME_CUTOFF) continue;
        if (p >= (uint32_t)SIEVE_I) break;  /* Primes >= SIEVE_I go to bucket sieve */
        if (p > max_sieve_prime) break;
        if (p == skip_q) continue;

        uint8_t logp = fb->entries[fi].logp;
        uint32_t r = fb->entries[fi].r;

        uint32_t di, i_start;
        if (!compute_sieve_params(p, r, e0, e1, &di, &i_start)) continue;

        /* Walk rows with recurrence */
        uint32_t is = i_start;
        for (uint32_t j = 0; j < (uint32_t)SIEVE_J; j++) {
            uint8_t *row = sieve + j * (uint32_t)SIEVE_I;
            for (uint32_t i = is; i < (uint32_t)SIEVE_I; i += p)
                row[i] += logp;
            is = is >= di ? is - di : is + p - di;
        }
    }
}

/* Bucket-fill phase for large primes: p >= BUCKET_THRESH.
 *
 * For p >= SIEVE_I, each row has at most 1 hit. We use flat-stride
 * enumeration: the flat sieve positions of all hits form a sequence
 * that can be computed using two alternating strides, avoiding the
 * per-row loop entirely.
 *
 * Key insight: in the flat array, hits step by s1 = SIEVE_I - di or
 * s2 = SIEVE_I - di + p depending on whether the column index wraps.
 * We precompute both strides and walk through hits directly.
 */
static void bucket_fill(bucket_t *buckets, const factor_base_t *fb,
                          const vec2_t *e0, const vec2_t *e1,
                          uint32_t skip_q) {
    for (uint32_t bi = 0; bi < (uint32_t)NUM_BLOCKS; bi++)
        bucket_reset(&buckets[bi]);

    for (uint32_t fi = 0; fi < fb->count; fi++) {
        uint32_t p = fb->entries[fi].p;
        if (p < (uint32_t)SIEVE_I) continue;  /* handled by line sieve */
        if (p > max_sieve_prime) break;  /* FB is sorted by p */
        if (p == skip_q) continue;

        uint8_t logp = fb->entries[fi].logp;
        uint32_t r = fb->entries[fi].r;

        uint32_t di, i_start;
        if (!compute_sieve_params(p, r, e0, e1, &di, &i_start)) continue;

        /* All primes here have p >= SIEVE_I, so at most 1 hit per row.
         * Walk rows with recurrence. */
        uint32_t is = i_start;
        for (uint32_t j = 0; j < (uint32_t)SIEVE_J; j++) {
            if (is < (uint32_t)SIEVE_I) {
                uint32_t flat = j * (uint32_t)SIEVE_I + is;
                bucket_add(&buckets[flat >> 15], (uint16_t)(flat & BLOCK_MASK), logp);
            }
            /* Branchless modular subtract */
            uint32_t a = is - di;
            uint32_t b = is + (p - di);
            is = (is >= di) ? a : b;
        }
    }
}

/* Apply bucket entries to the sieve array, one 32KB block at a time.
 * This is cache-friendly: each block fits in L1 cache. */
static void bucket_apply(uint8_t *sieve, const bucket_t *buckets) {
    for (uint32_t bi = 0; bi < (uint32_t)NUM_BLOCKS; bi++) {
        uint8_t *block = sieve + bi * (uint32_t)BLOCK_SIZE;
        const bucket_t *bk = &buckets[bi];
        uint32_t cnt = bk->count;
        const bucket_entry_t *ep = bk->entries;

        /* Unrolled: process 4 entries at a time for throughput */
        uint32_t k = 0;
        for (; k + 4 <= cnt; k += 4) {
            block[ep[k + 0].pos] += ep[k + 0].logp;
            block[ep[k + 1].pos] += ep[k + 1].logp;
            block[ep[k + 2].pos] += ep[k + 2].logp;
            block[ep[k + 3].pos] += ep[k + 3].logp;
        }
        for (; k < cnt; k++)
            block[ep[k].pos] += ep[k].logp;
    }
}

/* ============================================================
 * AVX512 sieve-array scan for candidates
 *
 * We scan both alg_sieve and rat_sieve simultaneously.  A position
 * is a candidate if alg_sieve[idx] >= alg_thresh AND rat_sieve[idx] >= rat_thresh.
 *
 * Using AVX512BW: compare 64 bytes at a time, AND the two masks,
 * then iterate over set bits.
 * ============================================================ */

/*
 * Optimized cofactorization.
 *
 * Key insight: for NFS, p | F(a,b) iff a ≡ r*b (mod p) for some root r.
 * Instead of trial-dividing by every FB prime (93K entries!), we:
 *   1. Check a ≡ r*b mod p using cheap integer arithmetic for each FB entry
 *   2. Only call expensive mpz_divexact when we know p divides the norm
 *   3. For primes p where p fits in 32 bits and a,b fit in 64 bits,
 *      the check (a mod p) == (r * (b mod p)) mod p is a single division
 *
 * But even checking all 93K entries with cheap arithmetic is still slow.
 * Better approach: trial divide by small primes (up to TRIAL_DIV_LIMIT)
 * since they're few and frequent, then check the cofactor size.
 * For medium/large primes, we rely on the sieve having identified the
 * right positions and just need to find which specific primes divide.
 *
 * Hybrid strategy:
 *   - Primes up to TRIAL_DIV_LIMIT: direct mpz_tdiv_ui (fast for small p)
 *   - Primes above TRIAL_DIV_LIMIT: use (a,b) mod p check against FB roots
 *     but only scan those FB entries in a precomputed "medium" range
 *   - Very large primes (near lpb): the cofactor IS the large prime
 */
#define TRIAL_DIV_LIMIT 4096

/* Indices into FB where primes cross thresholds */
static uint32_t alg_fb_medium_start;  /* first entry with p >= TRIAL_DIV_LIMIT */
static uint32_t rat_fb_medium_start;

static void compute_fb_split(void) {
    alg_fb_medium_start = alg_fb.count;
    for (uint32_t i = 0; i < alg_fb.count; i++) {
        if (alg_fb.entries[i].p >= TRIAL_DIV_LIMIT) {
            alg_fb_medium_start = i;
            break;
        }
    }
    rat_fb_medium_start = rat_fb.count;
    for (uint32_t i = 0; i < rat_fb.count; i++) {
        if (rat_fb.entries[i].p >= TRIAL_DIV_LIMIT) {
            rat_fb_medium_start = i;
            break;
        }
    }
    printf("FB split: alg small=%u medium=%u, rat small=%u medium=%u\n",
           alg_fb_medium_start, alg_fb.count - alg_fb_medium_start,
           rat_fb_medium_start, rat_fb.count - rat_fb_medium_start);
}

/* Trial-divide a cofactor by factor base, using root-check for medium/large.
 * fb: factor base, fb_med_start: index where medium primes start.
 * a, b: lattice coordinates (b > 0).
 * Returns 1 if cofactor is smooth (or has acceptable large prime). */
static int trial_divide_side(mpz_t cofactor, const mpz_t norm,
                               const factor_base_t *fb, uint32_t fb_med_start,
                               int64_t a, int64_t b,
                               uint32_t lpb_bits, uint32_t mfb_bits,
                               uint32_t skip_q,
                               uint32_t *primes_out, int *np_out) {
    int np = *np_out;
    mpz_set(cofactor, norm);

    /* Remove skip_q if given */
    if (skip_q > 1) {
        while (mpz_divisible_ui_p(cofactor, skip_q))
            mpz_divexact_ui(cofactor, cofactor, skip_q);
    }

    /* Phase 1: Small primes by direct mpz trial division.
     * These are few (< ~600 primes below 4096) and hit frequently. */
    for (uint32_t fi = 0; fi < fb_med_start; fi++) {
        uint32_t p = fb->entries[fi].p;
        if (p == skip_q) continue;
        /* Quick check: if cofactor < p, done with small primes */
        if (mpz_cmp_ui(cofactor, p) < 0) break;
        while (mpz_divisible_ui_p(cofactor, p)) {
            mpz_divexact_ui(cofactor, cofactor, p);
            if (np < MAX_REL_PRIMES) primes_out[np++] = p;
        }
    }

    if (mpz_cmp_ui(cofactor, 1) == 0) { *np_out = np; return 1; }

    /* Phase 2: Medium and large primes using root check.
     * For each FB entry (p, r), p divides F(a,b) iff a ≡ r*b (mod p).
     * We compute a_mod = ((a % p) + p) % p and check against (r * b_mod) % p.
     * This avoids expensive GMP division for non-divisors. */
    uint64_t b_pos = (uint64_t)b;  /* b > 0 guaranteed */

    for (uint32_t fi = fb_med_start; fi < fb->count; fi++) {
        uint32_t p = fb->entries[fi].p;
        if (p == skip_q) continue;
        /* Early exit: if cofactor < p^2, remaining factor (if any) is prime */
        {
            uint64_t p2 = (uint64_t)p * p;
            if (mpz_cmp_ui(cofactor, (unsigned long)p2) < 0) break;
        }

        uint32_t r = fb->entries[fi].r;

        /* Projective root sentinel */
        if (r == p) {
            /* p | b check */
            if (b_pos % p == 0) {
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    if (np < MAX_REL_PRIMES) primes_out[np++] = p;
                }
            }
            continue;
        }

        /* Root check: a ≡ r*b (mod p) */
        int64_t a_mod = a % (int64_t)p;
        if (a_mod < 0) a_mod += p;
        uint64_t b_mod = b_pos % p;
        uint64_t rb_mod = (uint64_t)r * b_mod % p;

        if ((uint64_t)a_mod != rb_mod) continue;

        /* This prime divides F(a,b) — extract all powers */
        while (mpz_divisible_ui_p(cofactor, p)) {
            mpz_divexact_ui(cofactor, cofactor, p);
            if (np < MAX_REL_PRIMES) primes_out[np++] = p;
        }
    }

    *np_out = np;

    /* Check cofactor */
    if (mpz_cmp_ui(cofactor, 1) == 0) return 1;

    size_t cof_bits = mpz_sizeinbase(cofactor, 2);
    if (cof_bits <= lpb_bits) {
        /* Single large prime */
        if (np < MAX_REL_PRIMES)
            primes_out[(*np_out)++] = (uint32_t)mpz_get_ui(cofactor);
        return 1;
    }
    if (mfb_bits > 0 && cof_bits <= mfb_bits) {
        /* Double large prime candidate: reject if cofactor is prime
         * (it would be a prime > 2^lpb, too large for single LP) */
        if (mpz_probab_prime_p(cofactor, 1)) return 0;
        /* Accept DLP — store cofactor value */
        if (np < MAX_REL_PRIMES)
            primes_out[(*np_out)++] = (uint32_t)mpz_get_ui(cofactor);
        return 1;
    }
    return 0;  /* cofactor too large */
}

/* Process candidate at flat index idx.  Returns 1 if a relation was found. */
static int process_candidate(uint32_t idx, const vec2_t *e0, const vec2_t *e1,
                               uint32_t q, mpz_t anorm, mpz_t rnorm,
                               mpz_t cofactor) {
    int32_t si = (int32_t)(idx % SIEVE_I) - SIEVE_I / 2;
    int32_t sj = (int32_t)(idx / SIEVE_I);

    int64_t a, b;
    ij_to_ab(&a, &b, si, sj, e0, e1);

    if (b == 0 || a == 0) return 0;
    if (b < 0) { a = -a; b = -b; }

    /* Evaluate norms */
    eval_alg_poly(anorm, &poly, a, b);
    mpz_abs(anorm, anorm);
    eval_rat_poly(rnorm, &poly, a, b);

    /* --- Algebraic side --- */
    uint32_t alg_primes[MAX_REL_PRIMES];
    int alg_np = 0;
    alg_primes[alg_np++] = q;

    if (!trial_divide_side(cofactor, anorm, &alg_fb, alg_fb_medium_start,
                            a, b, poly.lpba, poly.mfba, q,
                            alg_primes, &alg_np))
        return 0;

    /* --- Rational side --- */
    uint32_t rat_primes[MAX_REL_PRIMES];
    int rat_np = 0;

    if (!trial_divide_side(cofactor, rnorm, &rat_fb, rat_fb_medium_start,
                            a, b, poly.lpbr, poly.mfbr, 0,
                            rat_primes, &rat_np))
        return 0;

    /* Output relation */
    fprintf(outfp, "%ld,%ld:", (long)a, (long)b);
    for (int k = 0; k < alg_np; k++) {
        if (k > 0) fputc(',', outfp);
        fprintf(outfp, "%x", alg_primes[k]);
    }
    fputc(':', outfp);
    for (int k = 0; k < rat_np; k++) {
        if (k > 0) fputc(',', outfp);
        fprintf(outfp, "%x", rat_primes[k]);
    }
    fputc('\n', outfp);

    return 1;
}

static uint64_t scan_and_process(uint8_t *asieve, uint8_t *rsieve,
                                   uint8_t athresh, uint8_t rthresh,
                                   const vec2_t *e0, const vec2_t *e1,
                                   uint32_t q) {
    uint64_t rels = 0;
    mpz_t anorm, rnorm, cofactor;
    mpz_inits(anorm, rnorm, cofactor, NULL);

#ifdef __AVX512BW__
    /* AVX512BW path: compare 64 bytes at a time */
    __m512i va_thr = _mm512_set1_epi8((char)(athresh - 1));
    __m512i vr_thr = _mm512_set1_epi8((char)(rthresh - 1));

    uint32_t idx = 0;
    for (; idx + 64 <= SIEVE_AREA; idx += 64) {
        __m512i va = _mm512_loadu_si512((__m512i *)(asieve + idx));
        __m512i vr = _mm512_loadu_si512((__m512i *)(rsieve + idx));

        /* Compare unsigned: a > (thresh-1) means a >= thresh.
         * _mm512_cmpgt_epu8_mask does unsigned comparison. */
        __mmask64 ma = _mm512_cmpgt_epu8_mask(va, va_thr);
        __mmask64 mr = _mm512_cmpgt_epu8_mask(vr, vr_thr);
        __mmask64 both = ma & mr;

        while (both) {
            int bit = __builtin_ctzll(both);
            rels += process_candidate(idx + bit, e0, e1, q, anorm, rnorm, cofactor);
            both &= both - 1;  /* clear lowest set bit */
        }
    }
    /* Tail */
    for (; idx < SIEVE_AREA; idx++) {
        if (asieve[idx] >= athresh && rsieve[idx] >= rthresh)
            rels += process_candidate(idx, e0, e1, q, anorm, rnorm, cofactor);
    }
#else
    /* Scalar fallback: process 8 bytes at a time using uint64_t trick */
    for (uint32_t idx = 0; idx < SIEVE_AREA; idx++) {
        if (asieve[idx] >= athresh && rsieve[idx] >= rthresh)
            rels += process_candidate(idx, e0, e1, q, anorm, rnorm, cofactor);
    }
#endif

    mpz_clears(anorm, rnorm, cofactor, NULL);
    return rels;
}

/* ============================================================
 * Sieve one special-q
 * ============================================================ */

/* Accumulated timing for phases (in clock ticks) */
static clock_t time_clear = 0, time_line = 0, time_bucket = 0;
static clock_t time_apply = 0, time_scan = 0;

static void sieve_special_q(uint32_t q, uint32_t qroot) {
    vec2_t e0 = {q, 0};
    vec2_t e1 = {qroot, 1};
    reduce_lattice(&e0, &e1);

    clock_t t0, t1;

    /* Clear sieve arrays */
    t0 = clock();
    memset(alg_sieve, 0, SIEVE_AREA);
    memset(rat_sieve, 0, SIEVE_AREA);
    t1 = clock(); time_clear += t1 - t0;

    /* ---- Compute thresholds ---- */
    double log2_alg_norm = 0;
    {
        double max_a = (double)SIEVE_I * fmax(fabs((double)e0.x), fabs((double)e1.x));
        double max_b = (double)SIEVE_J * fmax(fabs((double)e0.y), fabs((double)e1.y));
        if (max_b < 1) max_b = 1;
        if (max_a < 1) max_a = 1;
        for (int i = 0; i <= poly.degree; i++) {
            double ci = fabs(mpz_get_d(poly.c[i]));
            if (ci > 0) {
                double term = log2(ci) + i * log2(max_a) + (poly.degree - i) * log2(max_b);
                if (term > log2_alg_norm) log2_alg_norm = term;
            }
        }
    }
    double log2_rat_norm = 0;
    {
        double max_a = (double)SIEVE_I * fmax(fabs((double)e0.x), fabs((double)e1.x));
        double max_b = (double)SIEVE_J * fmax(fabs((double)e0.y), fabs((double)e1.y));
        double y1 = fabs(mpz_get_d(poly.Y1));
        double y0 = fabs(mpz_get_d(poly.Y0));
        log2_rat_norm = log2(y1 * max_a + y0 * max_b + 1);
    }

    /* Threshold: we want the sieved log-sum to be close to the norm,
     * minus the contribution of unsieved small primes.
     * lambda controls how aggressively we filter.
     * thresh = norm_bits - norm_bits / (lambda + 1) = norm_bits * lambda / (lambda + 1) */
    uint8_t alg_thresh = (uint8_t)(log2_alg_norm * poly.alambda / (poly.alambda + 1));
    uint8_t rat_thresh = (uint8_t)(log2_rat_norm * poly.rlambda / (poly.rlambda + 1));

    /* Subtract special-q contribution from algebraic threshold */
    uint8_t logq = (uint8_t)(log2((double)q) + 0.5);
    if (alg_thresh > logq) alg_thresh -= logq;

    /* Subtract estimated unsieved prime contribution.
     *
     * Small primes (< SMALL_PRIME_CUTOFF): not sieved, handled by trial div.
     *   sum(log2(p) * degree/p) for p < 64 ≈ degree * 5 bits for algebraic.
     *
     * Large primes (> max_sieve_prime): not sieved, handled by root-check.
     *   sum(log2(p)/p) for p in (max_sieve_prime, alim] ≈
     *     integral of ln(x)/(x*ln(2)) dx = [ln(x)^2 / (2*ln(2))]
     *   from max_sieve_prime to alim.
     *   For degree-d polynomial, multiply by ~d (number of roots).
     */
    uint8_t small_est = (uint8_t)(poly.degree * 5);
    double large_est_alg = 0, large_est_rat = 0;
    if (max_sieve_prime < poly.alim) {
        double lnhi = log((double)poly.alim);
        double lnlo = log((double)max_sieve_prime);
        large_est_alg = poly.degree * (lnhi * lnhi - lnlo * lnlo) / (2.0 * log(2.0));
    }
    if (max_sieve_prime < poly.rlim) {
        double lnhi = log((double)poly.rlim);
        double lnlo = log((double)max_sieve_prime);
        large_est_rat = (lnhi * lnhi - lnlo * lnlo) / (2.0 * log(2.0));
    }

    uint8_t unsieved_alg = (uint8_t)(small_est + large_est_alg);
    uint8_t unsieved_rat = (uint8_t)(small_est / 2 + large_est_rat);

    static int debug_thresh = 1;
    if (debug_thresh) {
        printf("Threshold debug: log2_alg_norm=%.1f log2_rat_norm=%.1f\n",
               log2_alg_norm, log2_rat_norm);
        printf("  alg_thresh (before sub)=%u, unsieved_alg=%u (small=%u large=%.1f)\n",
               alg_thresh, unsieved_alg, small_est, large_est_alg);
        printf("  rat_thresh (before sub)=%u, unsieved_rat=%u (small=%u large=%.1f)\n",
               rat_thresh, unsieved_rat, small_est / 2, large_est_rat);
        debug_thresh = 0;
    }

    if (alg_thresh > unsieved_alg) alg_thresh -= unsieved_alg;
    else alg_thresh = 0;
    if (rat_thresh > unsieved_rat) rat_thresh -= unsieved_rat;
    else rat_thresh = 0;

    /* Safety floor: don't let thresholds go below a minimum */
    if (alg_thresh < 15) alg_thresh = 15;
    if (rat_thresh < 10) rat_thresh = 10;

    /* ---- Phase 1: Line sieve for medium primes ---- */
    t0 = clock();
    line_sieve(alg_sieve, &alg_fb, &e0, &e1, q);
    line_sieve(rat_sieve, &rat_fb, &e0, &e1, 0);
    t1 = clock(); time_line += t1 - t0;

    /* ---- Phase 2: Bucket sieve for large primes ---- */
    t0 = clock();
    bucket_fill(alg_buckets, &alg_fb, &e0, &e1, q);
    bucket_fill(rat_buckets, &rat_fb, &e0, &e1, 0);
    t1 = clock(); time_bucket += t1 - t0;

    t0 = clock();
    bucket_apply(alg_sieve, alg_buckets);
    bucket_apply(rat_sieve, rat_buckets);
    t1 = clock(); time_apply += t1 - t0;

    /* ---- Phase 3: Scan for candidates and trial-factor ---- */
    t0 = clock();
    uint64_t rels_this_q = scan_and_process(alg_sieve, rat_sieve,
                                             alg_thresh, rat_thresh,
                                             &e0, &e1, q);
    t1 = clock(); time_scan += t1 - t0;

    total_rels += rels_this_q;

    if (rels_this_q > 0) {
        printf("q=%u: %lu relations (total: %lu)\n",
               q, (unsigned long)rels_this_q, (unsigned long)total_rels);
        fflush(stdout);
    }
}

/* ============================================================
 * Main
 * ============================================================ */

int main(int argc, char *argv[]) {
    uint32_t start_q = 0, q_range = 0;
    char *job_file = NULL;
    char *out_file = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) start_q = atoi(argv[++i]);
        else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) q_range = atoi(argv[++i]);
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) out_file = argv[++i];
        else if (strcmp(argv[i], "-a") == 0 && i + 1 < argc) job_file = argv[++i];
        else if (strcmp(argv[i], "-S") == 0 && i + 1 < argc) max_sieve_prime = (uint32_t)atol(argv[++i]);
    }

    if (!job_file || !start_q || !q_range) {
        fprintf(stderr, "Usage: %s -f <startq> -c <qrange> -o <outfile> -a <jobfile> [-S max_sieve_prime]\n",
                argv[0]);
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

    /* Build prime sieve up to max(rlim, alim) */
    uint32_t max_lim = poly.rlim > poly.alim ? poly.rlim : poly.alim;
    printf("Building prime sieve up to %u...\n", max_lim);
    init_prime_sieve(max_lim);

    /* Build factor bases */
    build_factor_bases();

    /* Compute FB split points for optimized cofactorization */
    compute_fb_split();

    /* FB is built in order of p, so it's already sorted. */

    /* Allocate sieve arrays */
    alg_sieve = calloc(SIEVE_AREA, 1);
    rat_sieve = calloc(SIEVE_AREA, 1);
    if (!alg_sieve || !rat_sieve) {
        fprintf(stderr, "Failed to allocate sieve arrays (%u bytes each)\n", SIEVE_AREA);
        return 1;
    }

    /* Allocate bucket arrays */
    alg_buckets = malloc(NUM_BLOCKS * sizeof(bucket_t));
    rat_buckets = malloc(NUM_BLOCKS * sizeof(bucket_t));
    /* Estimate bucket fill: total large-prime hits ≈ SIEVE_AREA * sum(1/p for p>=BUCKET_THRESH)
     * For alim=1.2M: sum ≈ ln(1200000/256) ≈ 8.5, so total ≈ 8M/256*8.5 ≈ hmm.
     * A simpler heuristic: alloc proportional to SIEVE_AREA / BUCKET_THRESH per block.
     * With NUM_BLOCKS = 256, and ~degree * sum(1/p) per position ≈ 35:
     * total_hits ≈ 35 * 8M = 280M — no, that's too big.
     * Actually sum(logp/p) is what's sieved, and sum(1/p) for p in [256,1.2M]
     * is about 8.5, so there are about 8.5 * SIEVE_AREA / (avg_p) hits?
     * No — there is 1 hit per p per SIEVE_I entries (per row), and SIEVE_J rows.
     * hits_per_prime ≈ SIEVE_AREA / p.
     * total_hits = sum_{p>=256}(SIEVE_AREA/p) ≈ SIEVE_AREA * sum(1/p, p=256..1.2M)
     *            ≈ 8M * 8.5 ≈ 68M — too many to store all!
     * But per block: ≈ 68M / 256 ≈ 266K per block. Use 512K to be safe. */
    uint32_t bucket_alloc_per_block = (1 << 19);  /* 512K entries per block */
    for (uint32_t i = 0; i < (uint32_t)NUM_BLOCKS; i++) {
        bucket_init(&alg_buckets[i], bucket_alloc_per_block);
        bucket_init(&rat_buckets[i], bucket_alloc_per_block);
    }

    printf("Sieve area: I=%d, J=%d (%u bytes)\n", SIEVE_I, SIEVE_J, SIEVE_AREA);
    printf("Blocks: %d x %d bytes, bucket alloc: %u/block\n",
           NUM_BLOCKS, BLOCK_SIZE, bucket_alloc_per_block);
    printf("Sieving q from %u to %u...\n", start_q, start_q + q_range);

    clock_t start_time = clock();
    uint32_t nq = 0;

    for (uint32_t q = start_q; q < start_q + q_range; q++) {
        if (!is_prime(q)) continue;

        uint32_t roots[MAX_DEGREE];
        int nroots = find_alg_roots(roots, q);

        for (int ri = 0; ri < nroots; ri++) {
            sieve_special_q(q, roots[ri]);
        }

        nq++;
        if (nq % 100 == 0) {
            double elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
            printf("[progress] q=%u, %u special-q's, %lu rels, %.1f sec, %.1f rels/sec\n",
                   q, nq, (unsigned long)total_rels, elapsed,
                   total_rels / (elapsed + 0.001));
            fflush(stdout);
        }
    }

    double elapsed = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    printf("\nTotal: %lu relations from %u special-q's in %.1f seconds (%.1f rels/sec)\n",
           (unsigned long)total_rels, nq, elapsed, total_rels / (elapsed + 0.001));

    printf("Phase timing: clear=%.2fs line=%.2fs bucket=%.2fs apply=%.2fs scan=%.2fs\n",
           (double)time_clear / CLOCKS_PER_SEC,
           (double)time_line / CLOCKS_PER_SEC,
           (double)time_bucket / CLOCKS_PER_SEC,
           (double)time_apply / CLOCKS_PER_SEC,
           (double)time_scan / CLOCKS_PER_SEC);

    fclose(outfp);
    free(alg_sieve);
    free(rat_sieve);
    free(prime_sieve_bits);
    for (uint32_t i = 0; i < (uint32_t)NUM_BLOCKS; i++) {
        free(alg_buckets[i].entries);
        free(rat_buckets[i].entries);
    }
    free(alg_buckets);
    free(rat_buckets);
    free(rat_fb.entries);
    free(alg_fb.entries);

    return 0;
}
