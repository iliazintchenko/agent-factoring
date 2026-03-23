/*
 * batch_siqs.c - SIQS with Batch Smoothness Detection
 *
 * Instead of traditional sieving (O(M * B) per polynomial), this uses
 * Bernstein's batch smoothness test: generate Q(x) values from many
 * SIQS polynomials, then test them all at once using product/remainder
 * trees and iterated GCDs.
 *
 * The key insight: traditional sieve is fast because it uses byte-level
 * L1-cache operations. But it scales as O(FB_size) per polynomial.
 * Batch smoothness is O(1) per candidate (amortized over the batch),
 * independent of factor base size. For large FB, batch could win.
 *
 * Compile: gcc -O3 -march=native -o batch_siqs library/batch_siqs.c -lgmp -lm
 * Usage: ./batch_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_FB 50000
#define MAX_RELATIONS 200000
#define SEED 42

typedef struct {
    unsigned int p;
    int r;        /* sqrt(kN) mod p */
    double logp;
} fb_entry_t;

typedef struct {
    mpz_t qx;           /* Q(x) value = ((ax+b)^2 - kN) / a */
    mpz_t axpb;         /* ax + b (for x-side of congruence) */
    unsigned char *exponents; /* exponent vector mod 2 */
    unsigned int lp;     /* large prime (0 if fully smooth) */
} relation_t;

static mpz_t N, kN;
static int k_mult;     /* Knuth-Schroeppel multiplier */
static fb_entry_t fb[MAX_FB];
static int fb_count;
static relation_t relations[MAX_RELATIONS];
static int rel_count;

/* Simple sieve for primes up to limit */
static int *gen_primes(int limit, int *count) {
    char *sieve = calloc(limit + 1, 1);
    int *primes = malloc(sizeof(int) * (limit / 2 + 100));
    *count = 0;
    for (int i = 2; i <= limit; i++) {
        if (!sieve[i]) {
            primes[(*count)++] = i;
            for (long j = (long)i * i; j <= limit; j += i)
                sieve[j] = 1;
        }
    }
    free(sieve);
    return primes;
}

/* Tonelli-Shanks modular square root */
static int modsqrt(unsigned long n, unsigned int p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    /* Euler criterion */
    mpz_t tmp, mod;
    mpz_init_set_ui(tmp, n);
    mpz_init_set_ui(mod, p);
    mpz_powm_ui(tmp, tmp, (p - 1) / 2, mod);
    if (mpz_cmp_ui(tmp, 1) != 0) {
        mpz_clear(tmp); mpz_clear(mod);
        return -1;
    }

    if (p % 4 == 3) {
        mpz_set_ui(tmp, n);
        mpz_t exp;
        mpz_init_set_ui(exp, (p + 1) / 4);
        mpz_powm(tmp, tmp, exp, mod);
        int r = mpz_get_ui(tmp);
        mpz_clear(tmp); mpz_clear(mod); mpz_clear(exp);
        return r;
    }

    /* Full Tonelli-Shanks */
    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned int z = 2;
    while (1) {
        mpz_set_ui(tmp, z);
        mpz_powm_ui(tmp, tmp, (p - 1) / 2, mod);
        if (mpz_cmp_ui(tmp, p - 1) == 0) break;
        z++;
    }

    mpz_t M, c, t, R;
    mpz_init_set_ui(M, S);
    mpz_init(c); mpz_init(t); mpz_init(R);

    mpz_set_ui(c, z);
    mpz_t qq; mpz_init_set_ui(qq, Q);
    mpz_powm(c, c, qq, mod);

    mpz_set_ui(t, n);
    mpz_powm(t, t, qq, mod);

    mpz_set_ui(R, n);
    mpz_set_ui(qq, (Q + 1) / 2);
    mpz_powm(R, R, qq, mod);

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            int r = mpz_get_ui(R);
            mpz_clear(M); mpz_clear(c); mpz_clear(t); mpz_clear(R);
            mpz_clear(tmp); mpz_clear(mod); mpz_clear(qq);
            return r;
        }

        int i = 0;
        mpz_set(qq, t);
        while (mpz_cmp_ui(qq, 1) != 0) {
            mpz_mul(qq, qq, qq);
            mpz_mod(qq, qq, mod);
            i++;
        }

        int m = mpz_get_ui(M);
        mpz_t b; mpz_init_set(b, c);
        for (int j = 0; j < m - i - 1; j++) {
            mpz_mul(b, b, b);
            mpz_mod(b, b, mod);
        }

        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, mod);
        mpz_mul(t, t, c); mpz_mod(t, t, mod);
        mpz_mul(R, R, b); mpz_mod(R, R, mod);
        mpz_clear(b);
    }
}

/* Knuth-Schroeppel multiplier selection */
static int select_multiplier(mpz_t n) {
    static const int candidates[] = {1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 0};
    double best_score = -1e9;
    int best_k = 1;

    for (int i = 0; candidates[i]; i++) {
        int k = candidates[i];
        mpz_t kn;
        mpz_init(kn);
        mpz_mul_ui(kn, n, k);

        double score = -0.5 * log(k);
        unsigned long kn_mod8 = mpz_fdiv_ui(kn, 8);

        if (kn_mod8 == 1) score += 2 * log(2);
        else if (kn_mod8 == 5) score += log(2);
        else score += 0.5 * log(2);

        /* Check a few small primes */
        int pc;
        int *ps = gen_primes(200, &pc);
        for (int j = 0; j < pc && ps[j] < 200; j++) {
            int p = ps[j];
            if (p == 2) continue;
            if (k % p == 0) continue;
            unsigned long kn_mod_p = mpz_fdiv_ui(kn, p);
            mpz_t tmp, mod;
            mpz_init_set_ui(tmp, kn_mod_p);
            mpz_init_set_ui(mod, p);
            mpz_powm_ui(tmp, tmp, (p - 1) / 2, mod);
            if (mpz_cmp_ui(tmp, 1) == 0)
                score += 2.0 * log(p) / (p - 1);
            mpz_clear(tmp); mpz_clear(mod);
        }
        free(ps);

        if (score > best_score) {
            best_score = score;
            best_k = k;
        }
        mpz_clear(kn);
    }
    return best_k;
}

/* Build factor base for kN */
static int build_fb(int smooth_bound) {
    fb_count = 0;

    /* Sign */
    fb[fb_count].p = 0;
    fb[fb_count].r = 0;
    fb[fb_count].logp = 0;
    fb_count++;

    int pc;
    int *ps = gen_primes(smooth_bound, &pc);

    for (int i = 0; i < pc; i++) {
        int p = ps[i];
        unsigned long kn_mod = mpz_fdiv_ui(kN, p);
        int r = modsqrt(kn_mod, p);
        if (r < 0) continue;

        fb[fb_count].p = p;
        fb[fb_count].r = r;
        fb[fb_count].logp = log2(p);
        fb_count++;
        if (fb_count >= MAX_FB) break;
    }

    free(ps);
    return fb_count;
}

/*
 * SIQS polynomial generation:
 * Q(x) = ((a*x + b)^2 - kN) / a = a*x^2 + 2*b*x + c
 * where c = (b^2 - kN) / a
 *
 * a = product of s factor base primes (each ~sqrt(2*kN) / M)
 * b chosen such that b^2 ≡ kN (mod a)
 *
 * For batch smoothness, we generate ALL Q(x) values for x in [-M, M]
 * across many polynomials, then batch-test.
 */

typedef struct {
    mpz_t a, b, c;
    int a_primes[20]; /* indices into fb of primes dividing a */
    int s;            /* number of primes in a */
} siqs_poly_t;

/* Generate a SIQS polynomial */
static int gen_polynomial(siqs_poly_t *poly, gmp_randstate_t rstate,
                         int M, int digits) {
    /* Target: a ≈ sqrt(2*kN) / M */
    mpz_t target, two_kn;
    mpz_init(target);
    mpz_init(two_kn);
    mpz_mul_ui(two_kn, kN, 2);
    mpz_sqrt(target, two_kn);
    mpz_tdiv_q_ui(target, target, M);

    /* Choose s factor base primes whose product ≈ target */
    /* Use primes from middle of factor base */
    int s = 0;
    double log_target = mpz_sizeinbase(target, 2);
    double log_a = 0;

    mpz_init_set_ui(poly->a, 1);
    mpz_init(poly->b);
    mpz_init(poly->c);

    /* Start from random position in upper half of FB */
    int start_idx = fb_count / 3 + gmp_urandomm_ui(rstate, fb_count / 3);
    int used[MAX_FB] = {0};

    while (log_a < log_target - 5 && s < 15) {
        int idx = start_idx + gmp_urandomm_ui(rstate, fb_count / 3);
        if (idx >= fb_count || idx < 2) idx = 2 + gmp_urandomm_ui(rstate, fb_count - 2);
        if (used[idx]) continue;
        if (fb[idx].p < 100) continue;

        used[idx] = 1;
        poly->a_primes[s] = idx;
        mpz_mul_ui(poly->a, poly->a, fb[idx].p);
        log_a += fb[idx].logp;
        s++;
    }
    poly->s = s;

    if (s < 2) {
        mpz_clear(target); mpz_clear(two_kn);
        return 0;
    }

    /* Compute b via CRT: b^2 ≡ kN (mod a) */
    /* For each prime p_i | a, compute b_i = sqrt(kN) mod p_i */
    /* Then use CRT to combine */

    /* Simple approach for small s: use Garner's algorithm */
    mpz_t *bi = malloc(sizeof(mpz_t) * s);
    for (int i = 0; i < s; i++) {
        mpz_init(bi[i]);
        int idx = poly->a_primes[i];
        unsigned int p = fb[idx].p;
        int r = fb[idx].r;

        /* b_i = r * (a/p_i)^(-1) mod p_i * (a/p_i) */
        mpz_t a_over_p, inv;
        mpz_init(a_over_p);
        mpz_init(inv);
        mpz_divexact_ui(a_over_p, poly->a, p);

        /* inv = (a/p_i)^(-1) mod p_i */
        mpz_t ptmp;
        mpz_init_set_ui(ptmp, p);
        mpz_invert(inv, a_over_p, ptmp);

        /* bi = r * inv mod p_i * (a/p_i) */
        mpz_mul_ui(bi[i], inv, r);
        mpz_mod_ui(bi[i], bi[i], p);
        mpz_mul(bi[i], bi[i], a_over_p);

        mpz_clear(a_over_p); mpz_clear(inv); mpz_clear(ptmp);
    }

    /* b = sum of bi mod a, adjusted to be in [-a/2, a/2] */
    mpz_set_ui(poly->b, 0);
    for (int i = 0; i < s; i++) {
        mpz_add(poly->b, poly->b, bi[i]);
    }
    mpz_mod(poly->b, poly->b, poly->a);

    /* Ensure b^2 ≡ kN (mod a) */
    mpz_t check;
    mpz_init(check);
    mpz_mul(check, poly->b, poly->b);
    mpz_sub(check, check, kN);
    if (!mpz_divisible_p(check, poly->a)) {
        /* Try negating some bi components */
        /* CRT gives us b and a-b; try a-b */
        mpz_sub(poly->b, poly->a, poly->b);
        mpz_mul(check, poly->b, poly->b);
        mpz_sub(check, check, kN);
        if (!mpz_divisible_p(check, poly->a)) {
            /* Failed to find valid b */
            for (int i = 0; i < s; i++) mpz_clear(bi[i]);
            free(bi);
            mpz_clear(check); mpz_clear(target); mpz_clear(two_kn);
            return 0;
        }
    }

    /* c = (b^2 - kN) / a */
    mpz_mul(poly->c, poly->b, poly->b);
    mpz_sub(poly->c, poly->c, kN);
    mpz_divexact(poly->c, poly->c, poly->a);

    for (int i = 0; i < s; i++) mpz_clear(bi[i]);
    free(bi);
    mpz_clear(check); mpz_clear(target); mpz_clear(two_kn);
    return 1;
}

/* Compute Q(x) = a*x^2 + 2*b*x + c for SIQS polynomial */
static void compute_qx(mpz_t result, siqs_poly_t *poly, int x) {
    /* result = a*x*x + 2*b*x + c */
    mpz_mul_si(result, poly->a, x);
    mpz_add(result, result, poly->b);
    mpz_add(result, result, poly->b);
    mpz_mul_si(result, result, x);
    mpz_add(result, result, poly->c);
    /* Alternative: result = ((a*x + b)^2 - kN) / a
     * But we already have c = (b^2 - kN)/a, so Q(x) = a*x^2 + 2*b*x + c */
}

/* Compute ax + b for the x-side of congruence */
static void compute_axpb(mpz_t result, siqs_poly_t *poly, int x) {
    mpz_mul_si(result, poly->a, x);
    mpz_add(result, result, poly->b);
}

/*
 * Batch smoothness detection using product/remainder trees.
 *
 * Given values v[0..n-1], test each for B-smoothness.
 * Returns: number of smooth/SLP values found and added to relations[].
 */
static int batch_detect(mpz_t *vals, mpz_t *axpbs, int n,
                       unsigned long lp_bound, int *a_prime_indices,
                       int a_prime_count) {
    if (n == 0) return 0;

    /* Compute primorial = product of p^k for all FB primes where p^k fits */
    mpz_t primorial;
    mpz_init_set_ui(primorial, 1);

    for (int i = 1; i < fb_count; i++) {
        unsigned int p = fb[i].p;
        unsigned long pk = p;
        /* Use p^k where p^k < 2^60 (so product stays reasonable) */
        while (pk <= 1000000000000000UL / p) pk *= p;
        mpz_mul_ui(primorial, primorial, pk);
    }

    /* Take absolute values */
    mpz_t *abs_vals = malloc(sizeof(mpz_t) * n);
    for (int i = 0; i < n; i++) {
        mpz_init(abs_vals[i]);
        mpz_abs(abs_vals[i], vals[i]);
        if (mpz_sgn(abs_vals[i]) == 0) mpz_set_ui(abs_vals[i], 1);
    }

    /* Build product tree */
    /* Product tree: level 0 = leaves, level k = pairwise products */
    int levels = 1;
    { int sz = n; while (sz > 1) { sz = (sz + 1) / 2; levels++; } }

    int *lev_size = malloc(sizeof(int) * levels);
    lev_size[0] = n;
    for (int l = 1; l < levels; l++)
        lev_size[l] = (lev_size[l-1] + 1) / 2;

    mpz_t **tree = malloc(sizeof(mpz_t*) * levels);
    tree[0] = abs_vals; /* Level 0 is the input */

    for (int l = 1; l < levels; l++) {
        tree[l] = malloc(sizeof(mpz_t) * lev_size[l]);
        for (int i = 0; i < lev_size[l]; i++) {
            mpz_init(tree[l][i]);
            if (2*i+1 < lev_size[l-1])
                mpz_mul(tree[l][i], tree[l-1][2*i], tree[l-1][2*i+1]);
            else
                mpz_set(tree[l][i], tree[l-1][2*i]);
        }
    }

    /* Remainder tree: compute primorial mod each leaf */
    mpz_t **rem = malloc(sizeof(mpz_t*) * levels);
    for (int l = 0; l < levels; l++) {
        rem[l] = malloc(sizeof(mpz_t) * lev_size[l]);
        for (int i = 0; i < lev_size[l]; i++)
            mpz_init(rem[l][i]);
    }

    int top = levels - 1;
    mpz_mod(rem[top][0], primorial, tree[top][0]);

    for (int l = top - 1; l >= 0; l--) {
        for (int i = 0; i < lev_size[l]; i++) {
            mpz_mod(rem[l][i], rem[l+1][i/2], tree[l][i]);
        }
    }

    /* For each candidate, check smoothness */
    int found = 0;
    mpz_t cofactor, gcd_val;
    mpz_init(cofactor);
    mpz_init(gcd_val);
    unsigned char *exponents = calloc(fb_count, 1);

    for (int i = 0; i < n && rel_count < MAX_RELATIONS; i++) {
        /* Iterated GCD to extract smooth part */
        mpz_set(cofactor, abs_vals[i]);

        /* First pass: GCD with remainder */
        mpz_gcd(gcd_val, rem[0][i], cofactor);
        if (mpz_cmp_ui(gcd_val, 1) > 0)
            mpz_divexact(cofactor, cofactor, gcd_val);

        /* Iterate until stable */
        int iters = 0;
        while (mpz_cmp_ui(cofactor, 1) > 0 && iters < 20) {
            mpz_gcd(gcd_val, cofactor, primorial);
            if (mpz_cmp_ui(gcd_val, 1) <= 0) break;
            mpz_divexact(cofactor, cofactor, gcd_val);
            iters++;
        }

        /* Check result */
        int smooth_type = 0;
        unsigned int lp_val = 0;

        if (mpz_cmp_ui(cofactor, 1) == 0) {
            smooth_type = 1; /* Fully smooth */
        } else if (mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= lp_bound) {
            unsigned long cf = mpz_get_ui(cofactor);
            /* Verify cofactor is prime (or 1) */
            mpz_t cfz;
            mpz_init_set_ui(cfz, cf);
            if (mpz_probab_prime_p(cfz, 10)) {
                smooth_type = 2; /* SLP */
                lp_val = cf;
            }
            mpz_clear(cfz);
        }

        if (smooth_type) {
            /* Trial divide to get exact exponents */
            int sign = (mpz_sgn(vals[i]) < 0) ? 1 : 0;
            memset(exponents, 0, fb_count);
            exponents[0] = sign;

            mpz_t tdval;
            mpz_init(tdval);
            mpz_abs(tdval, vals[i]);

            for (int j = 1; j < fb_count; j++) {
                unsigned int p = fb[j].p;
                while (mpz_divisible_ui_p(tdval, p)) {
                    exponents[j] ^= 1;
                    mpz_divexact_ui(tdval, tdval, p);
                }
            }

            /* Mark a-primes in exponent vector (they divide a, which divides (ax+b)^2) */
            /* Actually for SIQS: Q(x) = ((ax+b)^2 - kN)/a, so a divides out */
            /* The a-primes appear in the factorization of a*Q(x) = (ax+b)^2 - kN */
            /* But we're tracking Q(x) = a*x^2 + 2bx + c directly */
            /* a-primes DON'T appear in Q(x) exponents (they were divided out) */

            /* Verify cofactor */
            if (smooth_type == 1 && mpz_cmp_ui(tdval, 1) != 0) {
                mpz_clear(tdval);
                continue; /* false positive */
            }
            if (smooth_type == 2) {
                unsigned long rem_val = mpz_get_ui(tdval);
                if (rem_val > lp_bound) {
                    mpz_clear(tdval);
                    continue;
                }
                lp_val = rem_val;
            }
            mpz_clear(tdval);

            /* Add relation */
            relation_t *r = &relations[rel_count];
            mpz_init_set(r->qx, vals[i]);
            mpz_init_set(r->axpb, axpbs[i]);
            r->exponents = malloc(fb_count);
            memcpy(r->exponents, exponents, fb_count);
            r->lp = lp_val;
            rel_count++;
            found++;
        }
    }

    /* Cleanup */
    free(exponents);
    mpz_clear(cofactor);
    mpz_clear(gcd_val);

    for (int l = 1; l < levels; l++) {
        for (int i = 0; i < lev_size[l]; i++)
            mpz_clear(tree[l][i]);
        free(tree[l]);
    }
    for (int l = 0; l < levels; l++) {
        for (int i = 0; i < lev_size[l]; i++)
            mpz_clear(rem[l][i]);
        free(rem[l]);
    }
    for (int i = 0; i < n; i++) mpz_clear(abs_vals[i]);
    free(abs_vals);
    free(tree); free(rem); free(lev_size);
    mpz_clear(primorial);

    return found;
}

/* Gaussian elimination over GF(2) */
static int gaussian_elim(unsigned char **matrix, int rows, int cols,
                        unsigned char **history) {
    for (int i = 0; i < rows; i++) {
        history[i] = calloc((rows + 7) / 8, 1);
        history[i][i / 8] |= (1 << (i % 8));
    }

    int rank = 0;
    for (int col = 0; col < cols && rank < rows; col++) {
        int piv = -1;
        for (int r = rank; r < rows; r++) {
            if (matrix[r][col]) { piv = r; break; }
        }
        if (piv < 0) continue;

        if (piv != rank) {
            unsigned char *tmp;
            tmp = matrix[piv]; matrix[piv] = matrix[rank]; matrix[rank] = tmp;
            tmp = history[piv]; history[piv] = history[rank]; history[rank] = tmp;
        }

        for (int r = 0; r < rows; r++) {
            if (r != rank && matrix[r][col]) {
                for (int c = 0; c < cols; c++)
                    matrix[r][c] ^= matrix[rank][c];
                int hb = (rows + 7) / 8;
                for (int b = 0; b < hb; b++)
                    history[r][b] ^= history[rank][b];
            }
        }
        rank++;
    }
    return rank;
}

/* Try to extract factor from a dependency */
static int try_factor(int *dep, int dep_count, mpz_t factor) {
    mpz_t x, y, tmp;
    mpz_init_set_ui(x, 1);
    mpz_init_set_ui(y, 1);
    mpz_init(tmp);

    /* x = product of (ax+b) mod N */
    /* y = sqrt(product of a*Q(x)) mod N */
    /* Note: (ax+b)^2 = a*Q(x) + kN, so (ax+b)^2 ≡ a*Q(x) (mod kN) */
    /* Actually: (ax+b)^2 ≡ a*Q(x) (mod N) since kN ≡ 0 would need k|N */
    /* Correct: (ax+b)^2 - kN = a*Q(x), so (ax+b)^2 ≡ a*Q(x) (mod N) only if k=1 */
    /* For k>1: (ax+b)^2 ≡ a*Q(x) (mod kN). We need mod N. */
    /* Solution: x = product of (ax+b) mod N, y² = product of a*Q(x) mod N */

    /* Track full exponents for y computation */
    int *total_exp = calloc(fb_count, sizeof(int));

    for (int i = 0; i < dep_count; i++) {
        relation_t *r = &relations[dep[i]];

        /* x *= (ax+b) mod N */
        mpz_mul(x, x, r->axpb);
        mpz_mod(x, x, N);

        /* Count exponents of |Q(x)| */
        mpz_abs(tmp, r->qx);
        for (int j = 1; j < fb_count; j++) {
            unsigned int p = fb[j].p;
            while (mpz_divisible_ui_p(tmp, p)) {
                total_exp[j]++;
                mpz_divexact_ui(tmp, tmp, p);
            }
        }
        if (mpz_sgn(r->qx) < 0) total_exp[0]++;
    }

    /* Check all exponents are even */
    for (int j = 0; j < fb_count; j++) {
        if (total_exp[j] % 2 != 0) {
            free(total_exp);
            mpz_clear(x); mpz_clear(y); mpz_clear(tmp);
            return 0;
        }
    }

    /* y = product of p^(e/2) mod N */
    for (int j = 1; j < fb_count; j++) {
        if (total_exp[j] > 0) {
            mpz_set_ui(tmp, fb[j].p);
            mpz_powm_ui(tmp, tmp, total_exp[j] / 2, N);
            mpz_mul(y, y, tmp);
            mpz_mod(y, y, N);
        }
    }
    free(total_exp);

    /* factor = gcd(x - y, N) or gcd(x + y, N) */
    mpz_sub(tmp, x, y);
    mpz_gcd(factor, tmp, N);
    int found = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);

    if (!found) {
        mpz_add(tmp, x, y);
        mpz_gcd(factor, tmp, N);
        found = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);
    }

    mpz_clear(x); mpz_clear(y); mpz_clear(tmp);
    return found;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_init(N);
    mpz_init(kN);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);
    fprintf(stderr, "Factoring %d-digit number\n", digits);

    /* Check small factors */
    for (int p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact_ui(cofactor, N, p);
            gmp_printf("%Zd = %d * %Zd\n", N, p, cofactor);
            mpz_clear(cofactor);
            return 0;
        }
    }

    /* Multiplier */
    k_mult = select_multiplier(N);
    mpz_mul_ui(kN, N, k_mult);
    fprintf(stderr, "Multiplier k = %d\n", k_mult);

    /* Smoothness bound */
    double ln_n = digits * log(10);
    double ln_ln_n = log(ln_n);
    int smooth_bound = (int)exp(0.5 * sqrt(ln_n * ln_ln_n));
    /* Tuning: slightly smaller FB to get more smooth numbers relative to FB size */
    smooth_bound = (int)(smooth_bound * 0.7);
    if (smooth_bound < 300) smooth_bound = 300;
    if (smooth_bound > 2000000) smooth_bound = 2000000;

    build_fb(smooth_bound);
    fprintf(stderr, "Factor base: %d primes (bound=%d), largest=%u\n",
            fb_count, smooth_bound, fb[fb_count-1].p);

    /* Sieve interval half-width M */
    int M = (int)exp(0.5 * sqrt(ln_n * ln_ln_n)) * 2;
    if (M < 5000) M = 5000;
    if (M > 500000) M = 500000;

    /* Large prime bound */
    unsigned long lp_bound = (unsigned long)fb[fb_count-1].p * 50;
    if (lp_bound > (1UL << 28)) lp_bound = (1UL << 28);

    int target = fb_count + 30;
    fprintf(stderr, "M=%d, lp_bound=%lu, target=%d rels\n", M, lp_bound, target);

    /* Initialize RNG */
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, SEED);

    /* Relation collection */
    rel_count = 0;
    int poly_count = 0;
    int total_tested = 0;

    /* Batch buffer */
    int batch_alloc = 2 * M;
    if (batch_alloc > 50000) batch_alloc = 50000;
    mpz_t *batch_vals = malloc(sizeof(mpz_t) * batch_alloc);
    mpz_t *batch_axpb = malloc(sizeof(mpz_t) * batch_alloc);
    for (int i = 0; i < batch_alloc; i++) {
        mpz_init(batch_vals[i]);
        mpz_init(batch_axpb[i]);
    }

    while (rel_count < target) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        if (elapsed > 290.0) {
            fprintf(stderr, "Timeout at %.1fs with %d relations\n", elapsed, rel_count);
            break;
        }

        /* Generate batch of candidates from multiple polynomials */
        int bcnt = 0;

        while (bcnt < batch_alloc) {
            siqs_poly_t poly;
            if (!gen_polynomial(&poly, rstate, M, digits)) continue;
            poly_count++;

            /* Generate Q(x) for x in [-M/2, M/2] */
            int half = M / 2;
            for (int x = -half; x <= half && bcnt < batch_alloc; x++) {
                compute_qx(batch_vals[bcnt], &poly, x);
                if (mpz_sgn(batch_vals[bcnt]) == 0) continue;
                compute_axpb(batch_axpb[bcnt], &poly, x);
                bcnt++;
            }
            mpz_clear(poly.a); mpz_clear(poly.b); mpz_clear(poly.c);
        }

        total_tested += bcnt;

        /* Batch smoothness test */
        int found = batch_detect(batch_vals, batch_axpb, bcnt, lp_bound,
                                NULL, 0);

        clock_gettime(CLOCK_MONOTONIC, &t1);
        elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        fprintf(stderr, "[%.1fs] polys=%d tested=%d found=%d rels=%d/%d (%.0f/s)\n",
                elapsed, poly_count, total_tested, found,
                rel_count, target, rel_count / (elapsed > 0 ? elapsed : 1));
    }

    /* SLP pair merging */
    int slp_count = 0;
    int *slp_idx = malloc(sizeof(int) * rel_count);
    for (int i = 0; i < rel_count; i++)
        if (relations[i].lp > 0) slp_idx[slp_count++] = i;

    /* Sort by LP */
    for (int i = 1; i < slp_count; i++) {
        int key = slp_idx[i];
        unsigned int kv = relations[key].lp;
        int j = i - 1;
        while (j >= 0 && relations[slp_idx[j]].lp > kv) {
            slp_idx[j+1] = slp_idx[j]; j--;
        }
        slp_idx[j+1] = key;
    }

    /* Merge pairs */
    int merged_count = 0;
    for (int i = 0; i < slp_count - 1; i++) {
        if (relations[slp_idx[i]].lp == relations[slp_idx[i+1]].lp) {
            int a = slp_idx[i], b = slp_idx[i+1];
            if (rel_count < MAX_RELATIONS) {
                relation_t *r = &relations[rel_count];
                mpz_init(r->qx);
                mpz_mul(r->qx, relations[a].qx, relations[b].qx);
                mpz_init(r->axpb);
                /* axpb for merged = product of both axpb values */
                mpz_mul(r->axpb, relations[a].axpb, relations[b].axpb);
                r->exponents = malloc(fb_count);
                for (int j = 0; j < fb_count; j++)
                    r->exponents[j] = relations[a].exponents[j] ^ relations[b].exponents[j];
                r->lp = 0;
                rel_count++;
                merged_count++;
            }
            i++; /* skip partner */
        }
    }
    free(slp_idx);

    /* Count usable (fully smooth + merged) relations */
    int usable = 0;
    int *usable_idx = malloc(sizeof(int) * rel_count);
    for (int i = 0; i < rel_count; i++) {
        if (relations[i].lp == 0) {
            usable_idx[usable++] = i;
        }
    }

    fprintf(stderr, "Usable: %d (smooth: %d, merged: %d), need: %d\n",
            usable, usable - merged_count, merged_count, fb_count + 1);

    if (usable <= fb_count) {
        fprintf(stderr, "Not enough relations\n");
        goto cleanup;
    }

    /* Linear algebra */
    {
        int rows = usable;
        int cols = fb_count;
        unsigned char **matrix = malloc(sizeof(unsigned char*) * rows);
        unsigned char **history = malloc(sizeof(unsigned char*) * rows);

        for (int i = 0; i < rows; i++) {
            matrix[i] = malloc(cols);
            memcpy(matrix[i], relations[usable_idx[i]].exponents, cols);
        }

        int rank = gaussian_elim(matrix, rows, cols, history);
        fprintf(stderr, "Matrix: %d x %d, rank=%d, null=%d\n", rows, cols, rank, rows - rank);

        int *dep = malloc(sizeof(int) * rows);
        int factored = 0;
        mpz_t factor;
        mpz_init(factor);

        for (int t = rank; t < rows && !factored; t++) {
            int dc = 0;
            for (int i = 0; i < rows; i++) {
                if (history[t][i/8] & (1 << (i%8)))
                    dep[dc++] = usable_idx[i];
            }
            if (dc > 0)
                factored = try_factor(dep, dc, factor);
        }

        if (factored) {
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, factor);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
            fprintf(stderr, "Time: %.3f seconds\n", elapsed);
            mpz_clear(cofactor);
        } else {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "FAILED after %.3f seconds\n", elapsed);
            gmp_printf("FAIL %Zd\n", N);
        }

        mpz_clear(factor);
        free(dep);
        for (int i = 0; i < rows; i++) { free(matrix[i]); free(history[i]); }
        free(matrix); free(history);
    }

cleanup:
    free(usable_idx);
    for (int i = 0; i < batch_alloc; i++) {
        mpz_clear(batch_vals[i]);
        mpz_clear(batch_axpb[i]);
    }
    free(batch_vals); free(batch_axpb);
    gmp_randclear(rstate);
    mpz_clear(N); mpz_clear(kN);
    return 0;
}
