/*
 * specialq_qs.c - Quadratic Sieve with Special-Q technique
 *
 * Novel approach: Apply the NFS special-Q lattice sieve idea to QS.
 *
 * Standard QS sieve: For each polynomial, sieve over x in [-M, M],
 * checking all positions. This is O(M * B/log(B)) per polynomial.
 *
 * Special-Q QS: For each "special" prime q in the factor base,
 * only evaluate Q(x) at positions where q | Q(x). These positions
 * form arithmetic progressions x ≡ r (mod q) and x ≡ r' (mod q).
 * Each such Q(x)/q has one less prime factor needed for smoothness.
 *
 * The advantage: we process M/q positions instead of M, but each
 * position is q times more likely to be smooth (it's already
 * divisible by q). Net effect: same number of smooth values
 * found but with less sieve work per polynomial.
 *
 * For large q (say q > B^(1/2)), each position has Q(x)/q with
 * size |Q(x)|/q, which is significantly smaller and thus more
 * likely to be B-smooth.
 *
 * Compile: gcc -O3 -march=native -o specialq_qs library/specialq_qs.c -lgmp -lm
 * Usage: ./specialq_qs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define SEED 42
#define MAX_FB 50000
#define MAX_RELS 200000

typedef struct {
    unsigned int p;
    unsigned int r1, r2; /* Roots of kN mod p */
    double logp;
} fb_entry_t;

typedef struct {
    mpz_t value;      /* Q(x) */
    mpz_t xplusS;     /* x + floor(sqrt(kN)) */
    unsigned char *exponents;
    unsigned int lp;   /* large prime (0 if fully smooth) */
} relation_t;

static mpz_t N, kN, sqrtKN;
static int k_mult;
static fb_entry_t fb[MAX_FB];
static int fb_count;
static relation_t rels[MAX_RELS];
static int rel_count;

/* Tonelli-Shanks */
static int modsqrt_p(unsigned long n, unsigned int p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    mpz_t tmp, mod, base, exp;
    mpz_init_set_ui(tmp, n);
    mpz_init_set_ui(mod, p);
    mpz_init(base);
    mpz_init(exp);

    /* Euler criterion */
    mpz_set_ui(exp, (p - 1) / 2);
    mpz_powm(tmp, tmp, exp, mod);
    if (mpz_cmp_ui(tmp, 1) != 0) {
        mpz_clear(tmp); mpz_clear(mod); mpz_clear(base); mpz_clear(exp);
        return -1;
    }

    if (p % 4 == 3) {
        mpz_set_ui(base, n);
        mpz_set_ui(exp, (p + 1) / 4);
        mpz_powm(tmp, base, exp, mod);
        int r = mpz_get_ui(tmp);
        mpz_clear(tmp); mpz_clear(mod); mpz_clear(base); mpz_clear(exp);
        return r;
    }

    /* Full Tonelli-Shanks */
    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned int z = 2;
    while (1) {
        mpz_set_ui(tmp, z);
        mpz_set_ui(exp, (p - 1) / 2);
        mpz_powm(tmp, tmp, exp, mod);
        if (mpz_cmp_ui(tmp, p - 1) == 0) break;
        z++;
    }

    mpz_t M, c, t, R;
    mpz_init_set_ui(M, S);
    mpz_init(c); mpz_init(t); mpz_init(R);

    mpz_set_ui(c, z);
    mpz_set_ui(exp, Q);
    mpz_powm(c, c, exp, mod);

    mpz_set_ui(t, n);
    mpz_powm(t, t, exp, mod);

    mpz_set_ui(R, n);
    mpz_set_ui(exp, (Q + 1) / 2);
    mpz_powm(R, R, exp, mod);

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            int r = mpz_get_ui(R);
            mpz_clear(M); mpz_clear(c); mpz_clear(t); mpz_clear(R);
            mpz_clear(tmp); mpz_clear(mod); mpz_clear(base); mpz_clear(exp);
            return r;
        }

        int i = 0;
        mpz_set(tmp, t);
        while (mpz_cmp_ui(tmp, 1) != 0) {
            mpz_mul(tmp, tmp, tmp);
            mpz_mod(tmp, tmp, mod);
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

/* Build factor base */
static int build_fb(int smooth_bound) {
    fb_count = 0;
    fb[fb_count].p = 0; /* sign */
    fb[fb_count].logp = 0;
    fb_count++;

    char *sieve = calloc(smooth_bound + 1, 1);
    for (int i = 2; i <= smooth_bound; i++) {
        if (!sieve[i]) {
            unsigned long kn_mod = mpz_fdiv_ui(kN, i);
            int r = modsqrt_p(kn_mod, i);
            if (r >= 0) {
                fb[fb_count].p = i;
                fb[fb_count].r1 = r;
                fb[fb_count].r2 = (i - r) % i;
                fb[fb_count].logp = log2(i);
                fb_count++;
                if (fb_count >= MAX_FB) break;
            }
            for (long j = (long)i * i; j <= smooth_bound; j += i)
                sieve[j] = 1;
        }
    }
    free(sieve);
    return fb_count;
}

/* Select Knuth-Schroeppel multiplier */
static int select_k(mpz_t n) {
    static const int ks[] = {1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 0};
    double best = -1e9;
    int best_k = 1;
    for (int i = 0; ks[i]; i++) {
        mpz_t kn; mpz_init(kn);
        mpz_mul_ui(kn, n, ks[i]);
        double score = -0.5 * log(ks[i]);
        unsigned long mod8 = mpz_fdiv_ui(kn, 8);
        if (mod8 == 1) score += 2 * log(2);
        else if (mod8 == 5) score += log(2);
        else score += 0.5 * log(2);
        if (score > best) { best = score; best_k = ks[i]; }
        mpz_clear(kn);
    }
    return best_k;
}

/* Trial divide Q(x) by factor base, return 1 if smooth, 2 if SLP */
static int trial_divide(mpz_t qx, unsigned char *exponents,
                       unsigned int *lp, unsigned long lp_bound) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_abs(tmp, qx);

    memset(exponents, 0, fb_count);
    *lp = 0;

    if (mpz_sgn(qx) < 0) exponents[0] = 1;

    for (int i = 1; i < fb_count; i++) {
        unsigned int p = fb[i].p;
        while (mpz_divisible_ui_p(tmp, p)) {
            exponents[i] ^= 1;
            mpz_divexact_ui(tmp, tmp, p);
        }
    }

    int result = 0;
    if (mpz_cmp_ui(tmp, 1) == 0) {
        result = 1;
    } else if (mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp) <= lp_bound) {
        *lp = mpz_get_ui(tmp);
        result = 2;
    }

    mpz_clear(tmp);
    return result;
}

/* GF(2) Gaussian elimination */
static int gauss_elim(unsigned char **matrix, int rows, int cols,
                     unsigned char **history) {
    for (int i = 0; i < rows; i++) {
        history[i] = calloc((rows + 7) / 8, 1);
        history[i][i / 8] |= (1 << (i % 8));
    }

    int rank = 0;
    for (int col = 0; col < cols && rank < rows; col++) {
        int piv = -1;
        for (int r = rank; r < rows; r++)
            if (matrix[r][col]) { piv = r; break; }
        if (piv < 0) continue;

        if (piv != rank) {
            unsigned char *t;
            t = matrix[piv]; matrix[piv] = matrix[rank]; matrix[rank] = t;
            t = history[piv]; history[piv] = history[rank]; history[rank] = t;
        }

        for (int r = 0; r < rows; r++) {
            if (r != rank && matrix[r][col]) {
                for (int c = 0; c < cols; c++) matrix[r][c] ^= matrix[rank][c];
                for (int b = 0; b < (rows + 7) / 8; b++) history[r][b] ^= history[rank][b];
            }
        }
        rank++;
    }
    return rank;
}

/* Try to find factor from dependency */
static int try_factor(int *dep, int dc, mpz_t factor) {
    mpz_t x, y, tmp;
    mpz_init_set_ui(x, 1);
    mpz_init_set_ui(y, 1);
    mpz_init(tmp);

    int *total_exp = calloc(fb_count, sizeof(int));

    for (int i = 0; i < dc; i++) {
        relation_t *r = &rels[dep[i]];
        mpz_mul(x, x, r->xplusS);
        mpz_mod(x, x, N);

        mpz_abs(tmp, r->value);
        for (int j = 1; j < fb_count; j++) {
            while (mpz_divisible_ui_p(tmp, fb[j].p)) {
                total_exp[j]++;
                mpz_divexact_ui(tmp, tmp, fb[j].p);
            }
        }
        if (mpz_sgn(r->value) < 0) total_exp[0]++;
    }

    for (int j = 0; j < fb_count; j++) {
        if (total_exp[j] % 2 != 0) {
            free(total_exp);
            mpz_clear(x); mpz_clear(y); mpz_clear(tmp);
            return 0;
        }
    }

    for (int j = 1; j < fb_count; j++) {
        if (total_exp[j] > 0) {
            mpz_set_ui(tmp, fb[j].p);
            mpz_powm_ui(tmp, tmp, total_exp[j] / 2, N);
            mpz_mul(y, y, tmp);
            mpz_mod(y, y, N);
        }
    }
    free(total_exp);

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

/*
 * Special-Q sieve: for a chosen special prime q, only evaluate
 * Q(x) at positions where q | Q(x).
 *
 * Q(x) = (x + sqrtKN)^2 - kN
 * q | Q(x) iff (x + sqrtKN) ≡ ±r (mod q) where r² ≡ kN (mod q)
 * So x ≡ r - sqrtKN (mod q) or x ≡ -r - sqrtKN (mod q)
 *
 * For each such x = x0 + k*q (k = 0, 1, 2, ...),
 * Q(x)/q is smaller by factor q, making it more likely smooth.
 *
 * We also sieve Q(x)/q with the remaining factor base primes
 * (excluding q, which is already accounted for).
 */

static int specialq_collect(int q_idx, int M, unsigned long lp_bound,
                           unsigned char *sieve_buf, unsigned char *exp_buf) {
    unsigned int q = fb[q_idx].p;
    unsigned int r1 = fb[q_idx].r1;
    unsigned int r2 = fb[q_idx].r2;
    unsigned long sqrtKN_mod_q = mpz_fdiv_ui(sqrtKN, q);

    /* Two starting positions */
    int x0_1 = ((int)r1 - (int)sqrtKN_mod_q % q + q) % q;
    int x0_2 = ((int)r2 - (int)sqrtKN_mod_q % q + q) % q;

    int found = 0;
    mpz_t qx, xpS, tmp;
    mpz_init(qx);
    mpz_init(xpS);
    mpz_init(tmp);

    /* For each starting position, iterate x = x0 + k*q */
    for (int start = 0; start < 2; start++) {
        int x0 = (start == 0) ? x0_1 : x0_2;
        if (x0 == 0 && start == 1 && x0_1 == x0_2) continue;

        for (int x = x0; x < M && rel_count < MAX_RELS; x += q) {
            /* Q(x) = (x + sqrtKN)^2 - kN */
            mpz_set(xpS, sqrtKN);
            mpz_add_ui(xpS, xpS, x);
            mpz_mul(qx, xpS, xpS);
            mpz_sub(qx, qx, kN);

            /* Q(x) should be divisible by q */
            if (!mpz_divisible_ui_p(qx, q)) continue;

            /* Divide out q for smaller cofactor */
            /* But keep qx as the full value for the relation */

            /* Trial divide the full Q(x) */
            unsigned int lp = 0;
            int smooth = trial_divide(qx, exp_buf, &lp, lp_bound);

            if (smooth) {
                relation_t *r = &rels[rel_count];
                mpz_init_set(r->value, qx);
                mpz_init_set(r->xplusS, xpS);
                r->exponents = malloc(fb_count);
                memcpy(r->exponents, exp_buf, fb_count);
                r->lp = lp;
                rel_count++;
                found++;
            }
        }
    }

    mpz_clear(qx);
    mpz_clear(xpS);
    mpz_clear(tmp);
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
    mpz_init(sqrtKN);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);
    fprintf(stderr, "Special-Q QS: factoring %d-digit number\n", digits);

    /* Small factor check */
    for (int p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t c; mpz_init(c); mpz_divexact_ui(c, N, p);
            gmp_printf("%Zd = %d * %Zd\n", N, p, c);
            mpz_clear(c); return 0;
        }
    }

    /* Multiplier and FB */
    k_mult = select_k(N);
    mpz_mul_ui(kN, N, k_mult);
    mpz_sqrt(sqrtKN, kN);
    fprintf(stderr, "k = %d\n", k_mult);

    double ln_n = digits * log(10);
    double ln_ln_n = log(ln_n);
    int smooth_bound = (int)(exp(0.5 * sqrt(ln_n * ln_ln_n)) * 0.8);
    if (smooth_bound < 300) smooth_bound = 300;

    build_fb(smooth_bound);
    fprintf(stderr, "FB: %d primes up to %u\n", fb_count, fb[fb_count-1].p);

    unsigned long lp_bound = (unsigned long)fb[fb_count-1].p * 50;
    int target = fb_count + 30;
    int M = (int)(exp(0.5 * sqrt(ln_n * ln_ln_n)) * 10);
    if (M < 100000) M = 100000;

    fprintf(stderr, "M=%d, LP bound=%lu, target=%d\n", M, lp_bound, target);

    /* Sieve buffer and exponent buffer */
    unsigned char *sieve_buf = malloc(M);
    unsigned char *exp_buf = calloc(fb_count, 1);

    rel_count = 0;

    /* Special-Q collection: iterate over special primes */
    /* Use large FB primes as special primes (they give the biggest size reduction) */
    int total_tested = 0;
    int q_start = fb_count / 2; /* Start from middle of FB */

    for (int qi = q_start; qi < fb_count && rel_count < target; qi++) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        if (elapsed > 290.0) break;

        int found = specialq_collect(qi, M, lp_bound, sieve_buf, exp_buf);

        if (qi % 100 == 0 || found > 0) {
            fprintf(stderr, "[%.1fs] q=%u (idx %d/%d): found %d, total %d/%d\n",
                    elapsed, fb[qi].p, qi, fb_count, found, rel_count, target);
        }
    }

    /* Also do a vanilla pass (no special-Q) for x close to 0 */
    if (rel_count < target) {
        fprintf(stderr, "Adding vanilla pass for x near 0...\n");
        mpz_t qx, xpS;
        mpz_init(qx);
        mpz_init(xpS);

        for (int x = 1; x < M && rel_count < target; x++) {
            mpz_set(xpS, sqrtKN);
            mpz_add_ui(xpS, xpS, x);
            mpz_mul(qx, xpS, xpS);
            mpz_sub(qx, qx, kN);

            unsigned int lp = 0;
            int smooth = trial_divide(qx, exp_buf, &lp, lp_bound);
            if (smooth) {
                relation_t *r = &rels[rel_count];
                mpz_init_set(r->value, qx);
                mpz_init_set(r->xplusS, xpS);
                r->exponents = malloc(fb_count);
                memcpy(r->exponents, exp_buf, fb_count);
                r->lp = lp;
                rel_count++;
            }
        }
        mpz_clear(qx);
        mpz_clear(xpS);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Collected %d relations in %.1fs\n", rel_count, elapsed);

    /* SLP merging */
    int slp_count = 0;
    int *slp_idx = malloc(sizeof(int) * rel_count);
    for (int i = 0; i < rel_count; i++)
        if (rels[i].lp > 0) slp_idx[slp_count++] = i;

    /* Sort by LP */
    for (int i = 1; i < slp_count; i++) {
        int k = slp_idx[i];
        unsigned int kv = rels[k].lp;
        int j = i - 1;
        while (j >= 0 && rels[slp_idx[j]].lp > kv) {
            slp_idx[j+1] = slp_idx[j]; j--;
        }
        slp_idx[j+1] = k;
    }

    int merged = 0;
    for (int i = 0; i < slp_count - 1; i++) {
        if (rels[slp_idx[i]].lp == rels[slp_idx[i+1]].lp && rel_count < MAX_RELS) {
            int a = slp_idx[i], b = slp_idx[i+1];
            relation_t *r = &rels[rel_count];
            mpz_init(r->value);
            mpz_mul(r->value, rels[a].value, rels[b].value);
            mpz_init(r->xplusS);
            mpz_mul(r->xplusS, rels[a].xplusS, rels[b].xplusS);
            r->exponents = malloc(fb_count);
            for (int j = 0; j < fb_count; j++)
                r->exponents[j] = rels[a].exponents[j] ^ rels[b].exponents[j];
            r->lp = 0;
            rel_count++;
            merged++;
            i++;
        }
    }
    free(slp_idx);

    /* Count usable */
    int usable = 0;
    int *usable_idx = malloc(sizeof(int) * rel_count);
    for (int i = 0; i < rel_count; i++)
        if (rels[i].lp == 0) usable_idx[usable++] = i;

    fprintf(stderr, "Usable: %d (merged: %d), need: %d\n", usable, merged, fb_count + 1);

    if (usable <= fb_count) {
        fprintf(stderr, "Not enough relations\n");
        gmp_printf("FAIL %Zd\n", N);
        goto cleanup;
    }

    /* Linear algebra */
    {
        int rows = usable, cols = fb_count;
        unsigned char **matrix = malloc(sizeof(unsigned char*) * rows);
        unsigned char **history = malloc(sizeof(unsigned char*) * rows);
        for (int i = 0; i < rows; i++) {
            matrix[i] = malloc(cols);
            memcpy(matrix[i], rels[usable_idx[i]].exponents, cols);
        }

        int rank = gauss_elim(matrix, rows, cols, history);
        fprintf(stderr, "Matrix %d×%d, rank=%d\n", rows, cols, rank);

        int *dep = malloc(sizeof(int) * rows);
        mpz_t factor; mpz_init(factor);
        int factored = 0;

        for (int t = rank; t < rows && !factored; t++) {
            int dc = 0;
            for (int i = 0; i < rows; i++)
                if (history[t][i/8] & (1 << (i%8)))
                    dep[dc++] = usable_idx[i];
            if (dc > 0) factored = try_factor(dep, dc, factor);
        }

        if (factored) {
            mpz_t cof; mpz_init(cof);
            mpz_divexact(cof, N, factor);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cof);
            fprintf(stderr, "Time: %.3f seconds\n", elapsed);
            mpz_clear(cof);
        } else {
            gmp_printf("FAIL %Zd\n", N);
        }

        mpz_clear(factor);
        free(dep);
        for (int i = 0; i < rows; i++) { free(matrix[i]); free(history[i]); }
        free(matrix); free(history);
    }

cleanup:
    free(usable_idx);
    free(sieve_buf);
    free(exp_buf);
    mpz_clear(N); mpz_clear(kN); mpz_clear(sqrtKN);
    return 0;
}
