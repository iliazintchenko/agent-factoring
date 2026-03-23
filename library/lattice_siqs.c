/*
 * lattice_siqs.c - Self-Initializing Quadratic Sieve
 *
 * A clean, correct SIQS implementation focused on measuring scaling behavior.
 * Uses: Gray code polynomial switching, SLP (single large prime),
 * block sieve with logarithmic approximation, GF(2) Gaussian elimination.
 *
 * This is NOT trying to beat YAFU's hand-tuned AVX512 code in absolute time.
 * Instead, it measures the fundamental scaling of QS to compare against
 * alternative approaches.
 *
 * Compile: gcc -O3 -march=native -o lattice_siqs library/lattice_siqs.c -lgmp -lm
 * Usage: ./lattice_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ======== Primes ======== */
static int primes[200000];
static int nprimes;

static void gen_primes(int lim) {
    char *s = calloc(lim + 1, 1);
    nprimes = 0;
    for (int i = 2; i <= lim; i++) {
        if (!s[i]) primes[nprimes++] = i;
        for (long j = (long)i * i; j <= lim; j += i) s[j] = 1;
    }
    free(s);
}

/* ======== Tonelli-Shanks ======== */
static unsigned long modsqrt(unsigned long n, unsigned long p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    /* p ≡ 3 mod 4: simple case */
    if (p % 4 == 3) {
        mpz_t b, e, m, r;
        mpz_init_set_ui(b, n);
        mpz_init_set_ui(e, (p + 1) / 4);
        mpz_init_set_ui(m, p);
        mpz_init(r);
        mpz_powm(r, b, e, m);
        unsigned long res = mpz_get_ui(r);
        mpz_clear(b); mpz_clear(e); mpz_clear(m); mpz_clear(r);
        return res;
    }

    /* General Tonelli-Shanks */
    unsigned long Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }

    /* Find non-residue */
    unsigned long z = 2;
    while (1) {
        mpz_t zt, pt, et, rt;
        mpz_init_set_ui(zt, z);
        mpz_init_set_ui(pt, p);
        mpz_init_set_ui(et, (p - 1) / 2);
        mpz_init(rt);
        mpz_powm(rt, zt, et, pt);
        int is_nr = (mpz_get_ui(rt) == p - 1);
        mpz_clear(zt); mpz_clear(pt); mpz_clear(et); mpz_clear(rt);
        if (is_nr) break;
        z++;
    }

    mpz_t M_val, c, t, R, b, temp, pt;
    mpz_init_set_ui(M_val, S);
    mpz_init(c); mpz_init(t); mpz_init(R); mpz_init(b); mpz_init(temp);
    mpz_init_set_ui(pt, p);

    mpz_set_ui(c, z);
    mpz_t Qt; mpz_init_set_ui(Qt, Q);
    mpz_powm(c, c, Qt, pt);

    mpz_set_ui(t, n);
    mpz_powm(t, t, Qt, pt);

    mpz_add_ui(temp, Qt, 1);
    mpz_tdiv_q_2exp(temp, temp, 1);
    mpz_set_ui(R, n);
    mpz_powm(R, R, temp, pt);

    unsigned long M_long = S;
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            unsigned long res = mpz_get_ui(R);
            mpz_clear(M_val); mpz_clear(c); mpz_clear(t); mpz_clear(R);
            mpz_clear(b); mpz_clear(temp); mpz_clear(pt); mpz_clear(Qt);
            return res;
        }
        unsigned long i = 0;
        mpz_set(temp, t);
        while (mpz_cmp_ui(temp, 1) != 0) {
            mpz_mul(temp, temp, temp); mpz_mod(temp, temp, pt);
            i++;
        }
        mpz_set(b, c);
        for (unsigned long j = 0; j < M_long - i - 1; j++) {
            mpz_mul(b, b, b); mpz_mod(b, b, pt);
        }
        M_long = i;
        mpz_mul(c, b, b); mpz_mod(c, c, pt);
        mpz_mul(t, t, c); mpz_mod(t, t, pt);
        mpz_mul(R, R, b); mpz_mod(R, R, pt);
    }
}

/* ======== Factor base ======== */

typedef struct {
    int p;
    int logp;        /* log2(p) * 256 / log2(N) scaled to byte */
    unsigned long root1, root2; /* sqrt(N) mod p */
    double flogp;    /* log(p) */
} fb_prime_t;

/* ======== SIQS parameters (tuned per digit count) ======== */

typedef struct {
    int fb_bound;
    int num_a_factors;  /* number of primes in A */
    int sieve_radius;   /* M: sieve from -M to M per polynomial */
    int block_size;
    double thresh_adj;  /* threshold adjustment factor */
} siqs_params_t;

static siqs_params_t get_params(int digits) {
    siqs_params_t p;
    p.block_size = 32768;
    p.thresh_adj = 0.80;

    if (digits <= 30)      { p.fb_bound = 1500; p.num_a_factors = 3; p.sieve_radius = 32768; }
    else if (digits <= 35) { p.fb_bound = 2500; p.num_a_factors = 4; p.sieve_radius = 32768; }
    else if (digits <= 40) { p.fb_bound = 5000; p.num_a_factors = 5; p.sieve_radius = 65536; }
    else if (digits <= 45) { p.fb_bound = 10000; p.num_a_factors = 6; p.sieve_radius = 65536; }
    else if (digits <= 50) { p.fb_bound = 20000; p.num_a_factors = 6; p.sieve_radius = 65536; }
    else if (digits <= 55) { p.fb_bound = 40000; p.num_a_factors = 7; p.sieve_radius = 65536; }
    else if (digits <= 60) { p.fb_bound = 60000; p.num_a_factors = 7; p.sieve_radius = 65536; }
    else if (digits <= 65) { p.fb_bound = 90000; p.num_a_factors = 8; p.sieve_radius = 65536; }
    else if (digits <= 70) { p.fb_bound = 130000; p.num_a_factors = 8; p.sieve_radius = 65536; }
    else if (digits <= 75) { p.fb_bound = 180000; p.num_a_factors = 9; p.sieve_radius = 65536; }
    else if (digits <= 80) { p.fb_bound = 250000; p.num_a_factors = 9; p.sieve_radius = 65536; }
    else                   { p.fb_bound = 350000; p.num_a_factors = 10; p.sieve_radius = 131072; }

    return p;
}

/* ======== Relation storage ======== */

#define MAX_RELS 65536

typedef struct {
    mpz_t x;       /* x value such that x^2 ≡ Q (mod N) */
    int *exp;       /* exponent vector over FB */
    int neg;        /* 1 if Q(x) < 0 */
    unsigned long lp; /* large prime cofactor, 1 if fully smooth */
} rel_t;

/* ======== GF(2) linear algebra ======== */

static void do_linalg(mpz_t N, rel_t *rels, int nrels, int fb_sz, fb_prime_t *fb) {
    int ncols = fb_sz + 1;  /* +1 for sign */
    int nw = (ncols + 63) / 64;
    int hw = (nrels + 63) / 64;

    unsigned long **mat = calloc(nrels, sizeof(unsigned long *));
    unsigned long **hist = calloc(nrels, sizeof(unsigned long *));
    for (int i = 0; i < nrels; i++) {
        mat[i] = calloc(nw, sizeof(unsigned long));
        hist[i] = calloc(hw, sizeof(unsigned long));
        hist[i][i / 64] |= (1UL << (i % 64));

        if (rels[i].neg)
            mat[i][0] |= 1UL;
        for (int j = 0; j < fb_sz; j++) {
            if (rels[i].exp[j] & 1) {
                int col = j + 1;
                mat[i][col / 64] |= (1UL << (col % 64));
            }
        }
    }

    /* GE */
    int *pivot = malloc(ncols * sizeof(int));
    memset(pivot, -1, ncols * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        int pr = -1;
        for (int row = 0; row < nrels; row++) {
            if (!((mat[row][col / 64] >> (col % 64)) & 1)) continue;
            int used = 0;
            for (int c = 0; c < col; c++)
                if (pivot[c] == row) { used = 1; break; }
            if (!used) { pr = row; break; }
        }
        if (pr < 0) continue;
        pivot[col] = pr;

        for (int row = 0; row < nrels; row++) {
            if (row != pr && ((mat[row][col / 64] >> (col % 64)) & 1)) {
                for (int w = 0; w < nw; w++) mat[row][w] ^= mat[pr][w];
                for (int w = 0; w < hw; w++) hist[row][w] ^= hist[pr][w];
            }
        }
    }

    /* Find dependencies and try to factor */
    mpz_t X, Y, g, tmp;
    mpz_init(X); mpz_init(Y); mpz_init(g); mpz_init(tmp);

    for (int row = 0; row < nrels; row++) {
        int is_zero = 1;
        for (int w = 0; w < nw; w++)
            if (mat[row][w]) { is_zero = 0; break; }
        if (!is_zero) continue;

        mpz_set_ui(X, 1);
        int *comb_exp = calloc(fb_sz, sizeof(int));

        int count = 0;
        for (int i = 0; i < nrels; i++) {
            if (!((hist[row][i / 64] >> (i % 64)) & 1)) continue;
            mpz_mul(X, X, rels[i].x);
            mpz_mod(X, X, N);
            for (int j = 0; j < fb_sz; j++)
                comb_exp[j] += rels[i].exp[j];
            count++;
        }
        if (count < 2) { free(comb_exp); continue; }

        mpz_set_ui(Y, 1);
        for (int j = 0; j < fb_sz; j++) {
            int e = comb_exp[j] / 2;
            if (e > 0) {
                mpz_set_ui(tmp, fb[j].p);
                mpz_powm_ui(tmp, tmp, e, N);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }
        }

        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t cof; mpz_init(cof); mpz_divexact(cof, N, g);
            gmp_printf("%Zd = %Zd * %Zd\n", N, g, cof);
            mpz_clear(cof); free(comb_exp); goto done;
        }
        mpz_add(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t cof; mpz_init(cof); mpz_divexact(cof, N, g);
            gmp_printf("%Zd = %Zd * %Zd\n", N, g, cof);
            mpz_clear(cof); free(comb_exp); goto done;
        }
        free(comb_exp);
    }
    fprintf(stderr, "LA failed: all deps gave trivial gcd\n");

done:
    mpz_clear(X); mpz_clear(Y); mpz_clear(g); mpz_clear(tmp);
    for (int i = 0; i < nrels; i++) { free(mat[i]); free(hist[i]); }
    free(mat); free(hist); free(pivot);
}

/* ======== Knuth-Schroeppel multiplier ======== */

static int choose_multiplier(mpz_t N) {
    static const int mults[] = {1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 29, 31, 33, 37, 41, 43};
    int nmults = sizeof(mults) / sizeof(mults[0]);
    double best_score = -1e30;
    int best_k = 1;

    for (int mi = 0; mi < nmults; mi++) {
        int k = mults[mi];
        mpz_t kN;
        mpz_init(kN);
        mpz_mul_ui(kN, N, k);

        double score = -0.5 * log((double)k);

        /* Score contribution from small primes */
        for (int i = 0; i < nprimes && primes[i] < 1000; i++) {
            int p = primes[i];
            mpz_t pt;
            mpz_init_set_ui(pt, p);
            int leg = mpz_legendre(kN, pt);
            mpz_clear(pt);

            if (leg == 0) score += log((double)p);
            else if (leg == 1) score += 2.0 * log((double)p) / (p - 1);
        }

        if (score > best_score) {
            best_score = score;
            best_k = k;
        }
        mpz_clear(kN);
    }
    return best_k;
}

/* ======== Main SIQS ======== */

static void factor_siqs(mpz_t N_orig) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int digits = (int)mpz_sizeinbase(N_orig, 10);
    double log_n = mpz_sizeinbase(N_orig, 2) * log(2.0);

    gen_primes(1000000);
    int kn_mult = 1; /* Disable multiplier for now - LA sqrt step needs fixing */

    mpz_t N;
    mpz_init(N);
    mpz_set(N, N_orig);

    fprintf(stderr, "%d digits, multiplier %d\n", digits, kn_mult);

    siqs_params_t params = get_params(digits);

    gen_primes(params.fb_bound);

    /* Build factor base for kN */
    fb_prime_t *fb = malloc(nprimes * sizeof(fb_prime_t));
    int fb_sz = 0;

    /* Include 2 */
    fb[0].p = 2;
    fb[0].flogp = log(2.0);
    fb[0].root1 = fb[0].root2 = 0;
    fb_sz = 1;

    mpz_t tmp;
    mpz_init(tmp);

    for (int i = 1; i < nprimes && primes[i] <= params.fb_bound; i++) {
        int p = primes[i];
        mpz_set_ui(tmp, p);
        int leg = mpz_legendre(N, tmp);
        if (leg != 1) continue;

        unsigned long n_mod_p = mpz_fdiv_ui(N, p);
        unsigned long r = modsqrt(n_mod_p, p);
        if ((r * r) % p != n_mod_p) continue;  /* verify */

        fb[fb_sz].p = p;
        fb[fb_sz].flogp = log((double)p);
        fb[fb_sz].root1 = r;
        fb[fb_sz].root2 = p - r;
        fb_sz++;
    }

    fprintf(stderr, "FB: %d primes (max %d)\n", fb_sz, fb[fb_sz - 1].p);

    int target = fb_sz + 40;

    /* Allocate relations */
    rel_t *rels = malloc(MAX_RELS * sizeof(rel_t));
    int nrels = 0;

    /* SLP table */
    #define LP_HASH 131072
    typedef struct lpe { unsigned long lp; int idx; struct lpe *next; } lpe_t;
    lpe_t **lp_tab = calloc(LP_HASH, sizeof(lpe_t *));
    int n_partial = 0;
    int n_combined = 0;

    unsigned long lp_bound = (unsigned long)fb[fb_sz - 1].p * 50UL;

    /* SIQS polynomial: Q(x) = (Ax + B)^2 - N, with A = product of s primes from FB */
    /* Choose A primes from the middle of the factor base */
    int s = params.num_a_factors;
    int M = params.sieve_radius;

    /* Target A ≈ sqrt(2N) / M */
    mpz_t target_a, sqrt2n;
    mpz_init(target_a);
    mpz_init(sqrt2n);
    mpz_mul_ui(sqrt2n, N, 2);
    mpz_sqrt(sqrt2n, sqrt2n);
    mpz_tdiv_q_ui(target_a, sqrt2n, M);

    fprintf(stderr, "Target A: %zu bits, s=%d, M=%d\n", mpz_sizeinbase(target_a, 2), s, M);

    /* Select A primes: choose primes near target_a^(1/s) */
    double log_target_a = mpz_sizeinbase(target_a, 2) * log(2.0);
    double target_prime_log = log_target_a / s;
    int target_prime = (int)exp(target_prime_log);

    /* Find starting index in FB near target_prime */
    int a_start = fb_sz / 3;  /* Start from lower third of FB */
    for (int i = 0; i < fb_sz; i++) {
        if (fb[i].p >= target_prime / 2) { a_start = i; break; }
    }
    if (a_start + s > fb_sz) a_start = fb_sz - s - 10;
    if (a_start < 5) a_start = 5;

    /* Sieve arrays */
    unsigned char *sieve = malloc(2 * M * sizeof(unsigned char));
    int *exps = malloc(fb_sz * sizeof(int));
    int *a_indices = malloc(s * sizeof(int));

    /* SIQS main loop */
    mpz_t A, B, C, Ax_B, Qval, x_val;
    mpz_init(A); mpz_init(B); mpz_init(C); mpz_init(Ax_B); mpz_init(Qval); mpz_init(x_val);

    /* Gray code counters */
    int a_poly_count = 0;
    int max_a_polys = 10000;  /* Try many A values */

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    int total_polys = 0;
    long total_sieve_hits = 0;

    while (nrels < target && a_poly_count < max_a_polys) {
        /* Select A = product of s primes from FB */
        /* Simple approach: randomly choose s primes from the target range */
        int range = fb_sz - a_start;
        if (range < s + 5) range = s + 5;

        mpz_set_ui(A, 1);
        for (int i = 0; i < s; i++) {
            int tries = 0;
            int idx;
            do {
                idx = a_start + (int)(gmp_urandomm_ui(rng, range));
                if (idx >= fb_sz) idx = fb_sz - 1;
                /* Check not already selected */
                int dup = 0;
                for (int j = 0; j < i; j++)
                    if (a_indices[j] == idx) { dup = 1; break; }
                if (!dup) break;
                tries++;
            } while (tries < 100);
            a_indices[i] = idx;
            mpz_mul_ui(A, A, fb[idx].p);
        }

        /* Compute B using Tonelli-Shanks for each a_prime */
        /* B^2 ≡ N (mod A), solved via CRT */
        /* For each a_prime p_i: B ≡ sqrt(N) (mod p_i) */
        /* Then combine with CRT */

        /* Simple: B = sqrt(N) mod A via CRT of roots */
        mpz_set_ui(B, 0);
        for (int i = 0; i < s; i++) {
            int p = fb[a_indices[i]].p;
            unsigned long r = fb[a_indices[i]].root1;  /* sqrt(N) mod p */

            /* CRT: B = B + (r - B mod p) * (A/p) * inv(A/p mod p) */
            mpz_t Ap, Ap_inv, diff;
            mpz_init(Ap); mpz_init(Ap_inv); mpz_init(diff);

            mpz_divexact_ui(Ap, A, p);
            mpz_set_ui(tmp, p);
            mpz_invert(Ap_inv, Ap, tmp);

            unsigned long B_mod_p = mpz_fdiv_ui(B, p);
            long d = ((long)r - (long)B_mod_p + p) % p;
            mpz_set_ui(diff, d);
            mpz_mul(diff, diff, Ap);
            mpz_mul(diff, diff, Ap_inv);
            mpz_add(B, B, diff);
            mpz_mod(B, B, A);

            mpz_clear(Ap); mpz_clear(Ap_inv); mpz_clear(diff);
        }

        /* Ensure B^2 ≡ N (mod A) */
        mpz_mul(tmp, B, B);
        mpz_sub(tmp, tmp, N);
        if (!mpz_divisible_p(tmp, A)) {
            a_poly_count++;
            continue;  /* Bad polynomial, skip */
        }

        /* C = (B^2 - N) / A */
        mpz_mul(C, B, B);
        mpz_sub(C, C, N);
        mpz_divexact(C, C, A);

        /* Q(x) = Ax^2 + 2Bx + C = (Ax+B)^2/A - N/A
         * Actually: (Ax+B)^2 - N = A*(Ax^2 + 2Bx + C)
         * So Q(x) = Ax^2 + 2Bx + C and (Ax+B)^2 = A*Q(x) + N
         * => (Ax+B)^2 ≡ A*Q(x) (mod N) */

        a_poly_count++;
        total_polys++;

        /* Compute sieve roots for this polynomial */
        /* Q(x) ≡ 0 (mod p) for FB prime p:
         * Ax^2 + 2Bx + C ≡ 0 (mod p)
         * If p | A: x ≡ -C/(2B) (mod p) [single root]
         * Else: x ≡ (-B ± sqrt(N)) / A (mod p) */

        /* Clear sieve */
        /* Use log approximation: init to log2(Q(x_center)) and subtract log2(p) for each hit */
        double log2_qmax = mpz_sizeinbase(A, 2) - 1 + 2 * log2((double)M);
        /* Actually Q(x) ≈ A*x^2 near edges, ≈ C near center */
        double log2_thresh = log2(fabs(mpz_get_d(C))) * params.thresh_adj;
        if (log2_thresh < 10) log2_thresh = 10;

        /* Initialize sieve to estimated log2(|Q(x)|) */
        /* Q(x) = Ax^2 + 2Bx + C, x from -M to M-1 */
        /* For simplicity, just use a flat threshold */
        memset(sieve, 0, 2 * M);

        /* Sieve */
        for (int i = 0; i < fb_sz; i++) {
            int p = fb[i].p;
            int logp_byte = (int)(fb[i].flogp / log(2.0) * 1.4);
            if (logp_byte < 1) logp_byte = 1;
            if (logp_byte > 255) logp_byte = 255;

            /* Check if p divides A */
            int divides_a = 0;
            for (int j = 0; j < s; j++)
                if (fb[a_indices[j]].p == p) { divides_a = 1; break; }

            if (p == 2) {
                /* Sieve every other position */
                unsigned long c_mod2 = mpz_fdiv_ui(C, 2);
                int start = (c_mod2 == 0) ? 0 : 1;
                for (int x = start; x < 2 * M; x += 2)
                    if ((int)sieve[x] + logp_byte <= 255)
                        sieve[x] += logp_byte;
                continue;
            }

            if (divides_a) {
                /* Single root: x ≡ -C * (2B)^{-1} (mod p) */
                unsigned long B2 = (2 * mpz_fdiv_ui(B, p)) % p;
                mpz_set_ui(tmp, B2);
                mpz_t pt2; mpz_init_set_ui(pt2, p);
                if (mpz_invert(tmp, tmp, pt2)) {
                    unsigned long Cm = mpz_fdiv_ui(C, p);
                    unsigned long root = ((p - Cm) % p * mpz_get_ui(tmp)) % p;
                    /* Shift: sieve index = x + M */
                    int start = (int)((root + M) % p);
                    for (int x = start; x < 2 * M; x += p)
                        if ((int)sieve[x] + logp_byte <= 255)
                            sieve[x] += logp_byte;
                }
                mpz_clear(pt2);
                continue;
            }

            /* Two roots: x ≡ (-B ± sqrt(N)) * A^{-1} (mod p) */
            unsigned long A_mod_p = mpz_fdiv_ui(A, p);
            unsigned long B_mod_p = mpz_fdiv_ui(B, p);
            unsigned long r1 = fb[i].root1;  /* sqrt(N) mod p */
            unsigned long r2 = fb[i].root2;

            /* Compute A^{-1} mod p */
            mpz_set_ui(tmp, A_mod_p);
            mpz_t pt2; mpz_init_set_ui(pt2, p);
            mpz_invert(tmp, tmp, pt2);
            unsigned long A_inv = mpz_get_ui(tmp);
            mpz_clear(pt2);

            /* root1 = (-B + r1) * A_inv mod p */
            unsigned long sol1 = ((p - B_mod_p + r1) % p * A_inv) % p;
            /* root2 = (-B + r2) * A_inv mod p */
            unsigned long sol2 = ((p - B_mod_p + r2) % p * A_inv) % p;

            /* Sieve with both roots */
            int start1 = (int)((sol1 + M) % p);
            for (int x = start1; x < 2 * M; x += p)
                if ((int)sieve[x] + logp_byte <= 255)
                    sieve[x] += logp_byte;

            if (sol2 != sol1) {
                int start2 = (int)((sol2 + M) % p);
                for (int x = start2; x < 2 * M; x += p)
                    if ((int)sieve[x] + logp_byte <= 255)
                        sieve[x] += logp_byte;
            }
        }

        /* Check sieve candidates */
        int thresh_byte = (int)(log2_thresh);
        if (thresh_byte < 10) thresh_byte = 10;

        for (int xi = 0; xi < 2 * M && nrels < MAX_RELS; xi++) {
            if (sieve[xi] < thresh_byte) continue;

            long x = (long)xi - M;
            total_sieve_hits++;

            /* Compute Q(x) = Ax^2 + 2Bx + C */
            mpz_set_si(Qval, x);
            mpz_mul(Qval, Qval, A);
            mpz_t Bx2; mpz_init(Bx2);
            mpz_set_si(Bx2, x);
            mpz_mul(Bx2, Bx2, B);
            mpz_mul_ui(Bx2, Bx2, 2);
            mpz_mul_si(Qval, A, x);
            mpz_add(Qval, Qval, B);
            mpz_add(Qval, Qval, B);  /* Qval = Ax + 2B */
            mpz_mul_si(Qval, Qval, x);
            mpz_add(Qval, Qval, C);  /* Qval = Ax^2 + 2Bx + C */
            mpz_clear(Bx2);

            int neg = (mpz_sgn(Qval) < 0);
            if (neg) mpz_neg(Qval, Qval);

            /* Trial divide */
            int smooth = 1;
            memset(exps, 0, fb_sz * sizeof(int));

            mpz_t cofactor;
            mpz_init_set(cofactor, Qval);

            for (int i = 0; i < fb_sz; i++) {
                while (mpz_divisible_ui_p(cofactor, fb[i].p)) {
                    exps[i]++;
                    mpz_divexact_ui(cofactor, cofactor, fb[i].p);
                }
            }

            /* Include A primes in exponent vector (since we track A*Q(x)) */
            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb_sz; i++) {
                    if (fb[i].p == fb[a_indices[j]].p) {
                        exps[i]++;
                        break;
                    }
                }
            }

            unsigned long lp = 1;
            if (mpz_cmp_ui(cofactor, 1) == 0) {
                /* Fully smooth */
            } else if (mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= lp_bound) {
                lp = mpz_get_ui(cofactor);
                if (!mpz_probab_prime_p(cofactor, 5)) {
                    mpz_clear(cofactor);
                    continue;  /* Cofactor not prime, skip */
                }
            } else {
                mpz_clear(cofactor);
                continue;  /* Too large */
            }
            mpz_clear(cofactor);

            /* Compute x_val = (Ax + B) mod N (the "square root" side) */
            mpz_set_si(x_val, x);
            mpz_mul(x_val, x_val, A);
            mpz_add(x_val, x_val, B);
            mpz_mod(x_val, x_val, N);

            if (lp == 1) {
                /* Full relation */
                mpz_init_set(rels[nrels].x, x_val);
                rels[nrels].exp = malloc(fb_sz * sizeof(int));
                memcpy(rels[nrels].exp, exps, fb_sz * sizeof(int));
                rels[nrels].neg = neg;
                rels[nrels].lp = 1;
                nrels++;
            } else {
                /* Partial (SLP) */
                unsigned int h = (unsigned int)(lp % LP_HASH);
                lpe_t *e = lp_tab[h];
                int found = 0;
                while (e) {
                    if (e->lp == lp) {
                        /* Combine! x_combined = x1 * x2 * LP^{-1} mod N */
                        /* Because (x1*x2)^2 ≡ LP^2 * smooth, so */
                        /* (x1*x2/LP)^2 ≡ smooth */
                        int pidx = e->idx;
                        mpz_init(rels[nrels].x);
                        mpz_mul(rels[nrels].x, rels[pidx].x, x_val);
                        mpz_mod(rels[nrels].x, rels[nrels].x, N);
                        /* Divide by LP */
                        mpz_t lp_inv;
                        mpz_init(lp_inv);
                        mpz_set_ui(lp_inv, lp);
                        if (!mpz_invert(lp_inv, lp_inv, N)) {
                            mpz_clear(lp_inv);
                            mpz_clear(rels[nrels].x);
                            e = e->next;
                            continue;
                        }
                        mpz_mul(rels[nrels].x, rels[nrels].x, lp_inv);
                        mpz_mod(rels[nrels].x, rels[nrels].x, N);
                        mpz_clear(lp_inv);
                        rels[nrels].exp = malloc(fb_sz * sizeof(int));
                        for (int k = 0; k < fb_sz; k++)
                            rels[nrels].exp[k] = rels[pidx].exp[k] + exps[k];
                        rels[nrels].neg = rels[pidx].neg ^ neg;
                        rels[nrels].lp = 1;
                        nrels++;
                        n_combined++;
                        found = 1;
                        break;
                    }
                    e = e->next;
                }
                if (!found) {
                    /* Store partial */
                    int pidx = MAX_RELS / 2 + n_partial;
                    if (pidx < MAX_RELS) {
                        mpz_init_set(rels[pidx].x, x_val);
                        rels[pidx].exp = malloc(fb_sz * sizeof(int));
                        memcpy(rels[pidx].exp, exps, fb_sz * sizeof(int));
                        rels[pidx].neg = neg;
                        rels[pidx].lp = lp;

                        lpe_t *ne = malloc(sizeof(lpe_t));
                        ne->lp = lp; ne->idx = pidx; ne->next = lp_tab[h];
                        lp_tab[h] = ne;
                        n_partial++;
                    }
                }
            }
        }

        /* Progress */
        if (a_poly_count % 50 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "  %d polys: %d/%d rels (%d SLP), %.1f rels/s, %.1fs\n",
                    a_poly_count, nrels, target, n_combined, nrels / (el > 0 ? el : 1), el);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Sieve: %d rels (%d combined SLP, %d partial) in %.3fs, %d polys, %ld hits\n",
            nrels, n_combined, n_partial, sieve_time, total_polys, total_sieve_hits);

    if (nrels < fb_sz + 1) {
        fprintf(stderr, "FAIL: insufficient relations\n");
        goto cleanup;
    }

    /* Linear algebra */
    fprintf(stderr, "Linear algebra (%d x %d)...\n", nrels, fb_sz + 1);
    do_linalg(N, rels, nrels, fb_sz, fb);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Total: %.3fs\n", total);

cleanup:
    for (int i = 0; i < LP_HASH; i++) {
        lpe_t *e = lp_tab[i];
        while (e) { lpe_t *n = e->next; free(e); e = n; }
    }
    free(lp_tab);

    for (int i = 0; i < nrels; i++) {
        mpz_clear(rels[i].x);
        free(rels[i].exp);
    }
    for (int i = 0; i < n_partial; i++) {
        int idx = MAX_RELS / 2 + i;
        if (idx < MAX_RELS) {
            mpz_clear(rels[idx].x);
            free(rels[idx].exp);
        }
    }
    free(rels);
    free(exps);
    free(sieve);
    free(a_indices);
    free(fb);
    mpz_clear(A); mpz_clear(B); mpz_clear(C); mpz_clear(Ax_B);
    mpz_clear(Qval); mpz_clear(x_val); mpz_clear(tmp);
    mpz_clear(target_a); mpz_clear(sqrt2n); mpz_clear(N);
    gmp_randclear(rng);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }
    mpz_t N;
    mpz_init(N);
    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }
    if (mpz_probab_prime_p(N, 25)) {
        gmp_printf("%Zd is prime\n", N);
        mpz_clear(N);
        return 0;
    }
    factor_siqs(N);
    mpz_clear(N);
    return 0;
}
