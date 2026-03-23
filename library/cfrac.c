/*
 * CFRAC - Continued Fraction Factoring Method
 *
 * Uses the continued fraction expansion of sqrt(kN) to find smooth
 * numbers. Each partial quotient gives a value Q_n where
 * P_n^2 ≡ (-1)^n * Q_n (mod N), and Q_n ≤ 2*sqrt(N).
 *
 * Novel feature: batch smooth-part extraction via GCD with primorial,
 * avoiding per-prime trial division.
 *
 * Compile: gcc -O3 -march=native -o cfrac library/cfrac.c -lgmp -lm
 * Usage: ./cfrac <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

#define MAX_FB     100000
#define MAX_RELS   200000
#define MAX_DEPS   128
#define SEED       42

/* ==================== Factor base ==================== */

static int fb_primes[MAX_FB];
static double fb_logp[MAX_FB];
static int fb_size, fb_bound;

static int tonelli_shanks(unsigned long n, int p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    /* Check if QR */
    long long r = 1, base = n, exp = (p - 1) / 2;
    while (exp > 0) {
        if (exp & 1) r = r * base % p;
        exp >>= 1; base = base * base % p;
    }
    if (r != 1) return -1;

    if (p % 4 == 3) {
        r = 1; base = n; exp = (p + 1) / 4;
        while (exp > 0) {
            if (exp & 1) r = r * base % p;
            exp >>= 1; base = base * base % p;
        }
        return (int)r;
    }

    /* Full Tonelli-Shanks for p ≡ 1 (mod 4) */
    int s = 0; long long q = p - 1;
    while (!(q & 1)) { s++; q >>= 1; }
    long long z = 2;
    while (1) { /* find non-residue */
        r = 1; base = z; exp = (p - 1) / 2;
        while (exp > 0) {
            if (exp & 1) r = r * base % p;
            exp >>= 1; base = base * base % p;
        }
        if (r == p - 1) break;
        z++;
    }
    long long m = s;
    long long c = 1; base = z; exp = q;
    while (exp > 0) { if (exp & 1) c = c * base % p; exp >>= 1; base = base * base % p; }
    long long t = 1; base = n; exp = q;
    while (exp > 0) { if (exp & 1) t = t * base % p; exp >>= 1; base = base * base % p; }
    r = 1; base = n; exp = (q + 1) / 2;
    while (exp > 0) { if (exp & 1) r = r * base % p; exp >>= 1; base = base * base % p; }

    while (1) {
        if (t == 1) return (int)r;
        long long i = 0, tmp = t;
        while (tmp != 1) { tmp = tmp * tmp % p; i++; }
        long long b = c;
        for (long long j = 0; j < m - i - 1; j++) b = b * b % p;
        m = i; c = b * b % p; t = t * c % p; r = r * b % p;
    }
}

/* Sieve of Eratosthenes */
static int *gen_primes(int limit, int *count) {
    char *s = calloc(limit + 1, 1);
    for (int i = 2; (long)i * i <= limit; i++)
        if (!s[i]) for (int j = i * i; j <= limit; j += i) s[j] = 1;
    int cnt = 0;
    for (int i = 2; i <= limit; i++) if (!s[i]) cnt++;
    int *p = malloc(cnt * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= limit; i++) if (!s[i]) p[idx++] = i;
    free(s);
    *count = cnt;
    return p;
}

static int select_multiplier(const mpz_t N) {
    static const int cands[] = {1,3,5,7,11,13,17,19,23,29,31,37,41,43};
    double best = -1e30; int best_k = 1;
    mpz_t tmp; mpz_init(tmp);
    for (int ci = 0; ci < 14; ci++) {
        int k = cands[ci];
        mpz_mul_ui(tmp, N, k);
        double sc = -0.5 * log(k);
        int r8 = mpz_fdiv_ui(tmp, 8);
        if (r8 == 1) sc += 2*log(2.0);
        else if (r8 == 5) sc += log(2.0);
        for (int p = 3; p < 200; p += 2) {
            /* Check if p is prime */
            int is_prime = 1;
            for (int d = 2; d * d <= p; d++) if (p % d == 0) { is_prime = 0; break; }
            if (!is_prime) continue;
            if (mpz_kronecker_si(tmp, p) == 1) sc += 2.0*log(p)/(p-1);
        }
        if (sc > best) { best = sc; best_k = k; }
    }
    mpz_clear(tmp);
    return best_k;
}

static void build_fb(const mpz_t kN, int target) {
    int plim = target * 20;
    if (plim < 10000) plim = 10000;
    int nprimes; int *primes = gen_primes(plim, &nprimes);

    fb_size = 0;
    fb_primes[0] = -1; fb_logp[0] = 0; fb_size = 1; /* sign */
    fb_primes[1] = 2; fb_logp[1] = log(2); fb_size = 2;

    for (int i = 0; i < nprimes && fb_size < target; i++) {
        int p = primes[i]; if (p == 2) continue;
        if (mpz_kronecker_si(kN, p) != 1) continue;
        fb_primes[fb_size] = p;
        fb_logp[fb_size] = log(p);
        fb_size++;
    }
    fb_bound = fb_primes[fb_size - 1];
    free(primes);
}

/* ==================== Relations ==================== */

typedef struct {
    mpz_t P;            /* P_n: P^2 ≡ (-1)^n * Q (mod N) */
    int *exponents;     /* exponent vector [sign, p0, p1, ...] */
    int large_prime;    /* 0 for full, >0 for SLP */
} rel_t;

static rel_t rels[MAX_RELS];
static int nrels = 0;

/* LP hash for SLP combining */
#define LP_HASH_SZ (1 << 18)
static int lp_hash[LP_HASH_SZ];
static int lp_next[MAX_RELS];
static int ncombined = 0;

static void init_lp(void) { memset(lp_hash, -1, sizeof(lp_hash)); }

/* ==================== Primorial for batch GCD ==================== */
static mpz_t primorial;

static void compute_primorial(void) {
    mpz_init_set_ui(primorial, 1);
    for (int i = 1; i < fb_size; i++) { /* skip sign at index 0 */
        int p = fb_primes[i];
        long long pk = p;
        while (pk <= (1LL << 28)) {
            mpz_mul_ui(primorial, primorial, p);
            if (pk > (1LL << 28) / p) break;
            pk *= p;
        }
    }
}

/* Extract smooth part via repeated GCD with primorial */
/* Returns: 1 if fully smooth or SLP, 0 otherwise */
/* Sets exponents[] and *lp */
static int check_smooth(int *exponents, int *lp,
                        const mpz_t Q_val, int sign, long lp_bound) {
    memset(exponents, 0, fb_size * sizeof(int));
    exponents[0] = sign; /* sign bit */

    mpz_t cofactor, g, smooth;
    mpz_init(cofactor); mpz_init(g); mpz_init(smooth);
    mpz_abs(cofactor, Q_val);

    /* Extract smooth part via primorial GCD */
    for (int iter = 0; iter < 30; iter++) {
        mpz_gcd(g, cofactor, primorial);
        if (mpz_cmp_ui(g, 1) == 0) break;
        mpz_divexact(cofactor, cofactor, g);
    }

    *lp = 0;
    int ok = 0;

    if (mpz_cmp_ui(cofactor, 1) == 0) {
        ok = 1; /* fully smooth */
    } else if (mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= (unsigned long)lp_bound) {
        if (mpz_probab_prime_p(cofactor, 2)) {
            *lp = mpz_get_ui(cofactor);
            ok = 1; /* SLP */
        }
    }

    if (ok) {
        /* Trial divide to get exponents */
        mpz_abs(smooth, Q_val);
        for (int i = 1; i < fb_size; i++) {
            int p = fb_primes[i];
            while (mpz_divisible_ui_p(smooth, p)) {
                mpz_divexact_ui(smooth, smooth, p);
                exponents[i]++;
            }
            if (mpz_cmp_ui(smooth, 1) == 0) break;
        }
    }

    mpz_clear(cofactor); mpz_clear(g); mpz_clear(smooth);
    return ok;
}

static void add_relation(const mpz_t P, const int *exp, int lp_val) {
    if (nrels >= MAX_RELS) return;
    rel_t *r = &rels[nrels];
    mpz_init_set(r->P, P);
    r->exponents = malloc(fb_size * sizeof(int));
    memcpy(r->exponents, exp, fb_size * sizeof(int));
    r->large_prime = lp_val;

    if (lp_val > 1) {
        int h = lp_val % LP_HASH_SZ;
        int idx = lp_hash[h];
        while (idx >= 0) {
            if (rels[idx].large_prime == lp_val) { ncombined++; break; }
            idx = lp_next[idx];
        }
        lp_next[nrels] = lp_hash[h];
        lp_hash[h] = nrels;
    }
    nrels++;
}

/* ==================== CFRAC: continued fraction expansion ==================== */

/*
 * Expand the continued fraction of sqrt(kN).
 * At step n:
 *   P_n = a_n * P_{n-1} + P_{n-2}  (mod N)
 *   Q_n satisfies P_n^2 ≡ (-1)^n * Q_n (mod kN)
 *   Q_n ≤ 2 * sqrt(kN)
 *
 * Morrison-Brillhart (1975):
 *   g_{-1} = 0, g_0 = floor(sqrt(kN))
 *   q_0 = g_0^2 - kN (always negative)
 *   Recurrence:
 *     g_{n+1} = floor((sqrt(kN) + P_n) / |Q_n|)
 *     P_{n+1} = g_{n+1} * |Q_n| - P_n
 *     Q_{n+1} = Q_n + g_{n+1} * (P_n - P_{n+1})   [actually: Q_{n-1} + g_n*(P_{n-1} - P_n)]
 *
 * Standard recurrence:
 *   P_0 = floor(sqrt(kN)), Q_0 = 1, A_{-1} = 1, A_0 = P_0
 *   Then: g_n = floor((P_0 + P_n) / Q_n)
 *          P_{n+1} = g_n * Q_n - P_n
 *          Q_{n+1} = (kN - P_{n+1}^2) / Q_n  [exact division]
 *          A_n = g_n * A_{n-1} + A_{n-2}  (mod N)
 *   And: A_n^2 ≡ (-1)^{n+1} * Q_{n+1} (mod N)
 */

static void cfrac_factor(const mpz_t N) {
    int digits = mpz_sizeinbase(N, 10);
    int bits = mpz_sizeinbase(N, 2);
    fprintf(stderr, "CFRAC: %d digits (%d bits)\n", digits, bits);

    /* Parameters */
    int target_fb;
    if (digits <= 30) target_fb = 100;
    else if (digits <= 35) target_fb = 150;
    else if (digits <= 40) target_fb = 300;
    else if (digits <= 45) target_fb = 600;
    else if (digits <= 50) target_fb = 1200;
    else if (digits <= 55) target_fb = 2200;
    else if (digits <= 60) target_fb = 4000;
    else if (digits <= 65) target_fb = 7000;
    else if (digits <= 70) target_fb = 12000;
    else target_fb = 20000;

    int mult = select_multiplier(N);
    mpz_t kN;
    mpz_init(kN);
    mpz_mul_ui(kN, N, mult);
    fprintf(stderr, "Multiplier: %d\n", mult);

    build_fb(kN, target_fb);
    fprintf(stderr, "FB: %d primes, bound=%d\n", fb_size, fb_bound);

    long lp_bound = 0; /* full relations only for now */

    compute_primorial();
    fprintf(stderr, "Primorial: %zu bits\n", mpz_sizeinbase(primorial, 2));

    init_lp();

    /* CF expansion */
    mpz_t sqrtKN, P_prev, P_curr, Q_prev, Q_curr, A_prev, A_curr;
    mpz_t g, temp, Qval;
    mpz_init(sqrtKN); mpz_init(P_prev); mpz_init(P_curr);
    mpz_init(Q_prev); mpz_init(Q_curr);
    mpz_init(A_prev); mpz_init(A_curr);
    mpz_init(g); mpz_init(temp); mpz_init(Qval);

    mpz_sqrt(sqrtKN, kN);

    /* Initial values */
    mpz_set(P_curr, sqrtKN);          /* P_0 = floor(sqrt(kN)) */
    mpz_set_ui(Q_prev, 1);            /* Q_{-1} = 1 */
    /* Q_0 = kN - P_0^2 */
    mpz_mul(temp, P_curr, P_curr);
    mpz_sub(Q_curr, kN, temp);
    if (mpz_sgn(Q_curr) == 0) {
        /* kN is a perfect square - trivial factoring */
        fprintf(stderr, "kN is a perfect square!\n");
        mpz_clear(kN); return;
    }

    mpz_set_ui(A_prev, 1);            /* A_{-1} = 1 */
    mpz_set(A_curr, sqrtKN);          /* A_0 = floor(sqrt(kN)) */

    int target_rels = fb_size + 20; /* full rels only, just need fb+small excess */
    int *exponents = malloc(fb_size * sizeof(int));

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    long step = 0;
    int sign = 1; /* (-1)^(n+1): starts at 1 for n=0 since Q_0 = kN - P_0^2 > 0 if kN not perfect square */

    /* Check Q_0 */
    int lp_val;
    if (check_smooth(exponents, &lp_val, Q_curr, 0, lp_bound)) {
        mpz_mod(temp, A_curr, N);
        add_relation(temp, exponents, lp_val);
    }
    step++;

    while (nrels + ncombined < target_rels) {
        /* g = floor((sqrtKN + P_curr) / Q_curr) */
        mpz_add(temp, sqrtKN, P_curr);
        mpz_tdiv_q(g, temp, Q_curr);

        /* P_next = g * Q_curr - P_curr */
        mpz_t P_next;
        mpz_init(P_next);
        mpz_mul(P_next, g, Q_curr);
        mpz_sub(P_next, P_next, P_curr);

        /* Q_next = Q_prev + g * (P_curr - P_next) */
        mpz_t Q_next;
        mpz_init(Q_next);
        mpz_sub(temp, P_curr, P_next);
        mpz_mul(temp, g, temp);
        mpz_add(Q_next, Q_prev, temp);

        /* Alternatively: Q_next = (kN - P_next^2) / Q_curr (exact) */
        /* This is more numerically stable */
        mpz_mul(temp, P_next, P_next);
        mpz_sub(temp, kN, temp);
        mpz_divexact(Q_next, temp, Q_curr);

        /* A_next = g * A_curr + A_prev (mod N) */
        mpz_t A_next;
        mpz_init(A_next);
        mpz_mul(A_next, g, A_curr);
        mpz_add(A_next, A_next, A_prev);
        mpz_mod(A_next, A_next, N);

        sign = 1 - sign; /* alternating sign */

        /* Check if Q_next is smooth */
        /* We have A_next^2 ≡ (-1)^(step) * Q_next (mod kN)
         * So A_next^2 ≡ (-1)^(step) * Q_next (mod N) since N | kN/mult
         * Actually A_next^2 ≡ (-1)^(step) * Q_next (mod kN), so also (mod N)
         * since kN = mult * N and gcd(mult, N) = 1
         */
        if (check_smooth(exponents, &lp_val, Q_next, sign, lp_bound)) {
            add_relation(A_next, exponents, lp_val);
        }

        /* Shift */
        mpz_set(P_prev, P_curr);
        mpz_set(P_curr, P_next);
        mpz_set(Q_prev, Q_curr);
        mpz_set(Q_curr, Q_next);
        mpz_set(A_prev, A_curr);
        mpz_set(A_curr, A_next);

        mpz_clear(P_next); mpz_clear(Q_next); mpz_clear(A_next);

        step++;

        if (step % 50000 == 0) {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double el = (now.tv_sec - t0.tv_sec) + (now.tv_nsec - t0.tv_nsec) / 1e9;
            int full = 0, partial = 0;
            for (int i = 0; i < nrels; i++) {
                if (rels[i].large_prime <= 1) full++; else partial++;
            }
            fprintf(stderr, "\r[%.1fs] step=%ld rels=%d(f=%d p=%d c=%d) "
                    "need=%d rate=%.0f steps/s  ",
                    el, step, nrels, full, partial, ncombined,
                    target_rels, step / el);
        }

        if (step > 100000000) {
            fprintf(stderr, "\nGiving up after %ld steps\n", step);
            break;
        }
    }

    struct timespec t1;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "\nCFRAC: %d rels in %.2fs (%ld steps, %.0f steps/s)\n",
            nrels, sieve_time, step, step / sieve_time);

    free(exponents);

    /* ==================== LA: GF(2) Gaussian Elimination ==================== */
    /* Build matrix including LP columns */
    int *distinct_lps = malloc(nrels * sizeof(int));
    int n_dlps = 0;
    for (int i = 0; i < nrels; i++) {
        if (rels[i].large_prime <= 1) continue;
        int lp = rels[i].large_prime;
        int found = 0;
        for (int j = 0; j < n_dlps; j++)
            if (distinct_lps[j] == lp) { found = 1; break; }
        if (!found) {
            int cnt = 0;
            for (int j = i; j < nrels; j++)
                if (rels[j].large_prime == lp) { cnt++; if (cnt >= 2) break; }
            if (cnt >= 2) distinct_lps[n_dlps++] = lp;
        }
    }

    /* Filter usable relations */
    int *usable = calloc(nrels, sizeof(int));
    int n_usable = 0;
    for (int i = 0; i < nrels; i++) {
        if (rels[i].large_prime <= 1) { usable[i] = 1; n_usable++; }
        else {
            for (int j = 0; j < n_dlps; j++)
                if (rels[i].large_prime == distinct_lps[j]) {
                    usable[i] = 1; n_usable++; break;
                }
        }
    }

    int ncols = fb_size + n_dlps;
    fprintf(stderr, "LA: %d rels, %d cols (fb=%d + %d LPs)\n",
            n_usable, ncols, fb_size, n_dlps);

    /* Build matrix */
    int rw = (ncols + 63) / 64;
    int aw = (n_usable + 63) / 64;
    uint64_t **mat = malloc(n_usable * sizeof(uint64_t *));
    uint64_t **aug = malloc(n_usable * sizeof(uint64_t *));
    int *rel_idx = malloc(n_usable * sizeof(int)); /* maps row -> original rel index */

    int row = 0;
    for (int i = 0; i < nrels; i++) {
        if (!usable[i]) continue;
        mat[row] = calloc(rw, sizeof(uint64_t));
        aug[row] = calloc(aw, sizeof(uint64_t));
        aug[row][row / 64] |= (1ULL << (row % 64));

        /* Set FB exponent parities */
        for (int j = 0; j < fb_size; j++) {
            if (rels[i].exponents[j] & 1) {
                mat[row][j / 64] |= (1ULL << (j % 64));
            }
        }
        /* Set LP column */
        if (rels[i].large_prime > 1) {
            for (int j = 0; j < n_dlps; j++) {
                if (distinct_lps[j] == rels[i].large_prime) {
                    int col = fb_size + j;
                    mat[row][col / 64] |= (1ULL << (col % 64));
                    break;
                }
            }
        }
        rel_idx[row] = i;
        row++;
    }

    /* Gaussian elimination */
    int *pcol = malloc(n_usable * sizeof(int));
    memset(pcol, -1, n_usable * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        int piv = -1;
        for (int r = 0; r < n_usable; r++) {
            if (pcol[r] >= 0) continue;
            if ((mat[r][col / 64] >> (col % 64)) & 1) { piv = r; break; }
        }
        if (piv < 0) continue;
        pcol[piv] = col;
        for (int r = 0; r < n_usable; r++) {
            if (r == piv) continue;
            if (!((mat[r][col / 64] >> (col % 64)) & 1)) continue;
            for (int w = 0; w < rw; w++) mat[r][w] ^= mat[piv][w];
            for (int w = 0; w < aw; w++) aug[r][w] ^= aug[piv][w];
        }
    }

    /* Extract and try dependencies */
    mpz_t X, Y, factor;
    mpz_init(X); mpz_init(Y); mpz_init(factor);

    int factored = 0;
    for (int r = 0; r < n_usable && !factored; r++) {
        /* Check if zero row */
        int zero = 1;
        for (int w = 0; w < rw && zero; w++) if (mat[r][w]) zero = 0;
        if (!zero) continue;

        /* Extract dependency */
        int *total_exp = calloc(fb_size, sizeof(int));
        mpz_set_ui(X, 1);

        int dep_cnt = 0;
        int *lps = malloc(n_usable * sizeof(int));
        int nlps = 0;

        for (int i = 0; i < n_usable; i++) {
            if (!((aug[r][i / 64] >> (i % 64)) & 1)) continue;
            dep_cnt++;
            int ri = rel_idx[i];

            mpz_mul(X, X, rels[ri].P);
            mpz_mod(X, X, N);

            for (int j = 0; j < fb_size; j++)
                total_exp[j] += rels[ri].exponents[j];

            if (rels[ri].large_prime > 1)
                lps[nlps++] = rels[ri].large_prime;
        }

        if (dep_cnt < 2) { free(total_exp); free(lps); continue; }

        /* Check all exponents even */
        int ok = 1;
        for (int j = 0; j < fb_size && ok; j++)
            if (total_exp[j] % 2) ok = 0;

        /* Check LP exponents */
        if (ok) {
            for (int i = 0; i < nlps - 1; i++)
                for (int j = i + 1; j < nlps; j++)
                    if (lps[i] > lps[j]) { int t = lps[i]; lps[i] = lps[j]; lps[j] = t; }
            for (int i = 0; i < nlps && ok; ) {
                int cnt = 1;
                while (i + cnt < nlps && lps[i + cnt] == lps[i]) cnt++;
                if (cnt % 2) { ok = 0; break; }
                i += cnt;
            }
        }

        if (!ok) { free(total_exp); free(lps); continue; }

        /* Compute Y = sqrt(product) mod N */
        mpz_set_ui(Y, 1);
        for (int j = 1; j < fb_size; j++) { /* skip sign at j=0 */
            if (total_exp[j] == 0) continue;
            mpz_set_ui(temp, fb_primes[j]);
            mpz_powm_ui(temp, temp, total_exp[j] / 2, N);
            mpz_mul(Y, Y, temp);
            mpz_mod(Y, Y, N);
        }
        /* Include LP contribution */
        for (int i = 0; i < nlps; ) {
            int cnt = 1;
            while (i + cnt < nlps && lps[i + cnt] == lps[i]) cnt++;
            mpz_set_ui(temp, lps[i]);
            mpz_powm_ui(temp, temp, cnt / 2, N);
            mpz_mul(Y, Y, temp);
            mpz_mod(Y, Y, N);
            i += cnt;
        }

        /* Handle sign: if total sign exponent is odd, multiply by sqrt(-1) mod N
         * Actually, total_exp[0] should be even (dependency guarantees it) */

        /* Try gcd(X ± Y, N) */
        mpz_sub(temp, X, Y);
        mpz_gcd(factor, temp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            factored = 1;
        } else {
            mpz_add(temp, X, Y);
            mpz_gcd(factor, temp, N);
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0)
                factored = 1;
        }

        free(total_exp); free(lps);
    }

    if (factored) {
        struct timespec t2;
        clock_gettime(CLOCK_MONOTONIC, &t2);
        double total = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec) / 1e9;

        mpz_t cofn; mpz_init(cofn);
        mpz_divexact(cofn, N, factor);
        if (mpz_cmp(factor, cofn) > 0) mpz_swap(factor, cofn);

        gmp_printf("Factor: %Zd\n", factor);
        gmp_printf("Cofactor: %Zd\n", cofn);
        fprintf(stderr, "Time: %.3fs (cfrac %.3fs, LA %.3fs)\n",
                total, sieve_time, total - sieve_time);
        mpz_clear(cofn);
    } else {
        fprintf(stderr, "FAILED: %d rels, fb=%d\n", nrels, fb_size);
    }

    /* Cleanup */
    for (int i = 0; i < n_usable; i++) { free(mat[i]); free(aug[i]); }
    free(mat); free(aug); free(rel_idx); free(pcol);
    free(distinct_lps); free(usable);
    mpz_clear(X); mpz_clear(Y); mpz_clear(factor);
    mpz_clear(kN); mpz_clear(sqrtKN);
    mpz_clear(P_prev); mpz_clear(P_curr);
    mpz_clear(Q_prev); mpz_clear(Q_curr);
    mpz_clear(A_prev); mpz_clear(A_curr);
    mpz_clear(g); mpz_clear(temp); mpz_clear(Qval);
    mpz_clear(primorial);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    /* Quick trial division */
    int nprimes; int *primes = gen_primes(1000000, &nprimes);
    for (int i = 0; i < nprimes; i++) {
        if (mpz_divisible_ui_p(N, primes[i])) {
            mpz_t cofn; mpz_init(cofn);
            mpz_divexact_ui(cofn, N, primes[i]);
            gmp_printf("Factor: %d\nCofactor: %Zd\n", primes[i], cofn);
            mpz_clear(cofn); free(primes); mpz_clear(N);
            return 0;
        }
    }
    free(primes);

    cfrac_factor(N);
    mpz_clear(N);
    return 0;
}
