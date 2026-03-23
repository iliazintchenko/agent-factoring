/*
 * balanced_nfs.c - NFS polynomial selection optimized for balanced semiprimes
 *
 * For N = p*q with p ≈ q ≈ √N, we can construct better NFS polynomials:
 *
 * Standard GNFS for degree d: choose m ≈ N^{1/(d+1)}, f(x) has coefficients ≈ N^{1/(d+1)}
 * Balanced semiprime trick: since p ≈ √N, we can use degree-2 polynomial
 * f(x) = x^2 + a*x + b where f(p) ≡ 0 (mod p), giving:
 *   p^2 + a*p + b ≡ 0 (mod N)  =>  b = -(p^2 + a*p) mod N
 *
 * But we don't know p! However, we know p ≈ √N. So we search near √N.
 *
 * Novel approach: "Near-SNFS" for balanced semiprimes
 * Choose m = ⌊√N⌋. Then N = m^2 + r where r = N - m^2 ≈ m (since p,q differ by O(m)).
 * This gives us: x^2 - N has root m, so f(x) = x^2 - N, g(x) = x - m.
 * ||f|| = N (huge), but f(m) = m^2 - N = -r (small: |r| ≈ m ≈ √N).
 *
 * Better: use degree 1 on rational side, degree 2 on algebraic side:
 *   f(x) = x^2 + c1*x + c0  (algebraic)
 *   g(x) = x - m              (rational)
 * where f(m) ≡ 0 (mod N). So c0 = N - m^2 - c1*m.
 * Choose c1 to minimize ||f||.
 *
 * Norms: for (a,b) coprime pair:
 *   F(a,b) = b^2 * f(a/b) = a^2 + c1*a*b + c0*b^2
 *   G(a,b) = a - m*b
 *
 * F(a,b) ≈ a^2 + c1*a*b + (N-m^2)*b^2
 * If c1 ≈ 0 and |a|,|b| ≤ A: |F| ≈ A^2 + N*A^2/m^2... no, that's wrong.
 *
 * Actually: c0 = N - m^2 - c1*m. If m ≈ √N, then N - m^2 ≈ (p+q-2m)*m + (p-m)*(q-m).
 * Since p+q ≈ 2m (balanced), N - m^2 ≈ small * m + small^2. Hmm, N - m^2 can be up to m.
 *
 * For a random balanced semiprime: |N - m^2| ≈ m = √N. So c0 ≈ √N.
 * With c1 = 0: f(x) = x^2 + c0, ||f||_∞ = max(1, c0) = √N.
 * Norms: |F(a,b)| ≈ A^2 + √N * A^2 = (1 + √N) * A^2 ≈ √N * A^2.
 *
 * Compare with standard degree-4 GNFS: norms ≈ N^{1/5} * A^4.
 * For our degree-2: norms ≈ N^{1/2} * A^2.
 *
 * The degree-2 norms are LARGER: N^{1/2} * A^2 > N^{1/5} * A^4 when
 * A < N^{3/10}. For A ≈ N^{2/(d+2)} = N^{1/2} (degree 2): norms ≈ N^{3/2}. Bad.
 *
 * Conclusion: degree-2 NFS for balanced semiprimes is WORSE than degree-4 GNFS.
 * The "balanced" structure doesn't help because the polynomial coefficients
 * are still O(√N), which is huge.
 *
 * BUT: what if we use a HIGHER degree polynomial that exploits the balanced structure?
 *
 * Alternative: "Centered polynomial" trick
 * Let m = ⌊√N⌋. Write N = m^2 + r, |r| < 2m+1.
 * Then (m+k)^2 = m^2 + 2mk + k^2 = N + 2mk + k^2 - r.
 * So (m+k)^2 ≡ 2mk + k^2 - r (mod N).
 * The right side is O(mk) ≈ O(√N * k) which is much smaller than N for small k.
 *
 * For k ≤ K: RHS ≈ 2m*K ≈ 2√N * K.
 * Smoothness probability: u = log(2√N * K) / log(B).
 *
 * This is essentially Fermat's method: try x = m+1, m+2, ... and check if x^2 - N is
 * a perfect square. But instead of checking for perfect square, check for smoothness.
 *
 * "Fermat sieve": sieve x^2 - N for x near √N, using standard QS-style sieve.
 * The values x^2 - N ≈ 2√N * (x - m) are smaller than QS values when x - m is small.
 *
 * This is actually a well-known variant called "Lehman's method" or "Fermat sieve".
 * It's a QS with a specific polynomial f(x) = (m+x)^2 - N = 2mx + x^2 + (m^2-N).
 *
 * Let's implement this and measure scaling!
 *
 * Compile: gcc -O3 -march=native -o balanced_nfs library/balanced_nfs.c -lgmp -lm
 * Usage: ./balanced_nfs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

#define MAX_FB     100000
#define MAX_RELS   300000
#define SEED       42

/* Factor base */
static int fb_primes[MAX_FB];
static int fb_size, fb_bound;

static int tonelli_shanks_fast(long long n, int p) {
    n = ((n % p) + p) % p;
    if (n == 0) return 0;
    if (p == 2) return n & 1;

    /* Check QR */
    long long r = 1, base = n, exp = (p-1)/2;
    while (exp > 0) {
        if (exp & 1) r = (__int128)r * base % p;
        exp >>= 1; base = (__int128)base * base % p;
    }
    if (r != 1) return -1;

    if (p % 4 == 3) {
        r = 1; base = n; exp = (p+1)/4;
        while (exp > 0) {
            if (exp & 1) r = (__int128)r * base % p;
            exp >>= 1; base = (__int128)base * base % p;
        }
        return (int)r;
    }

    /* Full T-S */
    int s = 0; long long q = p-1;
    while (!(q & 1)) { s++; q >>= 1; }
    long long z = 2;
    while (1) {
        r = 1; base = z; exp = (p-1)/2;
        while (exp > 0) {
            if (exp & 1) r = (__int128)r * base % p;
            exp >>= 1; base = (__int128)base * base % p;
        }
        if (r == p-1) break; z++;
    }
    long long m = s;
    long long c, t;
    r = 1; base = z; exp = q;
    while (exp > 0) { if (exp & 1) r = (__int128)r * base % p; exp >>= 1; base = (__int128)base * base % p; }
    c = r;
    r = 1; base = n; exp = q;
    while (exp > 0) { if (exp & 1) r = (__int128)r * base % p; exp >>= 1; base = (__int128)base * base % p; }
    t = r;
    r = 1; base = n; exp = (q+1)/2;
    while (exp > 0) { if (exp & 1) r = (__int128)r * base % p; exp >>= 1; base = (__int128)base * base % p; }

    while (1) {
        if (t == 1) return (int)r;
        long long i = 0, tmp = t;
        while (tmp != 1) { tmp = (__int128)tmp * tmp % p; i++; }
        long long b = c;
        for (long long j = 0; j < m-i-1; j++) b = (__int128)b * b % p;
        m = i; c = (__int128)b * b % p;
        t = (__int128)t * c % p; r = (__int128)r * b % p;
    }
}

static void build_fb(const mpz_t N, int target) {
    /* Sieve of Eratosthenes */
    int plim = target * 20;
    if (plim < 50000) plim = 50000;
    char *sieve = calloc(plim + 1, 1);
    for (int i = 2; (long)i*i <= plim; i++)
        if (!sieve[i]) for (int j = i*i; j <= plim; j += i) sieve[j] = 1;

    fb_size = 0;
    fb_primes[fb_size++] = -1; /* sign */
    fb_primes[fb_size++] = 2;

    for (int p = 3; p <= plim && fb_size < target; p += 2) {
        if (sieve[p]) continue;
        if (mpz_kronecker_si(N, p) != 1) continue;
        fb_primes[fb_size++] = p;
    }
    fb_bound = fb_primes[fb_size - 1];
    free(sieve);
}

/* ==================== Fermat Sieve ==================== */
/*
 * Polynomial: Q(x) = (m+x)^2 - N where m = floor(sqrt(N))
 * Q(x) = 2mx + x^2 + (m^2 - N)
 * Let r = m^2 - N (negative since m < sqrt(N))
 * Q(x) = x^2 + 2mx + r
 *
 * For small x: Q(x) ≈ 2mx + r ≈ 2m*x (linear in x!)
 * Size: |Q(x)| ≈ 2m*|x| for small x
 *
 * Compare with SIQS: |Q(x)| ≈ M * sqrt(N) ≈ M * m
 * Fermat sieve: |Q(x)| ≈ 2m * |x|
 *
 * For the same sieve interval M, values are similar. But the Fermat sieve
 * has the advantage that values grow LINEARLY (not quadratically) near x=0.
 *
 * Also: Q(x) = (m+x)^2 - N, so if Q(x) is smooth, we have
 * (m+x)^2 ≡ Q(x) (mod N) directly!
 * No need for self-initializing polynomials or A-factor tracking.
 */

typedef struct {
    mpz_t Y;          /* m + x */
    int *exponents;    /* exponents[0] = sign, exponents[1..] = FB primes */
    int large_prime;
} rel_t;

static rel_t rels[MAX_RELS];
static int nrels = 0;

/* LP hash */
#define LP_HASH_SZ (1 << 20)
static int lp_hash[LP_HASH_SZ];
static int lp_next[MAX_RELS];
static int ncombined = 0;

static void init_lp(void) { memset(lp_hash, -1, sizeof(lp_hash)); }

static void add_rel(const mpz_t Y, const int *exp, int lp) {
    if (nrels >= MAX_RELS) return;
    rel_t *r = &rels[nrels];
    mpz_init_set(r->Y, Y);
    r->exponents = malloc(fb_size * sizeof(int));
    memcpy(r->exponents, exp, fb_size * sizeof(int));
    r->large_prime = lp;

    if (lp > 1) {
        int h = lp % LP_HASH_SZ;
        int idx = lp_hash[h];
        while (idx >= 0) {
            if (rels[idx].large_prime == lp) { ncombined++; break; }
            idx = lp_next[idx];
        }
        lp_next[nrels] = lp_hash[h];
        lp_hash[h] = nrels;
    }
    nrels++;
}

/* Sieve Q(x) = x^2 + 2mx + r for smoothness */
static int sieve_block(unsigned char *sieve_arr, long block_start, int block_size,
                       const mpz_t m, const mpz_t r_val, const mpz_t N,
                       int *root1, int *root2,
                       long lp_bound, unsigned char thresh) {
    int found = 0;
    memset(sieve_arr, 0, block_size);

    /* Sieve: for each FB prime p, Q(x) ≡ 0 mod p when
     * x^2 + 2mx + r ≡ 0 mod p
     * x ≡ -m ± sqrt(m^2 - r) ≡ -m ± sqrt(N) mod p
     * So roots are: x ≡ (-m + sqrt(N)) mod p and x ≡ (-m - sqrt(N)) mod p
     */
    for (int i = 2; i < fb_size; i++) { /* skip sign and 2 */
        int p = fb_primes[i];
        int logp = (int)(log(p) * 1.44 + 0.5);
        if (logp < 1) logp = 1;
        if (logp > 127) logp = 127;

        int r1 = root1[i], r2 = root2[i];
        if (r1 < 0) continue;

        int off1 = ((int)((r1 - block_start % p + p) % p));
        for (int j = off1; j < block_size; j += p) {
            int v = sieve_arr[j] + logp;
            sieve_arr[j] = (v > 255) ? 255 : (unsigned char)v;
        }

        if (r2 >= 0 && r2 != r1) {
            int off2 = ((int)((r2 - block_start % p + p) % p));
            for (int j = off2; j < block_size; j += p) {
                int v = sieve_arr[j] + logp;
                sieve_arr[j] = (v > 255) ? 255 : (unsigned char)v;
            }
        }
    }

    /* Scan and trial divide candidates */
    mpz_t Qx, Y_val, cofactor;
    mpz_init(Qx); mpz_init(Y_val); mpz_init(cofactor);
    int *exponents = malloc(fb_size * sizeof(int));

    for (int j = 0; j < block_size; j++) {
        if (sieve_arr[j] < thresh) continue;

        long x = block_start + j;
        if (x <= 0) continue; /* only x > 0 for Fermat sieve */

        /* Q(x) = x^2 + 2mx + r = (m+x)^2 - N */
        mpz_set_si(Y_val, x);
        mpz_add(Y_val, Y_val, m);    /* Y = m + x */
        mpz_mul(Qx, Y_val, Y_val);
        mpz_sub(Qx, Qx, N);          /* Q = Y^2 - N */

        if (mpz_sgn(Qx) <= 0) continue;

        /* Trial divide */
        memset(exponents, 0, fb_size * sizeof(int));
        mpz_set(cofactor, Qx);

        /* Handle sign */
        if (mpz_sgn(cofactor) < 0) {
            exponents[0] = 1;
            mpz_neg(cofactor, cofactor);
        }

        /* Trial divide by FB primes */
        for (int i = 1; i < fb_size; i++) {
            int p = fb_primes[i];
            while (mpz_divisible_ui_p(cofactor, p)) {
                mpz_divexact_ui(cofactor, cofactor, p);
                exponents[i]++;
            }
            if (mpz_cmp_ui(cofactor, 1) == 0) break;
        }

        int lp = 0;
        if (mpz_cmp_ui(cofactor, 1) == 0) {
            /* Fully smooth */
        } else if (mpz_fits_ulong_p(cofactor) &&
                   mpz_get_ui(cofactor) <= (unsigned long)lp_bound &&
                   mpz_probab_prime_p(cofactor, 2)) {
            lp = mpz_get_ui(cofactor);
        } else {
            continue; /* Not smooth enough */
        }

        add_rel(Y_val, exponents, lp);
        found++;
    }

    free(exponents);
    mpz_clear(Qx); mpz_clear(Y_val); mpz_clear(cofactor);
    return found;
}

/* ==================== GF(2) solve + sqrt ==================== */

static int try_factor(mpz_t factor, const mpz_t N,
                      const int *dep, int dep_size) {
    mpz_t X, Y, temp;
    mpz_init_set_ui(X, 1);
    mpz_init_set_ui(Y, 1);
    mpz_init(temp);

    int *total_exp = calloc(fb_size, sizeof(int));
    int *lps = malloc(dep_size * sizeof(int));
    int nlps = 0;

    for (int i = 0; i < dep_size; i++) {
        int ri = dep[i];
        /* X = product of Y_i = product of (m + x_i) mod N */
        mpz_mul(X, X, rels[ri].Y);
        mpz_mod(X, X, N);

        for (int j = 0; j < fb_size; j++)
            total_exp[j] += rels[ri].exponents[j];

        if (rels[ri].large_prime > 1)
            lps[nlps++] = rels[ri].large_prime;
    }

    /* Check all exponents even */
    int ok = 1;
    for (int j = 0; j < fb_size && ok; j++)
        if (total_exp[j] % 2) ok = 0;

    /* Check LP exponents */
    if (ok && nlps > 0) {
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

    if (!ok) {
        free(total_exp); free(lps);
        mpz_clear(X); mpz_clear(Y); mpz_clear(temp);
        return 0;
    }

    /* Y = product of p^(exp/2) mod N */
    for (int j = 1; j < fb_size; j++) {
        if (total_exp[j] == 0) continue;
        mpz_set_ui(temp, fb_primes[j]);
        mpz_powm_ui(temp, temp, total_exp[j] / 2, N);
        mpz_mul(Y, Y, temp);
        mpz_mod(Y, Y, N);
    }

    /* Include LP contributions */
    for (int i = 0; i < nlps; ) {
        int cnt = 1;
        while (i + cnt < nlps && lps[i + cnt] == lps[i]) cnt++;
        mpz_set_ui(temp, lps[i]);
        mpz_powm_ui(temp, temp, cnt / 2, N);
        mpz_mul(Y, Y, temp);
        mpz_mod(Y, Y, N);
        i += cnt;
    }

    /* gcd(X - Y, N) or gcd(X + Y, N) */
    mpz_sub(temp, X, Y);
    mpz_gcd(factor, temp, N);
    int success = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);

    if (!success) {
        mpz_add(temp, X, Y);
        mpz_gcd(factor, temp, N);
        success = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);
    }

    free(total_exp); free(lps);
    mpz_clear(X); mpz_clear(Y); mpz_clear(temp);
    return success;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N, m, r_val, factor;
    mpz_init_set_str(N, argv[1], 10);
    mpz_init(m); mpz_init(r_val); mpz_init(factor);

    int digits = mpz_sizeinbase(N, 10);
    int bits = mpz_sizeinbase(N, 2);
    fprintf(stderr, "Balanced NFS (Fermat sieve): %d digits (%d bits)\n", digits, bits);

    /* m = floor(sqrt(N)), r = m^2 - N */
    mpz_sqrt(m, N);
    mpz_mul(r_val, m, m);
    mpz_sub(r_val, r_val, N); /* r = m^2 - N (negative) */
    fprintf(stderr, "m = floor(sqrt(N)): %zu bits, |r| = |m^2 - N|: %zu bits\n",
            mpz_sizeinbase(m, 2), mpz_sizeinbase(r_val, 2));

    /* Parameters */
    int target_fb;
    if (digits <= 30) target_fb = 100;
    else if (digits <= 35) target_fb = 200;
    else if (digits <= 40) target_fb = 400;
    else if (digits <= 45) target_fb = 800;
    else if (digits <= 50) target_fb = 1500;
    else if (digits <= 55) target_fb = 2800;
    else if (digits <= 60) target_fb = 5000;
    else target_fb = 10000;

    build_fb(N, target_fb);
    fprintf(stderr, "FB: %d primes, bound=%d\n", fb_size, fb_bound);

    long lp_bound = (long)fb_bound * 80;

    /* Compute sieve roots: for each p, Q(x) ≡ 0 mod p when
     * (m+x)^2 ≡ N mod p, i.e., x ≡ sqrt(N) - m mod p or x ≡ -sqrt(N) - m mod p */
    int *root1 = malloc(fb_size * sizeof(int));
    int *root2 = malloc(fb_size * sizeof(int));

    for (int i = 0; i < fb_size; i++) {
        int p = fb_primes[i];
        if (p <= 0) { root1[i] = root2[i] = -1; continue; }
        if (p == 2) { root1[i] = root2[i] = -1; continue; } /* handle 2 separately */

        int sqrtN_p = tonelli_shanks_fast(mpz_fdiv_ui(N, p), p);
        if (sqrtN_p < 0) { root1[i] = root2[i] = -1; continue; }

        int m_mod_p = mpz_fdiv_ui(m, p);
        root1[i] = ((sqrtN_p - m_mod_p) % p + p) % p;
        root2[i] = ((p - sqrtN_p - m_mod_p) % p + p) % p;
    }

    /* Sieve interval: x from 1 to M */
    int block_size = 32768;
    int num_blocks;
    if (digits <= 30) num_blocks = 4;
    else if (digits <= 35) num_blocks = 8;
    else if (digits <= 40) num_blocks = 16;
    else if (digits <= 45) num_blocks = 32;
    else if (digits <= 50) num_blocks = 64;
    else num_blocks = 128;

    int M = num_blocks * block_size;
    fprintf(stderr, "Sieve interval: [1, %d], threshold computing...\n", M);

    /* Threshold: log2(2*m*M) * scale */
    double log2_qmax = mpz_sizeinbase(m, 2) + 1 + log2(M);
    unsigned char thresh = (unsigned char)(log2_qmax * 0.62);
    if (thresh < 20) thresh = 20;
    if (thresh > 200) thresh = 200;
    fprintf(stderr, "Threshold: %d (est log2(Q_max) = %.1f)\n", thresh, log2_qmax);

    int target_rels = fb_size + 64;
    init_lp();
    unsigned char *sieve_arr = malloc(block_size);

    fprintf(stderr, "Need %d relations\n", target_rels);

    /* Sieve blocks */
    int total_blocks = 0;
    long sieve_start = 1;

    while (nrels + ncombined < target_rels) {
        for (int b = 0; b < num_blocks && nrels + ncombined < target_rels; b++) {
            sieve_block(sieve_arr, sieve_start + b * block_size, block_size,
                        m, r_val, N, root1, root2, lp_bound, thresh);
            total_blocks++;
        }
        sieve_start += num_blocks * block_size;

        /* Progress */
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        double el = (now.tv_sec - t0.tv_sec) + (now.tv_nsec - t0.tv_nsec) / 1e9;
        int full = 0, partial = 0;
        for (int i = 0; i < nrels; i++) {
            if (rels[i].large_prime <= 1) full++; else partial++;
        }
        fprintf(stderr, "\r[%.1fs] blocks=%d x_max=%ld rels=%d(f=%d p=%d c=%d) "
                "need=%d rate=%.0f/s  ",
                el, total_blocks, sieve_start,
                nrels, full, partial, ncombined,
                target_rels, nrels > 0 ? nrels / el : 0);

        if (sieve_start > 1000000000L) {
            fprintf(stderr, "\nGiving up: sieve interval exhausted\n");
            break;
        }
    }

    struct timespec t1;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_time = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "\nSieve: %d rels in %.2fs (%.1f/s)\n", nrels, sieve_time, nrels / sieve_time);

    free(sieve_arr);

    /* ==================== LA ==================== */
    /* Collect usable relations (full + paired SLP) */
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
    fprintf(stderr, "LA: %d usable rels, %d cols (%d FB + %d LP)\n",
            n_usable, ncols, fb_size, n_dlps);

    if (n_usable <= ncols) {
        fprintf(stderr, "FAILED: not enough usable rels (%d <= %d cols)\n", n_usable, ncols);
        free(distinct_lps); free(usable);
        goto cleanup;
    }

    /* Build and solve GF(2) matrix */
    int rw = (ncols + 63) / 64;
    int aw = (n_usable + 63) / 64;
    uint64_t **mat = malloc(n_usable * sizeof(uint64_t *));
    uint64_t **aug = malloc(n_usable * sizeof(uint64_t *));
    int *rel_idx = malloc(n_usable * sizeof(int));

    int row = 0;
    for (int i = 0; i < nrels; i++) {
        if (!usable[i]) continue;
        mat[row] = calloc(rw, sizeof(uint64_t));
        aug[row] = calloc(aw, sizeof(uint64_t));
        aug[row][row / 64] |= (1ULL << (row % 64));

        for (int j = 0; j < fb_size; j++)
            if (rels[i].exponents[j] & 1)
                mat[row][j / 64] |= (1ULL << (j % 64));

        if (rels[i].large_prime > 1)
            for (int j = 0; j < n_dlps; j++)
                if (distinct_lps[j] == rels[i].large_prime) {
                    int col = fb_size + j;
                    mat[row][col / 64] |= (1ULL << (col % 64));
                    break;
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

    /* Try dependencies */
    int factored = 0;
    for (int r = 0; r < n_usable && !factored; r++) {
        int zero = 1;
        for (int w = 0; w < rw && zero; w++) if (mat[r][w]) zero = 0;
        if (!zero) continue;

        /* Extract dependency */
        int *dep = malloc(n_usable * sizeof(int));
        int dep_size = 0;
        for (int i = 0; i < n_usable; i++)
            if ((aug[r][i / 64] >> (i % 64)) & 1)
                dep[dep_size++] = rel_idx[i];

        if (dep_size >= 2 && try_factor(factor, N, dep, dep_size)) {
            struct timespec t2;
            clock_gettime(CLOCK_MONOTONIC, &t2);
            double total = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec) / 1e9;

            mpz_t cofn; mpz_init(cofn);
            mpz_divexact(cofn, N, factor);
            if (mpz_cmp(factor, cofn) > 0) mpz_swap(factor, cofn);

            gmp_printf("Factor: %Zd\n", factor);
            gmp_printf("Cofactor: %Zd\n", cofn);
            fprintf(stderr, "Time: %.3fs (sieve %.3fs, LA %.3fs)\n",
                    total, sieve_time, total - sieve_time);
            factored = 1;
            mpz_clear(cofn);
        }
        free(dep);
    }

    if (!factored)
        fprintf(stderr, "FAILED to find factor from %d dependencies\n", n_usable);

    for (int i = 0; i < n_usable; i++) { free(mat[i]); free(aug[i]); }
    free(mat); free(aug); free(rel_idx); free(pcol);
    free(distinct_lps); free(usable);

cleanup:
    free(root1); free(root2);
    mpz_clear(N); mpz_clear(m); mpz_clear(r_val); mpz_clear(factor);
    return factored ? 0 : 1;
}
