/*
 * specialq_siqs.c - Special-Q SIQS: novel QS variant using lattice sieving
 *
 * Key idea: Apply NFS-style "special-Q" technique to QS.
 * For each large prime Q not in the factor base:
 *   1. Find positions x where Q | QS_poly(x)
 *   2. Divide out Q from QS_poly(x) at those positions
 *   3. The quotient QS_poly(x)/Q is smaller → more likely B-smooth
 *   4. Relations from special-Q: QS_poly(x) = Q * (smooth cofactor)
 *
 * This is different from standard large prime variation:
 * - LP variation: sieve normally, accept relations with 1-2 large prime factors
 * - Special-Q: CHOOSE a large prime Q upfront, sieve only positions divisible by Q
 *
 * The advantage: for each special-Q, the effective polynomial values are
 * QS_poly(x)/Q ≈ sqrt(N)/Q, which is Q times smaller than standard QS values.
 * This dramatically increases smoothness probability per candidate.
 *
 * The cost: only ~2/Q fraction of sieve positions are divisible by Q,
 * so we need Q different special-Q values to cover the same sieve space.
 *
 * Net effect: smoothness increases by Q^u (where u = ln(sqrt(N))/ln(B)),
 * but we process Q times fewer candidates. For u > 1, the smoothness
 * gain dominates, giving a net speedup.
 *
 * Compile: gcc -O3 -march=native -o specialq_siqs library/specialq_siqs.c -lgmp -lm
 * Usage: ./specialq_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_FB 50000
#define SIEVE_BLOCK 32768
#define SEED 42

static int primes[MAX_FB];
static int fb_roots[MAX_FB]; /* sqrt(N) mod p */
static int fb_size;
static int fb_logp[MAX_FB];

static void gen_primes(int limit) {
    char *s = calloc(limit + 1, 1);
    fb_size = 0;
    for (int i = 2; i <= limit; i++) s[i] = 1;
    for (int i = 2; (long)i * i <= limit; i++)
        if (s[i]) for (int j = i * i; j <= limit; j += i) s[j] = 0;
    for (int i = 2; i <= limit && fb_size < MAX_FB; i++)
        if (s[i]) primes[fb_size++] = i;
    free(s);
}

static uint32_t modsqrt(uint32_t a, uint32_t p) {
    if (a == 0) return 0;
    if (p == 2) return a & 1;
    /* Euler criterion */
    uint64_t test = 1, base = a % p;
    uint32_t exp = (p - 1) / 2;
    while (exp > 0) {
        if (exp & 1) test = test * base % p;
        base = base * base % p;
        exp >>= 1;
    }
    if (test != 1) return UINT32_MAX; /* not QR */

    /* Tonelli-Shanks */
    uint32_t Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }
    if (S == 1) {
        base = a; exp = (p + 1) / 4;
        uint64_t r = 1;
        while (exp > 0) {
            if (exp & 1) r = r * base % p;
            base = base * base % p;
            exp >>= 1;
        }
        return (uint32_t)r;
    }
    uint32_t z = 2;
    while (1) {
        base = z; exp = (p - 1) / 2; test = 1;
        while (exp > 0) {
            if (exp & 1) test = test * base % p;
            base = base * base % p;
            exp >>= 1;
        }
        if (test == p - 1) break;
        z++;
    }
    uint32_t M = S;
    uint64_t c = 1; base = z; exp = Q;
    while (exp > 0) { if (exp & 1) c = c * base % p; base = base * base % p; exp >>= 1; }
    uint64_t t = 1; base = a; exp = Q;
    while (exp > 0) { if (exp & 1) t = t * base % p; base = base * base % p; exp >>= 1; }
    uint64_t R = 1; base = a; exp = (Q + 1) / 2;
    while (exp > 0) { if (exp & 1) R = R * base % p; base = base * base % p; exp >>= 1; }

    while (t != 1) {
        uint32_t i = 1; uint64_t tmp = t * t % p;
        while (tmp != 1) { tmp = tmp * tmp % p; i++; }
        uint64_t b = c;
        for (uint32_t j = 0; j < M - i - 1; j++) b = b * b % p;
        M = i; c = b * b % p; t = t * c % p; R = R * b % p;
    }
    return (uint32_t)R;
}

static uint32_t modinv(uint32_t a, uint32_t m) {
    int64_t g = m, x = 0, y = 1, a1 = a;
    while (a1) { int64_t q = g / a1, t = g - q * a1; g = a1; a1 = t; t = x - q * y; x = y; y = t; }
    return (uint32_t)((x % (int64_t)m + m) % m);
}

/* Build factor base for kN */
static int build_fb(mpz_t kN, int target) {
    gen_primes(target * 15);
    int count = 0;
    int *new_primes = malloc(target * sizeof(int));
    int *new_roots = malloc(target * sizeof(int));
    int *new_logp = malloc(target * sizeof(int));

    for (int i = 0; i < fb_size && count < target; i++) {
        int p = primes[i];
        uint32_t kn_mod = mpz_fdiv_ui(kN, p);
        uint32_t s = modsqrt(kn_mod, p);
        if (s == UINT32_MAX) continue;
        new_primes[count] = p;
        new_roots[count] = s;
        new_logp[count] = (int)(log2(p) * 1.44 + 0.5);
        count++;
    }

    memcpy(primes, new_primes, count * sizeof(int));
    memcpy(fb_roots, new_roots, count * sizeof(int));
    memcpy(fb_logp, new_logp, count * sizeof(int));
    fb_size = count;
    free(new_primes); free(new_roots); free(new_logp);
    return count;
}

/*
 * Special-Q sieving:
 * For special prime Q with root r (QS_poly(r) ≡ 0 mod Q):
 *   Sieve positions: x = r + k*Q for k = 0, ±1, ±2, ...
 *   At each position: evaluate QS_poly(x) / Q, test for B-smoothness
 *
 * QS polynomial: Q(x) = (x + ceil(sqrt(kN)))^2 - kN
 * Q(r) ≡ 0 mod Q means (r + s)^2 ≡ kN mod Q where s = ceil(sqrt(kN))
 * So r ≡ ±sqrt(kN mod Q) - s mod Q
 */
static int specialq_sieve(mpz_t N, mpz_t kN, int special_q, int q_root,
                           int sieve_half, mpz_t sqrt_kN) {
    int rels_found = 0;
    uint8_t *sieve = calloc(SIEVE_BLOCK, 1);

    /* Map: sieve position j → x = q_root + (j - sieve_half) * special_q */
    /* But that spreads positions too far. Instead, sieve the reduced values. */

    /* For positions x where Q | (x+s)^2 - kN:
     * x = q_root + k*special_q for integer k
     * QS_val(x) = (x+s)^2 - kN, which is divisible by special_q
     * Reduced value: QS_val(x) / special_q
     */

    /* Sieve the reduced values for B-smoothness */
    /* For each FB prime p:
     *   Find sieve positions where p | reduced_val
     *   p | QS_val(x)/Q iff p | QS_val(x) iff (x+s)^2 ≡ kN (mod p)
     *   Same roots as normal QS, but only at x = q_root + k*Q
     *   In sieve coordinates: k such that q_root + k*Q ≡ root_p (mod p)
     *   k ≡ (root_p - q_root) * Q^(-1) (mod p)
     */

    uint32_t Q_inv[MAX_FB]; /* Q^(-1) mod p for each FB prime */
    for (int i = 0; i < fb_size; i++) {
        int p = primes[i];
        if (p == special_q) { Q_inv[i] = 0; continue; }
        Q_inv[i] = modinv(special_q % p, p);
    }

    /* Sieve threshold: log2(QS_val/Q) ≈ log2(sqrt(kN) * M) - log2(Q)
     * where M = sieve_half * Q */
    int bits_kN = mpz_sizeinbase(kN, 2);
    int bits_val = bits_kN / 2 + (int)log2(sieve_half) + (int)log2(special_q);
    int bits_reduced = bits_val - (int)log2(special_q);
    int threshold = (int)(bits_reduced * 0.7); /* 70% of expected size */

    /* Initialize sieve */
    memset(sieve, 0, SIEVE_BLOCK);
    int sieve_len = 2 * sieve_half;
    if (sieve_len > SIEVE_BLOCK) sieve_len = SIEVE_BLOCK;

    /* Sieve: for each FB prime, add logp at divisible positions */
    for (int i = 0; i < fb_size; i++) {
        int p = primes[i];
        if (p == special_q || p < 3) continue;
        int logp = fb_logp[i];
        int s = fb_roots[i];

        /* Two roots: x ≡ s - sqrt_kN_mod_p, x ≡ -s - sqrt_kN_mod_p (mod p) */
        int sqrt_kN_mod_p = mpz_fdiv_ui(sqrt_kN, p);
        int root1 = ((s - sqrt_kN_mod_p) % p + p) % p;
        int root2 = ((-s - sqrt_kN_mod_p) % p + p) % p;

        /* Convert to sieve coordinates: k such that q_root + k*Q ≡ root (mod p) */
        /* k ≡ (root - q_root) * Q_inv (mod p) */
        int k1 = ((int64_t)(root1 - q_root % p + p) % p * Q_inv[i]) % p;
        int k2 = ((int64_t)(root2 - q_root % p + p) % p * Q_inv[i]) % p;

        /* Add sieve_half offset */
        int pos1 = (k1 + sieve_half) % p;
        int pos2 = (k2 + sieve_half) % p;

        for (int j = pos1; j < sieve_len; j += p) sieve[j] += logp;
        if (root1 != root2) {
            for (int j = pos2; j < sieve_len; j += p) sieve[j] += logp;
        }
    }

    /* Scan for candidates */
    mpz_t qval, reduced, x_val, cofactor;
    mpz_inits(qval, reduced, x_val, cofactor, NULL);

    for (int j = 0; j < sieve_len; j++) {
        if (sieve[j] < threshold) continue;

        /* Compute x = q_root + (j - sieve_half) * special_q */
        int64_t k = j - sieve_half;
        mpz_set_si(x_val, k);
        mpz_mul_ui(x_val, x_val, special_q);
        mpz_add_ui(x_val, x_val, q_root);

        /* QS_val = (x + sqrt_kN)^2 - kN */
        mpz_add(qval, x_val, sqrt_kN);
        mpz_mul(qval, qval, qval);
        mpz_sub(qval, qval, kN);

        if (mpz_sgn(qval) <= 0) continue;

        /* Check Q divides QS_val */
        if (!mpz_divisible_ui_p(qval, special_q)) continue;

        /* Reduced value = QS_val / Q */
        mpz_divexact_ui(reduced, qval, special_q);

        /* Trial divide reduced value by factor base */
        mpz_set(cofactor, reduced);
        int smooth = 1;
        for (int i = 0; i < fb_size && mpz_cmp_ui(cofactor, 1) > 0; i++) {
            while (mpz_divisible_ui_p(cofactor, primes[i]))
                mpz_divexact_ui(cofactor, cofactor, primes[i]);
        }

        if (mpz_cmp_ui(cofactor, 1) == 0) {
            rels_found++;
        } else if (mpz_cmp_ui(cofactor, (uint64_t)primes[fb_size-1] * 30) < 0) {
            /* SLP relation */
            rels_found++;
        }
    }

    free(sieve);
    mpz_clears(qval, reduced, x_val, cofactor, NULL);
    return rels_found;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [fb_size]\n", argv[0]);
        return 1;
    }

    mpz_t N, kN, sqrt_kN;
    mpz_inits(N, kN, sqrt_kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);
    int bits = mpz_sizeinbase(N, 2);

    /* Parameters */
    int fb_target = argc > 2 ? atoi(argv[2]) :
        (digits <= 30 ? 200 : digits <= 40 ? 500 : digits <= 50 ? 1500 :
         digits <= 60 ? 4000 : digits <= 70 ? 10000 : 30000);

    /* Use multiplier k=1 for simplicity */
    mpz_set(kN, N);
    mpz_sqrt(sqrt_kN, kN);
    mpz_add_ui(sqrt_kN, sqrt_kN, 1);

    fprintf(stderr, "N: %d digits, %d bits\n", digits, bits);

    /* Build factor base */
    build_fb(kN, fb_target);
    fprintf(stderr, "Factor base: %d primes, max %d\n", fb_size, primes[fb_size-1]);

    int target_rels = fb_size + 30;
    int total_rels = 0;

    /* Special-Q range: use primes from fb_max to 10*fb_max */
    int q_min = primes[fb_size - 1] + 1;
    int q_max = q_min * 10;
    int sieve_half = SIEVE_BLOCK / 2;

    fprintf(stderr, "Special-Q range: %d to %d\n", q_min, q_max);
    fprintf(stderr, "Target: %d relations\n", target_rels);

    /* Generate special-Q primes */
    int sq_limit = q_max + 1000;
    char *is_prime = calloc(sq_limit + 1, 1);
    for (int i = 2; i <= sq_limit; i++) is_prime[i] = 1;
    for (int i = 2; (long)i * i <= sq_limit; i++)
        if (is_prime[i]) for (int j = i * i; j <= sq_limit; j += i) is_prime[j] = 0;

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int q_count = 0;
    for (int q = q_min; q <= q_max && total_rels < target_rels; q++) {
        if (!is_prime[q]) continue;

        /* Find roots of kN mod q */
        uint32_t kn_mod_q = mpz_fdiv_ui(kN, q);
        uint32_t sq_root = modsqrt(kn_mod_q, q);
        if (sq_root == UINT32_MAX) continue;

        /* sqrt_kN mod q */
        uint32_t sqrt_mod = mpz_fdiv_ui(sqrt_kN, q);

        /* Two roots: r1 = sq_root - sqrt_mod, r2 = q - sq_root - sqrt_mod */
        int r1 = ((int64_t)sq_root - sqrt_mod + q) % q;
        int r2 = ((int64_t)q - sq_root - sqrt_mod + q) % q;

        int found1 = specialq_sieve(N, kN, q, r1, sieve_half, sqrt_kN);
        int found2 = (r1 != r2) ? specialq_sieve(N, kN, q, r2, sieve_half, sqrt_kN) : 0;
        total_rels += found1 + found2;
        q_count++;

        if (q_count % 100 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
            fprintf(stderr, "Q=%d, %d special-Qs, %d rels, %.1fs, %.1f rels/sec\n",
                    q, q_count, total_rels, elapsed, total_rels / elapsed);
            if (elapsed > 280) break;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &now);
    double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "Total: %d rels from %d special-Qs in %.1fs (%.1f rels/sec)\n",
            total_rels, q_count, elapsed, total_rels / elapsed);

    if (total_rels >= target_rels) {
        fprintf(stderr, "Enough relations! (LA not implemented yet)\n");
    } else {
        fprintf(stderr, "Need %d more relations\n", target_rels - total_rels);
    }

    free(is_prime);
    mpz_clears(N, kN, sqrt_kN, NULL);
    return total_rels >= target_rels ? 0 : 1;
}
