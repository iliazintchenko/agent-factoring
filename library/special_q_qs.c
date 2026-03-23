/*
 * special_q_qs.c - Special-Q Enhanced Quadratic Sieve
 *
 * Novel approach: combines NFS-style special-Q lattice sieving with QS.
 * Instead of sieving Q(x) = (Ax+B)² - N and hoping for smooth values,
 * we FIX one large prime factor q and sieve Q(x)/q for smoothness.
 *
 * For each SIQS polynomial Q(x) = ((Ax+B)² - kN)/A:
 * - For a chosen "special Q" prime q, find x values where q | Q(x)
 * - These are x ≡ r (mod q) where (Ax+B)² ≡ kN (mod q)
 * - Then Q(x)/q is ~q times smaller, so MUCH more likely smooth
 * - This gives guaranteed DLP relations (one LP = q)
 *
 * The advantage: directed search for DLP relations is more efficient
 * than the standard sieve approach of hoping to find them.
 *
 * Compile: gcc -O3 -march=native -mavx512bw -o special_q_qs library/special_q_qs.c -lgmp -lm
 * Usage: ./special_q_qs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>

#define SIEVE_SIZE    32768
#define MAX_FB        50000
#define MAX_RELS      300000

/* ==================== Factor Base ==================== */
static int primes[MAX_FB + 1000]; // extra for special-Q primes
static int nprimes;
static int fb[MAX_FB];
static int fb_sqrtN[MAX_FB]; // sqrt(kN) mod fb[i]
static int fb_count;
static uint8_t fb_logp[MAX_FB];

static void sieve_primes(int limit) {
    char *is_composite = calloc(limit + 1, 1);
    nprimes = 0;
    for (int i = 2; i <= limit; i++) {
        if (!is_composite[i]) {
            primes[nprimes++] = i;
            for (long j = (long)i * i; j <= limit; j += i)
                is_composite[j] = 1;
        }
    }
    free(is_composite);
}

// Tonelli-Shanks for modular square root
static int mod_sqrt(int n, int p) {
    n = ((n % p) + p) % p;
    if (n == 0) return 0;
    if (p == 2) return n & 1;

    // Check QR
    long long pw = 1;
    int exp = (p - 1) / 2;
    long long base = n;
    while (exp > 0) {
        if (exp & 1) pw = pw * base % p;
        base = base * base % p;
        exp >>= 1;
    }
    if (pw != 1) return -1;

    if (p % 4 == 3) {
        pw = 1; exp = (p + 1) / 4; base = n;
        while (exp > 0) {
            if (exp & 1) pw = pw * base % p;
            base = base * base % p;
            exp >>= 1;
        }
        return (int)pw;
    }

    // Full Tonelli-Shanks
    int q = p - 1, s = 0;
    while (q % 2 == 0) { q /= 2; s++; }

    int z = 2;
    while (1) {
        pw = 1; exp = (p - 1) / 2; base = z;
        while (exp > 0) {
            if (exp & 1) pw = pw * base % p;
            base = base * base % p;
            exp >>= 1;
        }
        if (pw == p - 1) break;
        z++;
    }

    int M = s;
    long long c = 1; base = z; exp = q;
    while (exp > 0) { if (exp & 1) c = c * base % p; base = base * base % p; exp >>= 1; }
    long long t = 1; base = n; exp = q;
    while (exp > 0) { if (exp & 1) t = t * base % p; base = base * base % p; exp >>= 1; }
    long long r = 1; base = n; exp = (q + 1) / 2;
    while (exp > 0) { if (exp & 1) r = r * base % p; base = base * base % p; exp >>= 1; }

    while (t != 1) {
        int i = 0; long long tmp = t;
        while (tmp != 1) { tmp = tmp * tmp % p; i++; }
        long long b = c;
        for (int j = 0; j < M - i - 1; j++) b = b * b % p;
        M = i;
        c = b * b % p;
        t = t * c % p;
        r = r * b % p;
    }
    return (int)r;
}

static int knuth_schroeppel(mpz_t N) {
    static const int multipliers[] = {1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 29, 31, 33, 35, 37, 39, 41, 43};
    double best_score = -1e30;
    int best_k = 1;

    for (int mi = 0; mi < 19; mi++) {
        int k = multipliers[mi];
        mpz_t kN;
        mpz_init(kN);
        mpz_mul_ui(kN, N, k);

        double score = -0.5 * log(k);
        int kN_mod8 = mpz_fdiv_ui(kN, 8);
        if (kN_mod8 == 1) score += 2 * log(2);
        else if (kN_mod8 == 5) score += log(2);

        for (int j = 1; j < nprimes && primes[j] < 200; j++) {
            int p = primes[j];
            int kN_mod_p = mpz_fdiv_ui(kN, p);
            int r = mod_sqrt(kN_mod_p, p);
            if (r >= 0) {
                double contrib = 2.0 * log(p) / (p - 1);
                if (k % p == 0) contrib *= 0.5;
                score += contrib;
            }
        }

        if (score > best_score) { best_score = score; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Parameter Selection ==================== */
typedef struct {
    int fb_size;
    int sieve_radius; // M = half sieve interval
    int num_a_factors;
    int lp_bound_mult;
    int special_q_count; // number of special-Q primes to use
} params_t;

static params_t get_params(int digits) {
    if (digits <= 30) return (params_t){200,   32768,  4, 30, 100};
    if (digits <= 35) return (params_t){400,   65536,  5, 40, 200};
    if (digits <= 40) return (params_t){800,   65536,  6, 50, 400};
    if (digits <= 45) return (params_t){1500,  131072, 7, 60, 600};
    if (digits <= 50) return (params_t){3000,  131072, 7, 80, 1000};
    if (digits <= 55) return (params_t){5000,  196608, 8, 100, 1500};
    if (digits <= 60) return (params_t){8000,  262144, 9, 120, 2000};
    if (digits <= 65) return (params_t){14000, 327680, 10, 150, 3000};
    if (digits <= 70) return (params_t){25000, 393216, 10, 200, 5000};
    return (params_t){40000, 524288, 11, 250, 8000};
}

/* ==================== Relation Storage ==================== */
typedef struct {
    int *exponents;  // exponent vector over FB
    int nexp;
    mpz_t x_val;     // x such that x² ≡ product (mod N)
    uint32_t lp1, lp2; // large primes (0 if none)
    int sign;         // 1 if Q(x) < 0
} relation_t;

static relation_t rels[MAX_RELS];
static int nrels;

/* Hash table for partial relations (single LP) */
#define LP_HASH_SIZE (1 << 20)
#define LP_HASH_MASK (LP_HASH_SIZE - 1)
static int lp_hash[LP_HASH_SIZE];
static int lp_hash_count;

static void init_lp_hash(void) {
    memset(lp_hash, -1, sizeof(lp_hash));
    lp_hash_count = 0;
}

/* ==================== Special-Q Sieve Core ==================== */
/*
 * For a given SIQS polynomial Q(x) = ((Ax+B)² - kN) / A
 * and a special-Q prime q (not in FB), find x values where q | Q(x).
 *
 * q | Q(x) means (Ax+B)² ≡ kN (mod q)
 * So Ax+B ≡ ±sqrt(kN) (mod q)
 * x ≡ (±sqrt(kN) - B) * A^(-1) (mod q)
 *
 * For each such x₀, sieve Q(x)/q over x = x₀ + j*q for j in sieve range.
 * Q(x₀ + j*q)/q is a value of size ~M*sqrt(N)/(A*q).
 * This is q times smaller than Q(x₀), hence more likely smooth.
 */

static uint8_t sieve_arr[SIEVE_SIZE] __attribute__((aligned(64)));

/* Compute sieve solution for prime p given polynomial (A, B)
 * Q(x) = ((Ax+B)² - kN) / A
 * p | Q(x) when (Ax+B)² ≡ kN (mod p), i.e., Ax ≡ sqrt(kN) - B (mod p)
 */
static void compute_sieve_roots(int *root1, int *root2,
                                 mpz_t A, mpz_t B, mpz_t kN,
                                 int p, int sqrtN_mod_p) {
    int A_mod = (int)mpz_fdiv_ui(A, p);
    int B_mod = (int)mpz_fdiv_ui(B, p);

    if (A_mod == 0) { *root1 = *root2 = -1; return; }

    // Compute A^(-1) mod p
    long long A_inv = 1, base = A_mod, exp = p - 2;
    while (exp > 0) {
        if (exp & 1) A_inv = A_inv * base % p;
        base = base * base % p;
        exp >>= 1;
    }

    // root1 = (sqrtN - B) * A^(-1) mod p
    int r1 = (int)(((long long)(sqrtN_mod_p - B_mod + p) % p) * A_inv % p);
    // root2 = (-sqrtN - B) * A^(-1) mod p
    int r2 = (int)(((long long)(p - sqrtN_mod_p - B_mod + p) % p) * A_inv % p);

    *root1 = r1;
    *root2 = (r1 != r2) ? r2 : -1;
}

/* ==================== Main Sieve ==================== */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N, kN, A, B, Qval, x_val, factor, cofactor, temp;
    mpz_inits(N, kN, A, B, Qval, x_val, factor, cofactor, temp, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    fprintf(stderr, "Factoring %d-digit number\n", digits);

    // Check small factors
    for (int p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_set_ui(factor, p);
            mpz_divexact(cofactor, N, factor);
            gmp_printf("%Zd\n%Zd\n", factor, cofactor);
            return 0;
        }
    }

    sieve_primes(2000000);
    int k = knuth_schroeppel(N);
    mpz_mul_ui(kN, N, k);
    fprintf(stderr, "Multiplier: %d\n", k);

    params_t params = get_params(digits);

    // Build factor base
    fb_count = 0;
    fb[fb_count] = -1; fb_sqrtN[fb_count] = 0; fb_logp[fb_count] = 0; fb_count++;
    fb[fb_count] = 2; fb_sqrtN[fb_count] = (int)mpz_fdiv_ui(kN, 2); fb_logp[fb_count] = 1; fb_count++;

    for (int i = 1; i < nprimes && fb_count < params.fb_size; i++) {
        int p = primes[i];
        if (p == 2) continue;
        int kN_mod_p = (int)mpz_fdiv_ui(kN, p);
        int r = mod_sqrt(kN_mod_p, p);
        if (r >= 0) {
            fb[fb_count] = p;
            fb_sqrtN[fb_count] = r;
            fb_logp[fb_count] = (uint8_t)(log2(p) + 0.5);
            fb_count++;
        }
    }

    int largest_fb = fb[fb_count - 1];
    uint64_t lp_bound = (uint64_t)largest_fb * params.lp_bound_mult;
    int target = fb_count + params.fb_size / 10;

    fprintf(stderr, "FB: %d primes, largest=%d, LP_bound=%lu, target=%d rels\n",
            fb_count, largest_fb, lp_bound, target);

    // Select special-Q primes: primes just above the factor base
    int *sq_primes = malloc(params.special_q_count * sizeof(int));
    int *sq_sqrt = malloc(params.special_q_count * sizeof(int));
    int nsq = 0;

    for (int i = 0; i < nprimes && nsq < params.special_q_count; i++) {
        int p = primes[i];
        if (p <= largest_fb) continue;
        if ((uint64_t)p > lp_bound) break;

        int kN_mod_p = (int)mpz_fdiv_ui(kN, p);
        int r = mod_sqrt(kN_mod_p, p);
        if (r >= 0) {
            sq_primes[nsq] = p;
            sq_sqrt[nsq] = r;
            nsq++;
        }
    }

    fprintf(stderr, "Special-Q: %d primes in range [%d, %lu]\n",
            nsq, largest_fb + 1, lp_bound);

    // Initialize relation storage
    nrels = 0;
    init_lp_hash();

    int fb_sieve_start = 3; // skip -1, 2, and very small primes
    int M = params.sieve_radius;
    int num_blocks = (2 * M + SIEVE_SIZE - 1) / SIEVE_SIZE;

    // RNG for A coefficient selection
    uint32_t rng = 42;
    int poly_count = 0;
    int sq_idx = 0; // current special-Q index

    while (nrels < target && sq_idx < nsq) {
        int q = sq_primes[sq_idx];
        int q_sqrt = sq_sqrt[sq_idx];
        sq_idx++;

        // For this special-Q, generate multiple SIQS polynomials
        // and check the special positions where q | Q(x)

        // Select A coefficient
        rng = rng * 1103515245u + 12345u;

        // Simple A selection: product of s primes from FB
        int s = params.num_a_factors;
        int a_primes[16];
        mpz_set_ui(A, 1);

        double log_target = (mpz_sizeinbase(kN, 2) * log(2.0)) / 2.0 - log(M);
        double ideal_prime_log = log_target / s;
        double ideal_prime = exp(ideal_prime_log);

        int center = fb_sieve_start;
        for (int i = fb_sieve_start; i < fb_count; i++) {
            if (fb[i] >= ideal_prime) { center = i; break; }
        }

        int spread = (fb_count - fb_sieve_start) / 4;
        if (spread < s * 3) spread = s * 3;
        int lo = center - spread / 2;
        int hi = center + spread / 2;
        if (lo < fb_sieve_start) lo = fb_sieve_start;
        if (hi >= fb_count) hi = fb_count - 1;

        int selected = 0;
        for (int j = 0; j < s; j++) {
            rng = rng * 1103515245u + 12345u;
            int idx = lo + (int)(rng % (uint32_t)(hi - lo + 1));

            int ok = 1;
            for (int k2 = 0; k2 < selected; k2++)
                if (a_primes[k2] == idx) { ok = 0; break; }
            if (!ok) { j--; continue; }

            a_primes[selected++] = idx;
            mpz_mul_ui(A, A, fb[idx]);
        }

        // Compute B using Tonelli-Shanks (simplified)
        mpz_sqrt(B, kN);
        mpz_mod(B, B, A);
        // Adjust B so B² ≡ kN (mod A)
        // This is simplified - proper SIQS B computation is more complex
        // For now, just use B = floor(sqrt(kN)) mod A
        mpz_t B2;
        mpz_init(B2);
        mpz_mul(B2, B, B);
        mpz_sub(B2, B2, kN);
        mpz_mod(B2, B2, A);
        if (mpz_sgn(B2) != 0) {
            // B is not a valid square root, try sqrt(kN) mod A more carefully
            // Use CRT with individual a_primes
            mpz_set_ui(B, 0);
            for (int j = 0; j < s; j++) {
                int p = fb[a_primes[j]];
                int r = fb_sqrtN[a_primes[j]];
                mpz_t Ap, Ap_inv, mod_p;
                mpz_inits(Ap, Ap_inv, mod_p, NULL);
                mpz_divexact_ui(Ap, A, p);
                mpz_set_ui(mod_p, p);
                mpz_invert(Ap_inv, Ap, mod_p);
                long long gamma = ((long long)r * mpz_get_ui(Ap_inv)) % p;
                if (gamma > p / 2) gamma = p - gamma;
                mpz_addmul_ui(B, Ap, (unsigned long)gamma);
                mpz_clears(Ap, Ap_inv, mod_p, NULL);
            }
            // Verify B² ≡ kN (mod A)
            mpz_mul(B2, B, B);
            mpz_sub(B2, B2, kN);
            mpz_mod(B2, B2, A);
            if (mpz_sgn(B2) != 0) {
                mpz_sub(B, A, B); // try other sign
            }
        }
        mpz_clear(B2);

        // Find x positions where q | Q(x)
        // Q(x) = ((Ax+B)² - kN) / A
        // q | Q(x) ⟺ (Ax+B)² ≡ kN (mod q)
        // Ax+B ≡ ±sqrt(kN) mod q
        // x ≡ (±sqrt(kN) - B) * A^(-1) mod q

        int A_mod_q = (int)mpz_fdiv_ui(A, q);
        int B_mod_q = (int)mpz_fdiv_ui(B, q);

        if (A_mod_q == 0) continue; // A divisible by q, skip

        // A^(-1) mod q
        long long A_inv_q = 1, base_q = A_mod_q;
        int exp_q = q - 2;
        while (exp_q > 0) {
            if (exp_q & 1) A_inv_q = A_inv_q * base_q % q;
            base_q = base_q * base_q % q;
            exp_q >>= 1;
        }

        // Two roots
        int x_mod_q[2];
        x_mod_q[0] = (int)(((long long)(q_sqrt - B_mod_q + q) % q) * A_inv_q % q);
        x_mod_q[1] = (int)(((long long)(q - q_sqrt - B_mod_q + q) % q) * A_inv_q % q);

        // For each root, generate x values in the sieve interval [-M, M]
        for (int ri = 0; ri < 2; ri++) {
            int x0 = x_mod_q[ri];
            // x = x0 + j*q, offset so x is in [-M, M]
            int first_x = x0 - M;
            if (first_x < 0) {
                int shift = (-first_x + q - 1) / q;
                first_x += shift * q;
            }

            for (int x = first_x - M; x <= M; x += q) {
                if (x < -M || x > M) continue;

                // Compute Q(x) = ((Ax+B)² - kN) / A
                mpz_mul_si(Qval, A, x);
                mpz_add(Qval, Qval, B);
                mpz_mul(Qval, Qval, Qval);
                mpz_sub(Qval, Qval, kN);
                mpz_divexact(Qval, Qval, A);

                // Q(x) should be divisible by q
                if (!mpz_divisible_ui_p(Qval, q)) continue;

                // Divide out q to get cofactor
                mpz_divexact_ui(Qval, Qval, q);

                // Now trial-divide the cofactor by FB primes
                int sign = (mpz_sgn(Qval) < 0) ? 1 : 0;
                mpz_abs(Qval, Qval);

                int exponents[MAX_FB];
                int nexp = 0;

                if (sign) exponents[nexp++] = 0; // sign bit

                // Add q to exponent vector (as a large prime)
                // Actually, we handle q as a known LP, not in the exponent vector

                mpz_t remaining;
                mpz_init_set(remaining, Qval);

                for (int fi = 1; fi < fb_count; fi++) {
                    int p = fb[fi];
                    if (p <= 0) continue;
                    int exp2 = 0;
                    while (mpz_divisible_ui_p(remaining, p)) {
                        mpz_divexact_ui(remaining, remaining, p);
                        exp2++;
                    }
                    if (exp2 & 1) exponents[nexp++] = fi;
                }

                // Check cofactor
                uint64_t cofactor_val = 0;
                if (mpz_fits_ulong_p(remaining)) {
                    cofactor_val = mpz_get_ui(remaining);
                }

                if (cofactor_val == 1) {
                    // Full smooth (with q as guaranteed LP) = DLP relation!
                    // Store relation with LP = q
                    if (nrels < MAX_RELS) {
                        rels[nrels].exponents = malloc(nexp * sizeof(int));
                        memcpy(rels[nrels].exponents, exponents, nexp * sizeof(int));
                        rels[nrels].nexp = nexp;
                        mpz_init(rels[nrels].x_val);
                        mpz_mul_si(rels[nrels].x_val, A, x);
                        mpz_add(rels[nrels].x_val, rels[nrels].x_val, B);
                        rels[nrels].lp1 = q;
                        rels[nrels].lp2 = 0;
                        rels[nrels].sign = sign;
                        nrels++;
                        poly_count++;
                    }
                } else if (cofactor_val > 1 && cofactor_val <= lp_bound) {
                    // DLP relation: q and cofactor_val are two large primes
                    if (nrels < MAX_RELS) {
                        rels[nrels].exponents = malloc(nexp * sizeof(int));
                        memcpy(rels[nrels].exponents, exponents, nexp * sizeof(int));
                        rels[nrels].nexp = nexp;
                        mpz_init(rels[nrels].x_val);
                        mpz_mul_si(rels[nrels].x_val, A, x);
                        mpz_add(rels[nrels].x_val, rels[nrels].x_val, B);
                        rels[nrels].lp1 = q;
                        rels[nrels].lp2 = (uint32_t)cofactor_val;
                        rels[nrels].sign = sign;
                        nrels++;
                        poly_count++;
                    }
                }

                mpz_clear(remaining);
            }
        }

        if (sq_idx % 100 == 0) {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - t0.tv_sec) + (now.tv_nsec - t0.tv_nsec) * 1e-9;
            fprintf(stderr, "  sq=%d/%d: %d/%d rels, %.1f rels/s, %.1fs\n",
                    sq_idx, nsq, nrels, target,
                    elapsed > 0 ? nrels / elapsed : 0, elapsed);
        }
    }

    struct timespec t1;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;

    fprintf(stderr, "Special-Q sieve: %d rels in %.1fs (%d special-Q primes used)\n",
            nrels, total, sq_idx);

    if (nrels < fb_count + 1) {
        fprintf(stderr, "FAIL: insufficient relations (%d < %d)\n", nrels, fb_count + 1);
        return 1;
    }

    fprintf(stderr, "TODO: Linear algebra not yet implemented\n");
    fprintf(stderr, "Have %d relations over %d-element factor base\n", nrels, fb_count);

    // Cleanup
    free(sq_primes);
    free(sq_sqrt);
    mpz_clears(N, kN, A, B, Qval, x_val, factor, cofactor, temp, NULL);

    return 1; // Not yet complete
}
