/*
 * lattice_factor.c - Lattice-based factoring for balanced semiprimes
 *
 * For N = p*q where p ≈ q ≈ √N, we know that:
 *   p = floor(√N) + a,  q = floor(√N) + b  (approximately)
 * where a + b ≈ 0 and a*b ≈ N - floor(√N)²
 *
 * Strategy 1: Enhanced Fermat with lattice search
 *   N = ((p+q)/2)² - ((p-q)/2)²
 *   Search for s = (p+q)/2 near √N using multiple modular constraints
 *
 * Strategy 2: Coppersmith-inspired small root finding
 *   f(x) = (√N + x)² - N has a root at x = (p - √N) which is "small"
 *   relative to N. Use LLL to find this root.
 *
 * Strategy 3: Multivariate Coppersmith
 *   f(x,y) = (√N + x)(√N + y) - N = 0
 *   Has small roots (x,y) = (p-√N, q-√N)
 *
 * Strategy 4: Chinese Remainder + lattice
 *   For many small primes r, compute p mod r candidates
 *   Use CRT + lattice reduction to reconstruct p
 *
 * Compile: gcc -O3 -march=native -o lattice_factor library/lattice_factor.c -lgmp -lm
 * Usage: ./lattice_factor <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ==================== Strategy 1: Enhanced Fermat ==================== */
/*
 * Standard Fermat: try s = ceil(√N), ceil(√N)+1, ... until s²-N is a perfect square
 * Enhancement: use modular constraints to skip non-candidates
 * For each small prime p, s² ≡ N (mod p), so s must be in specific residue classes
 * Use sieve of these constraints to skip most candidates
 */

// Check if n is a perfect square, if so set root
static int is_perfect_square(mpz_t n, mpz_t root) {
    if (mpz_sgn(n) < 0) return 0;
    mpz_sqrt(root, n);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul(tmp, root, root);
    int result = (mpz_cmp(tmp, n) == 0);
    mpz_clear(tmp);
    return result;
}

// Small primes for sieve
static const int SMALL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113
};
#define NUM_SMALL_PRIMES 29

// Modular square root via Tonelli-Shanks for small primes
static int modsqrt(int n, int p) {
    n = ((n % p) + p) % p;
    if (n == 0) return 0;
    if (p == 2) return n;

    // Check if n is QR mod p
    int pow = 1;
    for (int i = 0; i < (p-1)/2; i++) pow = (int)((long long)pow * n % p);
    if (pow != 1) return -1; // not a QR

    if (p % 4 == 3) {
        pow = 1;
        int exp = (p + 1) / 4;
        int base = n;
        while (exp > 0) {
            if (exp & 1) pow = (int)((long long)pow * base % p);
            base = (int)((long long)base * base % p);
            exp >>= 1;
        }
        return pow;
    }

    // General Tonelli-Shanks
    int q = p - 1, s = 0;
    while (q % 2 == 0) { q /= 2; s++; }

    int z = 2;
    while (1) {
        pow = 1;
        for (int i = 0; i < (p-1)/2; i++) pow = (int)((long long)pow * z % p);
        if (pow == p - 1) break;
        z++;
    }

    int m = s;
    int c = 1;
    { int base = z, exp = q;
      while (exp > 0) {
        if (exp & 1) c = (int)((long long)c * base % p);
        base = (int)((long long)base * base % p);
        exp >>= 1;
      }
    }
    int t = 1;
    { int base = n, exp = q;
      while (exp > 0) {
        if (exp & 1) t = (int)((long long)t * base % p);
        base = (int)((long long)base * base % p);
        exp >>= 1;
      }
    }
    int r = 1;
    { int base = n, exp = (q + 1) / 2;
      while (exp > 0) {
        if (exp & 1) r = (int)((long long)r * base % p);
        base = (int)((long long)base * base % p);
        exp >>= 1;
      }
    }

    while (1) {
        if (t == 1) return r;
        int i = 0;
        int tmp = t;
        while (tmp != 1) {
            tmp = (int)((long long)tmp * tmp % p);
            i++;
        }
        int b = c;
        for (int j = 0; j < m - i - 1; j++) b = (int)((long long)b * b % p);
        m = i;
        c = (int)((long long)b * b % p);
        t = (int)((long long)t * c % p);
        r = (int)((long long)r * b % p);
    }
}

int fermat_enhanced(mpz_t N, mpz_t factor) {
    mpz_t s, s2, diff, root;
    mpz_inits(s, s2, diff, root, NULL);

    // s = ceil(√N)
    mpz_sqrt(s, N);
    mpz_mul(s2, s, s);
    if (mpz_cmp(s2, N) < 0) mpz_add_ui(s, s, 1);

    // Compute allowed residues for s modulo small primes
    // s² ≡ N (mod p) means s ≡ ±√(N mod p) (mod p)
    unsigned char *sieve = NULL;
    int sieve_mod = 1;
    int prime_list[NUM_SMALL_PRIMES];
    int roots1[NUM_SMALL_PRIMES], roots2[NUM_SMALL_PRIMES];
    int num_sieve_primes = 0;

    // Use product of first few primes as sieve modulus
    // We'll use CRT-style sieve
    for (int i = 0; i < NUM_SMALL_PRIMES && SMALL_PRIMES[i] <= 113; i++) {
        int p = SMALL_PRIMES[i];
        int n_mod_p = (int)mpz_fdiv_ui(N, p);
        int r = modsqrt(n_mod_p, p);
        if (r < 0) {
            // N is not QR mod p, so p | N or no solution
            // Actually this means s² can never be ≡ N mod p
            // So we can skip s values where s² ≡ N (mod p) fails
            // But since N is a semiprime, if p doesn't divide N,
            // then N is not QR mod p means s²-N is never 0 mod p.
            // This is fine - we just note no constraint from this prime.
            prime_list[num_sieve_primes] = p;
            roots1[num_sieve_primes] = -1; // no valid residue
            roots2[num_sieve_primes] = -1;
            num_sieve_primes++;
            continue;
        }
        prime_list[num_sieve_primes] = p;
        roots1[num_sieve_primes] = r;
        roots2[num_sieve_primes] = (p - r) % p;
        num_sieve_primes++;
    }

    // Combined sieve using product of first few primes
    // Use primes up to 31 for combined modulus: 3*5*7*11*13*17*19*23*29*31 = 100280245065
    // Too large. Use smaller set: 3*5*7*11*13*17*19*23 = 223092870
    long long combined_mod = 1;
    int num_combined = 0;
    for (int i = 0; i < num_sieve_primes; i++) {
        if (combined_mod * prime_list[i] > 200000000LL) break;
        combined_mod *= prime_list[i];
        num_combined = i + 1;
    }

    // Build sieve: for each residue 0..combined_mod-1, check if it's valid
    // A residue r is valid if for each prime p in our set:
    //   r ≡ roots1[i] or roots2[i] (mod prime_list[i])
    // Valid residues are candidates for s mod combined_mod

    // For a combined modulus of ~200M, we can't enumerate all residues
    // Instead, use CRT to generate valid residues
    // Number of valid residues: ∏(2 or 0 or 1 per prime) ≈ 2^8 = 256 per combined_mod

    // Actually, let's use a simpler approach: sieve with individual primes
    // Check each candidate s against all primes

    long long limit = 1000000000LL; // 10^9 candidates max
    long long s_start_ui = 0; // offset from ceil(√N)

    mpz_t s_base;
    mpz_init_set(s_base, s); // s_base = ceil(√N)

    int found = 0;
    for (long long offset = 0; offset < limit; offset++) {
        // Check modular constraints
        int valid = 1;
        // Quick check: s = s_base + offset
        // s mod p = (s_base mod p + offset mod p) mod p
        for (int i = 0; i < num_sieve_primes; i++) {
            int p = prime_list[i];
            if (roots1[i] == -1) continue; // N not QR mod p, skip
            int s_mod_p = (int)((mpz_fdiv_ui(s_base, p) + (offset % p)) % p);
            if (s_mod_p != roots1[i] && s_mod_p != roots2[i]) {
                valid = 0;
                break;
            }
        }
        if (!valid) continue;

        // This candidate passes modular screening
        mpz_add_ui(s, s_base, offset);
        mpz_mul(s2, s, s);
        mpz_sub(diff, s2, N);

        if (is_perfect_square(diff, root)) {
            // s² - N = root² => N = (s-root)(s+root)
            mpz_sub(factor, s, root);
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                found = 1;
                break;
            }
        }
    }

    mpz_clears(s, s2, diff, root, s_base, NULL);
    return found;
}

/* ==================== Strategy 2: Dixon/Lehman hybrid ==================== */
/*
 * Lehman's method: for balanced semiprimes with |p-q| < N^(1/3),
 * searches for (a, b) with 4aN = (2a*s + b)² for small a and b.
 *
 * We extend this by using modular information to reduce the search space.
 * For each k = 1, 2, ..., we look for s near √(kN) such that
 * 4kN - s² is a perfect square times something small.
 */

int lehman_enhanced(mpz_t N, mpz_t factor, int k_max) {
    mpz_t kN, s, s2, diff, root, four_kN, g;
    mpz_inits(kN, s, s2, diff, root, four_kN, g, NULL);

    for (int k = 1; k <= k_max; k++) {
        mpz_mul_ui(kN, N, k);
        mpz_mul_ui(four_kN, kN, 4);

        // s starts at ceil(2*√(kN))
        mpz_sqrt(s, four_kN);
        mpz_mul(s2, s, s);
        if (mpz_cmp(s2, four_kN) < 0) mpz_add_ui(s, s, 1);

        // Search range: from s to s + ceil(N^(1/6) / (4√k))
        // For balanced semiprimes this range should be small
        mpz_t range;
        mpz_init(range);
        mpz_root(range, N, 6); // N^(1/6)
        unsigned long range_ul = mpz_get_ui(range);
        if (range_ul > 10000000UL) range_ul = 10000000UL;
        range_ul = range_ul / (2 * (unsigned long)ceil(sqrt(k))) + 1;

        for (unsigned long j = 0; j <= range_ul; j++) {
            mpz_mul(s2, s, s);
            mpz_sub(diff, s2, four_kN);

            if (mpz_sgn(diff) >= 0 && is_perfect_square(diff, root)) {
                // 4kN = s² - root² = (s-root)(s+root)
                mpz_add(g, s, root);
                mpz_gcd(factor, g, N);
                if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                    mpz_clears(kN, s, s2, diff, root, four_kN, g, range, NULL);
                    return 1;
                }
                mpz_sub(g, s, root);
                mpz_gcd(factor, g, N);
                if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                    mpz_clears(kN, s, s2, diff, root, four_kN, g, range, NULL);
                    return 1;
                }
            }
            mpz_add_ui(s, s, 1);
        }
        mpz_clear(range);
    }

    mpz_clears(kN, s, s2, diff, root, four_kN, g, NULL);
    return 0;
}

/* ==================== Strategy 3: Smooth congruence via lattice ==================== */
/*
 * Novel approach: Generate congruences x² ≡ y (mod N) where y is small,
 * using lattice reduction.
 *
 * Build a lattice where short vectors correspond to (x, y) pairs with
 * x² - y ≡ 0 (mod N) and |y| is small.
 *
 * Lattice basis:
 *   B = [ N   0 ]
 *       [ a   1 ]
 * where a ≈ √N (mod N).
 *
 * Short vectors (v1, v2) of this lattice satisfy v1 ≡ a*v2 (mod N).
 * If v1 = x, v2 = 1, then x ≡ a (mod N) and x² ≡ a² ≡ N (mod N).
 *
 * But we need x² - kN = y where y is smooth. The lattice gives us
 * small |x| values where x² mod N is also small.
 *
 * Extended: use a higher-dimensional lattice with multiple "small root" vectors.
 */

/* Simple 2x2 lattice reduction (Gauss reduction) */
typedef struct {
    mpz_t a11, a12, a21, a22;
} lattice2;

void lattice2_init(lattice2 *L) {
    mpz_inits(L->a11, L->a12, L->a21, L->a22, NULL);
}

void lattice2_clear(lattice2 *L) {
    mpz_clears(L->a11, L->a12, L->a21, L->a22, NULL);
}

// Gauss lattice reduction for 2x2 lattice
void gauss_reduce(lattice2 *L) {
    mpz_t dot, norm1, norm2, q, tmp;
    mpz_inits(dot, norm1, norm2, q, tmp, NULL);

    for (int iter = 0; iter < 1000; iter++) {
        // norm1 = |v1|² = a11² + a12²
        mpz_mul(norm1, L->a11, L->a11);
        mpz_mul(tmp, L->a12, L->a12);
        mpz_add(norm1, norm1, tmp);

        // norm2 = |v2|² = a21² + a22²
        mpz_mul(norm2, L->a21, L->a21);
        mpz_mul(tmp, L->a22, L->a22);
        mpz_add(norm2, norm2, tmp);

        // Ensure v1 is shorter
        if (mpz_cmp(norm1, norm2) > 0) {
            mpz_swap(L->a11, L->a21);
            mpz_swap(L->a12, L->a22);
            mpz_swap(norm1, norm2);
        }

        // dot = v1 · v2
        mpz_mul(dot, L->a11, L->a21);
        mpz_mul(tmp, L->a12, L->a22);
        mpz_add(dot, dot, tmp);

        // q = round(dot / norm1)
        // q = (2*dot + norm1) / (2*norm1)  for rounding
        mpz_mul_ui(tmp, dot, 2);
        if (mpz_sgn(dot) >= 0)
            mpz_add(tmp, tmp, norm1);
        else
            mpz_sub(tmp, tmp, norm1);
        mpz_mul_ui(norm1, norm1, 2);
        mpz_tdiv_q(q, tmp, norm1);
        mpz_divexact_ui(norm1, norm1, 2);

        if (mpz_sgn(q) == 0) break;

        // v2 = v2 - q * v1
        mpz_mul(tmp, q, L->a11);
        mpz_sub(L->a21, L->a21, tmp);
        mpz_mul(tmp, q, L->a12);
        mpz_sub(L->a22, L->a22, tmp);
    }

    mpz_clears(dot, norm1, norm2, q, tmp, NULL);
}

/*
 * Lattice-based smooth congruence search:
 * For various multipliers k, build lattice:
 *   v1 = (N, 0)
 *   v2 = (floor(√(kN)), 1) * scaling
 * After reduction, short vectors give (x, c) with x ≈ c * √(kN) (mod N)
 * so x² ≈ c² * kN (mod N) and x² mod N ≈ c² * kN - q*N for small q
 */

int lattice_smooth_factor(mpz_t N, mpz_t factor, int time_limit_sec) {
    time_t start = time(NULL);

    mpz_t sqrtN, x, x2, residue, g;
    mpz_inits(sqrtN, x, x2, residue, g, NULL);
    mpz_sqrt(sqrtN, N);

    // Try various scaling factors and multipliers
    for (int k = 1; k <= 10000 && time(NULL) - start < time_limit_sec; k++) {
        lattice2 L;
        lattice2_init(&L);

        // Lattice: v1 = (N, 0), v2 = (floor(√(kN)) mod N, k)
        mpz_t kN, sqrtkN;
        mpz_inits(kN, sqrtkN, NULL);
        mpz_mul_ui(kN, N, k);
        mpz_sqrt(sqrtkN, kN);

        // Scale: we want both components comparable
        // v1 = (N, 0), v2 = (sqrtkN, 1)
        // After reduction, short vector (a, b) satisfies a ≡ b*sqrtkN (mod N)
        // So a² ≡ b²*kN (mod N)
        // If b is small, a² - b²*kN = qN for some q, and |a² - b²*kN| is small

        mpz_set(L.a11, N);
        mpz_set_ui(L.a12, 0);
        mpz_mod(L.a21, sqrtkN, N);
        mpz_set_ui(L.a22, 1);

        gauss_reduce(&L);

        // Check both short vectors
        for (int vec = 0; vec < 2; vec++) {
            mpz_t *ax = (vec == 0) ? &L.a11 : &L.a21;
            mpz_t *bx = (vec == 0) ? &L.a12 : &L.a22;

            // x = a, coefficient = b
            // x² - b²*kN ≡ 0 (mod N)
            // gcd(x - b*sqrtkN, N) or gcd(x + b*sqrtkN, N) might give factor

            if (mpz_sgn(*bx) == 0) continue;

            mpz_mul(x, *bx, sqrtkN);
            mpz_sub(g, *ax, x);
            mpz_gcd(g, g, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                lattice2_clear(&L);
                mpz_clears(kN, sqrtkN, sqrtN, x, x2, residue, g, NULL);
                return 1;
            }

            mpz_add(g, *ax, x);
            mpz_gcd(g, g, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                lattice2_clear(&L);
                mpz_clears(kN, sqrtkN, sqrtN, x, x2, residue, g, NULL);
                return 1;
            }
        }

        lattice2_clear(&L);
        mpz_clears(kN, sqrtkN, NULL);
    }

    mpz_clears(sqrtN, x, x2, residue, g, NULL);
    return 0;
}

/* ==================== Strategy 4: Smooth number accumulation ==================== */
/*
 * Novel idea: Instead of sieving, generate random x values and compute
 * x² mod N. Accumulate partial factorizations using a "factor tree":
 *
 * 1. For many x values, compute r = x² mod N
 * 2. For each pair (r_i, r_j), compute gcd(r_i, r_j)
 * 3. Common factors reveal shared prime factors
 * 4. Build factor base dynamically from discovered primes
 * 5. Once enough relations found, do linear algebra
 *
 * This avoids the traditional sieve entirely and may discover
 * structure in the residues.
 */

/* ==================== Main ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N, factor, cofactor;
    mpz_inits(N, factor, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    if (mpz_cmp_ui(N, 4) < 0) {
        gmp_printf("N=%Zd is too small\n", N);
        mpz_clears(N, factor, cofactor, NULL);
        return 1;
    }

    // Check for small factors first
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_set_ui(factor, p);
            mpz_divexact(cofactor, N, factor);
            gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
            mpz_clears(N, factor, cofactor, NULL);
            return 0;
        }
    }

    int digits = (int)(mpz_sizeinbase(N, 10));
    printf("Factoring %d-digit number\n", digits);

    time_t start = time(NULL);

    // Strategy 1: Enhanced Fermat (fast for close factors)
    printf("Trying enhanced Fermat...\n");
    if (fermat_enhanced(N, factor)) {
        mpz_divexact(cofactor, N, factor);
        double elapsed = difftime(time(NULL), start);
        gmp_printf("%Zd = %Zd * %Zd (Fermat, %.1fs)\n", N, factor, cofactor, elapsed);
        mpz_clears(N, factor, cofactor, NULL);
        return 0;
    }
    printf("Fermat: no close factors found (%.0fs)\n", difftime(time(NULL), start));

    // Strategy 2: Lehman enhanced
    printf("Trying enhanced Lehman...\n");
    int k_max = 1000000; // search up to k=10^6
    if (lehman_enhanced(N, factor, k_max)) {
        mpz_divexact(cofactor, N, factor);
        double elapsed = difftime(time(NULL), start);
        gmp_printf("%Zd = %Zd * %Zd (Lehman k<=%d, %.1fs)\n", N, factor, cofactor, k_max, elapsed);
        mpz_clears(N, factor, cofactor, NULL);
        return 0;
    }
    printf("Lehman: no factors found with k<=%d (%.0fs)\n", k_max, difftime(time(NULL), start));

    // Strategy 3: Lattice smooth congruence
    printf("Trying lattice-based search...\n");
    if (lattice_smooth_factor(N, factor, 60)) {
        mpz_divexact(cofactor, N, factor);
        double elapsed = difftime(time(NULL), start);
        gmp_printf("%Zd = %Zd * %Zd (Lattice, %.1fs)\n", N, factor, cofactor, elapsed);
        mpz_clears(N, factor, cofactor, NULL);
        return 0;
    }
    printf("Lattice: no factors found (%.0fs)\n", difftime(time(NULL), start));

    printf("All strategies failed after %.0fs\n", difftime(time(NULL), start));
    mpz_clears(N, factor, cofactor, NULL);
    return 1;
}
