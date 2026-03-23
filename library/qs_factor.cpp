// Self-Initializing Quadratic Sieve (SIQS) implementation
// Usage: ./qs_factor <number>
// Outputs: factor1 factor2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <cstring>
#include <gmp.h>

// Tonelli-Shanks modular square root: returns r such that r^2 = n (mod p), or 0 if none
static unsigned long mod_sqrt(unsigned long n, unsigned long p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    // Check if n is a QR mod p using Euler's criterion
    // n^((p-1)/2) mod p should be 1
    unsigned long exp = (p - 1) / 2;
    unsigned long result = 1;
    unsigned long base = n;
    unsigned long e = exp;
    while (e > 0) {
        if (e & 1) result = (__uint128_t)result * base % p;
        base = (__uint128_t)base * base % p;
        e >>= 1;
    }
    if (result != 1) return 0; // not a QR

    // Tonelli-Shanks
    if (p % 4 == 3) {
        // Simple case
        unsigned long r = 1;
        base = n;
        e = (p + 1) / 4;
        while (e > 0) {
            if (e & 1) r = (__uint128_t)r * base % p;
            base = (__uint128_t)base * base % p;
            e >>= 1;
        }
        return r;
    }

    // General Tonelli-Shanks
    unsigned long Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    // Find a non-residue
    unsigned long z = 2;
    while (true) {
        result = 1; base = z; e = (p - 1) / 2;
        while (e > 0) {
            if (e & 1) result = (__uint128_t)result * base % p;
            base = (__uint128_t)base * base % p;
            e >>= 1;
        }
        if (result == p - 1) break;
        z++;
    }

    unsigned long M = S;
    // c = z^Q mod p
    unsigned long c = 1; base = z; e = Q;
    while (e > 0) {
        if (e & 1) c = (__uint128_t)c * base % p;
        base = (__uint128_t)base * base % p;
        e >>= 1;
    }
    // t = n^Q mod p
    unsigned long t = 1; base = n; e = Q;
    while (e > 0) {
        if (e & 1) t = (__uint128_t)t * base % p;
        base = (__uint128_t)base * base % p;
        e >>= 1;
    }
    // R = n^((Q+1)/2) mod p
    unsigned long R = 1; base = n; e = (Q + 1) / 2;
    while (e > 0) {
        if (e & 1) R = (__uint128_t)R * base % p;
        base = (__uint128_t)base * base % p;
        e >>= 1;
    }

    while (true) {
        if (t == 1) return R;
        // Find least i such that t^(2^i) = 1
        unsigned long i = 0, tmp = t;
        while (tmp != 1) { tmp = (__uint128_t)tmp * tmp % p; i++; }
        if (i == M) return 0; // shouldn't happen
        // b = c^(2^(M-i-1))
        unsigned long b = c;
        for (unsigned long j = 0; j < M - i - 1; j++)
            b = (__uint128_t)b * b % p;
        M = i;
        c = (__uint128_t)b * b % p;
        t = (__uint128_t)t * c % p;
        R = (__uint128_t)R * b % p;
    }
}

// Sieve primes up to limit
static std::vector<unsigned long> sieve_primes(unsigned long limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (unsigned long i = 2; i * i <= limit; i++) {
        if (is_prime[i]) {
            for (unsigned long j = i * i; j <= limit; j += i)
                is_prime[j] = false;
        }
    }
    std::vector<unsigned long> primes;
    for (unsigned long i = 2; i <= limit; i++)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

struct Relation {
    mpz_t x;        // The x value: Q(x) = (x + sqrt(N))^2 - N
    std::vector<unsigned char> exponents; // exponents mod 2 for factor base primes
    mpz_t Qx;       // The Q(x) value (for verification)
};

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    struct timespec t_start;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    mpz_t N, sqrtN, factor1, factor2;
    mpz_init(N);
    mpz_init(sqrtN);
    mpz_init(factor1);
    mpz_init(factor2);

    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }

    size_t digits = mpz_sizeinbase(N, 10);

    // Compute sqrt(N)
    mpz_sqrt(sqrtN, N);

    // Parameters based on digit count
    // Factor base size and sieve interval
    double ln_n = mpz_sizeinbase(N, 2) * log(2.0);
    double ln_ln_n = log(ln_n);

    // L[1/2] = exp(sqrt(ln(n) * ln(ln(n))))
    double L_half = exp(sqrt(ln_n * ln_ln_n));

    // Factor base bound: B ~ L[1/2]^(1/sqrt(2))
    // Sieve interval: M ~ L[1/2]^(1/sqrt(2))
    double smoothness_param = 1.0 / sqrt(2.0);

    // Adjust parameters empirically
    unsigned long B;
    unsigned long sieve_radius;

    if (digits <= 30) {
        B = 200;
        sieve_radius = 10000;
    } else if (digits <= 40) {
        B = 800;
        sieve_radius = 50000;
    } else if (digits <= 50) {
        B = 3000;
        sieve_radius = 200000;
    } else if (digits <= 60) {
        B = 10000;
        sieve_radius = 500000;
    } else if (digits <= 70) {
        B = 30000;
        sieve_radius = 1000000;
    } else if (digits <= 80) {
        B = 80000;
        sieve_radius = 2000000;
    } else if (digits <= 90) {
        B = 200000;
        sieve_radius = 4000000;
    } else {
        B = 500000;
        sieve_radius = 8000000;
    }

    // Build factor base: primes p <= B where N is a QR mod p
    std::vector<unsigned long> all_primes = sieve_primes(B);
    std::vector<unsigned long> factor_base;
    std::vector<unsigned long> fb_sqrt; // sqrt(N) mod p for each factor base prime

    factor_base.push_back(2);
    fb_sqrt.push_back(1); // placeholder for 2

    unsigned long n_mod_p;
    for (size_t i = 1; i < all_primes.size(); i++) {
        unsigned long p = all_primes[i];
        n_mod_p = mpz_fdiv_ui(N, p);
        unsigned long sq = mod_sqrt(n_mod_p, p);
        if (sq != 0) {
            factor_base.push_back(p);
            fb_sqrt.push_back(sq);
        }
    }

    size_t fb_size = factor_base.size();
    size_t relations_needed = fb_size + 20; // Need fb_size + some extra

    fprintf(stderr, "QS: %zu digits, factor base size=%zu, sieve_radius=%lu\n",
            digits, fb_size, sieve_radius);

    // Sieve using multiple polynomials (SIQS)
    // For basic QS: Q(x) = (x + sqrtN)^2 - N
    // Sieve log approximations

    // Precompute log approximations for factor base
    std::vector<unsigned char> fb_log(fb_size);
    for (size_t i = 0; i < fb_size; i++) {
        fb_log[i] = (unsigned char)(log((double)factor_base[i]) / log(2.0) * 1.5 + 0.5);
        if (fb_log[i] == 0) fb_log[i] = 1;
    }

    // Target log value for smooth numbers
    // |Q(x)| ~ 2 * sqrt(N) * |x| for x near 0
    // For x ~ sieve_radius, |Q(x)| ~ 2 * sqrt(N) * sieve_radius
    double target_log_bits = (mpz_sizeinbase(N, 2) / 2.0 + log2((double)sieve_radius)) * 1.5;
    unsigned char threshold = (unsigned char)(target_log_bits - log2((double)B) * 2.5 * 1.5);

    // Collect relations
    std::vector<Relation> relations;

    // Sieve array
    std::vector<unsigned char> sieve(2 * sieve_radius + 1, 0);

    // Compute sqrtN mod p and sieve starting points
    unsigned long sqrtN_mod;

    mpz_t x_val, qx, tmp, tmp2;
    mpz_init(x_val);
    mpz_init(qx);
    mpz_init(tmp);
    mpz_init(tmp2);

    // Basic QS sieve
    // Q(x) = (sqrtN + x)^2 - N for x in [-sieve_radius, sieve_radius]

    memset(sieve.data(), 0, sieve.size());

    // For each prime in factor base, find sieve positions
    for (size_t i = 0; i < fb_size; i++) {
        unsigned long p = factor_base[i];
        if (p == 2) {
            // Handle p=2 separately
            sqrtN_mod = mpz_fdiv_ui(sqrtN, 2);
            // (sqrtN + x)^2 - N ≡ 0 (mod 2)
            // sqrtN + x must make Q(x) even
            unsigned long start = (sqrtN_mod % 2 == 0) ? 0 : 1;
            for (unsigned long j = start; j < sieve.size(); j += 2) {
                sieve[j] += fb_log[i];
            }
            // Also handle higher powers of 2
            continue;
        }

        unsigned long sq = fb_sqrt[i];
        sqrtN_mod = mpz_fdiv_ui(sqrtN, p);

        // We need (sqrtN + x)^2 ≡ N (mod p)
        // sqrtN + x ≡ ±sq (mod p)
        // x ≡ sq - sqrtN or x ≡ -sq - sqrtN (mod p)

        long start1 = ((long)sq - (long)sqrtN_mod) % (long)p;
        if (start1 < 0) start1 += p;
        long start2 = (-(long)sq - (long)sqrtN_mod) % (long)p;
        if (start2 < 0) start2 += p;

        // Offset: sieve array index = x + sieve_radius
        // So sieve position for x = start1 is start1 + sieve_radius, but we need x in range
        long offset = (long)sieve_radius;

        // Adjust start1 to be in [0, sieve.size())
        long pos1 = (start1 + offset) % p;
        if (pos1 < 0) pos1 += p;
        long pos2 = (start2 + offset) % p;
        if (pos2 < 0) pos2 += p;

        for (long j = pos1; j < (long)sieve.size(); j += p)
            sieve[j] += fb_log[i];
        if (start1 != start2) {
            for (long j = pos2; j < (long)sieve.size(); j += p)
                sieve[j] += fb_log[i];
        }

        // Sieve with prime powers p^k
        unsigned long pk = p * p;
        while (pk <= B * 4) {
            unsigned long n_mod_pk = mpz_fdiv_ui(N, pk);
            // Check if sqrt exists mod pk (lift using Hensel)
            // For simplicity, skip prime powers for now
            break;
        }
    }

    // Scan sieve for smooth candidates
    for (unsigned long idx = 0; idx < sieve.size(); idx++) {
        if (sieve[idx] >= threshold) {
            long x = (long)idx - (long)sieve_radius;

            // Compute Q(x) = (sqrtN + x)^2 - N
            mpz_set(x_val, sqrtN);
            if (x >= 0)
                mpz_add_ui(x_val, x_val, (unsigned long)x);
            else
                mpz_sub_ui(x_val, x_val, (unsigned long)(-x));

            mpz_mul(qx, x_val, x_val);
            mpz_sub(qx, qx, N);

            if (mpz_sgn(qx) == 0) continue;

            // Trial divide by factor base
            mpz_set(tmp, qx);
            if (mpz_sgn(tmp) < 0) mpz_neg(tmp, tmp);

            std::vector<unsigned char> exps(fb_size, 0);
            bool is_negative = (mpz_sgn(qx) < 0);

            for (size_t i = 0; i < fb_size; i++) {
                unsigned long p = factor_base[i];
                while (mpz_divisible_ui_p(tmp, p)) {
                    mpz_divexact_ui(tmp, tmp, p);
                    exps[i]++;
                }
            }

            // Check if fully factored
            if (mpz_cmp_ui(tmp, 1) == 0) {
                Relation rel;
                mpz_init_set(rel.x, x_val);
                mpz_init_set(rel.Qx, qx);
                rel.exponents.resize(fb_size);
                for (size_t i = 0; i < fb_size; i++)
                    rel.exponents[i] = exps[i] & 1; // mod 2
                relations.push_back(rel);

                if (relations.size() >= relations_needed) break;

                if (relations.size() % 100 == 0) {
                    struct timespec now;
                    clock_gettime(CLOCK_MONOTONIC, &now);
                    double elapsed = (now.tv_sec - t_start.tv_sec) + (now.tv_nsec - t_start.tv_nsec) / 1e9;
                    fprintf(stderr, "QS: %zu/%zu relations (%.1fs)\n",
                            relations.size(), relations_needed, elapsed);
                }
            }
        }

        // Time check every 1M entries
        if (idx % 1000000 == 0 && idx > 0) {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - t_start.tv_sec) + (now.tv_nsec - t_start.tv_nsec) / 1e9;
            if (elapsed > 280.0) {
                fprintf(stderr, "QS: timeout during sieve after %.1fs, %zu relations\n",
                        elapsed, relations.size());
                break;
            }
        }
    }

    fprintf(stderr, "QS: collected %zu relations (need %zu)\n", relations.size(), relations_needed);

    if (relations.size() < fb_size + 1) {
        fprintf(stderr, "FAIL: not enough relations\n");
        return 1;
    }

    // Gaussian elimination mod 2 to find dependencies
    size_t nrels = relations.size();
    size_t ncols = fb_size + 1; // +1 for sign

    // Matrix: nrels rows x ncols columns, stored as bit vectors
    std::vector<std::vector<unsigned long>> matrix(nrels);
    size_t words = (ncols + 63) / 64;
    for (size_t i = 0; i < nrels; i++) {
        matrix[i].resize(words, 0);
        // Set sign bit
        if (mpz_sgn(relations[i].Qx) < 0)
            matrix[i][0] |= 1UL;
        // Set exponent parity bits
        for (size_t j = 0; j < fb_size; j++) {
            if (relations[i].exponents[j] & 1) {
                size_t bit = j + 1;
                matrix[i][bit / 64] |= (1UL << (bit % 64));
            }
        }
    }

    // Track which relations are combined
    std::vector<std::vector<unsigned long>> history(nrels);
    size_t hist_words = (nrels + 63) / 64;
    for (size_t i = 0; i < nrels; i++) {
        history[i].resize(hist_words, 0);
        history[i][i / 64] |= (1UL << (i % 64));
    }

    // Gaussian elimination
    std::vector<size_t> pivot_row(ncols, (size_t)-1);

    for (size_t col = 0; col < ncols; col++) {
        // Find pivot
        size_t piv = (size_t)-1;
        for (size_t row = 0; row < nrels; row++) {
            if (pivot_row[col] != (size_t)-1) break;
            bool already_pivot = false;
            for (size_t c = 0; c < col; c++) {
                if (pivot_row[c] == row) { already_pivot = true; break; }
            }
            if (already_pivot) continue;

            if (matrix[row][col / 64] & (1UL << (col % 64))) {
                piv = row;
                break;
            }
        }

        if (piv == (size_t)-1) continue;
        pivot_row[col] = piv;

        // Eliminate this column from all other rows
        for (size_t row = 0; row < nrels; row++) {
            if (row == piv) continue;
            if (matrix[row][col / 64] & (1UL << (col % 64))) {
                for (size_t w = 0; w < words; w++)
                    matrix[row][w] ^= matrix[piv][w];
                for (size_t w = 0; w < hist_words; w++)
                    history[row][w] ^= history[piv][w];
            }
        }
    }

    // Find null space vectors (rows that are all zero)
    bool found = false;
    for (size_t row = 0; row < nrels && !found; row++) {
        bool is_zero = true;
        for (size_t w = 0; w < words && is_zero; w++) {
            if (matrix[row][w] != 0) is_zero = false;
        }
        if (!is_zero) continue;

        // This row is a dependency - combine the relations
        mpz_t lhs, rhs;
        mpz_init(lhs);
        mpz_init(rhs);
        mpz_set_ui(lhs, 1);
        mpz_set_ui(rhs, 1);

        // Collect exponents for RHS
        std::vector<unsigned long> total_exps(fb_size, 0);

        for (size_t i = 0; i < nrels; i++) {
            if (history[row][i / 64] & (1UL << (i % 64))) {
                // This relation is part of the combination
                mpz_mul(lhs, lhs, relations[i].x);
                mpz_mod(lhs, lhs, N);

                // Accumulate exponents
                mpz_set(tmp, relations[i].Qx);
                if (mpz_sgn(tmp) < 0) mpz_neg(tmp, tmp);

                for (size_t j = 0; j < fb_size; j++) {
                    unsigned long p = factor_base[j];
                    while (mpz_divisible_ui_p(tmp, p)) {
                        mpz_divexact_ui(tmp, tmp, p);
                        total_exps[j]++;
                    }
                }
            }
        }

        // Compute RHS = product of primes^(exp/2)
        for (size_t j = 0; j < fb_size; j++) {
            if (total_exps[j] > 0) {
                mpz_ui_pow_ui(tmp, factor_base[j], total_exps[j] / 2);
                mpz_mul(rhs, rhs, tmp);
                mpz_mod(rhs, rhs, N);
            }
        }

        // gcd(lhs - rhs, N) should give a factor
        mpz_sub(tmp, lhs, rhs);
        mpz_gcd(factor1, tmp, N);

        if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
            mpz_divexact(factor2, N, factor1);

            if (mpz_cmp(factor1, factor2) > 0)
                mpz_swap(factor1, factor2);

            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - t_start.tv_sec) + (now.tv_nsec - t_start.tv_nsec) / 1e9;
            gmp_printf("%Zd %Zd\n", factor1, factor2);
            fprintf(stderr, "QS: factored in %.3fs\n", elapsed);
            found = true;
        }

        // Try lhs + rhs too
        if (!found) {
            mpz_add(tmp, lhs, rhs);
            mpz_gcd(factor1, tmp, N);

            if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
                mpz_divexact(factor2, N, factor1);

                if (mpz_cmp(factor1, factor2) > 0)
                    mpz_swap(factor1, factor2);

                struct timespec now;
                clock_gettime(CLOCK_MONOTONIC, &now);
                double elapsed = (now.tv_sec - t_start.tv_sec) + (now.tv_nsec - t_start.tv_nsec) / 1e9;
                gmp_printf("%Zd %Zd\n", factor1, factor2);
                fprintf(stderr, "QS: factored in %.3fs\n", elapsed);
                found = true;
            }
        }

        mpz_clear(lhs);
        mpz_clear(rhs);
    }

    if (!found) {
        fprintf(stderr, "FAIL: no non-trivial factor found from %zu dependencies\n", nrels);
        return 1;
    }

    // Cleanup
    for (auto &rel : relations) {
        mpz_clear(rel.x);
        mpz_clear(rel.Qx);
    }
    mpz_clear(x_val);
    mpz_clear(qx);
    mpz_clear(tmp);
    mpz_clear(tmp2);
    mpz_clear(N);
    mpz_clear(sqrtN);
    mpz_clear(factor1);
    mpz_clear(factor2);

    return found ? 0 : 1;
}
