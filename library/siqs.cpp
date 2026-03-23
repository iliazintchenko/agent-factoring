/*
 * Quadratic Sieve (basic QS, not self-initializing) implementation.
 * Usage: ./siqs <number>
 * Single-threaded, seed=42.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <gmp.h>

static long power_mod(long base, long exp, long mod) {
    long r = 1;
    base %= mod;
    if (base < 0) base += mod;
    while (exp > 0) {
        if (exp & 1) r = (__int128)r * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return r;
}

static int legendre(long a, long p) {
    long r = power_mod(a, (p - 1) / 2, p);
    return r == p - 1 ? -1 : (int)r;
}

static long tonelli_shanks(long n, long p) {
    if (p == 2) return n & 1;
    n = ((n % p) + p) % p;
    if (n == 0) return 0;

    long Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    if (S == 1) return power_mod(n, (p + 1) / 4, p);

    long z = 2;
    while (power_mod(z, (p - 1) / 2, p) != p - 1) z++;

    long M = S;
    long c = power_mod(z, Q, p);
    long t = power_mod(n, Q, p);
    long R = power_mod(n, (Q + 1) / 2, p);

    while (t != 1) {
        long i = 0, tmp = t;
        while (tmp != 1) { tmp = (__int128)tmp * tmp % p; i++; }
        long b = c;
        for (long j = 0; j < M - i - 1; j++) b = (__int128)b * b % p;
        M = i;
        c = (__int128)b * b % p;
        t = (__int128)t * c % p;
        R = (__int128)R * b % p;
    }
    return R;
}

static std::vector<long> sieve_primes(long limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (long i = 2; i * i <= limit; i++)
        if (is_prime[i])
            for (long j = i * i; j <= limit; j += i)
                is_prime[j] = false;
    std::vector<long> primes;
    for (long i = 2; i <= limit; i++)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

/* GF(2) row-reduction with tracking */
struct BitRow {
    std::vector<uint64_t> bits;
    std::vector<uint64_t> history; /* which relations are combined */
    int words_bits, words_hist;

    void init(int nbits, int nhist) {
        words_bits = (nbits + 63) / 64;
        words_hist = (nhist + 63) / 64;
        bits.assign(words_bits, 0);
        history.assign(words_hist, 0);
    }

    void set_bit(int i) { bits[i / 64] |= (1ULL << (i % 64)); }
    int get_bit(int i) const { return (bits[i / 64] >> (i % 64)) & 1; }
    void set_hist(int i) { history[i / 64] |= (1ULL << (i % 64)); }
    int get_hist(int i) const { return (history[i / 64] >> (i % 64)) & 1; }

    void xor_with(const BitRow &o) {
        for (int i = 0; i < words_bits; i++) bits[i] ^= o.bits[i];
        for (int i = 0; i < words_hist; i++) history[i] ^= o.history[i];
    }

    bool is_zero_bits() const {
        for (int i = 0; i < words_bits; i++) if (bits[i]) return false;
        return true;
    }
};

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    mpz_t N, sqrtN, tmp, g;
    mpz_init(N); mpz_init(sqrtN); mpz_init(tmp); mpz_init(g);

    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    struct timespec tstart;
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    size_t digits = mpz_sizeinbase(N, 10);
    double ln_n = log(2.0) * mpz_sizeinbase(N, 2);
    double ln_ln_n = log(ln_n);

    /* L[1/2] optimal smoothness bound: exp(0.5 * sqrt(ln(n) * ln(ln(n)))) */
    double L_half = exp(0.5 * sqrt(ln_n * ln_ln_n));

    /* Factor base bound - scale with L */
    long B = (long)(L_half * 1.2);
    if (B < 3000) B = 3000;
    if (B > 10000000) B = 10000000;

    /* Sieve interval - needs to be large enough to find enough smooth values */
    long M = B * 20;
    if (M < 100000) M = 100000;
    if (M > 50000000) M = 50000000;

    /* Build factor base */
    std::vector<long> all_primes = sieve_primes(B);
    std::vector<long> fb;  /* factor base primes (index 0 = -1 for sign) */
    fb.push_back(-1);

    mpz_t nmodp;
    mpz_init(nmodp);
    for (long p : all_primes) {
        mpz_mod_ui(nmodp, N, p);
        long nmod = mpz_get_ui(nmodp);
        if (p == 2 || legendre(nmod, p) == 1) {
            fb.push_back(p);
        }
    }
    mpz_clear(nmodp);

    int fb_size = fb.size();
    int target_relations = fb_size + 20;

    fprintf(stderr, "QS: %zu digits, B=%ld, fb=%d, M=%ld, target=%d relations\n",
            digits, B, fb_size, M, target_relations);

    mpz_sqrt(sqrtN, N);
    /* Ensure sqrtN^2 >= N (ceiling sqrt) */
    mpz_mul(tmp, sqrtN, sqrtN);
    if (mpz_cmp(tmp, N) < 0) mpz_add_ui(sqrtN, sqrtN, 1);

    /* Precompute sqrt(N) mod p and sieve start positions */
    std::vector<long> sqrt_nmodp(fb_size);
    std::vector<long> sqrtN_modp(fb_size);
    for (int i = 1; i < fb_size; i++) {
        long p = fb[i];
        mpz_mod_ui(tmp, N, p);
        sqrt_nmodp[i] = tonelli_shanks(mpz_get_ui(tmp), p);
        mpz_mod_ui(tmp, sqrtN, p);
        sqrtN_modp[i] = mpz_get_ui(tmp);
    }

    /* Log approximation sieving */
    std::vector<float> sieve_arr(2 * M + 1, 0.0f);

    /* Target: log2(Q(x)) for x near 0 is approximately log2(2*sqrt(N)*|x|)
     * For the threshold, we want: sum of logs >= log2(Q) - slack */

    fprintf(stderr, "Sieving...\n");
    for (int i = 1; i < fb_size; i++) {
        long p = fb[i];
        if (p > M) break;
        float logp = log2f((float)p);

        if (p == 2) {
            /* Every position where Q(x) is even */
            for (long j = 0; j <= 2 * M; j++) {
                /* Q(x) = (sqrtN+x)^2 - N, check mod 2 */
                /* Just add log2 for all even positions as approximation */
                sieve_arr[j] += logp;
            }
            continue;
        }

        /* Two roots: x = sqrt(N) - sqrtN mod p, and x = -sqrt(N) - sqrtN mod p */
        long r1 = ((sqrt_nmodp[i] - sqrtN_modp[i]) % p + p) % p;
        long r2 = ((p - sqrt_nmodp[i] - sqrtN_modp[i]) % p + p) % p;

        /* Sieve positive direction from r1 */
        for (long j = M + r1; j <= 2 * M; j += p) sieve_arr[j] += logp;
        /* Sieve negative direction from r1 */
        for (long j = M + r1 - p; j >= 0; j -= p) sieve_arr[j] += logp;

        if (r2 != r1) {
            for (long j = M + r2; j <= 2 * M; j += p) sieve_arr[j] += logp;
            for (long j = M + r2 - p; j >= 0; j -= p) sieve_arr[j] += logp;
        }

        /* Also sieve prime powers */
        long pk = p * p;
        while (pk <= B && pk > 0) {
            /* For prime powers, there may be 0 or 2 roots mod p^k */
            /* Simplified: just lift the roots */
            /* Skip for now - basic QS */
            break;
        }
    }

    fprintf(stderr, "Collecting smooth relations...\n");

    struct Relation {
        mpz_t x_plus_sqrt; /* sqrtN + x */
        std::vector<int> exponents; /* full exponents for each fb element */
    };
    std::vector<Relation> relations;

    mpz_t Q_val, x_plus_sqrt;
    mpz_init(Q_val); mpz_init(x_plus_sqrt);

    /* Compute threshold: Q(x) ≈ 2*sqrt(N)*|x| near center
     * log2(Q(0)) ≈ bits(N)/2
     * We want sum_logs >= target - slack where slack allows for large prime variation */
    float base_log2_Q = (float)mpz_sizeinbase(N, 2) / 2.0f;
    float slack = 28.0f; /* Allow ~28 bits of slack for large prime and approximation errors */

    /* Scan from center outward for better Q values */
    /* Build a list of candidate positions sorted by sieve value */
    for (long idx = 0; idx <= 2 * M && (int)relations.size() < target_relations; idx++) {
        long xval = idx - M;

        /* Approximate log2(Q(x)): for x near 0, Q ≈ 2*sqrtN*|x|.
         * More precisely, log2(|Q|) varies. Use a rough estimate. */
        float approx_log2_Q;
        if (xval == 0) {
            /* Q(0) = sqrtN^2 - N, which is small */
            approx_log2_Q = 10.0f; /* small */
        } else {
            long ax = xval < 0 ? -xval : xval;
            approx_log2_Q = base_log2_Q + 1.0f + log2f((float)ax);
        }

        if (sieve_arr[idx] < approx_log2_Q - slack) continue;

        /* Trial division to confirm smoothness */
        mpz_set(x_plus_sqrt, sqrtN);
        if (xval >= 0) mpz_add_ui(x_plus_sqrt, x_plus_sqrt, xval);
        else mpz_sub_ui(x_plus_sqrt, x_plus_sqrt, -xval);

        mpz_mul(Q_val, x_plus_sqrt, x_plus_sqrt);
        mpz_sub(Q_val, Q_val, N);

        if (mpz_sgn(Q_val) == 0) {
            /* Perfect square! Direct factor */
            mpz_gcd(g, x_plus_sqrt, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                gmp_printf("FACTOR: %Zd\n", g);
                mpz_t cofactor; mpz_init(cofactor);
                mpz_divexact(cofactor, N, g);
                gmp_printf("COFACTOR: %Zd\n", cofactor);
                fprintf(stderr, "QS found perfect square factor in %.3fs\n", elapsed);
                mpz_clear(cofactor);
                return 0;
            }
        }

        std::vector<int> exps(fb_size, 0);
        int negative = (mpz_sgn(Q_val) < 0);
        if (negative) { exps[0] = 1; mpz_neg(Q_val, Q_val); }

        for (int j = 1; j < fb_size; j++) {
            long p = fb[j];
            while (mpz_divisible_ui_p(Q_val, p)) {
                mpz_divexact_ui(Q_val, Q_val, p);
                exps[j]++;
            }
        }

        if (mpz_cmp_ui(Q_val, 1) == 0) {
            /* Smooth! */
            Relation rel;
            mpz_init_set(rel.x_plus_sqrt, x_plus_sqrt);
            rel.exponents = exps;
            relations.push_back(rel);

            if (relations.size() % 200 == 0 || relations.size() == 1) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                fprintf(stderr, "  %zu/%d relations (%.1fs)\n", relations.size(), target_relations, elapsed);
            }
        }
    }

    mpz_clear(Q_val); mpz_clear(x_plus_sqrt);

    fprintf(stderr, "Found %zu relations (target %d)\n", relations.size(), target_relations);

    if ((int)relations.size() <= fb_size) {
        fprintf(stderr, "Not enough relations. Increase M or B.\n");
        printf("FAIL\n");
        return 1;
    }

    /* Linear algebra: Gaussian elimination over GF(2) */
    fprintf(stderr, "Linear algebra...\n");
    int nrels = relations.size();
    std::vector<BitRow> matrix(nrels);
    for (int i = 0; i < nrels; i++) {
        matrix[i].init(fb_size, nrels);
        matrix[i].set_hist(i);
        for (int j = 0; j < fb_size; j++) {
            if (relations[i].exponents[j] & 1)
                matrix[i].set_bit(j);
        }
    }

    /* Row reduce */
    std::vector<int> pivots(fb_size, -1);
    int rank = 0;
    for (int col = 0; col < fb_size; col++) {
        int pivot = -1;
        for (int row = rank; row < nrels; row++) {
            if (matrix[row].get_bit(col)) { pivot = row; break; }
        }
        if (pivot == -1) continue;
        if (pivot != rank) std::swap(matrix[pivot], matrix[rank]);
        pivots[col] = rank;
        for (int row = 0; row < nrels; row++) {
            if (row != rank && matrix[row].get_bit(col))
                matrix[row].xor_with(matrix[rank]);
        }
        rank++;
    }

    fprintf(stderr, "Rank = %d, null space dim = %d\n", rank, nrels - rank);

    /* Try each null space vector */
    int found = 0;
    mpz_t prod_x, prod_y, sq;
    mpz_init(prod_x); mpz_init(prod_y); mpz_init(sq);

    for (int row = rank; row < nrels && !found; row++) {
        if (matrix[row].is_zero_bits()) {
            /* This is a dependency */
            mpz_set_ui(prod_x, 1);
            std::vector<long> sum_exp(fb_size, 0);

            for (int i = 0; i < nrels; i++) {
                if (!matrix[row].get_hist(i)) continue;

                mpz_mul(prod_x, prod_x, relations[i].x_plus_sqrt);
                mpz_mod(prod_x, prod_x, N);

                for (int j = 0; j < fb_size; j++)
                    sum_exp[j] += relations[i].exponents[j];
            }

            /* Verify all exponents are even */
            bool all_even = true;
            for (int j = 0; j < fb_size; j++) {
                if (sum_exp[j] & 1) { all_even = false; break; }
            }
            if (!all_even) continue;

            /* y = product of fb[j]^(sum_exp[j]/2) mod N */
            mpz_set_ui(prod_y, 1);
            for (int j = 1; j < fb_size; j++) {
                if (sum_exp[j] == 0) continue;
                mpz_set_ui(tmp, fb[j]);
                mpz_powm_ui(tmp, tmp, sum_exp[j] / 2, N);
                mpz_mul(prod_y, prod_y, tmp);
                mpz_mod(prod_y, prod_y, N);
            }
            /* Handle sign: if sum_exp[0] is odd (it shouldn't be now), skip */

            /* gcd(x-y, N) */
            mpz_sub(tmp, prod_x, prod_y);
            mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                gmp_printf("FACTOR: %Zd\n", g);
                mpz_t cofactor; mpz_init(cofactor);
                mpz_divexact(cofactor, N, g);
                gmp_printf("COFACTOR: %Zd\n", cofactor);
                fprintf(stderr, "QS factored in %.3fs\n", elapsed);
                mpz_clear(cofactor);
                found = 1;
                break;
            }

            /* Try x+y */
            mpz_add(tmp, prod_x, prod_y);
            mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                gmp_printf("FACTOR: %Zd\n", g);
                mpz_t cofactor; mpz_init(cofactor);
                mpz_divexact(cofactor, N, g);
                gmp_printf("COFACTOR: %Zd\n", cofactor);
                fprintf(stderr, "QS factored in %.3fs\n", elapsed);
                mpz_clear(cofactor);
                found = 1;
                break;
            }
        }
    }

    if (!found) {
        /* Also try null vectors that weren't at the bottom */
        for (int row = 0; row < rank && !found; row++) {
            if (!matrix[row].is_zero_bits()) continue;
            /* Same logic as above */
            mpz_set_ui(prod_x, 1);
            std::vector<long> sum_exp(fb_size, 0);
            for (int i = 0; i < nrels; i++) {
                if (!matrix[row].get_hist(i)) continue;
                mpz_mul(prod_x, prod_x, relations[i].x_plus_sqrt);
                mpz_mod(prod_x, prod_x, N);
                for (int j = 0; j < fb_size; j++)
                    sum_exp[j] += relations[i].exponents[j];
            }
            bool all_even = true;
            for (int j = 0; j < fb_size; j++) {
                if (sum_exp[j] & 1) { all_even = false; break; }
            }
            if (!all_even) continue;
            mpz_set_ui(prod_y, 1);
            for (int j = 1; j < fb_size; j++) {
                if (sum_exp[j] == 0) continue;
                mpz_set_ui(tmp, fb[j]);
                mpz_powm_ui(tmp, tmp, sum_exp[j] / 2, N);
                mpz_mul(prod_y, prod_y, tmp);
                mpz_mod(prod_y, prod_y, N);
            }
            mpz_sub(tmp, prod_x, prod_y);
            mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                gmp_printf("FACTOR: %Zd\n", g);
                mpz_t cofactor; mpz_init(cofactor);
                mpz_divexact(cofactor, N, g);
                gmp_printf("COFACTOR: %Zd\n", cofactor);
                fprintf(stderr, "QS factored in %.3fs\n", elapsed);
                mpz_clear(cofactor);
                found = 1;
            }
            mpz_add(tmp, prod_x, prod_y);
            mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                gmp_printf("FACTOR: %Zd\n", g);
                mpz_t cofactor; mpz_init(cofactor);
                mpz_divexact(cofactor, N, g);
                gmp_printf("COFACTOR: %Zd\n", cofactor);
                fprintf(stderr, "QS factored in %.3fs\n", elapsed);
                mpz_clear(cofactor);
                found = 1;
            }
        }
    }

    if (!found) {
        fprintf(stderr, "No non-trivial factor found from %d dependencies.\n", nrels - rank);
        printf("FAIL\n");
    }

    /* Cleanup */
    mpz_clear(prod_x); mpz_clear(prod_y); mpz_clear(sq);
    for (auto &rel : relations) mpz_clear(rel.x_plus_sqrt);
    mpz_clear(N); mpz_clear(sqrtN); mpz_clear(tmp); mpz_clear(g);

    return found ? 0 : 1;
}
