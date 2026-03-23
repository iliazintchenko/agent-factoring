/*
 * Self-Initializing Quadratic Sieve (SIQS)
 *
 * Usage: ./siqs <number>
 * Outputs: factor on stdout, timing on stderr
 *
 * Single-threaded, uses GMP for big integers.
 * Seed is always 42.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cassert>
#include <cstdint>
#include <gmp.h>

using namespace std;

// ============================================================
// Parameters tuned by digit count
// ============================================================
struct SIQSParams {
    int factor_base_size;
    int sieve_radius;
    int large_prime_mult;  // LP bound = largest_fb_prime * this
    int max_polynomials;
    double threshold_offset; // bits subtracted from log2(sqrt(N)*M)
};

static SIQSParams get_params(int digits) {
    if (digits <= 30) return {150, 30000, 40, 100000, 3.5};
    if (digits <= 35) return {300, 50000, 50, 150000, 3.5};
    if (digits <= 40) return {500, 65536, 60, 200000, 4.0};
    if (digits <= 45) return {800, 100000, 80, 300000, 4.5};
    if (digits <= 50) return {1200, 150000, 100, 500000, 5.0};
    if (digits <= 55) return {1800, 200000, 120, 800000, 5.5};
    if (digits <= 60) return {2500, 300000, 150, 1500000, 6.0};
    if (digits <= 65) return {3500, 400000, 200, 2500000, 6.5};
    if (digits <= 70) return {5000, 500000, 250, 4000000, 7.0};
    if (digits <= 75) return {7000, 700000, 300, 6000000, 7.5};
    if (digits <= 80) return {10000, 900000, 400, 10000000, 8.0};
    if (digits <= 85) return {14000, 1200000, 500, 15000000, 8.5};
    if (digits <= 90) return {20000, 1500000, 600, 25000000, 9.0};
    if (digits <= 95) return {28000, 2000000, 800, 40000000, 9.5};
    return {40000, 2500000, 1000, 60000000, 10.0};
}

// ============================================================
// Sieve of Eratosthenes
// ============================================================
static vector<int> sieve_primes(int limit) {
    vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; (long long)i * i <= limit; i++)
        if (is_prime[i])
            for (int j = i * i; j <= limit; j += i)
                is_prime[j] = false;
    vector<int> primes;
    for (int i = 2; i <= limit; i++)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

// ============================================================
// Tonelli-Shanks: find r such that r^2 ≡ n (mod p)
// Returns r or -1 if n is not a QR mod p
// ============================================================
static long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = (__int128)result * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return result;
}

static long long tonelli_shanks(long long n, long long p) {
    if (p == 2) return n & 1;
    n = ((n % p) + p) % p;
    if (n == 0) return 0;
    if (mod_pow(n, (p - 1) / 2, p) != 1) return -1;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);

    long long Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    long long z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) z++;

    long long M = S;
    long long c = mod_pow(z, Q, p);
    long long t = mod_pow(n, Q, p);
    long long R = mod_pow(n, (Q + 1) / 2, p);

    while (true) {
        if (t == 1) return R;
        long long i = 0, tmp = t;
        while (tmp != 1) { tmp = (__int128)tmp * tmp % p; i++; }
        long long b = c;
        for (long long j = 0; j < M - i - 1; j++)
            b = (__int128)b * b % p;
        M = i;
        c = (__int128)b * b % p;
        t = (__int128)t * c % p;
        R = (__int128)R * b % p;
    }
}

// ============================================================
// GF(2) Gaussian elimination using 64-bit words
// ============================================================
struct BitMatrix {
    int rows, cols;
    int words_per_row;
    vector<vector<uint64_t>> data;

    BitMatrix() : rows(0), cols(0), words_per_row(0) {}
    BitMatrix(int r, int c) : rows(r), cols(c), words_per_row((c + 63) / 64) {
        data.assign(r, vector<uint64_t>(words_per_row, 0));
    }
    void set(int r, int c) { data[r][c >> 6] |= (1ULL << (c & 63)); }
    bool get(int r, int c) const { return (data[r][c >> 6] >> (c & 63)) & 1; }
    void xor_row(int dst, int src) {
        for (int i = 0; i < words_per_row; i++)
            data[dst][i] ^= data[src][i];
    }
};

// Find null space vectors of a matrix over GF(2)
// Returns list of vectors in null space (each vector = set of row indices)
static vector<vector<int>> find_null_space(int num_rows, int num_cols,
                                            const vector<vector<int>>& exp_vectors) {
    // Build augmented matrix [M | I] where M has rows=relations, cols=fb primes
    int total_cols = num_cols + num_rows;
    BitMatrix mat(num_rows, total_cols);

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            if (exp_vectors[i][j] & 1)
                mat.set(i, j);
        }
        mat.set(i, num_cols + i); // identity part
    }

    // Row echelon form on the M part
    vector<int> pivot_col(num_rows, -1);
    int cur_row = 0;
    for (int col = 0; col < num_cols && cur_row < num_rows; col++) {
        // Find pivot in this column
        int piv = -1;
        for (int row = cur_row; row < num_rows; row++) {
            if (mat.get(row, col)) { piv = row; break; }
        }
        if (piv < 0) continue;

        // Swap rows
        if (piv != cur_row) swap(mat.data[piv], mat.data[cur_row]);
        pivot_col[cur_row] = col;

        // Eliminate all other rows
        for (int row = 0; row < num_rows; row++) {
            if (row != cur_row && mat.get(row, col))
                mat.xor_row(row, cur_row);
        }
        cur_row++;
    }

    // Rows from cur_row..num_rows-1 should have all-zero M part => null space
    vector<vector<int>> null_vectors;
    for (int row = cur_row; row < num_rows; row++) {
        // Verify M part is zero
        bool all_zero = true;
        for (int col = 0; col < num_cols; col++) {
            if (mat.get(row, col)) { all_zero = false; break; }
        }
        if (!all_zero) continue;

        // Extract which original rows are involved
        vector<int> involved;
        for (int i = 0; i < num_rows; i++) {
            if (mat.get(row, num_cols + i))
                involved.push_back(i);
        }
        if (!involved.empty())
            null_vectors.push_back(involved);
    }

    // Also check rows before cur_row that might have all-zero M part
    // (shouldn't happen if elimination is correct, but just in case)

    return null_vectors;
}

// ============================================================
// Main SIQS
// ============================================================

struct FBEntry {
    int p;
    int logp;       // floor(log2(p) * 256) / 16  i.e. scaled log
    int sqrt_n;     // sqrt(N) mod p
    int soln1, soln2;
};

struct Relation {
    mpz_t ax_plus_b; // (ax+b) such that (ax+b)^2 - N = a * Q(x)
    vector<int> exponents; // full exponent vector over factor base (including -1)
};

// Simple deterministic PRNG
static uint64_t rng_state = 42;
static uint64_t next_rand() {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return rng_state >> 16;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    auto elapsed_sec = [&]() -> double {
        struct timespec t1;
        clock_gettime(CLOCK_MONOTONIC, &t1);
        return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    };

    mpz_t N, sqrtN, factor, tmp, tmp2;
    mpz_init_set_str(N, argv[1], 10);
    mpz_init(sqrtN);
    mpz_init(factor);
    mpz_init(tmp);
    mpz_init(tmp2);

    // Quick checks
    if (mpz_even_p(N)) { printf("2\n"); return 0; }
    if (mpz_perfect_square_p(N)) {
        mpz_sqrt(sqrtN, N);
        gmp_printf("%Zd\n", sqrtN);
        return 0;
    }

    // Small trial division up to 1M
    {
        vector<int> small_primes = sieve_primes(1000000);
        for (int p : small_primes) {
            if (mpz_divisible_ui_p(N, p)) {
                printf("%d\n", p);
                fprintf(stderr, "SIQS: found small factor %d in %.3fs\n", p, elapsed_sec());
                return 0;
            }
        }
    }

    int digits = (int)mpz_sizeinbase(N, 10);
    SIQSParams params = get_params(digits);

    fprintf(stderr, "SIQS: factoring %d-digit number, FB=%d, M=%d\n",
            digits, params.factor_base_size, params.sieve_radius);

    mpz_sqrt(sqrtN, N);

    // ---- Build factor base ----
    int prime_limit = 200000;
    while (true) {
        vector<int> ap = sieve_primes(prime_limit);
        if ((int)ap.size() >= params.factor_base_size * 3) break;
        prime_limit *= 2;
    }
    vector<int> all_primes = sieve_primes(prime_limit);

    vector<FBEntry> fb;
    // Index 0: -1 (sign)
    fb.push_back({-1, 0, 0, 0, 0});

    for (size_t i = 0; i < all_primes.size() && (int)fb.size() <= params.factor_base_size; i++) {
        int p = all_primes[i];
        int n_mod_p = (int)mpz_fdiv_ui(N, p);
        long long sq = tonelli_shanks(n_mod_p, p);
        if (sq >= 0) {
            int logp = (int)(log2((double)p) * 16.0 + 0.5);
            fb.push_back({p, logp, (int)sq, 0, 0});
        }
    }

    int fb_size = (int)fb.size();
    fprintf(stderr, "Factor base: %d primes, largest = %d\n", fb_size, fb.back().p);

    long long lp_bound = (long long)fb.back().p * params.large_prime_mult;

    int M = params.sieve_radius;
    int needed = fb_size + 30;

    // Store relations
    vector<Relation> relations;
    // Large prime partials: lp_value -> index in partials_list
    unordered_map<int, int> partial_map;
    struct Partial { Relation rel; int lp; };
    vector<Partial> partials_list;

    int total_full = 0, total_combined = 0, total_partials = 0;
    int poly_count = 0;

    // target_a ≈ sqrt(2N) / M
    mpz_t target_a, a_val, b_val, c_val;
    mpz_init(target_a);
    mpz_init(a_val);
    mpz_init(b_val);
    mpz_init(c_val);

    mpz_mul_ui(target_a, N, 2);
    mpz_sqrt(target_a, target_a);
    mpz_tdiv_q_ui(target_a, target_a, M);

    // Number of primes composing 'a'
    int num_a_factors;
    if (digits <= 35) num_a_factors = 3;
    else if (digits <= 45) num_a_factors = 4;
    else if (digits <= 55) num_a_factors = 5;
    else if (digits <= 65) num_a_factors = 6;
    else if (digits <= 75) num_a_factors = 7;
    else if (digits <= 85) num_a_factors = 8;
    else num_a_factors = 9;

    // Range of FB indices to pick 'a' factors from
    int a_lo = max(2, fb_size / 5);
    int a_hi = min(fb_size - 1, fb_size * 4 / 5);

    // Sieve array - int16 to hold values up to ~1500
    vector<int16_t> sieve_arr(2 * M + 2);

    // Precompute log2 threshold
    // At x=M, Q(x)/a ≈ aM^2 + 2bM + c ≈ aM^2 ≈ sqrt(2N)*M
    // log2(sqrt(2N)*M) = log2(2N)/2 + log2(M)
    double log2_qmax = mpz_sizeinbase(N, 2) / 2.0 + log2((double)M);
    // We want sieve values that accumulate close to this
    int threshold = (int)((log2_qmax - params.threshold_offset) * 16.0);
    if (threshold < 0) threshold = 0;

    fprintf(stderr, "Threshold: %d (log2_qmax*16=%.0f, offset=%.1f)\n",
            threshold, log2_qmax * 16, params.threshold_offset);

    // Main sieving loop
    while (total_full + total_combined < needed && poly_count < params.max_polynomials) {
        // ---- Generate new polynomial ----
        // Select random primes for 'a'
        vector<int> a_idx;
        mpz_set_ui(a_val, 1);

        vector<int> pool;
        for (int i = a_lo; i <= a_hi; i++) pool.push_back(i);

        for (int i = 0; i < num_a_factors && !pool.empty(); i++) {
            int r = next_rand() % pool.size();
            int idx = pool[r];
            a_idx.push_back(idx);
            mpz_mul_ui(a_val, a_val, fb[idx].p);
            pool.erase(pool.begin() + r);
        }
        sort(a_idx.begin(), a_idx.end());

        // Compute b via CRT: b^2 ≡ N (mod a), b ≡ sqrt(N) (mod q_i) for each q_i | a
        // Using Garner's algorithm / CRT
        mpz_set_ui(b_val, 0);
        for (int i = 0; i < num_a_factors; i++) {
            int qi = fb[a_idx[i]].p;
            int si = fb[a_idx[i]].sqrt_n;

            // Compute a/qi
            mpz_t a_div_qi, qi_mpz, inv_val, contrib;
            mpz_init(a_div_qi);
            mpz_init_set_ui(qi_mpz, qi);
            mpz_init(inv_val);
            mpz_init(contrib);

            mpz_tdiv_q_ui(a_div_qi, a_val, qi);

            // inverse of (a/qi) mod qi
            mpz_invert(inv_val, a_div_qi, qi_mpz);

            // b += si * inv * (a/qi)
            mpz_mul_ui(contrib, inv_val, si);
            mpz_mul(contrib, contrib, a_div_qi);
            mpz_add(b_val, b_val, contrib);

            mpz_clear(a_div_qi);
            mpz_clear(qi_mpz);
            mpz_clear(inv_val);
            mpz_clear(contrib);
        }
        mpz_mod(b_val, b_val, a_val);

        // Ensure b < a/2 (pick smaller root)
        mpz_tdiv_q_ui(tmp, a_val, 2);
        if (mpz_cmp(b_val, tmp) > 0)
            mpz_sub(b_val, a_val, b_val);

        // Verify b^2 ≡ N (mod a)
        mpz_mul(tmp, b_val, b_val);
        mpz_sub(tmp, tmp, N);
        mpz_mod(tmp, tmp, a_val);
        if (mpz_sgn(tmp) != 0) {
            poly_count++;
            continue;
        }

        // c = (b^2 - N) / a
        mpz_mul(c_val, b_val, b_val);
        mpz_sub(c_val, c_val, N);
        mpz_tdiv_q(c_val, c_val, a_val);
        // Note: c should be negative since b^2 < N typically
        // Q(x) = a*x^2 + 2*b*x + c, and a*Q(x) = (a*x+b)^2 - N

        // Update sieve roots for each FB prime
        for (int i = 1; i < fb_size; i++) {
            int p = fb[i].p;

            // Check if p divides a
            bool divides_a = false;
            for (int j : a_idx) {
                if (j == i) { divides_a = true; break; }
            }
            if (divides_a) {
                fb[i].soln1 = fb[i].soln2 = -1;
                continue;
            }

            // Roots of Q(x) mod p: x ≡ (±sqrt(N) - b) * a^(-1) (mod p)
            long long a_mod_p = mpz_fdiv_ui(a_val, p);
            long long b_mod_p = mpz_fdiv_ui(b_val, p);
            long long ainv = mod_pow(a_mod_p, p - 2, p);
            int sq = fb[i].sqrt_n;

            long long s1 = ((__int128)(sq - b_mod_p + p) % p * ainv) % p;
            long long s2 = ((__int128)(p - sq - b_mod_p + p) % p * ainv) % p;

            fb[i].soln1 = (int)s1;
            fb[i].soln2 = (int)s2;
        }

        // ---- Sieve ----
        memset(sieve_arr.data(), 0, (2 * M + 2) * sizeof(int16_t));

        for (int i = 1; i < fb_size; i++) {
            int p = fb[i].p;
            int logp = fb[i].logp;

            if (p < 3) {
                // Handle p=2 specially (every other position)
                if (fb[i].soln1 >= 0) {
                    int off = ((fb[i].soln1 % 2) + M) % 2;
                    for (int j = off; j <= 2 * M; j += 2)
                        sieve_arr[j] += logp;
                }
                continue;
            }

            if (fb[i].soln1 < 0) continue; // divides a

            // Sieve soln1
            {
                int start = (int)(((long long)fb[i].soln1 + M) % p);
                for (int j = start; j <= 2 * M; j += p)
                    sieve_arr[j] += logp;
            }

            // Sieve soln2 (if different from soln1)
            if (fb[i].soln1 != fb[i].soln2) {
                int start = (int)(((long long)fb[i].soln2 + M) % p);
                for (int j = start; j <= 2 * M; j += p)
                    sieve_arr[j] += logp;
            }
        }

        // Handle primes dividing a (single root sieving)
        for (int ai : a_idx) {
            int p = fb[ai].p;
            int logp = fb[ai].logp;
            long long b_mod_p = mpz_fdiv_ui(b_val, p);
            long long c_mod_p = mpz_fdiv_ui(c_val, p); // this gives positive remainder

            // For q | a: Q(x) ≡ 2bx + c (mod q) since a ≡ 0 (mod q)
            // Root: x ≡ -c * (2b)^(-1) (mod q)
            long long twob = (2 * b_mod_p) % p;
            if (twob == 0) continue;
            long long twob_inv = mod_pow(twob, p - 2, p);

            // c_val might be negative; mpz_fdiv_ui always returns non-negative
            // x ≡ (p - c_mod_p) * twob_inv (mod p)
            long long root;
            // Actually: c_val could be negative, mpz_fdiv_ui returns the positive remainder
            // Q(x) = ax^2 + 2bx + c, and c = (b^2 - N)/a
            // We need 2bx + c ≡ 0 (mod p)  =>  x ≡ -c/(2b) (mod p)
            // mpz_fdiv_ui gives c mod p (non-negative)
            mpz_t c_mod;
            mpz_init(c_mod);
            mpz_fdiv_r_ui(c_mod, c_val, p);
            long long c_pos = mpz_get_ui(c_mod);
            mpz_clear(c_mod);

            root = ((p - c_pos) % p * twob_inv) % p;

            int off = (int)(((long long)root + M) % p);
            for (int j = off; j <= 2 * M; j += p)
                sieve_arr[j] += logp;
        }

        // ---- Check smooth candidates ----
        for (int idx = 0; idx <= 2 * M; idx++) {
            if (sieve_arr[idx] < threshold) continue;

            int x = idx - M;

            // Compute Q(x) = a*x^2 + 2*b*x + c
            mpz_t Qval;
            mpz_init(Qval);
            // Q(x) = (a*x + 2b)*x + c
            mpz_mul_si(Qval, a_val, x);    // a*x
            mpz_addmul_ui(Qval, b_val, 2); // a*x + 2b
            mpz_mul_si(Qval, Qval, x);     // (a*x + 2b)*x
            mpz_add(Qval, Qval, c_val);    // a*x^2 + 2*b*x + c

            // Trial divide Q(x) over factor base
            mpz_t rem;
            mpz_init(rem);
            mpz_abs(rem, Qval);
            bool negative = (mpz_sgn(Qval) < 0);

            vector<int> exponents(fb_size, 0);
            if (negative) exponents[0] = 1; // sign

            bool smooth = true;
            for (int i = 1; i < fb_size; i++) {
                int p = fb[i].p;
                while (mpz_divisible_ui_p(rem, p)) {
                    mpz_divexact_ui(rem, rem, p);
                    exponents[i]++;
                }
            }

            int lp = 0;
            if (mpz_cmp_ui(rem, 1) == 0) {
                // Fully smooth
            } else if (mpz_fits_slong_p(rem) && mpz_get_si(rem) <= lp_bound && mpz_get_si(rem) > 1) {
                lp = (int)mpz_get_si(rem);
            } else {
                smooth = false;
            }

            if (smooth) {
                // Build relation: (a*x + b)^2 ≡ a * Q(x) (mod N)
                // The RHS factorization is: a * Q(x) = product of FB primes (and sign)
                // So we need the exponents of Q(x) plus the exponents of a
                Relation rel;
                mpz_init(rel.ax_plus_b);
                mpz_mul_si(rel.ax_plus_b, a_val, x);
                mpz_add(rel.ax_plus_b, rel.ax_plus_b, b_val);
                rel.exponents = exponents;

                // Add exponents for prime factors of a
                for (int ai : a_idx) {
                    rel.exponents[ai]++;
                }

                if (lp == 0) {
                    // Full relation
                    relations.push_back(rel);
                    total_full++;
                } else {
                    // Partial relation with large prime
                    auto it = partial_map.find(lp);
                    if (it != partial_map.end()) {
                        // Combine with existing partial
                        Partial& other = partials_list[it->second];
                        Relation combined;
                        mpz_init(combined.ax_plus_b);
                        mpz_mul(combined.ax_plus_b, rel.ax_plus_b, other.rel.ax_plus_b);
                        // Divide by lp to cancel the lp^2 on the RHS
                        // (ax1+b1)(ax2+b2)/lp squared = (a1*Q1*a2*Q2)/lp^2
                        mpz_t lp_mpz, lp_inv;
                        mpz_init_set_ui(lp_mpz, lp);
                        mpz_init(lp_inv);
                        mpz_invert(lp_inv, lp_mpz, N);
                        mpz_mul(combined.ax_plus_b, combined.ax_plus_b, lp_inv);
                        mpz_mod(combined.ax_plus_b, combined.ax_plus_b, N);
                        mpz_clear(lp_mpz);
                        mpz_clear(lp_inv);

                        combined.exponents.resize(fb_size);
                        for (int k = 0; k < fb_size; k++)
                            combined.exponents[k] = rel.exponents[k] + other.rel.exponents[k];
                        // lp^2 is canceled by dividing X by lp

                        relations.push_back(combined);
                        total_combined++;

                        // Remove used partial
                        partial_map.erase(it);
                    } else {
                        partial_map[lp] = partials_list.size();
                        partials_list.push_back({rel, lp});
                        total_partials++;
                    }
                }
            }

            mpz_clear(Qval);
            mpz_clear(rem);

            if (total_full + total_combined >= needed) break;
        }

        poly_count++;

        if (poly_count % 1000 == 0) {
            fprintf(stderr, "  [%.1fs] poly %d: %d full + %d combined = %d/%d (partials: %d)\n",
                    elapsed_sec(), poly_count, total_full, total_combined,
                    total_full + total_combined, needed, total_partials);
        }
    }

    fprintf(stderr, "Collected %d relations (%d full + %d combined) in %d polys, %.1fs\n",
            total_full + total_combined, total_full, total_combined, poly_count, elapsed_sec());

    if (total_full + total_combined < fb_size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n",
                total_full + total_combined, fb_size + 1);
        return 1;
    }

    // ---- Linear algebra ----
    int nrels = (int)relations.size();
    fprintf(stderr, "Linear algebra: %d relations x %d columns\n", nrels, fb_size);

    // Extract exponent vectors mod 2
    vector<vector<int>> exp_vecs(nrels, vector<int>(fb_size));
    for (int i = 0; i < nrels; i++) {
        for (int j = 0; j < fb_size; j++)
            exp_vecs[i][j] = relations[i].exponents[j];
    }

    vector<vector<int>> null_vecs = find_null_space(nrels, fb_size, exp_vecs);
    fprintf(stderr, "Found %d null space vectors\n", (int)null_vecs.size());

    bool factored = false;
    mpz_t X, Y;
    mpz_init(X);
    mpz_init(Y);

    for (auto& nv : null_vecs) {
        if (factored) break;

        // X = product of (ax+b) mod N
        mpz_set_ui(X, 1);
        vector<long long> total_exp(fb_size, 0);

        for (int idx : nv) {
            mpz_mul(X, X, relations[idx].ax_plus_b);
            mpz_mod(X, X, N);
            for (int j = 0; j < fb_size; j++)
                total_exp[j] += relations[idx].exponents[j];
        }

        // Check all exponents are even
        bool all_even = true;
        for (int j = 0; j < fb_size; j++) {
            if (total_exp[j] & 1) { all_even = false; break; }
        }
        if (!all_even) continue;

        // Y = product of p^(exp/2) mod N
        mpz_set_ui(Y, 1);
        for (int j = 1; j < fb_size; j++) {
            if (total_exp[j] == 0) continue;
            long long half = total_exp[j] / 2;
            mpz_set_ui(tmp, fb[j].p);
            mpz_powm_ui(tmp, tmp, half, N);
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, N);
        }

        // Try gcd(X - Y, N) and gcd(X + Y, N)
        mpz_sub(tmp, X, Y);
        mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            factored = true;
            break;
        }

        mpz_add(tmp, X, Y);
        mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            factored = true;
            break;
        }
    }

    double total_time = elapsed_sec();

    if (factored) {
        gmp_printf("%Zd\n", factor);
        fprintf(stderr, "SIQS: factored in %.3f seconds\n", total_time);
    } else {
        fprintf(stderr, "SIQS: FAILED after %.3f seconds (%d null vectors tried)\n",
                total_time, (int)null_vecs.size());
        return 1;
    }

    // Cleanup
    mpz_clear(N); mpz_clear(sqrtN); mpz_clear(factor);
    mpz_clear(tmp); mpz_clear(tmp2);
    mpz_clear(target_a); mpz_clear(a_val); mpz_clear(b_val); mpz_clear(c_val);
    mpz_clear(X); mpz_clear(Y);
    for (auto& r : relations) mpz_clear(r.ax_plus_b);
    for (auto& p : partials_list) mpz_clear(p.rel.ax_plus_b);

    return 0;
}
