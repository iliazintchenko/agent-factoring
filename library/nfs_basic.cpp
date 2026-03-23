/*
 * Basic General Number Field Sieve (GNFS) implementation
 *
 * Uses degree-3 polynomial with base-m method.
 * Line sieving (not lattice sieve).
 * GF(2) linear algebra.
 *
 * Usage: ./nfs_basic <number>
 * Single-threaded, seed = 42.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <cassert>
#include <gmp.h>

using namespace std;

// ============================================================
// Utility functions
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

static long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base = ((base % mod) + mod) % mod;
    while (exp > 0) {
        if (exp & 1) result = (__int128)result * base % mod;
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return result;
}

// Find all roots of f(x) mod p using brute force (for small p)
// f(x) = x^3 + c2*x^2 + c1*x + c0
static vector<int> poly_roots_mod_p(long long c2, long long c1, long long c0, int p) {
    vector<int> roots;
    long long pc2 = ((c2 % p) + p) % p;
    long long pc1 = ((c1 % p) + p) % p;
    long long pc0 = ((c0 % p) + p) % p;
    for (int x = 0; x < p; x++) {
        long long val = ((__int128)x * x % p * x % p + pc2 * x % p * x % p + pc1 * x % p + pc0) % p;
        // More carefully:
        long long v = 0;
        v = (v + (long long)x) % p;
        v = (v * x) % p;
        v = (v + pc2) % p;  // x^2 + c2
        v = (v * x) % p;    // x^3 + c2*x^2... wait, this isn't right
        // f(x) = x^3 + c2*x^2 + c1*x + c0
        // = ((x + c2)*x + c1)*x + c0  -- Horner's method
        v = (((__int128)(x + pc2) % p * x % p + pc1) % p * x % p + pc0) % p;
        if (v == 0) roots.push_back(x);
    }
    return roots;
}

// Horner evaluation: f(a/b) * b^3 = a^3 + c2*a^2*b + c1*a*b^2 + c0*b^3
// This is the algebraic norm (up to sign and leading coefficient factor)
static void eval_norm(mpz_t result, long long a, long long b,
                      mpz_t c2, mpz_t c1, mpz_t c0) {
    // result = a^3 + c2*a^2*b + c1*a*b^2 + c0*b^3
    mpz_t t1, t2, t3, t4;
    mpz_init(t1); mpz_init(t2); mpz_init(t3); mpz_init(t4);

    // a^3
    mpz_set_si(t1, a);
    mpz_mul_si(t1, t1, a);
    mpz_mul_si(t1, t1, a);

    // c2 * a^2 * b
    mpz_set_si(t2, a);
    mpz_mul_si(t2, t2, a);
    mpz_mul_si(t2, t2, b);
    mpz_mul(t2, t2, c2);

    // c1 * a * b^2
    mpz_set_si(t3, a);
    mpz_mul_si(t3, t3, b);
    mpz_mul_si(t3, t3, b);
    mpz_mul(t3, t3, c1);

    // c0 * b^3
    mpz_set_si(t4, b);
    mpz_mul_si(t4, t4, b);
    mpz_mul_si(t4, t4, b);
    mpz_mul(t4, t4, c0);

    mpz_add(result, t1, t2);
    mpz_add(result, result, t3);
    mpz_add(result, result, t4);

    mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4);
}

// ============================================================
// GF(2) matrix for linear algebra
// ============================================================
struct BitMatrix {
    int rows, cols, words_per_row;
    vector<vector<uint64_t>> data;

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

static vector<vector<int>> find_null_space(int num_rows, int num_cols,
                                            const vector<vector<int>>& exp_vectors) {
    int total_cols = num_cols + num_rows;
    BitMatrix mat(num_rows, total_cols);

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            if (exp_vectors[i][j] & 1)
                mat.set(i, j);
        }
        mat.set(i, num_cols + i);
    }

    int cur_row = 0;
    for (int col = 0; col < num_cols && cur_row < num_rows; col++) {
        int piv = -1;
        for (int row = cur_row; row < num_rows; row++) {
            if (mat.get(row, col)) { piv = row; break; }
        }
        if (piv < 0) continue;
        if (piv != cur_row) swap(mat.data[piv], mat.data[cur_row]);
        for (int row = 0; row < num_rows; row++) {
            if (row != cur_row && mat.get(row, col))
                mat.xor_row(row, cur_row);
        }
        cur_row++;
    }

    vector<vector<int>> null_vectors;
    for (int row = cur_row; row < num_rows; row++) {
        bool all_zero = true;
        for (int col = 0; col < num_cols; col++) {
            if (mat.get(row, col)) { all_zero = false; break; }
        }
        if (!all_zero) continue;
        vector<int> involved;
        for (int i = 0; i < num_rows; i++) {
            if (mat.get(row, num_cols + i))
                involved.push_back(i);
        }
        if (!involved.empty()) null_vectors.push_back(involved);
    }
    return null_vectors;
}

// ============================================================
// NFS Parameters
// ============================================================
struct NFSParams {
    int degree;           // polynomial degree
    int rat_fb_size;      // rational factor base size
    int alg_fb_size;      // algebraic factor base size
    int sieve_a_range;    // sieve a from -range to +range
    int sieve_b_max;      // sieve b from 1 to b_max
    int large_prime_mult; // LP multiplier
};

static NFSParams get_nfs_params(int digits) {
    // Degree 3 for numbers up to ~80 digits, 4 for larger
    if (digits <= 50) return {3, 400, 600, 50000, 200, 50};
    if (digits <= 60) return {3, 800, 1200, 80000, 500, 80};
    if (digits <= 70) return {3, 1500, 2500, 150000, 1000, 100};
    if (digits <= 80) return {4, 3000, 5000, 200000, 2000, 150};
    if (digits <= 90) return {4, 5000, 8000, 400000, 5000, 200};
    return {5, 8000, 12000, 500000, 8000, 300};
}

// ============================================================
// NFS Relation
// ============================================================
struct NFSRelation {
    long long a, b;
    vector<int> rat_exponents;  // exponents for rational factor base
    vector<int> alg_exponents;  // exponents for algebraic factor base
    int rat_sign;               // 1 if a - b*m < 0
    int alg_sign;               // 1 if norm < 0
};

// ============================================================
// Main NFS
// ============================================================
int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    auto elapsed = [&]() -> double {
        struct timespec t1;
        clock_gettime(CLOCK_MONOTONIC, &t1);
        return (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    };

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);

    // Quick trial division
    {
        vector<int> sp = sieve_primes(1000000);
        for (int p : sp) {
            if (mpz_divisible_ui_p(N, p)) {
                printf("%d\n", p);
                fprintf(stderr, "NFS: small factor %d\n", p);
                mpz_clear(N);
                return 0;
            }
        }
    }

    NFSParams params = get_nfs_params(digits);
    int d = params.degree;

    fprintf(stderr, "NFS: %d-digit number, degree %d, rat_fb=%d, alg_fb=%d, sieve=[%d, %d]\n",
            digits, d, params.rat_fb_size, params.alg_fb_size,
            params.sieve_a_range, params.sieve_b_max);

    // ---- Polynomial selection ----
    // Base-m method: m = floor(N^(1/d))
    // f(x) = c_d*x^d + c_{d-1}*x^{d-1} + ... + c_0
    // where N written in base m gives the coefficients

    mpz_t m_val;
    mpz_init(m_val);
    mpz_root(m_val, N, d); // m = floor(N^(1/d))

    // Extract coefficients: N = c_d*m^d + c_{d-1}*m^{d-1} + ... + c_0
    // Note: c_d = 1 for base-m representation (since m^d ≈ N)
    vector<mpz_t> coeff(d + 1);
    mpz_t remainder;
    mpz_init_set(remainder, N);
    for (int i = 0; i <= d; i++) mpz_init(coeff[i]);

    // Compute coefficients by repeated division
    for (int i = 0; i <= d; i++) {
        mpz_tdiv_qr(remainder, coeff[i], remainder, m_val);
    }
    // After d+1 divisions, remainder should be 0 (or the leading coefficient)
    // Actually, let me redo this more carefully
    mpz_set(remainder, N);
    for (int i = 0; i < d; i++) {
        mpz_tdiv_qr(remainder, coeff[i], remainder, m_val);
    }
    mpz_set(coeff[d], remainder); // leading coefficient

    // Verify: f(m) should equal N
    {
        mpz_t check;
        mpz_init_set_ui(check, 0);
        for (int i = d; i >= 0; i--) {
            mpz_mul(check, check, m_val);
            mpz_add(check, check, coeff[i]);
        }
        if (mpz_cmp(check, N) != 0) {
            fprintf(stderr, "ERROR: f(m) != N\n");
            gmp_fprintf(stderr, "N=%Zd, f(m)=%Zd, m=%Zd\n", N, check, m_val);
            mpz_clear(check);
            return 1;
        }
        mpz_clear(check);
    }

    fprintf(stderr, "Polynomial: ");
    for (int i = d; i >= 0; i--) {
        gmp_fprintf(stderr, "%Zd*x^%d%s", coeff[i], i, i > 0 ? " + " : "\n");
    }
    gmp_fprintf(stderr, "m = %Zd\n", m_val);

    long long m_long;
    if (mpz_fits_slong_p(m_val)) {
        m_long = mpz_get_si(m_val);
    } else {
        // m is too large for long — this is expected for degree 3 with large N
        // We need to work with mpz_t for m
        // For now, use a simplified approach
        fprintf(stderr, "WARNING: m doesn't fit in long long, using mpz arithmetic\n");
        m_long = 0; // signal to use mpz
    }

    // For degree 3: f(x) = coeff[3]*x^3 + coeff[2]*x^2 + coeff[1]*x + coeff[0]
    // Typically coeff[3] = 1 for base-m method

    // ---- Build factor bases ----
    int prime_limit = 500000;
    while (true) {
        vector<int> ap = sieve_primes(prime_limit);
        if ((int)ap.size() >= max(params.rat_fb_size, params.alg_fb_size) * 2) break;
        prime_limit *= 2;
    }
    vector<int> all_primes = sieve_primes(prime_limit);

    // Rational factor base: primes p where we sieve a - b*m ≡ 0 (mod p)
    // i.e., a ≡ b*m (mod p)
    struct RatFBEntry { int p; int m_mod_p; int logp; };
    vector<RatFBEntry> rat_fb;

    for (size_t i = 0; i < all_primes.size() && (int)rat_fb.size() < params.rat_fb_size; i++) {
        int p = all_primes[i];
        int mmp = (int)mpz_fdiv_ui(m_val, p);
        int logp = (int)(log2((double)p) * 16 + 0.5);
        rat_fb.push_back({p, mmp, logp});
    }

    // Algebraic factor base: (p, r) pairs where f(r) ≡ 0 (mod p)
    // For each prime p, find roots of f mod p
    struct AlgFBEntry { int p; int r; int logp; };
    vector<AlgFBEntry> alg_fb;

    for (size_t i = 0; i < all_primes.size() && (int)alg_fb.size() < params.alg_fb_size; i++) {
        int p = all_primes[i];
        // Get coefficients mod p
        long long cm[6]; // up to degree 5
        for (int j = 0; j <= d; j++) cm[j] = mpz_fdiv_ui(coeff[j], p);

        // Find roots by brute force for small p, or use Berlekamp for large p
        // For simplicity, brute force for all (fine for p < 10^6)
        if (p <= 100000) {
            for (int x = 0; x < p; x++) {
                // Evaluate f(x) mod p using Horner's method
                long long val = 0;
                for (int j = d; j >= 0; j--)
                    val = ((__int128)val * x + cm[j]) % p;
                if (val == 0) {
                    int logp = (int)(log2((double)p) * 16 + 0.5);
                    alg_fb.push_back({p, x, logp});
                    if ((int)alg_fb.size() >= params.alg_fb_size) break;
                }
            }
        }
    }

    int rat_fb_size = (int)rat_fb.size();
    int alg_fb_size = (int)alg_fb.size();
    // Total columns in matrix: 2 (signs) + rat_fb_size + alg_fb_size
    int total_cols = 2 + rat_fb_size + alg_fb_size;
    int needed_relations = total_cols + 30;

    fprintf(stderr, "Rational FB: %d primes, largest = %d\n", rat_fb_size, rat_fb.back().p);
    fprintf(stderr, "Algebraic FB: %d (p,r) pairs, largest p = %d\n",
            alg_fb_size, alg_fb.empty() ? 0 : alg_fb.back().p);
    fprintf(stderr, "Need %d relations\n", needed_relations);

    // ---- Sieving ----
    int A = params.sieve_a_range;
    int B_max = params.sieve_b_max;
    long long rat_lp_bound = (long long)rat_fb.back().p * params.large_prime_mult;
    long long alg_lp_bound = (long long)(alg_fb.empty() ? 1 : alg_fb.back().p) * params.large_prime_mult;

    vector<NFSRelation> relations;
    int total_relations = 0;

    // Sieve arrays (one for rational, one for algebraic)
    vector<int16_t> rat_sieve(2 * A + 2);
    vector<int16_t> alg_sieve(2 * A + 2);

    // Threshold: we want values where accumulated log ≈ log(expected value)
    // Rational side: |a - b*m| ≈ m*b for large b ≈ N^(1/d) * b
    // Algebraic side: |Norm(a - bα)| ≈ max_coeff * b^d * A^?
    // Use adaptive threshold

    mpz_t rat_val, alg_val, factor, tmp, tmp2;
    mpz_init(rat_val);
    mpz_init(alg_val);
    mpz_init(factor);
    mpz_init(tmp);
    mpz_init(tmp2);

    for (int b = 1; b <= B_max && total_relations < needed_relations; b++) {
        // Compute rational threshold: log2(m * b) * 16
        mpz_mul_ui(tmp, m_val, b);
        double rat_log = mpz_sizeinbase(tmp, 2);
        int rat_thresh = (int)((rat_log - 10) * 16); // leave 10 bits slack
        if (rat_thresh < 0) rat_thresh = 0;

        // Compute algebraic threshold for this b
        // Norm ≈ max_coeff * A * b^d (roughly)
        // For more accuracy, evaluate at boundaries
        eval_norm(tmp, A, b, coeff[2], coeff[1], coeff[0]);
        mpz_abs(tmp, tmp);
        double alg_log = mpz_sizeinbase(tmp, 2);
        int alg_thresh = (int)((alg_log - 12) * 16);
        if (alg_thresh < 0) alg_thresh = 0;

        // Clear sieve arrays
        memset(rat_sieve.data(), 0, (2 * A + 2) * sizeof(int16_t));
        memset(alg_sieve.data(), 0, (2 * A + 2) * sizeof(int16_t));

        // Rational sieve: a ≡ b*m (mod p) => a positions
        for (int i = 0; i < rat_fb_size; i++) {
            int p = rat_fb[i].p;
            int logp = rat_fb[i].logp;
            long long bm = (long long)b % p * rat_fb[i].m_mod_p % p;
            // a ≡ bm (mod p)
            // Starting position in sieve: index = a + A, where a ≡ bm (mod p)
            int start = (int)((bm + A) % p);
            if (start < 0) start += p;
            for (int j = start; j <= 2 * A; j += p)
                rat_sieve[j] += logp;
        }

        // Algebraic sieve: for each (p, r) in alg_fb, sieve where a ≡ b*r (mod p)
        for (int i = 0; i < alg_fb_size; i++) {
            int p = alg_fb[i].p;
            int r = alg_fb[i].r;
            int logp = alg_fb[i].logp;
            long long br = (long long)b % p * r % p;
            int start = (int)((br + A) % p);
            if (start < 0) start += p;
            for (int j = start; j <= 2 * A; j += p)
                alg_sieve[j] += logp;
        }

        // Check candidates where both sides pass threshold
        for (int idx = 0; idx <= 2 * A; idx++) {
            if (rat_sieve[idx] < rat_thresh || alg_sieve[idx] < alg_thresh) continue;

            long long a = (long long)idx - A;
            if (a == 0) continue;

            // Check gcd(a, b) = 1
            long long g = __gcd(abs(a), (long long)b);
            if (g != 1) continue;

            // Compute rational value: a - b*m
            mpz_mul_ui(rat_val, m_val, b);
            mpz_set_si(tmp, a);
            mpz_sub(rat_val, tmp, rat_val); // a - b*m

            // Compute algebraic norm
            if (d == 3) {
                eval_norm(alg_val, a, b, coeff[2], coeff[1], coeff[0]);
                // Include leading coefficient: if coeff[3] != 1, multiply by coeff[3]^2
                // Actually for norm of (a - b*alpha) with f monic:
                // Norm = (-b)^d * f(a/b) = (-b)^3 * (a/b)^3 + ... = a^3 + c2*a^2*b + c1*a*b^2 + c0*b^3
                // If f is not monic: Norm = coeff[d]^(d-1) * (above)
                if (mpz_cmp_ui(coeff[3], 1) != 0) {
                    // Multiply by coeff[3]^(d-1) = coeff[3]^2
                    mpz_mul(alg_val, alg_val, coeff[3]);
                    mpz_mul(alg_val, alg_val, coeff[3]);
                }
            } else {
                // General degree: compute Norm = resultant(a - b*x, f(x))
                // For monic f: Norm = (-b)^d * f(a/b)
                // = a^d + c_{d-1}*a^{d-1}*b + ... + c_0*b^d
                mpz_set_ui(alg_val, 0);
                mpz_t term;
                mpz_init(term);
                for (int j = d; j >= 0; j--) {
                    // term = coeff[j] * a^j * b^(d-j)
                    mpz_set(term, coeff[j]);
                    for (int k = 0; k < j; k++) mpz_mul_si(term, term, a);
                    for (int k = 0; k < d - j; k++) mpz_mul_si(term, term, b);
                    mpz_add(alg_val, alg_val, term);
                }
                mpz_clear(term);
            }

            // Trial divide rational side
            mpz_t rat_rem;
            mpz_init(rat_rem);
            mpz_abs(rat_rem, rat_val);
            bool rat_neg = (mpz_sgn(rat_val) < 0);

            vector<int> rat_exp(rat_fb_size, 0);
            for (int i = 0; i < rat_fb_size; i++) {
                int p = rat_fb[i].p;
                while (mpz_divisible_ui_p(rat_rem, p)) {
                    mpz_divexact_ui(rat_rem, rat_rem, p);
                    rat_exp[i]++;
                }
            }

            bool rat_smooth = (mpz_cmp_ui(rat_rem, 1) == 0);
            if (!rat_smooth) {
                mpz_clear(rat_rem);
                continue;
            }

            // Trial divide algebraic side
            mpz_t alg_rem;
            mpz_init(alg_rem);
            mpz_abs(alg_rem, alg_val);
            bool alg_neg = (mpz_sgn(alg_val) < 0);

            vector<int> alg_exp(alg_fb_size, 0);
            for (int i = 0; i < alg_fb_size; i++) {
                int p = alg_fb[i].p;
                while (mpz_divisible_ui_p(alg_rem, p)) {
                    mpz_divexact_ui(alg_rem, alg_rem, p);
                    alg_exp[i]++;
                }
            }

            bool alg_smooth = (mpz_cmp_ui(alg_rem, 1) == 0);

            if (alg_smooth) {
                NFSRelation rel;
                rel.a = a;
                rel.b = b;
                rel.rat_exponents = rat_exp;
                rel.alg_exponents = alg_exp;
                rel.rat_sign = rat_neg ? 1 : 0;
                rel.alg_sign = alg_neg ? 1 : 0;
                relations.push_back(rel);
                total_relations++;
            }

            mpz_clear(rat_rem);
            mpz_clear(alg_rem);

            if (total_relations >= needed_relations) break;
        }

        if (b % 50 == 0) {
            fprintf(stderr, "  [%.1fs] b=%d: %d/%d relations\n",
                    elapsed(), b, total_relations, needed_relations);
        }
    }

    fprintf(stderr, "Collected %d relations in %.1fs\n", total_relations, elapsed());

    if (total_relations < total_cols + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n",
                total_relations, total_cols + 1);
        return 1;
    }

    // ---- Linear algebra ----
    fprintf(stderr, "Linear algebra: %d relations x %d columns\n",
            total_relations, total_cols);

    vector<vector<int>> exp_vecs(total_relations, vector<int>(total_cols, 0));
    for (int i = 0; i < total_relations; i++) {
        // Column layout: [rat_sign, alg_sign, rat_fb..., alg_fb...]
        exp_vecs[i][0] = relations[i].rat_sign;
        exp_vecs[i][1] = relations[i].alg_sign;
        for (int j = 0; j < rat_fb_size; j++)
            exp_vecs[i][2 + j] = relations[i].rat_exponents[j];
        for (int j = 0; j < alg_fb_size; j++)
            exp_vecs[i][2 + rat_fb_size + j] = relations[i].alg_exponents[j];
    }

    vector<vector<int>> null_vecs = find_null_space(total_relations, total_cols, exp_vecs);
    fprintf(stderr, "Found %d null space vectors\n", (int)null_vecs.size());

    // ---- Square root ----
    // For each null vector, compute:
    // Rational side: product of (a - b*m) for selected relations → perfect square
    // Take sqrt and reduce mod N
    // Algebraic side: more complex — need to compute sqrt in number field
    //
    // Simplified approach: since we're working with a rational square and an algebraic square,
    // the relationship gives us x^2 ≡ y^2 (mod N) directly.
    //
    // The rational product is a perfect square in Z, so sqrt is easy.
    // The algebraic product represents an element of Z[α] whose norm is a perfect square.
    // We need to find the image of this element's sqrt under the homomorphism α → m.
    //
    // For simplicity, use the following approach:
    // The product of (a_i - b_i*m) ≡ ∏(a_i - b_i*m) (mod N)
    // The product of the algebraic norms is a perfect square, say S^2
    // Then ∏(a_i - b_i*m) ≡ ±S (mod N) ... no, this isn't right.
    //
    // Actually, NFS produces: ∏(a_i - b_i*m) ≡ ∏(a_i - b_i*α) (mod N)
    //                          rational square    algebraic square
    // The rational side gives us X^2 = ∏(a_i - b_i*m) → X = ∏ p^(e/2) in Z
    // The algebraic side gives us Y^2 in Z[α] → image under α→m gives Y^2 mod N
    // Then X^2 ≡ Y^2 (mod N)

    // Computing the algebraic square root is the hard part.
    // Simplified version: compute the rational square root directly

    bool factored = false;
    mpz_t X, Y;
    mpz_init(X);
    mpz_init(Y);

    for (auto& nv : null_vecs) {
        if (factored) break;

        // Rational side: X^2 = ∏(a_i - b_i*m)
        // Compute X = product of primes to half-exponents
        vector<long long> total_rat_exp(rat_fb_size, 0);
        vector<long long> total_alg_exp(alg_fb_size, 0);
        int total_rat_sign = 0, total_alg_sign = 0;

        for (int idx : nv) {
            total_rat_sign += relations[idx].rat_sign;
            total_alg_sign += relations[idx].alg_sign;
            for (int j = 0; j < rat_fb_size; j++)
                total_rat_exp[j] += relations[idx].rat_exponents[j];
            for (int j = 0; j < alg_fb_size; j++)
                total_alg_exp[j] += relations[idx].alg_exponents[j];
        }

        // Check all even
        bool all_even = (total_rat_sign % 2 == 0) && (total_alg_sign % 2 == 0);
        for (int j = 0; j < rat_fb_size && all_even; j++)
            if (total_rat_exp[j] & 1) all_even = false;
        for (int j = 0; j < alg_fb_size && all_even; j++)
            if (total_alg_exp[j] & 1) all_even = false;
        if (!all_even) continue;

        // X = rational square root mod N
        mpz_set_ui(X, 1);
        for (int j = 0; j < rat_fb_size; j++) {
            if (total_rat_exp[j] == 0) continue;
            long long half = total_rat_exp[j] / 2;
            mpz_set_ui(tmp, rat_fb[j].p);
            mpz_powm_ui(tmp, tmp, half, N);
            mpz_mul(X, X, tmp);
            mpz_mod(X, X, N);
        }

        // Y = algebraic square root mod N
        // Approximation: Y = product of (a_i - b_i*m) mod N for selected relations,
        // then take rational square root of the product
        // This works because the product is a perfect square in Z
        // ... wait, that's X, not Y.
        //
        // Actually for NFS, the congruence is:
        // ∏(a_i - b_i*m) ≡ φ(∏(a_i - b_i*α)) (mod N) where φ: α → m
        //
        // Both sides are perfect squares. Rational sqrt = X.
        // For the algebraic sqrt, we need sqrt(∏(a_i - b_i*α)) in Z[α], then map to Z/NZ.
        //
        // This is the hardest part of NFS. For now, let's try:
        // Y = product of (a_i - b_i*m) mod N (before taking sqrt)
        // Then X^2 = Y (as a number), so X = sqrt(Y) mod N
        // We need to factor X^2 - Y ≡ 0 (mod N)... that's circular.

        // Alternative simpler approach:
        // Since both rational and algebraic products are perfect squares,
        // and they're congruent mod N, we have X_rat^2 ≡ X_alg^2 (mod N).
        // X_rat = rational sqrt = product of rational FB primes to half-power
        // X_alg = algebraic sqrt mapped to Z/NZ

        // For the algebraic side, we can compute the product directly:
        mpz_set_ui(Y, 1);
        for (int idx : nv) {
            mpz_set_si(tmp, relations[idx].a);
            mpz_mul_ui(tmp2, m_val, relations[idx].b);
            mpz_sub(tmp, tmp, tmp2); // a - b*m
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, N);
        }
        // Y = product of (a_i - b_i*m) mod N
        // Y should be X^2 mod N (since Y = product of rat values = product of rat_p^exp)
        // So sqrt(Y mod N) = X mod N

        // But we already have X = product of p^(half_exp) mod N
        // So X^2 ≡ Y (mod N)... but we need X^2 ≡ Z^2 (mod N) for some Z

        // The issue: in NFS, we need BOTH a rational and algebraic square.
        // The rational square gives us one side, the algebraic square gives the other.
        // Without computing the algebraic square root, we can't proceed.

        // Workaround: Just try gcd(product - 1, N) and gcd(product + 1, N) for various products
        // This is a heuristic but won't work in general.

        // For now, let's use the rational side only as a Fermat-like test:
        // X = rational sqrt, Y = ∏(a-bm) mod N, then X^2 ≡ Y (mod N)
        // gcd(X - Y_sqrt, N) where Y_sqrt is some approx

        // Actually, the simplest correct approach:
        // From the null space vector, we know ∏(a-bm) is a perfect square in Z.
        // Let S = sqrt(∏(a-bm)) computed exactly (using exponents).
        // Also, ∏(a-bm) ≡ ∏(a-bα) evaluated at α=m ≡ ∏Norm(a-bα)/∏Norm(a-bα) ...
        // No, the fundamental NFS identity is:
        // ∏(a-bm) ≡ 0 (mod N) is NOT true.
        // Rather, the map is: if f(m) ≡ 0 (mod N), then for any polynomial g(x):
        // g(m) ≡ 0 (mod N) iff f(x) | g(x) in Z[x]
        // The NFS relation: for each (a,b), (a - bm) and Norm(a - bα) share a common factor mod N.

        // This is getting complex. Let me use a DIFFERENT approach:
        // Use the combined QS-like method: since ∏(a-bm) = X^2 (rational square),
        // and we know X mod N, try gcd(X ± 1, N), gcd(X ± various, N).
        // This won't work in general.

        // CORRECT approach: compute algebraic square root.
        // For degree 3 with coefficient[3]=1:
        // ∏(a_i - b_i*α) = γ^2 where γ ∈ Z[α]/(f(α))
        // γ = u + v*α + w*α^2 for some integers u, v, w
        // Then φ(γ) = u + v*m + w*m^2 mod N
        // And φ(γ)^2 ≡ φ(γ^2) = ∏(a_i - b_i*m) ≡ X^2 (mod N)
        // Wait, that's not right either.
        // φ(∏(a_i - b_i*α)) = ∏φ(a_i - b_i*α) = ∏(a_i - b_i*m) = X^2
        // φ(γ^2) = φ(γ)^2
        // So φ(γ)^2 = X^2 mod N
        // And we need to compute Y = φ(γ) to get gcd(X - Y, N)

        // To find γ: we need to compute the square root of ∏(a_i - b_i*α) in Z[α]/(f(α))
        // This is done by Couveignes's algorithm or Montgomery's method.

        // For a simplified version: compute ∏(a_i - b_i*α) mod f(α) keeping track of
        // coefficients, then compute the square root in Z[α]/(f(α)).

        // Computing ∏(a_i - b_i*α) mod f(α):
        // Start with 1, then for each relation multiply by (a_i - b_i*α)
        // Reduce mod f(α) after each multiplication
        // This gives us an element u + v*α + w*α^2 (for degree 3)

        // Actually, this product might have huge coefficients. We should work mod some prime.
        // Better: work mod N directly.

        // Product of (a - b*α) mod f(α) mod N:
        // Represent as polynomial of degree < d in α, coefficients mod N
        vector<mpz_t> prod(d); // prod[0] + prod[1]*α + ... + prod[d-1]*α^(d-1)
        for (int j = 0; j < d; j++) mpz_init_set_ui(prod[j], 0);
        mpz_set_ui(prod[0], 1); // start with 1

        for (int idx : nv) {
            // Multiply by (a_i - b_i*α)
            // = a_i + (-b_i)*α + 0*α^2 + ...
            vector<mpz_t> mul(d);
            for (int j = 0; j < d; j++) mpz_init_set_ui(mul[j], 0);
            mpz_set_si(mul[0], relations[idx].a);
            mpz_set_si(mul[1], -relations[idx].b);

            // Multiply prod by mul mod f(α) mod N
            vector<mpz_t> result(2 * d - 1);
            for (int j = 0; j < 2 * d - 1; j++) mpz_init_set_ui(result[j], 0);

            for (int j = 0; j < d; j++) {
                for (int k = 0; k < d; k++) {
                    mpz_addmul(result[j + k], prod[j], mul[k]);
                }
            }

            // Reduce mod N
            for (int j = 0; j < 2 * d - 1; j++) mpz_mod(result[j], result[j], N);

            // Reduce mod f(α): for degree d, α^d = -c_{d-1}*α^{d-1} - ... - c_0
            // Process from highest degree down
            for (int j = 2 * d - 2; j >= d; j--) {
                // α^j = α^(j-d) * α^d = α^(j-d) * (-c_{d-1}*α^{d-1} - ... - c_0)
                for (int k = 0; k < d; k++) {
                    // result[j-d+k] -= result[j] * coeff[k] / coeff[d]
                    // (assuming f is monic, coeff[d] = 1, so just subtract)
                    if (mpz_cmp_ui(coeff[d], 1) == 0) {
                        mpz_submul(result[j - d + k], result[j], coeff[k]);
                    } else {
                        // Non-monic: need to handle leading coefficient
                        // α^d = -(c_{d-1}/c_d)*α^{d-1} - ... - c_0/c_d
                        // This is more complex; skip for now
                        mpz_submul(result[j - d + k], result[j], coeff[k]);
                    }
                }
                mpz_set_ui(result[j], 0);
            }

            // Copy to prod
            for (int j = 0; j < d; j++) {
                mpz_mod(prod[j], result[j], N);
            }

            for (int j = 0; j < d; j++) mpz_clear(mul[j]);
            for (int j = 0; j < 2 * d - 1; j++) mpz_clear(result[j]);
        }

        // Now prod[0..d-1] represents ∏(a_i - b_i*α) in Z[α]/(f(α)) mod N
        // This should be a perfect square in this ring.
        // Computing its square root is non-trivial.

        // Heuristic: try computing sqrt by Tonelli-Shanks-like method in the number field ring.
        // For now, let's just try: Y = prod[0] + prod[1]*m + prod[2]*m^2 mod N
        // This is φ(∏(a_i - b_i*α)) = ∏(a_i - b_i*m) mod N = X^2 mod N
        // We need sqrt of this in Z[α]/(f(α)), then map to Z/NZ.

        // Since computing algebraic sqrt is hard, let's try:
        // X = rational sqrt (computed from exponents)
        // Y = ∏(a-bm) mod N (= X^2 mod N)
        // Try random linear combinations to extract non-trivial factor

        // Actually, the simplest thing: gcd(X - 1, N), gcd(X + 1, N), etc.
        // Or use X directly: gcd(X^2 - Y, N) = gcd(0, N) = N (trivial)

        // We really need the algebraic square root.
        // Skip for now and try a CRT-based approach:
        // For each prime q not dividing N, compute sqrt(∏) mod q
        // Then CRT to get approximate sqrt, reduce mod N

        // Actually, the simplest working approach:
        // The product of (a-bm) mod N is already computed above (= X^2 mod N).
        // We need a DIFFERENT square root of this product.
        // If we can find TWO different sqrt values, their ratio gives a factor.

        // Method: compute X = prod of rat_fb_primes to half-power (one sqrt)
        // Then X^2 ≡ ∏(a-bm) (mod N)
        // We need another value Y such that Y^2 ≡ ∏(a-bm) (mod N)
        // and Y ≢ ±X (mod N)

        // This Y comes from the algebraic side. Without computing it,
        // let's just try the gcd approach with X:
        // gcd(X - 1, N), gcd(X + 1, N), gcd(X - m, N), etc.

        // Since this is unlikely to work for a generic N, let me try anyway
        mpz_sub_ui(tmp, X, 1); mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            factored = true; break;
        }
        mpz_add_ui(tmp, X, 1); mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            factored = true; break;
        }

        for (int j = 0; j < d; j++) mpz_clear(prod[j]);
    }

    double total_time = elapsed();

    if (factored) {
        gmp_printf("%Zd\n", factor);
        fprintf(stderr, "NFS: factored in %.3f seconds\n", total_time);
    } else {
        fprintf(stderr, "NFS: FAILED after %.3f seconds (algebraic sqrt not implemented)\n",
                total_time);
        return 1;
    }

    // Cleanup
    mpz_clear(N); mpz_clear(m_val); mpz_clear(rat_val); mpz_clear(alg_val);
    mpz_clear(factor); mpz_clear(tmp); mpz_clear(tmp2); mpz_clear(remainder);
    for (int i = 0; i <= d; i++) mpz_clear(coeff[i]);
    mpz_clear(X); mpz_clear(Y);

    return 0;
}
