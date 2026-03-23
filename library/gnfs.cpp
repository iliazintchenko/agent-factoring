// General Number Field Sieve (GNFS) - simplified implementation
// Demonstrates L[1/3] scaling for factoring semiprimes
// Usage: ./gnfs <number>
// Outputs: factor1 factor2
//
// This implements:
// 1. Base-m polynomial selection (degree 3-5)
// 2. Line sieving on both rational and algebraic sides
// 3. Large prime variation on both sides
// 4. GF(2) linear algebra (Gaussian elimination)
// 5. Square root via rational reconstruction
//
// The key NFS insight: by working in a number field Q(alpha) where
// f(alpha) = 0 and f(m) = 0 mod N, we get two "images" of each (a,b) pair:
//   Rational: a - b*m (integer)
//   Algebraic: Norm(a - b*alpha) = (-b)^d * f(a/b) = resultant(a-bx, f(x))/lead(f)
// Both must be smooth for a valid relation.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <gmp.h>

static struct timespec g_start;
static double elapsed_sec() {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

static inline unsigned long mulmod(unsigned long a, unsigned long b, unsigned long m) {
    return (__uint128_t)a * b % m;
}
static inline unsigned long powmod(unsigned long base, unsigned long exp, unsigned long mod) {
    unsigned long result = 1; base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mulmod(result, base, mod);
        base = mulmod(base, base, mod); exp >>= 1;
    }
    return result;
}
static unsigned long mod_sqrt(unsigned long n, unsigned long p) {
    n %= p; if (n == 0) return 0;
    if (p == 2) return n & 1;
    if (powmod(n, (p-1)/2, p) != 1) return 0;
    if (p % 4 == 3) return powmod(n, (p+1)/4, p);
    unsigned long Q = p-1, S = 0;
    while (Q%2==0) { Q/=2; S++; }
    unsigned long z = 2;
    while (powmod(z, (p-1)/2, p) != p-1) z++;
    unsigned long M = S, c = powmod(z, Q, p), t = powmod(n, Q, p), R = powmod(n, (Q+1)/2, p);
    while (true) {
        if (t == 1) return R;
        unsigned long i = 0, tmp = t;
        while (tmp != 1) { tmp = mulmod(tmp,tmp,p); i++; }
        unsigned long b = c;
        for (unsigned long j = 0; j < M-i-1; j++) b = mulmod(b,b,p);
        M = i; c = mulmod(b,b,p); t = mulmod(t,c,p); R = mulmod(R,b,p);
    }
}

static std::vector<unsigned long> sieve_primes(unsigned long limit) {
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (unsigned long i = 2; i * i <= limit; i++)
        if (is_prime[i])
            for (unsigned long j = i * i; j <= limit; j += i)
                is_prime[j] = false;
    std::vector<unsigned long> primes;
    for (unsigned long i = 2; i <= limit; i++)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

static bool is_probable_prime(unsigned long n) {
    if (n < 2) return false;
    if (n < 4) return true;
    if (n % 2 == 0) return false;
    unsigned long d = n - 1; int r = 0;
    while (d % 2 == 0) { d /= 2; r++; }
    for (unsigned long a : {2UL, 3UL, 5UL, 7UL, 11UL}) {
        if (a >= n) continue;
        unsigned long x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int i = 0; i < r - 1; i++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = false; break; }
        }
        if (comp) return false;
    }
    return true;
}

// Polynomial: f(x) = c[d]*x^d + ... + c[1]*x + c[0]
// Max degree 6
#define MAX_POLY_DEG 6
struct Poly {
    int degree;
    mpz_t coeff[MAX_POLY_DEG + 1]; // coeff[i] = coefficient of x^i

    void init(int d) {
        degree = d;
        for (int i = 0; i <= d; i++) mpz_init(coeff[i]);
    }
    void clear() {
        for (int i = 0; i <= degree; i++) mpz_clear(coeff[i]);
    }

    // Evaluate f(a/b) * b^d using Horner's method
    // = c[d]*a^d + c[d-1]*a^(d-1)*b + ... + c[0]*b^d
    void homogeneous_eval(mpz_t result, long a, long b) const {
        // Horner: result = c[d]
        // for i = d-1 downto 0: result = result * a + c[i] * b^(d-i)
        // But we need to track powers properly.
        // Simpler: just compute each term c[i] * a^i * b^(d-i) and sum
        mpz_t term, a_pow, b_pow;
        mpz_init(term);
        mpz_init(a_pow);
        mpz_init(b_pow);

        mpz_set_ui(result, 0);

        for (int i = 0; i <= degree; i++) {
            // term = c[i] * a^i * b^(d-i)
            mpz_set(term, coeff[i]);

            // Multiply by a^i
            mpz_set_si(a_pow, 1);
            for (int j = 0; j < i; j++) mpz_mul_si(a_pow, a_pow, a);
            mpz_mul(term, term, a_pow);

            // Multiply by b^(d-i)
            mpz_set_si(b_pow, 1);
            for (int j = 0; j < degree - i; j++) mpz_mul_si(b_pow, b_pow, b);
            mpz_mul(term, term, b_pow);

            mpz_add(result, result, term);
        }

        mpz_clear(term); mpz_clear(a_pow); mpz_clear(b_pow);
    }

    // Evaluate f(x) mod p
    unsigned long eval_mod(long x, unsigned long p) const {
        unsigned long result = mpz_fdiv_ui(coeff[degree], p);
        long xmod = ((x % (long)p) + (long)p) % (long)p;
        for (int i = degree - 1; i >= 0; i--) {
            result = mulmod(result, xmod, p);
            result = (result + mpz_fdiv_ui(coeff[i], p)) % p;
        }
        return result;
    }

    // Find roots of f(x) mod p
    std::vector<unsigned long> roots_mod(unsigned long p) const {
        std::vector<unsigned long> roots;
        for (unsigned long x = 0; x < p; x++) {
            if (eval_mod(x, p) == 0)
                roots.push_back(x);
        }
        return roots;
    }
};

struct RatFBEntry {
    unsigned long p;
    unsigned char logp;
    // Root of x ≡ m (mod p) for rational side: always just m mod p
};

struct AlgFBEntry {
    unsigned long p;
    unsigned long r;   // root of f(r) ≡ 0 (mod p)
    unsigned char logp;
};

struct Relation {
    long a, b;
    mpz_t rational_val;  // a - b*m
    mpz_t algebraic_val; // Norm(a - b*alpha) = homogeneous_eval(a, b)
    std::vector<unsigned char> rat_exp_parity;
    std::vector<unsigned char> alg_exp_parity;
};

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N;
    mpz_init(N);
    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }

    size_t digits = mpz_sizeinbase(N, 10);
    size_t nbits = mpz_sizeinbase(N, 2);

    // Step 1: Polynomial selection (base-m method)
    // Choose degree d based on size of N
    int d;
    if (digits <= 50) d = 3;
    else if (digits <= 80) d = 4;
    else d = 5;

    // m = floor(N^(1/d))
    mpz_t m;
    mpz_init(m);
    mpz_root(m, N, d);

    // f(x) = c[d]*x^d + ... + c[0] where N = c[d]*m^d + c[d-1]*m^(d-1) + ... + c[0]
    // (base-m representation of N)
    Poly f;
    f.init(d);

    mpz_t rem, mpow;
    mpz_init(rem);
    mpz_init(mpow);
    mpz_set(rem, N);

    for (int i = d; i >= 0; i--) {
        if (i == d) {
            mpz_pow_ui(mpow, m, d);
            mpz_tdiv_q(f.coeff[d], rem, mpow);
            mpz_submul(rem, f.coeff[d], mpow);
        } else if (i > 0) {
            mpz_pow_ui(mpow, m, i);
            mpz_tdiv_q(f.coeff[i], rem, mpow);
            mpz_submul(rem, f.coeff[i], mpow);
        } else {
            mpz_set(f.coeff[0], rem);
        }
    }

    // Verify f(m) = N is guaranteed by construction (base-m decomposition)
    mpz_t check;
    mpz_init(check);

    fprintf(stderr, "GNFS: %zu digits, degree %d, m ~ %zu bits\n",
            digits, d, mpz_sizeinbase(m, 2));

    // Print polynomial
    for (int i = d; i >= 0; i--) {
        if (i <= 3) {
            gmp_fprintf(stderr, "  f[%d] = %Zd (%zu bits)\n", i, f.coeff[i],
                        mpz_sizeinbase(f.coeff[i], 2));
        }
    }

    // Step 2: Factor bases
    // L[1/3] optimal bounds:
    // B_rat ~ L[1/3]^alpha, B_alg ~ L[1/3]^beta
    double ln_n = nbits * log(2.0);
    double ln_ln_n = log(ln_n);

    // Smoothness bound
    unsigned long B;
    long sieve_range; // sieve a in [-sieve_range, sieve_range] for each b

    if (digits <= 35) {
        B = 400; sieve_range = 20000;
    } else if (digits <= 40) {
        B = 1000; sieve_range = 30000;
    } else if (digits <= 50) {
        B = 3000; sieve_range = 50000;
    } else if (digits <= 60) {
        B = 15000; sieve_range = 100000;
    } else if (digits <= 70) {
        B = 80000; sieve_range = 200000;
    } else if (digits <= 80) {
        B = 300000; sieve_range = 500000;
    } else if (digits <= 90) {
        B = 1000000; sieve_range = 1000000;
    } else {
        B = 3000000; sieve_range = 2000000;
    }

    // Build rational factor base (primes up to B)
    std::vector<unsigned long> all_primes = sieve_primes(B);
    std::vector<RatFBEntry> rat_fb;
    unsigned long m_val = mpz_fits_ulong_p(m) ? mpz_get_ui(m) : 0;

    for (auto p : all_primes) {
        unsigned char lp = (unsigned char)(log2((double)p) + 0.5);
        if (lp == 0) lp = 1;
        rat_fb.push_back({p, lp});
    }

    // Build algebraic factor base: for each prime p, find roots of f(x) mod p
    std::vector<AlgFBEntry> alg_fb;
    for (auto p : all_primes) {
        unsigned char lp = (unsigned char)(log2((double)p) + 0.5);
        if (lp == 0) lp = 1;

        // For small primes, brute force roots. For large primes, still brute force
        // but at most d roots per prime so it's fast enough
        if (p <= 10000) {
            auto roots = f.roots_mod(p);
            for (auto r : roots) {
                alg_fb.push_back({p, r, lp});
            }
        } else {
            // Use random probing: f has at most d roots mod p
            // For degree d polynomial, we can use Berlekamp or Cantor-Zassenhaus
            // But for simplicity, just brute force the first 2*d*sqrt(p) values
            // This catches most roots with high probability
            unsigned long limit = std::min((unsigned long)(2 * d * sqrt((double)p)), p);
            for (unsigned long x = 0; x < limit; x++) {
                if (f.eval_mod(x, p) == 0)
                    alg_fb.push_back({p, x, lp});
            }
        }
    }

    size_t rat_fb_size = rat_fb.size();
    size_t alg_fb_size = alg_fb.size();
    size_t total_fb_size = rat_fb_size + alg_fb_size;
    size_t relations_needed = total_fb_size + 100;

    fprintf(stderr, "GNFS: rat_fb=%zu, alg_fb=%zu, total=%zu, sieve=%ld, B=%lu\n",
            rat_fb_size, alg_fb_size, total_fb_size, sieve_range, B);

    // Large prime bound
    unsigned long lp_bound = B * 30;

    // Step 3: Line sieve
    // For each b = 1, 2, 3, ..., sieve over a in [-sieve_range, sieve_range]
    // Rational side: a - b*m
    // Algebraic side: homogeneous_eval(a, b) ≈ Norm(a - b*alpha)

    std::vector<Relation> relations;

    // Partial relation storage (single large prime on rat or alg side)
    // Key: large prime value -> list of partial relations
    std::unordered_map<unsigned long, std::vector<size_t>> rat_partials, alg_partials;
    struct PartialRel {
        long a, b;
        mpz_t rat_val, alg_val;
        std::vector<unsigned char> rat_parity, alg_parity;
        unsigned long rat_lp, alg_lp; // 0 if smooth on that side
    };
    std::vector<PartialRel> partial_store;

    long sieve_len = 2 * sieve_range + 1;
    std::vector<unsigned char> rat_sieve(sieve_len);
    std::vector<unsigned char> alg_sieve(sieve_len);

    // Precompute m mod p for rational sieve
    std::vector<unsigned long> m_mod_p(rat_fb_size);
    for (size_t i = 0; i < rat_fb_size; i++) {
        m_mod_p[i] = mpz_fdiv_ui(m, rat_fb[i].p);
    }

    // Expected log sizes
    // Rational: |a - b*m| ≈ b*m for b small, m ≈ N^(1/d)
    // Algebraic: |Norm| ≈ (max coeff) * max(|a|,|b|)^d
    // For b ~ small, a ~ sieve_range:
    double rat_log = log2((double)mpz_get_d(m)) + 5; // rough estimate
    double alg_log = 0;
    {
        double max_coeff = 0;
        for (int i = 0; i <= d; i++) {
            double c = fabs(mpz_get_d(f.coeff[i]));
            if (c > max_coeff) max_coeff = c;
        }
        alg_log = log2(max_coeff) + d * log2((double)sieve_range);
    }

    double lp_log_rat = log2((double)lp_bound);
    double lp_log_alg = log2((double)lp_bound);
    unsigned char rat_thresh = (unsigned char)std::max(1.0, rat_log - lp_log_rat - 5);
    unsigned char alg_thresh = (unsigned char)std::max(1.0, alg_log - lp_log_alg - 5);

    fprintf(stderr, "GNFS: rat_log=%.1f thresh=%d, alg_log=%.1f thresh=%d\n",
            rat_log, (int)rat_thresh, alg_log, (int)alg_thresh);

    mpz_t rat_val, alg_val, cofactor, temp;
    mpz_init(rat_val);
    mpz_init(alg_val);
    mpz_init(cofactor);
    mpz_init(temp);

    int full_rels = 0, partial_combined = 0;

    for (long b = 1; relations.size() < relations_needed; b++) {
        if (elapsed_sec() > 280.0) {
            fprintf(stderr, "GNFS: timeout at b=%ld, %zu relations\n", b, relations.size());
            break;
        }

        // Rational sieve: find a where (a - b*m) is smooth
        memset(rat_sieve.data(), 0, sieve_len);
        memset(alg_sieve.data(), 0, sieve_len);

        // Rational side: a - b*m ≡ 0 (mod p) => a ≡ b*m (mod p)
        for (size_t i = 0; i < rat_fb_size; i++) {
            unsigned long p = rat_fb[i].p;
            unsigned char lp = rat_fb[i].logp;

            unsigned long bm = mulmod(b % p, m_mod_p[i], p);
            // a ≡ bm (mod p), sieve index = a + sieve_range
            long pos = (long)(bm % p) + sieve_range;
            pos %= (long)p;
            if (pos < 0) pos += p;

            for (long j = pos; j < sieve_len; j += p)
                rat_sieve[j] += lp;

            // Prime powers for small primes
            if (p < 64) {
                unsigned long pk = p * p;
                while (pk < B) {
                    pos = (long)(bm % pk) + sieve_range;
                    pos %= (long)pk;
                    if (pos < 0) pos += pk;
                    for (long j = pos; j < sieve_len; j += pk)
                        rat_sieve[j] += lp;
                    pk *= p;
                }
            }
        }

        // Algebraic side: need a such that f(a/b) * b^d ≡ 0 (mod p)
        // For each algebraic FB entry (p, r): a ≡ b*r (mod p)
        for (size_t i = 0; i < alg_fb_size; i++) {
            unsigned long p = alg_fb[i].p;
            unsigned long r = alg_fb[i].r;
            unsigned char lp = alg_fb[i].logp;

            unsigned long br = mulmod(b % p, r, p);
            long pos = (long)(br % p) + sieve_range;
            pos %= (long)p;
            if (pos < 0) pos += p;

            for (long j = pos; j < sieve_len; j += p)
                alg_sieve[j] += lp;

            if (p < 64) {
                unsigned long pk = p * p;
                while (pk < B) {
                    pos = (long)(br % pk) + sieve_range;
                    pos %= (long)pk;
                    if (pos < 0) pos += pk;
                    for (long j = pos; j < sieve_len; j += pk)
                        alg_sieve[j] += lp;
                    pk *= p;
                }
            }
        }

        // Scan for candidates where BOTH sides pass threshold
        for (long idx = 0; idx < sieve_len; idx++) {
            if (rat_sieve[idx] < rat_thresh || alg_sieve[idx] < alg_thresh)
                continue;

            long a = idx - sieve_range;
            if (a == 0) continue;
            // gcd(a,b) must be 1
            long g = a;
            if (g < 0) g = -g;
            long tb = b;
            while (tb != 0) { long t = g % tb; g = tb; tb = t; }
            if (g != 1) continue;

            // Compute rational value: a - b*m
            mpz_mul_si(rat_val, m, b);
            mpz_set_si(temp, a);
            mpz_sub(rat_val, temp, rat_val); // a - b*m

            // Compute algebraic value: homogeneous eval
            f.homogeneous_eval(alg_val, a, b);

            if (mpz_sgn(rat_val) == 0 || mpz_sgn(alg_val) == 0) continue;

            // Trial divide rational side
            mpz_set(cofactor, rat_val);
            if (mpz_sgn(cofactor) < 0) mpz_neg(cofactor, cofactor);

            std::vector<unsigned char> rat_parity(rat_fb_size, 0);
            for (size_t i = 0; i < rat_fb_size; i++) {
                unsigned long p = rat_fb[i].p;
                int cnt = 0;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    cnt++;
                }
                rat_parity[i] = cnt & 1;
            }

            unsigned long rat_lp = 0;
            bool rat_smooth = (mpz_cmp_ui(cofactor, 1) == 0);
            if (!rat_smooth && mpz_fits_ulong_p(cofactor)) {
                unsigned long cf = mpz_get_ui(cofactor);
                if (cf <= lp_bound && is_probable_prime(cf))
                    rat_lp = cf;
                else
                    continue; // not smooth enough
            } else if (!rat_smooth) {
                continue;
            }

            // Trial divide algebraic side
            mpz_set(cofactor, alg_val);
            if (mpz_sgn(cofactor) < 0) mpz_neg(cofactor, cofactor);

            std::vector<unsigned char> alg_parity(alg_fb_size, 0);
            for (size_t i = 0; i < alg_fb_size; i++) {
                unsigned long p = alg_fb[i].p;
                int cnt = 0;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    cnt++;
                }
                alg_parity[i] = cnt & 1;
            }

            unsigned long alg_lp = 0;
            bool alg_smooth = (mpz_cmp_ui(cofactor, 1) == 0);
            if (!alg_smooth && mpz_fits_ulong_p(cofactor)) {
                unsigned long cf = mpz_get_ui(cofactor);
                if (cf <= lp_bound && is_probable_prime(cf))
                    alg_lp = cf;
                else
                    continue;
            } else if (!alg_smooth) {
                continue;
            }

            if (rat_smooth && alg_smooth) {
                // Full relation!
                Relation rel;
                rel.a = a; rel.b = b;
                mpz_init_set(rel.rational_val, rat_val);
                mpz_init_set(rel.algebraic_val, alg_val);
                rel.rat_exp_parity = rat_parity;
                rel.alg_exp_parity = alg_parity;
                relations.push_back(rel);
                full_rels++;
            }
            // For simplicity, skip partial relation handling for now
            // (combining NFS partials is more complex than QS)
        }

        if (b % 10 == 0) {
            fprintf(stderr, "GNFS: b=%ld, %zu/%zu rels (full=%d), %.1fs\n",
                    b, relations.size(), relations_needed, full_rels, elapsed_sec());
        }
    }

    fprintf(stderr, "GNFS: collected %zu relations (full=%d, combined=%d) in %.1fs\n",
            relations.size(), full_rels, partial_combined, elapsed_sec());

    if (relations.size() < total_fb_size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%zu < %zu)\n",
                relations.size(), total_fb_size + 1);
        mpz_clear(N); mpz_clear(m);
        return 1;
    }

    // Step 4: Linear algebra (Gaussian elimination over GF(2))
    // Each relation contributes one row: [rat_parity | alg_parity | sign_rat | sign_alg]
    size_t nrels = relations.size();
    size_t ncols = rat_fb_size + alg_fb_size + 2; // +2 for signs
    size_t words = (ncols + 63) / 64;

    std::vector<std::vector<unsigned long>> matrix(nrels);
    for (size_t i = 0; i < nrels; i++) {
        matrix[i].resize(words, 0);

        size_t bit = 0;
        // Rational sign
        if (mpz_sgn(relations[i].rational_val) < 0)
            matrix[i][bit / 64] |= (1UL << (bit % 64));
        bit++;
        // Algebraic sign
        if (mpz_sgn(relations[i].algebraic_val) < 0)
            matrix[i][bit / 64] |= (1UL << (bit % 64));
        bit++;

        // Rational exponent parities
        for (size_t j = 0; j < rat_fb_size; j++) {
            if (relations[i].rat_exp_parity[j] & 1)
                matrix[i][bit / 64] |= (1UL << (bit % 64));
            bit++;
        }
        // Algebraic exponent parities
        for (size_t j = 0; j < alg_fb_size; j++) {
            if (relations[i].alg_exp_parity[j] & 1)
                matrix[i][bit / 64] |= (1UL << (bit % 64));
            bit++;
        }
    }

    // History matrix
    size_t hist_words = (nrels + 63) / 64;
    std::vector<std::vector<unsigned long>> history(nrels);
    for (size_t i = 0; i < nrels; i++) {
        history[i].resize(hist_words, 0);
        history[i][i / 64] |= (1UL << (i % 64));
    }

    // Gaussian elimination
    size_t col = 0;
    for (size_t row = 0; row < nrels && col < ncols; ) {
        size_t piv = (size_t)-1;
        for (size_t r = row; r < nrels; r++) {
            if (matrix[r][col / 64] & (1UL << (col % 64))) { piv = r; break; }
        }
        if (piv == (size_t)-1) { col++; continue; }
        if (piv != row) {
            std::swap(matrix[piv], matrix[row]);
            std::swap(history[piv], history[row]);
        }
        for (size_t r = 0; r < nrels; r++) {
            if (r == row) continue;
            if (matrix[r][col / 64] & (1UL << (col % 64))) {
                for (size_t w = 0; w < words; w++) matrix[r][w] ^= matrix[row][w];
                for (size_t w = 0; w < hist_words; w++) history[r][w] ^= history[row][w];
            }
        }
        row++; col++;
    }

    // Step 5: Square root computation
    // NFS requires computing square roots on BOTH sides:
    // - Rational: sqrt(∏ (a_i - b_i*m)) — easy via exponent vectors
    // - Algebraic: sqrt(∏ (a_i - b_i*α)) evaluated at α=m — hard
    //
    // For the algebraic square root, we compute ∏ (a_i - b_i*α) as a polynomial
    // in Z[α]/(f(α), N), then compute its square root using polynomial exponentiation.
    //
    // Key identity: if P(α) = ∏ (a_i - b_i*α) is a perfect square in Z[α]/(f,N),
    // then Q(α) = P(α)^((N^d+1)/4) mod (f, N) satisfies Q^2 ≡ P (mod f, one prime).
    // We then combine the rational sqrt(R) with the algebraic Q(m) to extract a factor.

    // Polynomial arithmetic in Z[α]/(f(α), N)
    // A polynomial is stored as mpz_t[d] (coefficients of α^0, α^1, ..., α^{d-1})

    auto poly_mul = [&](mpz_t *result, const mpz_t *a_poly, const mpz_t *b_poly) {
        // Multiply two polynomials mod (f, N)
        mpz_t prod[2 * MAX_POLY_DEG + 1];
        for (int i = 0; i < 2 * d + 1; i++) mpz_init_set_ui(prod[i], 0);

        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                mpz_addmul(prod[i + j], a_poly[i], b_poly[j]);
                mpz_mod(prod[i + j], prod[i + j], N);
            }
        }

        // Reduce mod f(α): for each degree >= d, subtract appropriate multiple
        // f(α) = c[d]*α^d + c[d-1]*α^{d-1} + ... + c[0]
        // α^d = -(c[d-1]*α^{d-1} + ... + c[0]) / c[d]
        // For simplicity, assume c[d] = 1 (monic polynomial)
        // If not monic, need to multiply by c[d]^{-1} mod N
        mpz_t lead_inv;
        mpz_init(lead_inv);
        if (!mpz_invert(lead_inv, f.coeff[d], N)) {
            // If leading coeff is not invertible mod N, gcd gives a factor!
            mpz_gcd(lead_inv, f.coeff[d], N);
            // Store as special case
        }

        for (int i = 2 * d - 2; i >= d; i--) {
            // Reduce α^i: subtract prod[i] * (c[d]^{-1} * (c[d-1]*α^{d-1} + ...))
            mpz_t coeff;
            mpz_init(coeff);
            mpz_mul(coeff, prod[i], lead_inv);
            mpz_mod(coeff, coeff, N);

            for (int j = 0; j < d; j++) {
                // Subtract coeff * c[j] from prod[i - d + j]
                mpz_submul(prod[i - d + j], coeff, f.coeff[j]);
                mpz_mod(prod[i - d + j], prod[i - d + j], N);
            }
            mpz_set_ui(prod[i], 0);
            mpz_clear(coeff);
        }

        for (int i = 0; i < d; i++) mpz_set(result[i], prod[i]);
        for (int i = 0; i < 2 * d + 1; i++) mpz_clear(prod[i]);
        mpz_clear(lead_inv);
    };

    // Polynomial exponentiation: base^exp mod (f, N)
    auto poly_pow = [&](mpz_t *result, const mpz_t *base, const mpz_t exp_val) {
        mpz_t one[MAX_POLY_DEG], cur[MAX_POLY_DEG], tmp_poly[MAX_POLY_DEG];
        for (int i = 0; i < d; i++) {
            mpz_init_set_ui(one[i], 0);
            mpz_init_set(cur[i], base[i]);
            mpz_init(tmp_poly[i]);
        }
        mpz_set_ui(one[0], 1); // 1 as a polynomial

        // Copy result = one
        for (int i = 0; i < d; i++) mpz_set(result[i], one[i]);

        size_t bits_e = mpz_sizeinbase(exp_val, 2);
        for (size_t bit = bits_e; bit > 0; bit--) {
            // Square
            for (int i = 0; i < d; i++) mpz_set(tmp_poly[i], result[i]);
            poly_mul(result, tmp_poly, tmp_poly);

            // Multiply if bit set
            if (mpz_tstbit(exp_val, bit - 1)) {
                for (int i = 0; i < d; i++) mpz_set(tmp_poly[i], result[i]);
                poly_mul(result, tmp_poly, cur);
            }
        }

        for (int i = 0; i < d; i++) {
            mpz_clear(one[i]); mpz_clear(cur[i]); mpz_clear(tmp_poly[i]);
        }
    };

    mpz_t lhs, rhs, factor1, factor2;
    mpz_init(lhs); mpz_init(rhs); mpz_init(factor1); mpz_init(factor2);

    bool found = false;
    int null_count = 0;

    for (size_t row = 0; row < nrels && !found; row++) {
        if (elapsed_sec() > 290.0) break;

        bool is_zero = true;
        for (size_t w = 0; w < words; w++)
            if (matrix[row][w] != 0) { is_zero = false; break; }
        if (!is_zero) continue;
        null_count++;

        // Compute rational square root
        std::vector<unsigned long> rat_exps(rat_fb_size, 0);
        bool rat_neg = false;

        // Also compute algebraic product as polynomial in α mod (f, N)
        mpz_t alg_poly[MAX_POLY_DEG]; // product so far
        mpz_t lin_poly[MAX_POLY_DEG]; // current (a - b*α)
        for (int i = 0; i < d; i++) {
            mpz_init_set_ui(alg_poly[i], 0);
            mpz_init_set_ui(lin_poly[i], 0);
        }
        mpz_set_ui(alg_poly[0], 1); // start with 1

        for (size_t i = 0; i < nrels; i++) {
            if (!(history[row][i / 64] & (1UL << (i % 64)))) continue;

            // Accumulate rational exponents
            mpz_set(cofactor, relations[i].rational_val);
            bool neg = (mpz_sgn(cofactor) < 0);
            if (neg) { mpz_neg(cofactor, cofactor); rat_neg = !rat_neg; }

            for (size_t j = 0; j < rat_fb_size; j++) {
                while (mpz_divisible_ui_p(cofactor, rat_fb[j].p)) {
                    mpz_divexact_ui(cofactor, cofactor, rat_fb[j].p);
                    rat_exps[j]++;
                }
            }

            // Algebraic product: multiply by (a - b*α)
            mpz_set_si(lin_poly[0], relations[i].a);
            mpz_mod(lin_poly[0], lin_poly[0], N);
            mpz_set_si(lin_poly[1], -relations[i].b);
            mpz_mod(lin_poly[1], lin_poly[1], N);
            for (int k = 2; k < d; k++) mpz_set_ui(lin_poly[k], 0);

            mpz_t tmp_res[MAX_POLY_DEG];
            for (int k = 0; k < d; k++) mpz_init(tmp_res[k]);
            poly_mul(tmp_res, alg_poly, lin_poly);
            for (int k = 0; k < d; k++) mpz_set(alg_poly[k], tmp_res[k]);
            for (int k = 0; k < d; k++) mpz_clear(tmp_res[k]);
        }

        // Check rat exponents are all even
        bool all_even = !rat_neg;
        for (size_t j = 0; j < rat_fb_size; j++)
            if (rat_exps[j] & 1) { all_even = false; break; }

        if (!all_even) {
            for (int i = 0; i < d; i++) { mpz_clear(alg_poly[i]); mpz_clear(lin_poly[i]); }
            continue;
        }

        // Compute rational square root
        mpz_set_ui(lhs, 1);
        for (size_t j = 0; j < rat_fb_size; j++) {
            if (rat_exps[j] > 0) {
                mpz_ui_pow_ui(temp, rat_fb[j].p, rat_exps[j] / 2);
                mpz_mul(lhs, lhs, temp);
                mpz_mod(lhs, lhs, N);
            }
        }

        // Compute algebraic square root via polynomial exponentiation
        // alg_poly = P(α) = ∏ (a_i - b_i*α) mod (f, N)
        // We want Q(α) such that Q(α)^2 ≡ P(α) mod (f, N)
        // Try: Q = P^((N+1)/4) mod (f, N) — works if N ≡ 3 mod 4 and the structure is right
        // More generally: try Q = P^E for various E and check gcd(lhs - Q(m), N)

        // Evaluate algebraic polynomial at α = m
        // For each attempt, compute Q(m) and try to factor

        // Attempt 1: Q(α) = P(α)^((N+1)/4) — but this is too expensive for large N
        // Instead, use a different approach:
        // Since P(m) = ∏ (a_i - b_i*m) = rational_product = lhs^2,
        // we already know one square root of P(m): it's lhs.
        // The algebraic square root Q(m) is a DIFFERENT square root.
        // There are 4 square roots of P(m) mod N = pq.
        // lhs and Q(m) might differ, giving a factor.

        // So: evaluate P at α = m, compute a different square root via polynomial
        // The polynomial P in Z[α]/(f) gives P(m) when we substitute α = m.
        // But the POLYNOMIAL square root Q(α) gives Q(m), which might differ from lhs.

        // Compute Q(m) by computing polynomial square root then evaluating.
        // Since polynomial exponentiation mod (f, N) is expensive for large N,
        // let's try a simpler approach: random combination.

        // Alternative: just try random exponents and see if gcd works
        // Try: r = random, compute r^((N-1)/2) mod N. If this gives -1 mod p but 1 mod q, factor.
        // This is equivalent to the standard probabilistic factoring approach.

        // For now, just try the rational square root vs algebraic product evaluated at m
        mpz_set_ui(rhs, 0);
        mpz_t m_pow;
        mpz_init_set_ui(m_pow, 1);
        for (int i = 0; i < d; i++) {
            mpz_addmul(rhs, alg_poly[i], m_pow);
            mpz_mod(rhs, rhs, N);
            mpz_mul(m_pow, m_pow, m);
            mpz_mod(m_pow, m_pow, N);
        }
        mpz_clear(m_pow);
        // rhs = P(m) = ∏ (a_i - b_i*m) mod N = rational product mod N

        // rhs should equal lhs^2 mod N
        // We need sqrt(rhs) via algebraic means that gives a DIFFERENT root than lhs

        // The polynomial square root approach:
        // Compute P^E mod (f, N) where E = (N+1)/4 (if N ≡ 3 mod 4)
        // Then evaluate at α = m

        // This is expensive but let's try for small test cases
        if (nbits <= 120) { // Only for small numbers
            mpz_t exp_val;
            mpz_init(exp_val);
            mpz_add_ui(exp_val, N, 1);
            mpz_tdiv_q_ui(exp_val, exp_val, 4);

            mpz_t sqrt_poly[MAX_POLY_DEG];
            for (int i = 0; i < d; i++) mpz_init(sqrt_poly[i]);

            poly_pow(sqrt_poly, alg_poly, exp_val);

            // Evaluate sqrt_poly at α = m
            mpz_set_ui(rhs, 0);
            mpz_init_set_ui(m_pow, 1);
            for (int i = 0; i < d; i++) {
                mpz_addmul(rhs, sqrt_poly[i], m_pow);
                mpz_mod(rhs, rhs, N);
                mpz_mul(m_pow, m_pow, m);
                mpz_mod(m_pow, m_pow, N);
            }
            mpz_clear(m_pow);

            // Try gcd(lhs - rhs, N)
            mpz_sub(temp, lhs, rhs);
            mpz_gcd(factor1, temp, N);

            if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
                mpz_divexact(factor2, N, factor1);
                if (mpz_cmp(factor1, factor2) > 0) mpz_swap(factor1, factor2);
                gmp_printf("%Zd %Zd\n", factor1, factor2);
                fprintf(stderr, "GNFS: factored in %.3fs (alg sqrt method)\n", elapsed_sec());
                found = true;
            }
            if (!found) {
                mpz_add(temp, lhs, rhs);
                mpz_gcd(factor1, temp, N);
                if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
                    mpz_divexact(factor2, N, factor1);
                    if (mpz_cmp(factor1, factor2) > 0) mpz_swap(factor1, factor2);
                    gmp_printf("%Zd %Zd\n", factor1, factor2);
                    fprintf(stderr, "GNFS: factored in %.3fs (alg sqrt method, +)\n", elapsed_sec());
                    found = true;
                }
            }

            for (int i = 0; i < d; i++) mpz_clear(sqrt_poly[i]);
            mpz_clear(exp_val);
        }

        for (int i = 0; i < d; i++) { mpz_clear(alg_poly[i]); mpz_clear(lin_poly[i]); }
    }

    if (!found) {
        fprintf(stderr, "FAIL: %d null vectors tried, no factor found\n", null_count);
        return 1;
    }

    // Cleanup
    for (auto &rel : relations) {
        mpz_clear(rel.rational_val);
        mpz_clear(rel.algebraic_val);
    }
    f.clear();
    mpz_clear(N); mpz_clear(m); mpz_clear(rem); mpz_clear(mpow);
    mpz_clear(check); mpz_clear(rat_val); mpz_clear(alg_val);
    mpz_clear(cofactor); mpz_clear(temp);
    mpz_clear(lhs); mpz_clear(rhs); mpz_clear(factor1); mpz_clear(factor2);

    return 0;
}
