// Smooth Number Chain Factoring
// Novel approach: instead of finding individual smooth numbers, find CHAINS
// where consecutive ratios are smooth. This reduces the smoothness requirement
// because we only need each step to be partially smooth.
//
// Key idea: Given x_0, x_1, ..., x_k with x_i^2 ≡ q_i (mod N),
// if q_i/q_{i-1} is smooth for each i, then ∏ (q_i/q_{i-1}) = q_k/q_0
// is the product of smooth numbers. Combined with the x_i values,
// this gives quadratic relations.
//
// More specifically: use continued fraction expansion of sqrt(N).
// The convergents p_k/q_k satisfy p_k^2 - N*q_k^2 = (-1)^k * A_k
// where A_k is small (~2*sqrt(N)). If A_k is smooth, we have a relation.
// But even if individual A_k are not smooth, RATIOS of consecutive A_k
// might be smooth because they share factors.
//
// Usage: ./chain_factor <number>

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

struct Relation {
    mpz_t x;    // x value (x^2 ≡ Qx mod N)
    mpz_t Qx;   // the smooth value
    std::vector<unsigned char> exp_parity;
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
    size_t bits = mpz_sizeinbase(N, 2);

    // Factor base
    double ln_n = bits * log(2.0);
    double ln_ln_n = log(ln_n);
    double L_half = exp(sqrt(ln_n * ln_ln_n));

    // CFRAC: the A_k values are ~2*sqrt(N), same as QS Q(x) values
    // Use similar factor base to QS
    unsigned long B = (unsigned long)(pow(L_half, 0.8 / sqrt(2.0)));
    if (B < 500) B = 500;
    if (B > 5000000) B = 5000000;

    std::vector<unsigned long> all_primes = sieve_primes(B);
    std::vector<unsigned long> fb;
    fb.push_back(2);

    for (size_t i = 1; i < all_primes.size(); i++) {
        unsigned long p = all_primes[i];
        unsigned long n_mod = mpz_fdiv_ui(N, p);
        if (mod_sqrt(n_mod, p) != 0) {
            fb.push_back(p);
        }
    }

    size_t fb_size = fb.size();
    size_t relations_needed = fb_size + 50;

    fprintf(stderr, "CFRAC+Chain: %zu digits, B=%lu, FB=%zu primes\n",
            digits, B, fb_size);

    // Large prime bound
    unsigned long lp_bound = B * 100;

    // Continued fraction expansion of sqrt(N)
    // Using the recurrence:
    // sqrt(N) = a_0 + 1/(a_1 + 1/(a_2 + ...))
    // P_{-1} = 1, P_0 = a_0, P_k = a_k * P_{k-1} + P_{k-2}
    // Q_{-1} = 0, Q_0 = 1, Q_k = a_k * Q_{k-1} + Q_{k-2}
    // (-1)^k * (P_k^2 - N*Q_k^2) = A_k (small positive value)

    mpz_t sqrtN, P_prev, P_curr, Q_prev, Q_curr, a_k;
    mpz_t g_val, r_val, m_val, A_val;
    mpz_init(sqrtN);
    mpz_init(P_prev);
    mpz_init(P_curr);
    mpz_init(Q_prev);
    mpz_init(Q_curr);
    mpz_init(a_k);
    mpz_init(g_val);
    mpz_init(r_val);
    mpz_init(m_val);
    mpz_init(A_val);

    mpz_sqrt(sqrtN, N);

    // Initialize CF expansion
    // g_0 = sqrtN, a_0 = floor(sqrt(N))
    // r_0 = sqrtN, m_0 = 0, d_0 = 1
    mpz_t d_val;
    mpz_init(d_val);

    mpz_set(a_k, sqrtN);
    mpz_set_ui(m_val, 0);
    mpz_set_ui(d_val, 1);

    // P_{-1} = 1, P_0 = a_0
    mpz_set_ui(P_prev, 1);
    mpz_set(P_curr, a_k);

    // Q_{-1} = 0, Q_0 = 1
    mpz_set_ui(Q_prev, 0);
    mpz_set_ui(Q_curr, 1);

    std::vector<Relation> relations;

    // Partial relations (single large prime)
    std::unordered_map<unsigned long, std::vector<size_t>> partials;
    struct PartialRel {
        mpz_t x, Qx;
        std::vector<unsigned char> exp_parity;
        unsigned long lp;
    };
    std::vector<PartialRel> partial_store;

    mpz_t temp, cofactor;
    mpz_init(temp);
    mpz_init(cofactor);

    int full_rels = 0, combined_rels = 0;
    long k = 0;

    while (relations.size() < relations_needed && elapsed_sec() < 280.0) {
        k++;

        // Compute next CF partial quotient using the standard algorithm:
        // m_{k+1} = d_k * a_k - m_k
        // d_{k+1} = (N - m_{k+1}^2) / d_k
        // a_{k+1} = floor((sqrtN + m_{k+1}) / d_{k+1})

        // m = d * a - m
        mpz_mul(temp, d_val, a_k);
        mpz_sub(temp, temp, m_val);
        mpz_set(m_val, temp);

        // d = (N - m^2) / d
        mpz_mul(temp, m_val, m_val);
        mpz_sub(temp, N, temp);
        mpz_divexact(temp, temp, d_val);
        mpz_set(d_val, temp);

        // a = floor((sqrtN + m) / d)
        mpz_add(temp, sqrtN, m_val);
        mpz_tdiv_q(a_k, temp, d_val);

        // Update convergents
        // P_{k+1} = a_{k+1} * P_k + P_{k-1}
        mpz_mul(temp, a_k, P_curr);
        mpz_add(temp, temp, P_prev);
        mpz_set(P_prev, P_curr);
        mpz_set(P_curr, temp);

        // Q_{k+1} = a_{k+1} * Q_k + Q_{k-1}
        mpz_mul(temp, a_k, Q_curr);
        mpz_add(temp, temp, Q_prev);
        mpz_set(Q_prev, Q_curr);
        mpz_set(Q_curr, temp);

        // A_k = (-1)^k * (P_k^2 - N * Q_k^2)
        // But it's easier to compute A_k = d_val (the denominator in the CF)
        // Actually: the key identity for CFRAC is:
        // P_k^2 ≡ (-1)^{k+1} * A_k (mod N)
        // where A_k = d_val at step k

        // P_k^2 mod N
        mpz_t x_val;
        mpz_init(x_val);
        mpz_mod(x_val, P_curr, N);

        // Q(x) = P_k^2 mod N, which should equal ±d_val
        // Actually: P_k^2 - N*Q_k^2 = (-1)^k * product of d_1...d_k-related terms
        // The value we test for smoothness is d_val (which is small, ~2*sqrt(N))

        // Actually in CFRAC, the relation is:
        // P_k^2 ≡ (-1)^{k+1} * (d_1 * d_2 * ... * d_k) (mod N)
        // No, that's not right either. Let me use the standard formulation:
        // P_k^2 - N * Q_k^2 = (-1)^{k+1} * A_{k+1}
        // where A_{k+1} is the (k+1)-th "complete quotient denominator"

        // Compute A = P^2 mod N (as signed value in [-N/2, N/2])
        mpz_mul(A_val, P_curr, P_curr);
        mpz_mod(A_val, A_val, N);
        // Make it the small residue
        mpz_t half_N;
        mpz_init(half_N);
        mpz_tdiv_q_ui(half_N, N, 2);
        if (mpz_cmp(A_val, half_N) > 0)
            mpz_sub(A_val, A_val, N);
        mpz_clear(half_N);

        // A_val should be small (~2*sqrt(N) at most)
        // Trial divide for smoothness
        mpz_set(cofactor, A_val);
        bool is_neg = (mpz_sgn(cofactor) < 0);
        if (is_neg) mpz_neg(cofactor, cofactor);

        // Skip if too large (shouldn't happen for CF but check)
        if (mpz_sizeinbase(cofactor, 2) > bits / 2 + 20) {
            mpz_clear(x_val);
            continue;
        }

        std::vector<unsigned char> exp_par(fb_size, 0);
        for (size_t i = 0; i < fb_size; i++) {
            unsigned long p = fb[i];
            int cnt = 0;
            while (mpz_divisible_ui_p(cofactor, p)) {
                mpz_divexact_ui(cofactor, cofactor, p);
                cnt++;
            }
            exp_par[i] = cnt & 1;
        }

        if (mpz_cmp_ui(cofactor, 1) == 0) {
            // Full relation: P_k^2 ≡ A_val (mod N)
            Relation rel;
            mpz_init_set(rel.x, x_val);
            mpz_init_set(rel.Qx, A_val);
            rel.exp_parity = exp_par;
            relations.push_back(rel);
            full_rels++;
        } else if (mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= lp_bound) {
            // Single large prime partial
            unsigned long lp = mpz_get_ui(cofactor);
            auto &plist = partials[lp];
            if (!plist.empty()) {
                size_t other_idx = plist[0];
                auto &other = partial_store[other_idx];

                Relation rel;
                mpz_init(rel.x);
                mpz_mul(rel.x, x_val, other.x);
                mpz_mod(rel.x, rel.x, N);

                mpz_init(rel.Qx);
                mpz_mul(rel.Qx, A_val, other.Qx);

                rel.exp_parity.resize(fb_size);
                for (size_t i = 0; i < fb_size; i++)
                    rel.exp_parity[i] = (exp_par[i] + other.exp_parity[i]) & 1;

                relations.push_back(rel);
                combined_rels++;
            } else {
                PartialRel pr;
                mpz_init_set(pr.x, x_val);
                mpz_init_set(pr.Qx, A_val);
                pr.exp_parity = exp_par;
                pr.lp = lp;
                plist.push_back(partial_store.size());
                partial_store.push_back(pr);
            }
        }

        mpz_clear(x_val);

        if (k % 10000 == 0) {
            fprintf(stderr, "CFRAC: k=%ld, %zu/%zu rels (full=%d, comb=%d), %.1fs\n",
                    k, relations.size(), relations_needed, full_rels, combined_rels, elapsed_sec());
        }
    }

    fprintf(stderr, "CFRAC: k=%ld, %zu rels (full=%d, combined=%d) in %.1fs\n",
            k, relations.size(), full_rels, combined_rels, elapsed_sec());

    if (relations.size() < fb_size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%zu < %zu)\n", relations.size(), fb_size + 1);
        return 1;
    }

    // Gaussian elimination (same as SIQS)
    size_t nrels = relations.size();
    size_t ncols = fb_size + 1; // +1 for sign
    size_t words = (ncols + 63) / 64;

    std::vector<std::vector<unsigned long>> matrix(nrels);
    for (size_t i = 0; i < nrels; i++) {
        matrix[i].resize(words, 0);
        if (mpz_sgn(relations[i].Qx) < 0)
            matrix[i][0] |= 1UL;
        for (size_t j = 0; j < fb_size; j++) {
            if (relations[i].exp_parity[j] & 1) {
                size_t bit = j + 1;
                matrix[i][bit / 64] |= (1UL << (bit % 64));
            }
        }
    }

    size_t hist_words = (nrels + 63) / 64;
    std::vector<std::vector<unsigned long>> history(nrels);
    for (size_t i = 0; i < nrels; i++) {
        history[i].resize(hist_words, 0);
        history[i][i / 64] |= (1UL << (i % 64));
    }

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

    // Extract factors
    mpz_t lhs, rhs, factor1, factor2, lp_product;
    mpz_init(lhs); mpz_init(rhs); mpz_init(factor1); mpz_init(factor2);
    mpz_init(lp_product);

    bool found = false;
    for (size_t row = 0; row < nrels && !found; row++) {
        bool is_zero = true;
        for (size_t w = 0; w < words; w++)
            if (matrix[row][w] != 0) { is_zero = false; break; }
        if (!is_zero) continue;

        mpz_set_ui(lhs, 1);
        mpz_set_ui(lp_product, 1);
        std::vector<unsigned long> total_exps(fb_size, 0);
        bool sign_neg = false;

        for (size_t i = 0; i < nrels; i++) {
            if (!(history[row][i / 64] & (1UL << (i % 64)))) continue;

            mpz_mul(lhs, lhs, relations[i].x);
            mpz_mod(lhs, lhs, N);

            mpz_set(temp, relations[i].Qx);
            bool neg = (mpz_sgn(temp) < 0);
            if (neg) { mpz_neg(temp, temp); sign_neg = !sign_neg; }

            for (size_t j = 0; j < fb_size; j++) {
                while (mpz_divisible_ui_p(temp, fb[j])) {
                    mpz_divexact_ui(temp, temp, fb[j]);
                    total_exps[j]++;
                }
            }
            if (mpz_cmp_ui(temp, 1) > 0)
                mpz_mul(lp_product, lp_product, temp);
        }

        bool lp_ok = (mpz_cmp_ui(lp_product, 1) == 0) || mpz_perfect_square_p(lp_product);
        int odd_count = 0;
        for (size_t j = 0; j < fb_size; j++)
            if (total_exps[j] & 1) odd_count++;

        if (sign_neg || odd_count > 0 || !lp_ok) continue;

        mpz_set_ui(rhs, 1);
        for (size_t j = 0; j < fb_size; j++) {
            if (total_exps[j] > 0) {
                mpz_ui_pow_ui(temp, fb[j], total_exps[j] / 2);
                mpz_mul(rhs, rhs, temp);
                mpz_mod(rhs, rhs, N);
            }
        }
        if (mpz_cmp_ui(lp_product, 1) > 0) {
            mpz_sqrt(temp, lp_product);
            mpz_mul(rhs, rhs, temp);
            mpz_mod(rhs, rhs, N);
        }

        mpz_sub(temp, lhs, rhs);
        mpz_gcd(factor1, temp, N);
        if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
            mpz_divexact(factor2, N, factor1);
            if (mpz_cmp(factor1, factor2) > 0) mpz_swap(factor1, factor2);
            gmp_printf("%Zd %Zd\n", factor1, factor2);
            fprintf(stderr, "CFRAC: factored in %.3fs\n", elapsed_sec());
            found = true;
        }
        if (!found) {
            mpz_add(temp, lhs, rhs);
            mpz_gcd(factor1, temp, N);
            if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
                mpz_divexact(factor2, N, factor1);
                if (mpz_cmp(factor1, factor2) > 0) mpz_swap(factor1, factor2);
                gmp_printf("%Zd %Zd\n", factor1, factor2);
                fprintf(stderr, "CFRAC: factored in %.3fs\n", elapsed_sec());
                found = true;
            }
        }
    }

    if (!found) {
        fprintf(stderr, "FAIL: no factor found\n");
        return 1;
    }

    // Cleanup
    for (auto &rel : relations) { mpz_clear(rel.x); mpz_clear(rel.Qx); }
    for (auto &pr : partial_store) { mpz_clear(pr.x); mpz_clear(pr.Qx); }
    mpz_clear(N); mpz_clear(sqrtN);
    mpz_clear(P_prev); mpz_clear(P_curr); mpz_clear(Q_prev); mpz_clear(Q_curr);
    mpz_clear(a_k); mpz_clear(g_val); mpz_clear(r_val); mpz_clear(m_val);
    mpz_clear(A_val); mpz_clear(d_val);
    mpz_clear(temp); mpz_clear(cofactor);
    mpz_clear(lhs); mpz_clear(rhs); mpz_clear(factor1); mpz_clear(factor2);
    mpz_clear(lp_product);

    return 0;
}
