// Self-Initializing Quadratic Sieve (SIQS) implementation
// Usage: ./siqs_factor <number>
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

static struct timespec g_start;

static double elapsed_sec() {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

// Modular arithmetic using 128-bit intermediates
static inline unsigned long mulmod(unsigned long a, unsigned long b, unsigned long m) {
    return (__uint128_t)a * b % m;
}

static inline unsigned long powmod(unsigned long base, unsigned long exp, unsigned long mod) {
    unsigned long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mulmod(result, base, mod);
        base = mulmod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

// Tonelli-Shanks modular square root
static unsigned long mod_sqrt(unsigned long n, unsigned long p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;
    if (powmod(n, (p - 1) / 2, p) != 1) return 0;
    if (p % 4 == 3) return powmod(n, (p + 1) / 4, p);

    unsigned long Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    unsigned long z = 2;
    while (powmod(z, (p - 1) / 2, p) != p - 1) z++;

    unsigned long M = S;
    unsigned long c = powmod(z, Q, p);
    unsigned long t = powmod(n, Q, p);
    unsigned long R = powmod(n, (Q + 1) / 2, p);

    while (true) {
        if (t == 1) return R;
        unsigned long i = 0, tmp = t;
        while (tmp != 1) { tmp = mulmod(tmp, tmp, p); i++; }
        unsigned long b = c;
        for (unsigned long j = 0; j < M - i - 1; j++)
            b = mulmod(b, b, p);
        M = i;
        c = mulmod(b, b, p);
        t = mulmod(t, c, p);
        R = mulmod(R, b, p);
    }
}

// Sieve primes
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

// GMP-based modular square root for computing sqrt(N) mod p
static unsigned long mpz_mod_sqrt(const mpz_t N, unsigned long p) {
    unsigned long n_mod_p = mpz_fdiv_ui(N, p);
    return mod_sqrt(n_mod_p, p);
}

struct FactorBase {
    unsigned long p;
    unsigned long sqrt_n; // sqrt(N) mod p
    unsigned char logp;   // log2(p) scaled
};

struct Relation {
    mpz_t x;        // a*x + b
    mpz_t Qx;       // Q(x) = (ax+b)^2 - N
    std::vector<unsigned char> exponents; // full exponents for factor base
};

// SIQS polynomial: g(x) = a*(x + b/a)^2 - (b^2 - N)/a
// where a is a product of factor base primes, b^2 ≡ N (mod a)
// g(x) = a*x^2 + 2*b*x + c, where c = (b^2 - N)/a

struct SIQSPoly {
    mpz_t a, b, c;  // g(x) = a*x^2 + 2*b*x + c
    // Precomputed sieve roots: for each fb prime p, the two values x such that g(x) ≡ 0 (mod p)
    std::vector<long> root1, root2;
};

static void init_poly(SIQSPoly &poly) {
    mpz_init(poly.a);
    mpz_init(poly.b);
    mpz_init(poly.c);
}

static void clear_poly(SIQSPoly &poly) {
    mpz_clear(poly.a);
    mpz_clear(poly.b);
    mpz_clear(poly.c);
}

// Compute SIQS polynomial
// a = product of selected factor base primes (each ~sqrt(2N/M^2)^(1/s) for s primes)
// b satisfies b^2 ≡ N (mod a)
static bool compute_poly(SIQSPoly &poly, const mpz_t N,
                         const std::vector<FactorBase> &fb,
                         const std::vector<size_t> &a_indices,
                         long sieve_radius) {
    size_t fb_size = fb.size();

    // Compute a = product of selected primes
    mpz_set_ui(poly.a, 1);
    for (size_t idx : a_indices) {
        mpz_mul_ui(poly.a, poly.a, fb[idx].p);
    }

    // Compute b using CRT: b^2 ≡ N (mod p_i) for each selected prime
    // Then combine using CRT to get b^2 ≡ N (mod a)

    // For each prime p_i in a, we know sqrt(N) mod p_i = fb[idx].sqrt_n
    // Use CRT to combine

    mpz_t b_crt, mod_so_far, temp, inv, diff;
    mpz_init(b_crt);
    mpz_init(mod_so_far);
    mpz_init(temp);
    mpz_init(inv);
    mpz_init(diff);

    // Start with first prime
    mpz_set_ui(b_crt, fb[a_indices[0]].sqrt_n);
    mpz_set_ui(mod_so_far, fb[a_indices[0]].p);

    for (size_t i = 1; i < a_indices.size(); i++) {
        unsigned long p = fb[a_indices[i]].p;
        unsigned long r = fb[a_indices[i]].sqrt_n;

        // Need b ≡ r (mod p) and b ≡ b_crt (mod mod_so_far)
        // b = b_crt + mod_so_far * t where t = (r - b_crt) * mod_so_far^(-1) (mod p)
        long diff_val = ((long)r - (long)mpz_fdiv_ui(b_crt, p)) % (long)p;
        if (diff_val < 0) diff_val += p;

        unsigned long msf_mod_p = mpz_fdiv_ui(mod_so_far, p);
        unsigned long inv_val = powmod(msf_mod_p, p - 2, p);
        unsigned long t = mulmod((unsigned long)diff_val, inv_val, p);

        mpz_mul_ui(temp, mod_so_far, t);
        mpz_add(b_crt, b_crt, temp);
        mpz_mul_ui(mod_so_far, mod_so_far, p);
        mpz_mod(b_crt, b_crt, mod_so_far);
    }

    mpz_set(poly.b, b_crt);

    // Make b odd if a is odd (doesn't matter much, just convention)
    // Ensure b < a/2 for smaller Q values
    mpz_t half_a;
    mpz_init(half_a);
    mpz_tdiv_q_ui(half_a, poly.a, 2);
    if (mpz_cmp(poly.b, half_a) > 0) {
        mpz_sub(poly.b, poly.a, poly.b);
    }
    mpz_clear(half_a);

    // Verify: b^2 ≡ N (mod a)
    mpz_mul(temp, poly.b, poly.b);
    mpz_sub(temp, temp, N);
    if (!mpz_divisible_p(temp, poly.a)) {
        // Try negating b
        mpz_sub(poly.b, poly.a, poly.b);
        mpz_mul(temp, poly.b, poly.b);
        mpz_sub(temp, temp, N);
        if (!mpz_divisible_p(temp, poly.a)) {
            mpz_clear(b_crt); mpz_clear(mod_so_far); mpz_clear(temp);
            mpz_clear(inv); mpz_clear(diff);
            return false;
        }
    }

    // c = (b^2 - N) / a
    mpz_mul(poly.c, poly.b, poly.b);
    mpz_sub(poly.c, poly.c, N);
    mpz_divexact(poly.c, poly.c, poly.a);

    // Compute sieve roots for each factor base prime
    poly.root1.resize(fb_size);
    poly.root2.resize(fb_size);

    // g(x) = a*x^2 + 2*b*x + c
    // g(x) ≡ 0 (mod p)
    // For primes not in 'a': x ≡ (-b ± sqrt(N)) * a^(-1) (mod p)
    // For primes in 'a': need special handling

    for (size_t i = 0; i < fb_size; i++) {
        unsigned long p = fb[i].p;
        if (p == 2) {
            // Handle p=2: just check parity
            // g(x) mod 2: need to figure out which x values give even g(x)
            unsigned long a2 = mpz_fdiv_ui(poly.a, 2);
            unsigned long b2 = mpz_fdiv_ui(poly.b, 2);
            unsigned long c2 = mpz_fdiv_ui(poly.c, 2);
            // g(x) = a2*x^2 + 2*b2*x + c2 mod 2 = a2*x^2 + c2 mod 2
            if (c2 == 0) {
                poly.root1[i] = 0;
                poly.root2[i] = 1;
            } else if (a2 == 1) {
                poly.root1[i] = 1;
                poly.root2[i] = -1; // only one root
            } else {
                poly.root1[i] = -1;
                poly.root2[i] = -1;
            }
            continue;
        }

        // Check if p divides a
        bool p_divides_a = false;
        for (size_t idx : a_indices) {
            if (fb[idx].p == p) { p_divides_a = true; break; }
        }

        if (p_divides_a) {
            // g(x) = a*x^2 + 2*b*x + c
            // Since p | a: g(x) ≡ 2*b*x + c (mod p)
            // x ≡ -c/(2b) (mod p)
            unsigned long b_mod = mpz_fdiv_ui(poly.b, p);
            unsigned long c_mod = mpz_fdiv_ui(poly.c, p);
            unsigned long two_b = (2 * b_mod) % p;
            if (two_b == 0) {
                poly.root1[i] = -1;
                poly.root2[i] = -1;
                continue;
            }
            unsigned long inv_2b = powmod(two_b, p - 2, p);
            unsigned long x0 = mulmod(p - c_mod, inv_2b, p);
            poly.root1[i] = (long)x0;
            poly.root2[i] = -1; // only one root
            continue;
        }

        // x = (-b ± sqrt_n) * a_inv (mod p)
        unsigned long a_mod = mpz_fdiv_ui(poly.a, p);
        unsigned long b_mod = mpz_fdiv_ui(poly.b, p);
        unsigned long a_inv = powmod(a_mod, p - 2, p);
        unsigned long sq = fb[i].sqrt_n;

        // root1 = (-b + sq) * a_inv mod p
        long val1 = ((long)sq - (long)b_mod) % (long)p;
        if (val1 < 0) val1 += p;
        poly.root1[i] = (long)mulmod((unsigned long)val1, a_inv, p);

        // root2 = (-b - sq) * a_inv mod p
        long val2 = (-(long)sq - (long)b_mod) % (long)p;
        if (val2 < 0) val2 += p;
        poly.root2[i] = (long)mulmod((unsigned long)val2, a_inv, p);
    }

    mpz_clear(b_crt); mpz_clear(mod_so_far); mpz_clear(temp);
    mpz_clear(inv); mpz_clear(diff);
    return true;
}

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

    // SIQS parameters tuned by digit count
    unsigned long B;         // factor base bound
    long sieve_radius;       // sieve half-interval M
    int num_a_primes;        // number of primes in 'a' polynomial coefficient

    double ln_n = bits * log(2.0);
    double ln_ln_n = log(ln_n);
    double L_half = exp(sqrt(ln_n * ln_ln_n));

    // Optimal: B ~ L^(1/2 * sqrt(2)), M ~ L^(1/2 * sqrt(2))
    // In practice, tune empirically
    // Parameters: B (smoothness bound), M (sieve half-interval)
    // num_a_primes chosen dynamically below based on target_a
    if (digits <= 30) {
        B = 1000; sieve_radius = 32768;
    } else if (digits <= 35) {
        B = 2000; sieve_radius = 32768;
    } else if (digits <= 40) {
        B = 5000; sieve_radius = 65536;
    } else if (digits <= 45) {
        B = 12000; sieve_radius = 65536;
    } else if (digits <= 50) {
        B = 25000; sieve_radius = 65536;
    } else if (digits <= 55) {
        B = 50000; sieve_radius = 65536;
    } else if (digits <= 60) {
        B = 100000; sieve_radius = 65536*2;
    } else if (digits <= 65) {
        B = 200000; sieve_radius = 65536*2;
    } else if (digits <= 70) {
        B = 400000; sieve_radius = 65536*2;
    } else if (digits <= 75) {
        B = 700000; sieve_radius = 65536*4;
    } else if (digits <= 80) {
        B = 1200000; sieve_radius = 65536*4;
    } else if (digits <= 85) {
        B = 2000000; sieve_radius = 65536*4;
    } else if (digits <= 90) {
        B = 3500000; sieve_radius = 65536*4;
    } else {
        B = 6000000; sieve_radius = 65536*8;
    }

    // Build factor base
    std::vector<unsigned long> all_primes = sieve_primes(B);
    std::vector<FactorBase> fb;

    // Always include 2
    fb.push_back({2, 1, 1});

    for (size_t i = 1; i < all_primes.size(); i++) {
        unsigned long p = all_primes[i];
        unsigned long sq = mpz_mod_sqrt(N, p);
        if (sq != 0) {
            unsigned char lp = (unsigned char)(log2((double)p) + 0.5);
            if (lp == 0) lp = 1;
            fb.push_back({p, sq, lp});
        }
    }

    size_t fb_size = fb.size();
    size_t relations_needed = fb_size + 50;

    fprintf(stderr, "SIQS: %zu digits, FB=%zu primes, M=%ld, a_primes=%d\n",
            digits, fb_size, sieve_radius, num_a_primes);

    // Determine target 'a' value: a ≈ sqrt(2*N) / M
    mpz_t target_a, two_n;
    mpz_init(target_a);
    mpz_init(two_n);
    mpz_mul_2exp(two_n, N, 1);
    mpz_sqrt(target_a, two_n);
    mpz_tdiv_q_ui(target_a, target_a, sieve_radius);

    double log_target_a = mpz_sizeinbase(target_a, 2) * log(2.0);

    // Dynamically choose num_a_primes: pick primes from upper 60% of factor base
    // Such that their product ≈ target_a
    // Use primes from index fb_size*0.4 to fb_size-1
    size_t range_start = std::max((size_t)2, fb_size * 2 / 5);
    size_t range_end = fb_size - 1;

    // Average log of primes in range
    double avg_log = 0;
    for (size_t i = range_start; i <= range_end; i++)
        avg_log += log((double)fb[i].p);
    avg_log /= (range_end - range_start + 1);

    num_a_primes = (int)(log_target_a / avg_log + 0.5);
    if (num_a_primes < 2) num_a_primes = 2;
    if (num_a_primes > (int)(range_end - range_start + 1) / 2)
        num_a_primes = (int)(range_end - range_start + 1) / 2;

    // Verify: if product of num_a_primes median primes is too far from target, adjust
    mpz_t target_prime;
    mpz_init(target_prime);

    fprintf(stderr, "SIQS: target_a ~ %zu bits, avg_log=%.1f, num_a_primes=%d, range=[%zu..%zu]\n",
            mpz_sizeinbase(target_a, 2), avg_log, num_a_primes, range_start, range_end);

    // Sieve array
    std::vector<unsigned char> sieve(2 * sieve_radius);
    long M = sieve_radius;

    // Collect relations
    std::vector<Relation> relations;

    // Seeded RNG
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    // Large prime bound for single large prime variation
    unsigned long large_prime_bound = fb.back().p * 100;

    // Track partial relations (single large prime)
    // Map from large prime -> index of first partial relation
    #include <unordered_map>

    struct PartialRelation {
        mpz_t x, Qx;
        std::vector<unsigned char> exponents;
        unsigned long large_prime;
    };
    std::unordered_map<unsigned long, std::vector<PartialRelation>> partials;

    SIQSPoly poly;
    init_poly(poly);

    int poly_count = 0;
    int max_polys = 100000;

    // Threshold: Q(x) for |x| ~ M has ~(bits/2 + log2(M)) bits.
    // Accept if sieve explains all but the large prime allowance.
    double expected_log2_Q = (double)bits / 2.0 + log2((double)sieve_radius) + 0.5;
    double lp_log = log2((double)large_prime_bound);
    unsigned char threshold = (unsigned char)(expected_log2_Q - lp_log - 4);
    if (threshold < 10) threshold = 10;

    fprintf(stderr, "SIQS: expected_log2_Q=%.1f, threshold=%d\n", expected_log2_Q, (int)threshold);

    mpz_t ax_plus_b, qval, cofactor, g1, g2, temp;
    mpz_init(ax_plus_b);
    mpz_init(qval);
    mpz_init(cofactor);
    mpz_init(g1);
    mpz_init(g2);
    mpz_init(temp);

    while (relations.size() < relations_needed && poly_count < max_polys) {
        if (elapsed_sec() > 280.0) {
            fprintf(stderr, "SIQS: timeout after %.1fs, %zu relations\n",
                    elapsed_sec(), relations.size());
            break;
        }

        // Pick random combination of primes for 'a'
        std::vector<size_t> a_indices;
        // Randomly select num_a_primes distinct indices from [range_start, range_end]
        std::vector<size_t> candidates;
        for (size_t i = range_start; i <= range_end; i++)
            candidates.push_back(i);

        // Fisher-Yates shuffle (seeded)
        for (size_t i = candidates.size() - 1; i > 0; i--) {
            size_t j = gmp_urandomm_ui(rng, i + 1);
            std::swap(candidates[i], candidates[j]);
        }
        for (int i = 0; i < num_a_primes && i < (int)candidates.size(); i++)
            a_indices.push_back(candidates[i]);

        std::sort(a_indices.begin(), a_indices.end());

        if (!compute_poly(poly, N, fb, a_indices, sieve_radius)) {
            poly_count++;
            continue;
        }

        poly_count++;

        // Sieve
        memset(sieve.data(), 0, sieve.size());

        for (size_t i = 0; i < fb_size; i++) {
            unsigned long p = fb[i].p;
            unsigned char lp = fb[i].logp;

            if (p == 2) {
                // Sieve every other position
                long r1 = poly.root1[i];
                if (r1 >= 0) {
                    for (long j = r1; j < 2 * M; j += 2)
                        sieve[j] += lp;
                }
                continue;
            }

            long r1 = poly.root1[i];
            long r2 = poly.root2[i];

            // Convert roots to sieve positions: root is in [0, p)
            // Sieve array index = x + M where x in [-M, M)
            // So position = root + M, adjusted modulo p

            if (r1 >= 0) {
                long pos = ((r1 - (-M)) % (long)p);
                if (pos < 0) pos += p;
                // Actually: x ≡ r1 (mod p), sieve index = x + M
                // x = r1, r1+p, r1+2p, ... and x = r1-p, r1-2p, ...
                // sieve index = x + M
                // Start from x = r1 - M_rounded (smallest x >= -M with x ≡ r1 mod p)
                pos = r1 + M;
                pos %= (long)p;
                if (pos < 0) pos += p;
                for (long j = pos; j < 2 * M; j += p)
                    sieve[j] += lp;
            }
            if (r2 >= 0 && r2 != r1) {
                long pos = r2 + M;
                pos %= (long)p;
                if (pos < 0) pos += p;
                for (long j = pos; j < 2 * M; j += p)
                    sieve[j] += lp;
            }

            // Sieve with prime powers for small primes
            if (p < 256) {
                unsigned long pk = p * p;
                while (pk < (unsigned long)sieve_radius) {
                    // Lift roots mod pk using Hensel's lemma
                    // For simplicity, just sieve with p^2 at existing roots
                    if (r1 >= 0) {
                        long pos = r1 + M;
                        pos %= (long)pk;
                        if (pos < 0) pos += pk;
                        for (long j = pos; j < 2 * M; j += pk)
                            sieve[j] += lp;
                    }
                    if (r2 >= 0 && r2 != r1) {
                        long pos = r2 + M;
                        pos %= (long)pk;
                        if (pos < 0) pos += pk;
                        for (long j = pos; j < 2 * M; j += pk)
                            sieve[j] += lp;
                    }
                    pk *= p;
                }
            }
        }

        // Scan sieve for smooth candidates
        for (long idx = 0; idx < 2 * M; idx++) {
            if (sieve[idx] < threshold) continue;

            long x = idx - M;

            // Compute Q(x) = a*x^2 + 2*b*x + c
            mpz_set_si(qval, x);
            mpz_mul_si(qval, qval, x);
            mpz_mul(qval, qval, poly.a);  // a*x^2

            mpz_set_si(temp, x);
            mpz_mul(temp, temp, poly.b);
            mpz_mul_2exp(temp, temp, 1);   // 2*b*x
            mpz_add(qval, qval, temp);
            mpz_add(qval, qval, poly.c);   // + c

            if (mpz_sgn(qval) == 0) continue;

            // The actual x value for the relation: ax + b
            mpz_mul_si(ax_plus_b, poly.a, x);
            mpz_add(ax_plus_b, ax_plus_b, poly.b);

            // Trial divide
            mpz_set(cofactor, qval);
            if (mpz_sgn(cofactor) < 0) mpz_neg(cofactor, cofactor);

            std::vector<unsigned char> exps(fb_size, 0);
            bool sign_negative = (mpz_sgn(qval) < 0);

            for (size_t i = 0; i < fb_size; i++) {
                unsigned long p = fb[i].p;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    exps[i]++;
                }
            }

            if (mpz_cmp_ui(cofactor, 1) == 0) {
                // Fully smooth - add relation
                // The relation is (ax+b)^2 ≡ a*g(x) (mod N)
                // So we need exponents of a*g(x), not just g(x)
                // Add exponents for the primes in 'a'
                for (size_t ai : a_indices)
                    exps[ai]++;

                Relation rel;
                mpz_init_set(rel.x, ax_plus_b);
                // Store a*g(x) as the Q value
                mpz_init(rel.Qx);
                mpz_mul(rel.Qx, qval, poly.a);
                rel.exponents = exps;

                relations.push_back(rel);

                if (relations.size() % 200 == 0) {
                    fprintf(stderr, "SIQS: %zu/%zu relations, %d polys, %.1fs\n",
                            relations.size(), relations_needed, poly_count, elapsed_sec());
                }
            } else if (mpz_cmp_ui(cofactor, large_prime_bound) <= 0 &&
                       mpz_fits_ulong_p(cofactor)) {
                // Single large prime variation
                unsigned long lp = mpz_get_ui(cofactor);
                // Add 'a' exponents (same as full relations)
                for (size_t ai : a_indices)
                    exps[ai]++;

                auto &plist = partials[lp];
                if (!plist.empty()) {
                    // Combine with existing partial to make full relation
                    auto &prev = plist[0];

                    // Combined: (x1 * x2)^2 ≡ Q1 * Q2 (mod N)
                    // The large prime appears once in each, squared cancels in mod 2
                    Relation rel;
                    mpz_init(rel.x);
                    mpz_mul(rel.x, ax_plus_b, prev.x);
                    mpz_mod(rel.x, rel.x, N);

                    mpz_init(rel.Qx);
                    mpz_mul(rel.Qx, qval, poly.a); // a*g(x)
                    mpz_mul(rel.Qx, rel.Qx, prev.Qx);

                    rel.exponents.resize(fb_size);
                    for (size_t i = 0; i < fb_size; i++)
                        rel.exponents[i] = (exps[i] + prev.exponents[i]) & 1;

                    relations.push_back(rel);
                } else {
                    PartialRelation pr;
                    mpz_init_set(pr.x, ax_plus_b);
                    mpz_init(pr.Qx);
                    mpz_mul(pr.Qx, qval, poly.a); // a*g(x)
                    pr.exponents = exps;
                    pr.large_prime = lp;
                    plist.push_back(pr);
                }
            }
        }
    }

    fprintf(stderr, "SIQS: collected %zu full relations from %d polynomials (%.1fs)\n",
            relations.size(), poly_count, elapsed_sec());

    if (relations.size() < fb_size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%zu < %zu)\n", relations.size(), fb_size + 1);
        return 1;
    }

    // Gaussian elimination mod 2
    size_t nrels = relations.size();
    size_t ncols = fb_size + 1; // +1 for sign
    size_t words = (ncols + 63) / 64;

    std::vector<std::vector<unsigned long>> matrix(nrels);
    for (size_t i = 0; i < nrels; i++) {
        matrix[i].resize(words, 0);
        if (mpz_sgn(relations[i].Qx) < 0)
            matrix[i][0] |= 1UL;
        for (size_t j = 0; j < fb_size; j++) {
            if (relations[i].exponents[j] & 1) {
                size_t bit = j + 1;
                matrix[i][bit / 64] |= (1UL << (bit % 64));
            }
        }
    }

    // History matrix for tracking combinations
    size_t hist_words = (nrels + 63) / 64;
    std::vector<std::vector<unsigned long>> history(nrels);
    for (size_t i = 0; i < nrels; i++) {
        history[i].resize(hist_words, 0);
        history[i][i / 64] |= (1UL << (i % 64));
    }

    // Row-reduce
    std::vector<size_t> pivot_col(nrels, (size_t)-1);
    size_t col = 0;
    for (size_t row = 0; row < nrels && col < ncols; ) {
        // Find pivot in this column
        size_t piv = (size_t)-1;
        for (size_t r = row; r < nrels; r++) {
            if (matrix[r][col / 64] & (1UL << (col % 64))) {
                piv = r;
                break;
            }
        }
        if (piv == (size_t)-1) {
            col++;
            continue;
        }

        // Swap rows
        if (piv != row) {
            std::swap(matrix[piv], matrix[row]);
            std::swap(history[piv], history[row]);
        }
        pivot_col[row] = col;

        // Eliminate
        for (size_t r = 0; r < nrels; r++) {
            if (r == row) continue;
            if (matrix[r][col / 64] & (1UL << (col % 64))) {
                for (size_t w = 0; w < words; w++)
                    matrix[r][w] ^= matrix[row][w];
                for (size_t w = 0; w < hist_words; w++)
                    history[r][w] ^= history[row][w];
            }
        }
        row++;
        col++;
    }

    // Find null vectors and try to extract factors
    mpz_t lhs, rhs, factor1, factor2;
    mpz_init(lhs);
    mpz_init(rhs);
    mpz_init(factor1);
    mpz_init(factor2);

    bool found = false;
    int null_count = 0, tried = 0;
    for (size_t row = 0; row < nrels && !found; row++) {
        // Check if this row is zero
        bool is_zero = true;
        for (size_t w = 0; w < words; w++) {
            if (matrix[row][w] != 0) { is_zero = false; break; }
        }
        if (!is_zero) continue;
        null_count++;

        // Combine relations in this null vector
        mpz_set_ui(lhs, 1);

        std::vector<unsigned long> total_exps(fb_size, 0);
        bool sign_neg = false;

        // Accumulate large prime cofactors (from combined partial relations)
        mpz_t lp_product;
        mpz_init(lp_product);
        mpz_set_ui(lp_product, 1);

        for (size_t i = 0; i < nrels; i++) {
            if (!(history[row][i / 64] & (1UL << (i % 64)))) continue;

            mpz_mul(lhs, lhs, relations[i].x);
            mpz_mod(lhs, lhs, N);

            // Accumulate exponents from Q values
            mpz_set(temp, relations[i].Qx);
            bool neg = (mpz_sgn(temp) < 0);
            if (neg) { mpz_neg(temp, temp); sign_neg = !sign_neg; }

            for (size_t j = 0; j < fb_size; j++) {
                unsigned long p = fb[j].p;
                while (mpz_divisible_ui_p(temp, p)) {
                    mpz_divexact_ui(temp, temp, p);
                    total_exps[j]++;
                }
            }
            // Remaining cofactor from large prime relations
            if (mpz_cmp_ui(temp, 1) > 0) {
                mpz_mul(lp_product, lp_product, temp);
            }
        }

        // The lp_product should be a perfect square
        bool lp_ok = true;
        if (mpz_cmp_ui(lp_product, 1) > 0) {
            if (!mpz_perfect_square_p(lp_product)) {
                lp_ok = false;
            }
        }

        // Check all exponents are even
        int odd_count = 0;
        for (size_t j = 0; j < fb_size; j++)
            if (total_exps[j] & 1) odd_count++;

        if (sign_neg || odd_count > 0 || !lp_ok) {
            tried++;
            mpz_clear(lp_product);
            continue;
        }

        // Compute rhs = product of p^(exp/2) * sqrt(lp_product) mod N
        mpz_set_ui(rhs, 1);
        for (size_t j = 0; j < fb_size; j++) {
            if (total_exps[j] > 0) {
                mpz_ui_pow_ui(temp, fb[j].p, total_exps[j] / 2);
                mpz_mul(rhs, rhs, temp);
                mpz_mod(rhs, rhs, N);
            }
        }
        // Include large prime square root
        if (mpz_cmp_ui(lp_product, 1) > 0) {
            mpz_sqrt(temp, lp_product);
            mpz_mul(rhs, rhs, temp);
            mpz_mod(rhs, rhs, N);
        }
        mpz_clear(lp_product);

        // Try gcd(lhs - rhs, N) and gcd(lhs + rhs, N)
        mpz_sub(temp, lhs, rhs);
        mpz_gcd(factor1, temp, N);

        if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
            mpz_divexact(factor2, N, factor1);
            if (mpz_cmp(factor1, factor2) > 0) mpz_swap(factor1, factor2);
            gmp_printf("%Zd %Zd\n", factor1, factor2);
            fprintf(stderr, "SIQS: factored in %.3fs\n", elapsed_sec());
            found = true;
            break;
        }

        mpz_add(temp, lhs, rhs);
        mpz_gcd(factor1, temp, N);

        if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
            mpz_divexact(factor2, N, factor1);
            if (mpz_cmp(factor1, factor2) > 0) mpz_swap(factor1, factor2);
            gmp_printf("%Zd %Zd\n", factor1, factor2);
            fprintf(stderr, "SIQS: factored in %.3fs\n", elapsed_sec());
            found = true;
            break;
        }
    }

    if (!found) {
        fprintf(stderr, "FAIL: %d null space vectors, %d had odd exps. No factor found.\n",
                null_count, tried);
        return 1;
    }

    // Cleanup
    for (auto &rel : relations) {
        mpz_clear(rel.x);
        mpz_clear(rel.Qx);
    }
    for (auto &[k, plist] : partials) {
        for (auto &pr : plist) {
            mpz_clear(pr.x);
            mpz_clear(pr.Qx);
        }
    }

    clear_poly(poly);
    gmp_randclear(rng);
    mpz_clear(N);
    mpz_clear(target_a);
    mpz_clear(two_n);
    mpz_clear(target_prime);
    mpz_clear(ax_plus_b);
    mpz_clear(qval);
    mpz_clear(cofactor);
    mpz_clear(g1);
    mpz_clear(g2);
    mpz_clear(temp);
    mpz_clear(lhs);
    mpz_clear(rhs);
    mpz_clear(factor1);
    mpz_clear(factor2);

    return 0;
}
