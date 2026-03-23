// Optimized SIQS with block sieving and double large prime variation
// Usage: ./siqs_opt <number>
// Outputs: factor1 factor2

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

struct FactorBase {
    unsigned long p;
    unsigned long sqrt_n;  // sqrt(N) mod p
    unsigned char logp;    // scaled log2(p)
    // For small primes, store prime powers
};

struct Relation {
    mpz_t x;        // the x value: (ax+b)
    mpz_t Qx;       // the smoothed value: a * g(x)
    std::vector<unsigned char> exp_parity; // parity of exponents
};

struct PartialRelation {
    mpz_t x, Qx;
    std::vector<unsigned char> exp_parity;
    unsigned long lp;  // large prime
};

// Double large prime relation
struct DLP_Relation {
    mpz_t x, Qx;
    std::vector<unsigned char> exp_parity;
    unsigned long lp1, lp2;  // two large primes (lp1 < lp2)
};

// SIQS polynomial: g(x) = a*x^2 + 2*b*x + c
struct SIQSPoly {
    mpz_t a, b, c;
    std::vector<unsigned long> root1, root2; // sieve roots mod p
    std::vector<size_t> a_indices;  // which FB primes are in 'a'
};

static bool compute_poly(SIQSPoly &poly, const mpz_t N,
                         const std::vector<FactorBase> &fb,
                         size_t fb_size) {
    // Compute b via CRT
    mpz_t b_crt, mod_so_far, temp;
    mpz_init(b_crt);
    mpz_init(mod_so_far);
    mpz_init(temp);

    mpz_set_ui(b_crt, fb[poly.a_indices[0]].sqrt_n);
    mpz_set_ui(mod_so_far, fb[poly.a_indices[0]].p);

    for (size_t i = 1; i < poly.a_indices.size(); i++) {
        unsigned long p = fb[poly.a_indices[i]].p;
        unsigned long r = fb[poly.a_indices[i]].sqrt_n;
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

    // Make b < a/2
    mpz_t half_a;
    mpz_init(half_a);
    mpz_tdiv_q_ui(half_a, poly.a, 2);
    if (mpz_cmp(poly.b, half_a) > 0)
        mpz_sub(poly.b, poly.a, poly.b);
    mpz_clear(half_a);

    // Verify and fix b
    mpz_mul(temp, poly.b, poly.b);
    mpz_sub(temp, temp, N);
    if (!mpz_divisible_p(temp, poly.a)) {
        mpz_sub(poly.b, poly.a, poly.b);
        mpz_mul(temp, poly.b, poly.b);
        mpz_sub(temp, temp, N);
        if (!mpz_divisible_p(temp, poly.a)) {
            mpz_clear(b_crt); mpz_clear(mod_so_far); mpz_clear(temp);
            return false;
        }
    }

    // c = (b^2 - N) / a
    mpz_mul(poly.c, poly.b, poly.b);
    mpz_sub(poly.c, poly.c, N);
    mpz_divexact(poly.c, poly.c, poly.a);

    // Compute sieve roots
    poly.root1.resize(fb_size);
    poly.root2.resize(fb_size);

    for (size_t i = 0; i < fb_size; i++) {
        unsigned long p = fb[i].p;

        if (p == 2) {
            unsigned long c2 = mpz_fdiv_ui(poly.c, 2);
            unsigned long a2 = mpz_fdiv_ui(poly.a, 2);
            if (c2 == 0) { poly.root1[i] = 0; poly.root2[i] = 1; }
            else if (a2 == 1) { poly.root1[i] = 1; poly.root2[i] = (unsigned long)-1; }
            else { poly.root1[i] = (unsigned long)-1; poly.root2[i] = (unsigned long)-1; }
            continue;
        }

        // Check if p divides a
        bool p_divides_a = false;
        for (size_t idx : poly.a_indices)
            if (fb[idx].p == p) { p_divides_a = true; break; }

        if (p_divides_a) {
            unsigned long b_mod = mpz_fdiv_ui(poly.b, p);
            unsigned long c_mod = mpz_fdiv_ui(poly.c, p);
            unsigned long two_b = (2 * b_mod) % p;
            if (two_b == 0) {
                poly.root1[i] = (unsigned long)-1;
                poly.root2[i] = (unsigned long)-1;
            } else {
                unsigned long inv_2b = powmod(two_b, p-2, p);
                poly.root1[i] = mulmod(p - c_mod, inv_2b, p);
                poly.root2[i] = (unsigned long)-1;
            }
            continue;
        }

        unsigned long a_mod = mpz_fdiv_ui(poly.a, p);
        unsigned long b_mod = mpz_fdiv_ui(poly.b, p);
        unsigned long a_inv = powmod(a_mod, p-2, p);
        unsigned long sq = fb[i].sqrt_n;

        long v1 = ((long)sq - (long)b_mod) % (long)p;
        if (v1 < 0) v1 += p;
        poly.root1[i] = mulmod((unsigned long)v1, a_inv, p);

        long v2 = (-(long)sq - (long)b_mod) % (long)p;
        if (v2 < 0) v2 += p;
        poly.root2[i] = mulmod((unsigned long)v2, a_inv, p);
    }

    mpz_clear(b_crt); mpz_clear(mod_so_far); mpz_clear(temp);
    return true;
}

// Miller-Rabin primality test
static bool is_probable_prime(unsigned long n) {
    if (n < 2) return false;
    if (n < 4) return true;
    if (n % 2 == 0) return false;

    unsigned long d = n - 1;
    int r = 0;
    while (d % 2 == 0) { d /= 2; r++; }

    // Test with bases 2, 3, 5, 7, 11
    unsigned long witnesses[] = {2, 3, 5, 7, 11};
    for (unsigned long a : witnesses) {
        if (a >= n) continue;
        unsigned long x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int i = 0; i < r - 1; i++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { composite = false; break; }
        }
        if (composite) return false;
    }
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

    // Parameters
    unsigned long B;
    long sieve_radius;
    static const int BLOCK_SIZE = 32768; // 32KB for cache-friendly sieving

    if (digits <= 30) { B = 1200; sieve_radius = 32768; }
    else if (digits <= 35) { B = 2500; sieve_radius = 32768; }
    else if (digits <= 40) { B = 6000; sieve_radius = 65536; }
    else if (digits <= 45) { B = 15000; sieve_radius = 65536; }
    else if (digits <= 50) { B = 30000; sieve_radius = 65536; }
    else if (digits <= 55) { B = 60000; sieve_radius = 65536; }
    else if (digits <= 60) { B = 120000; sieve_radius = 131072; }
    else if (digits <= 65) { B = 250000; sieve_radius = 131072; }
    else if (digits <= 70) { B = 500000; sieve_radius = 131072; }
    else if (digits <= 75) { B = 900000; sieve_radius = 262144; }
    else if (digits <= 80) { B = 1500000; sieve_radius = 262144; }
    else if (digits <= 85) { B = 2500000; sieve_radius = 262144; }
    else if (digits <= 90) { B = 4000000; sieve_radius = 262144; }
    else { B = 7000000; sieve_radius = 524288; }

    // Build factor base
    std::vector<unsigned long> all_primes = sieve_primes(B);
    std::vector<FactorBase> fb;
    fb.push_back({2, 1, 1});

    for (size_t i = 1; i < all_primes.size(); i++) {
        unsigned long p = all_primes[i];
        unsigned long n_mod_p = mpz_fdiv_ui(N, p);
        unsigned long sq = mod_sqrt(n_mod_p, p);
        if (sq != 0) {
            unsigned char lp = (unsigned char)(log2((double)p) + 0.5);
            if (lp == 0) lp = 1;
            fb.push_back({p, sq, lp});
        }
    }
    size_t fb_size = fb.size();
    size_t relations_needed = fb_size + 64;

    // Large prime bounds
    unsigned long lp1_bound = fb.back().p * 50;   // single large prime
    unsigned long lp2_bound = fb.back().p * 20;    // for double large prime

    // Dynamic num_a_primes
    mpz_t target_a, two_n;
    mpz_init(target_a);
    mpz_init(two_n);
    mpz_mul_2exp(two_n, N, 1);
    mpz_sqrt(target_a, two_n);
    mpz_tdiv_q_ui(target_a, target_a, sieve_radius);

    double log_target_a = mpz_sizeinbase(target_a, 2) * log(2.0);
    size_t range_start = std::max((size_t)2, fb_size * 2 / 5);
    size_t range_end = fb_size - 1;

    double avg_log = 0;
    for (size_t i = range_start; i <= range_end; i++)
        avg_log += log((double)fb[i].p);
    avg_log /= (range_end - range_start + 1);

    int num_a_primes = (int)(log_target_a / avg_log + 0.5);
    if (num_a_primes < 2) num_a_primes = 2;
    if (num_a_primes > (int)(range_end - range_start + 1) / 2)
        num_a_primes = (int)(range_end - range_start + 1) / 2;

    long M = sieve_radius;

    // Threshold
    double expected_log2_Q = (double)bits / 2.0 + log2((double)sieve_radius) + 0.5;
    double lp_log = log2((double)lp1_bound);
    unsigned char threshold = (unsigned char)(expected_log2_Q - lp_log - 3);
    if (threshold < 10) threshold = 10;

    fprintf(stderr, "SIQS-opt: %zu digits, FB=%zu, M=%ld, a_primes=%d, threshold=%d\n",
            digits, fb_size, M, num_a_primes, (int)threshold);

    // Storage for relations
    std::vector<Relation> relations;

    // Single large prime partials
    std::unordered_map<unsigned long, std::vector<PartialRelation>> slp_partials;

    // Double large prime graph
    // Map from large prime to list of DLP relations containing it
    std::unordered_map<unsigned long, std::vector<size_t>> dlp_graph;
    std::vector<DLP_Relation> dlp_rels;

    // Sieve array
    std::vector<unsigned char> sieve(2 * M);

    // RNG
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    SIQSPoly poly;
    mpz_init(poly.a);
    mpz_init(poly.b);
    mpz_init(poly.c);

    mpz_t ax_plus_b, qval, cofactor, temp;
    mpz_init(ax_plus_b);
    mpz_init(qval);
    mpz_init(cofactor);
    mpz_init(temp);

    int poly_count = 0;
    int full_rels = 0, slp_combined = 0, dlp_combined = 0;

    while (relations.size() < relations_needed) {
        if (elapsed_sec() > 280.0) {
            fprintf(stderr, "SIQS-opt: timeout, %zu relations\n", relations.size());
            break;
        }

        // Generate new polynomial
        std::vector<size_t> candidates;
        for (size_t i = range_start; i <= range_end; i++)
            candidates.push_back(i);

        for (size_t i = candidates.size() - 1; i > 0; i--) {
            size_t j = gmp_urandomm_ui(rng, i + 1);
            std::swap(candidates[i], candidates[j]);
        }

        poly.a_indices.clear();
        for (int i = 0; i < num_a_primes && i < (int)candidates.size(); i++)
            poly.a_indices.push_back(candidates[i]);
        std::sort(poly.a_indices.begin(), poly.a_indices.end());

        mpz_set_ui(poly.a, 1);
        for (size_t idx : poly.a_indices)
            mpz_mul_ui(poly.a, poly.a, fb[idx].p);

        if (!compute_poly(poly, N, fb, fb_size)) {
            poly_count++;
            continue;
        }
        poly_count++;

        // Block sieve
        for (long block_start = 0; block_start < 2 * M; block_start += BLOCK_SIZE) {
            long block_end = std::min(block_start + (long)BLOCK_SIZE, 2 * M);
            long block_len = block_end - block_start;

            memset(sieve.data() + block_start, 0, block_len);

            // Sieve this block
            for (size_t i = 0; i < fb_size; i++) {
                unsigned long p = fb[i].p;
                unsigned char lp = fb[i].logp;

                if (p == 2) {
                    unsigned long r1 = poly.root1[i];
                    if (r1 != (unsigned long)-1) {
                        long start = r1;
                        if (start < block_start) {
                            long skip = (block_start - start + 1) / 2 * 2;
                            start += skip;
                        }
                        for (long j = start; j < block_end; j += 2)
                            sieve[j] += lp;
                    }
                    continue;
                }

                // Root 1
                unsigned long r1 = poly.root1[i];
                if (r1 != (unsigned long)-1) {
                    // First position in block: ceil((block_start - (r1 + M)) / p) * p + (r1 + M) but simpler:
                    long pos = (long)((r1 + M) % p);
                    if (pos < 0) pos += p;
                    // Find first pos >= block_start
                    if (pos < block_start) {
                        long skip = ((block_start - pos + p - 1) / p) * p;
                        pos += skip;
                    }
                    for (long j = pos; j < block_end; j += p)
                        sieve[j] += lp;
                }

                // Root 2
                unsigned long r2 = poly.root2[i];
                if (r2 != (unsigned long)-1 && r2 != r1) {
                    long pos = (long)((r2 + M) % p);
                    if (pos < 0) pos += p;
                    if (pos < block_start) {
                        long skip = ((block_start - pos + p - 1) / p) * p;
                        pos += skip;
                    }
                    for (long j = pos; j < block_end; j += p)
                        sieve[j] += lp;
                }

                // Sieve with prime powers for small primes
                if (p < 128) {
                    unsigned long pk = p * p;
                    while (pk < (unsigned long)M) {
                        if (r1 != (unsigned long)-1) {
                            long pos = (long)((r1 + M) % pk);
                            if (pos < 0) pos += pk;
                            if (pos < block_start) {
                                long skip = ((block_start - pos + pk - 1) / pk) * pk;
                                pos += skip;
                            }
                            for (long j = pos; j < block_end; j += pk)
                                sieve[j] += lp;
                        }
                        if (r2 != (unsigned long)-1 && r2 != r1) {
                            long pos = (long)((r2 + M) % pk);
                            if (pos < 0) pos += pk;
                            if (pos < block_start) {
                                long skip = ((block_start - pos + pk - 1) / pk) * pk;
                                pos += skip;
                            }
                            for (long j = pos; j < block_end; j += pk)
                                sieve[j] += lp;
                        }
                        pk *= p;
                    }
                }
            }

            // Scan block for smooth candidates
            for (long idx = block_start; idx < block_end; idx++) {
                if (sieve[idx] < threshold) continue;

                long x = idx - M;

                // Compute Q(x) = a*x^2 + 2*b*x + c
                mpz_set_si(qval, x);
                mpz_mul_si(qval, qval, x);
                mpz_mul(qval, qval, poly.a);
                mpz_set_si(temp, x);
                mpz_mul(temp, temp, poly.b);
                mpz_mul_2exp(temp, temp, 1);
                mpz_add(qval, qval, temp);
                mpz_add(qval, qval, poly.c);

                if (mpz_sgn(qval) == 0) continue;

                mpz_mul_si(ax_plus_b, poly.a, x);
                mpz_add(ax_plus_b, ax_plus_b, poly.b);

                // Trial divide
                mpz_set(cofactor, qval);
                if (mpz_sgn(cofactor) < 0) mpz_neg(cofactor, cofactor);

                std::vector<unsigned char> exps(fb_size, 0);

                for (size_t i = 0; i < fb_size; i++) {
                    unsigned long p = fb[i].p;
                    while (mpz_divisible_ui_p(cofactor, p)) {
                        mpz_divexact_ui(cofactor, cofactor, p);
                        exps[i]++;
                    }
                }

                // Add 'a' exponents
                for (size_t ai : poly.a_indices)
                    exps[ai]++;

                // Convert to parity
                std::vector<unsigned char> parity(fb_size);
                for (size_t i = 0; i < fb_size; i++)
                    parity[i] = exps[i] & 1;

                if (mpz_cmp_ui(cofactor, 1) == 0) {
                    // Full relation
                    Relation rel;
                    mpz_init_set(rel.x, ax_plus_b);
                    mpz_init(rel.Qx);
                    mpz_mul(rel.Qx, qval, poly.a);
                    rel.exp_parity = parity;
                    relations.push_back(rel);
                    full_rels++;
                } else if (mpz_fits_ulong_p(cofactor)) {
                    unsigned long cf = mpz_get_ui(cofactor);

                    if (cf <= lp1_bound && is_probable_prime(cf)) {
                        // Single large prime
                        auto &plist = slp_partials[cf];
                        if (!plist.empty()) {
                            auto &prev = plist[0];
                            Relation rel;
                            mpz_init(rel.x);
                            mpz_mul(rel.x, ax_plus_b, prev.x);
                            mpz_mod(rel.x, rel.x, N);
                            mpz_init(rel.Qx);
                            mpz_mul(rel.Qx, qval, poly.a);
                            mpz_mul(rel.Qx, rel.Qx, prev.Qx);
                            rel.exp_parity.resize(fb_size);
                            for (size_t i = 0; i < fb_size; i++)
                                rel.exp_parity[i] = (parity[i] + prev.exp_parity[i]) & 1;
                            relations.push_back(rel);
                            slp_combined++;
                        } else {
                            PartialRelation pr;
                            mpz_init_set(pr.x, ax_plus_b);
                            mpz_init(pr.Qx);
                            mpz_mul(pr.Qx, qval, poly.a);
                            pr.exp_parity = parity;
                            pr.lp = cf;
                            plist.push_back(pr);
                        }
                    } else if (cf <= lp2_bound * lp2_bound) {
                        // Potential double large prime
                        // Try to factor the cofactor
                        // If cf = p1 * p2 with both prime and <= lp2_bound
                        unsigned long p1 = 0, p2 = 0;

                        // Quick trial divide by small primes (shouldn't happen but check)
                        for (unsigned long sp = 2; sp * sp <= cf && sp < 1000; sp++) {
                            if (cf % sp == 0) {
                                p1 = sp;
                                p2 = cf / sp;
                                break;
                            }
                        }

                        if (p1 == 0) {
                            // Try Pollard rho on the cofactor
                            if (cf > lp2_bound) {
                                // cf might be prime (skip) or product of two primes
                                if (is_probable_prime(cf)) continue; // single prime too large

                                // Simple Pollard rho
                                unsigned long x0 = 2, y0 = 2, d0 = 1;
                                auto f = [cf](unsigned long v) -> unsigned long {
                                    return ((__uint128_t)v * v + 1) % cf;
                                };
                                x0 = f(x0); y0 = f(f(y0));
                                for (int iter = 0; iter < 10000 && d0 == 1; iter++) {
                                    long diff = (long)x0 - (long)y0;
                                    if (diff < 0) diff = -diff;
                                    d0 = std::__gcd((unsigned long)diff, cf);
                                    x0 = f(x0); y0 = f(f(y0));
                                }
                                if (d0 > 1 && d0 < cf) {
                                    p1 = d0;
                                    p2 = cf / d0;
                                }
                            }
                        }

                        if (p1 > 0 && p2 > 0 && p1 <= lp2_bound && p2 <= lp2_bound
                            && is_probable_prime(p1) && is_probable_prime(p2)) {
                            if (p1 > p2) std::swap(p1, p2);

                            // Store DLP relation and check graph for cycles
                            DLP_Relation dlp;
                            mpz_init_set(dlp.x, ax_plus_b);
                            mpz_init(dlp.Qx);
                            mpz_mul(dlp.Qx, qval, poly.a);
                            dlp.exp_parity = parity;
                            dlp.lp1 = p1;
                            dlp.lp2 = p2;
                            size_t dlp_idx = dlp_rels.size();
                            dlp_rels.push_back(dlp);

                            // Add to graph
                            dlp_graph[p1].push_back(dlp_idx);
                            dlp_graph[p2].push_back(dlp_idx);

                            // Check if we can combine (find path in graph)
                            // Simple: if both p1 and p2 have other edges, try to combine
                            if (dlp_graph[p1].size() >= 2 && p1 == p2) {
                                // Self-loop: p1 == p2, combine the two relations
                                size_t other = dlp_graph[p1][0];
                                if (other == dlp_idx && dlp_graph[p1].size() > 1)
                                    other = dlp_graph[p1][1];
                                if (other != dlp_idx) {
                                    auto &o = dlp_rels[other];
                                    Relation rel;
                                    mpz_init(rel.x);
                                    mpz_mul(rel.x, dlp.x, o.x);
                                    mpz_mod(rel.x, rel.x, N);
                                    mpz_init(rel.Qx);
                                    mpz_mul(rel.Qx, dlp.Qx, o.Qx);
                                    rel.exp_parity.resize(fb_size);
                                    for (size_t i = 0; i < fb_size; i++)
                                        rel.exp_parity[i] = (dlp.exp_parity[i] + o.exp_parity[i]) & 1;
                                    relations.push_back(rel);
                                    dlp_combined++;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (poly_count % 500 == 0) {
            fprintf(stderr, "SIQS-opt: %zu/%zu rels (full=%d, slp=%d, dlp=%d), %d polys, %.1fs\n",
                    relations.size(), relations_needed, full_rels, slp_combined, dlp_combined,
                    poly_count, elapsed_sec());
        }
    }

    fprintf(stderr, "SIQS-opt: %zu relations from %d polys (full=%d, slp=%d, dlp=%d) in %.1fs\n",
            relations.size(), poly_count, full_rels, slp_combined, dlp_combined, elapsed_sec());

    if (relations.size() < fb_size + 1) {
        fprintf(stderr, "FAIL: not enough relations\n");
        return 1;
    }

    // Gaussian elimination (same as siqs_factor.cpp)
    size_t nrels = relations.size();
    size_t ncols = fb_size + 1;
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
                while (mpz_divisible_ui_p(temp, fb[j].p)) {
                    mpz_divexact_ui(temp, temp, fb[j].p);
                    total_exps[j]++;
                }
            }
            if (mpz_cmp_ui(temp, 1) > 0) {
                mpz_mul(lp_product, lp_product, temp);
            }
        }

        bool lp_ok = (mpz_cmp_ui(lp_product, 1) == 0) || mpz_perfect_square_p(lp_product);
        int odd_count = 0;
        for (size_t j = 0; j < fb_size; j++)
            if (total_exps[j] & 1) odd_count++;

        if (sign_neg || odd_count > 0 || !lp_ok) continue;

        mpz_set_ui(rhs, 1);
        for (size_t j = 0; j < fb_size; j++) {
            if (total_exps[j] > 0) {
                mpz_ui_pow_ui(temp, fb[j].p, total_exps[j] / 2);
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
            fprintf(stderr, "SIQS-opt: factored in %.3fs\n", elapsed_sec());
            found = true;
        }
        if (!found) {
            mpz_add(temp, lhs, rhs);
            mpz_gcd(factor1, temp, N);
            if (mpz_cmp_ui(factor1, 1) > 0 && mpz_cmp(factor1, N) < 0) {
                mpz_divexact(factor2, N, factor1);
                if (mpz_cmp(factor1, factor2) > 0) mpz_swap(factor1, factor2);
                gmp_printf("%Zd %Zd\n", factor1, factor2);
                fprintf(stderr, "SIQS-opt: factored in %.3fs\n", elapsed_sec());
                found = true;
            }
        }
    }

    if (!found) {
        fprintf(stderr, "FAIL: no factor found from linear algebra\n");
        return 1;
    }

    // Cleanup
    for (auto &rel : relations) { mpz_clear(rel.x); mpz_clear(rel.Qx); }
    for (auto &[k, v] : slp_partials) {
        for (auto &pr : v) { mpz_clear(pr.x); mpz_clear(pr.Qx); }
    }
    for (auto &dlp : dlp_rels) { mpz_clear(dlp.x); mpz_clear(dlp.Qx); }

    mpz_clear(poly.a); mpz_clear(poly.b); mpz_clear(poly.c);
    gmp_randclear(rng);
    mpz_clear(N); mpz_clear(target_a); mpz_clear(two_n);
    mpz_clear(ax_plus_b); mpz_clear(qval); mpz_clear(cofactor); mpz_clear(temp);
    mpz_clear(lhs); mpz_clear(rhs); mpz_clear(factor1); mpz_clear(factor2);
    mpz_clear(lp_product);

    return 0;
}
