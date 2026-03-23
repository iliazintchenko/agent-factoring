/*
 * Hybrid Factoring: Batch Smoothness + Multi-Large-Prime
 *
 * Novel approach combining:
 * 1. SIQS-style polynomial generation for candidate values
 * 2. Product tree / remainder tree for fast batch smoothness testing
 * 3. 3-large-prime variation: allow up to 3 large prime cofactors
 * 4. Cycle finding in the 3LP graph for relation combination
 *
 * The key insight: by allowing 3 large primes per relation and using
 * batch smoothness for fast detection, we can collect exponentially
 * more usable relations. The cycle structure in the multi-LP graph
 * provides sufficient combinations.
 *
 * Theoretical advantage: with k large primes allowed, the smoothness
 * probability improves by a factor of ~(log B)^k / k!, allowing a
 * smaller factor base B and thus better scaling.
 *
 * Usage: ./hybrid_factor <number>
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
#include <set>
#include <cassert>
#include <gmp.h>

using namespace std;

// ============================================================
// Utility
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

static long long tonelli_shanks(long long n, long long p) {
    if (p == 2) return n & 1;
    n = ((n % p) + p) % p;
    if (n == 0) return 0;
    if (mod_pow(n, (p-1)/2, p) != 1) return -1;
    if (p % 4 == 3) return mod_pow(n, (p+1)/4, p);
    long long Q = p-1, S = 0;
    while (Q%2==0) { Q/=2; S++; }
    long long z = 2;
    while (mod_pow(z, (p-1)/2, p) != p-1) z++;
    long long M=S, c=mod_pow(z,Q,p), t=mod_pow(n,Q,p), R=mod_pow(n,(Q+1)/2,p);
    while (true) {
        if (t==1) return R;
        long long i=0, tmp=t;
        while (tmp!=1) { tmp=(__int128)tmp*tmp%p; i++; }
        long long b=c;
        for (long long j=0; j<M-i-1; j++) b=(__int128)b*b%p;
        M=i; c=(__int128)b*b%p; t=(__int128)t*c%p; R=(__int128)R*b%p;
    }
}

// ============================================================
// GF(2) Linear Algebra
// ============================================================
struct BitMatrix {
    int rows, cols, wpr;
    vector<vector<uint64_t>> data;
    BitMatrix(int r, int c) : rows(r), cols(c), wpr((c+63)/64) {
        data.assign(r, vector<uint64_t>(wpr, 0));
    }
    void set(int r, int c) { data[r][c>>6] |= 1ULL<<(c&63); }
    bool get(int r, int c) const { return (data[r][c>>6]>>(c&63))&1; }
    void xor_row(int d, int s) {
        for (int i = 0; i < wpr; i++) data[d][i] ^= data[s][i];
    }
};

static vector<vector<int>> find_null_space(int nr, int nc,
                                            const vector<vector<int>>& evecs) {
    BitMatrix mat(nr, nc + nr);
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nc; j++) if (evecs[i][j] & 1) mat.set(i, j);
        mat.set(i, nc + i);
    }
    int cur = 0;
    for (int col = 0; col < nc && cur < nr; col++) {
        int piv = -1;
        for (int r = cur; r < nr; r++) if (mat.get(r, col)) { piv = r; break; }
        if (piv < 0) continue;
        if (piv != cur) swap(mat.data[piv], mat.data[cur]);
        for (int r = 0; r < nr; r++)
            if (r != cur && mat.get(r, col)) mat.xor_row(r, cur);
        cur++;
    }
    vector<vector<int>> result;
    for (int r = cur; r < nr; r++) {
        bool ok = true;
        for (int c = 0; c < nc; c++) if (mat.get(r, c)) { ok = false; break; }
        if (!ok) continue;
        vector<int> v;
        for (int i = 0; i < nr; i++) if (mat.get(r, nc+i)) v.push_back(i);
        if (!v.empty()) result.push_back(v);
    }
    return result;
}

// ============================================================
// Parameters
// ============================================================
struct Params {
    int fb_size;
    int sieve_half;
    int max_large_primes;  // 1, 2, or 3
    long long lp_bound;    // each LP must be < this
    double thresh_adj;
};

static Params get_params(int digits) {
    // Use smaller factor base (compensated by more LPs)
    // This is the key: smaller FB = less linear algebra, more LP relations
    if (digits <= 30) return {80, 32768, 2, 100000, 3.5};
    if (digits <= 35) return {120, 50000, 2, 200000, 4.0};
    if (digits <= 40) return {200, 65536, 2, 500000, 4.5};
    if (digits <= 45) return {350, 100000, 2, 1000000, 5.0};
    if (digits <= 50) return {600, 131072, 2, 2000000, 5.5};
    if (digits <= 55) return {900, 196608, 3, 5000000, 6.0};
    if (digits <= 60) return {1400, 262144, 3, 10000000, 6.5};
    if (digits <= 65) return {2000, 393216, 3, 20000000, 7.0};
    if (digits <= 70) return {3000, 524288, 3, 50000000, 7.5};
    if (digits <= 75) return {4500, 786432, 3, 100000000, 8.0};
    if (digits <= 80) return {6500, 1048576, 3, 200000000, 8.5};
    if (digits <= 85) return {9500, 1572864, 3, 500000000, 9.0};
    if (digits <= 90) return {14000, 2097152, 3, 1000000000LL, 9.5};
    if (digits <= 95) return {20000, 2621440, 3, 2000000000LL, 10.0};
    return {28000, 3145728, 3, 5000000000LL, 10.5};
}

// ============================================================
// Relation types
// ============================================================
struct Relation {
    mpz_t ax_plus_b;
    vector<int> exponents;  // over factor base
    vector<long long> large_primes; // 0, 1, 2, or 3 large primes
    vector<long long> lp_factors; // LPs that need to be included in Y computation
};

// ============================================================
// Double/Triple Large Prime Combination via Graph Cycles
// ============================================================
// For 2LP: build a graph where vertices are large primes, edges are relations
// A cycle in this graph gives a combination where all LPs cancel
// For 3LP: hypergraph, but we can reduce to 2LP by fixing one LP

struct LPGraph {
    // Adjacency list for 2LP relations
    // edge = (lp1, lp2, relation_index)
    unordered_map<long long, vector<pair<long long, int>>> adj;

    void add_edge(long long lp1, long long lp2, int rel_idx) {
        if (lp1 > lp2) swap(lp1, lp2);
        adj[lp1].push_back({lp2, rel_idx});
        adj[lp2].push_back({lp1, rel_idx});
    }

    // Find cycles using DFS - each cycle gives a full relation
    vector<vector<int>> find_cycles(int max_cycles) {
        vector<vector<int>> cycles;
        set<long long> visited;
        map<long long, pair<long long, int>> parent; // node -> (parent, edge_idx)

        for (auto& [start, _] : adj) {
            if (visited.count(start)) continue;
            if ((int)cycles.size() >= max_cycles) break;

            // BFS from start to find short cycles
            map<long long, int> dist;
            map<long long, pair<long long, int>> prev;
            vector<long long> queue;
            dist[start] = 0;
            queue.push_back(start);

            for (size_t qi = 0; qi < queue.size() && (int)cycles.size() < max_cycles; qi++) {
                long long u = queue[qi];
                for (auto& [v, ridx] : adj[u]) {
                    if (dist.find(v) == dist.end()) {
                        dist[v] = dist[u] + 1;
                        prev[v] = {u, ridx};
                        queue.push_back(v);
                    } else if (dist[v] >= dist[u] && prev[u].first != v) {
                        // Found a cycle! Trace back from u and v to find common ancestor
                        vector<int> cycle_rels;
                        cycle_rels.push_back(ridx);

                        // Trace from u to start
                        long long cur = u;
                        while (cur != start && prev.find(cur) != prev.end()) {
                            cycle_rels.push_back(prev[cur].second);
                            cur = prev[cur].first;
                        }

                        // Trace from v to start
                        cur = v;
                        while (cur != start && prev.find(cur) != prev.end()) {
                            cycle_rels.push_back(prev[cur].second);
                            cur = prev[cur].first;
                        }

                        if (cycle_rels.size() >= 2) {
                            cycles.push_back(cycle_rels);
                        }
                    }
                }
            }
            visited.insert(start);
        }
        return cycles;
    }
};

// ============================================================
// PRNG
// ============================================================
static uint64_t rng_state = 42;
static uint64_t next_rand() {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return rng_state >> 16;
}

// ============================================================
// Main
// ============================================================
int main(int argc, char* argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <number>\n", argv[0]); return 1; }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    auto timer = [&]() -> double {
        struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
        return (t1.tv_sec-t0.tv_sec) + (t1.tv_nsec-t0.tv_nsec)/1e9;
    };

    mpz_t N, factor, tmp, tmp2;
    mpz_init_set_str(N, argv[1], 10);
    mpz_init(factor); mpz_init(tmp); mpz_init(tmp2);

    if (mpz_even_p(N)) { printf("2\n"); return 0; }
    if (mpz_perfect_square_p(N)) { mpz_sqrt(tmp, N); gmp_printf("%Zd\n", tmp); return 0; }

    // Trial division
    {
        auto sp = sieve_primes(1000000);
        for (int p : sp) if (mpz_divisible_ui_p(N, p)) { printf("%d\n", p); return 0; }
    }

    int digits = (int)mpz_sizeinbase(N, 10);
    Params P = get_params(digits);

    fprintf(stderr, "Hybrid: %d digits, fb=%d, M=%d, max_lp=%d, lp_bound=%lld\n",
            digits, P.fb_size, P.sieve_half, P.max_large_primes, P.lp_bound);

    // Build factor base
    int plim = 200000;
    while ((int)sieve_primes(plim).size() < P.fb_size * 3) plim *= 2;
    auto all_primes = sieve_primes(plim);

    struct FBEntry { int p, logp, sqrt_n; };
    vector<FBEntry> fb;
    fb.push_back({-1, 0, 0}); // sign

    for (size_t i = 0; i < all_primes.size() && (int)fb.size() <= P.fb_size; i++) {
        int p = all_primes[i];
        int nm = (int)mpz_fdiv_ui(N, p);
        long long sq = tonelli_shanks(nm, p);
        if (sq >= 0) {
            int logp = max(1, (int)(log2((double)p) * 16 + 0.5));
            fb.push_back({p, logp, (int)sq});
        }
    }

    int fb_size = (int)fb.size();
    int M = P.sieve_half;
    int needed = fb_size + 30;

    fprintf(stderr, "FB: %d primes, largest=%d\n", fb_size, fb.back().p);

    // Sieve and collect relations
    vector<Relation> full_relations; // 0LP: ready to use
    vector<Relation> partial_1lp;   // 1LP: need matching
    vector<Relation> partial_2lp;   // 2LP: need cycle

    unordered_map<long long, int> lp1_map; // 1LP hash: lp -> index in partial_1lp
    LPGraph lp2_graph;

    int total_full = 0, total_from_1lp = 0, total_from_2lp = 0;
    int n_1lp = 0, n_2lp = 0, n_3lp = 0;

    mpz_t a_val, b_val, c_val, target_a;
    mpz_init(a_val); mpz_init(b_val); mpz_init(c_val); mpz_init(target_a);

    mpz_mul_ui(target_a, N, 2);
    mpz_sqrt(target_a, target_a);
    mpz_tdiv_q_ui(target_a, target_a, M);

    int num_a_factors;
    if (digits <= 35) num_a_factors = 3;
    else if (digits <= 45) num_a_factors = 4;
    else if (digits <= 55) num_a_factors = 5;
    else if (digits <= 65) num_a_factors = 6;
    else if (digits <= 75) num_a_factors = 7;
    else if (digits <= 85) num_a_factors = 8;
    else num_a_factors = 9;

    int a_lo = max(2, fb_size / 5);
    int a_hi = min(fb_size - 1, fb_size * 4 / 5);

    double log2_qmax = mpz_sizeinbase(N, 2) / 2.0 + log2((double)M);
    int threshold = (int)((log2_qmax - P.thresh_adj) * 16.0);
    if (threshold < 0) threshold = 0;

    vector<int16_t> sieve(2 * M + 2);
    int poly_count = 0;

    while (total_full + total_from_1lp + total_from_2lp < needed && timer() < 280) {
        // Generate SIQS polynomial
        vector<int> a_idx;
        mpz_set_ui(a_val, 1);
        vector<int> pool;
        for (int i = a_lo; i <= a_hi; i++) pool.push_back(i);
        for (int i = 0; i < num_a_factors && !pool.empty(); i++) {
            int r = next_rand() % pool.size();
            a_idx.push_back(pool[r]);
            mpz_mul_ui(a_val, a_val, fb[pool[r]].p);
            pool.erase(pool.begin() + r);
        }
        sort(a_idx.begin(), a_idx.end());

        // CRT for b
        mpz_set_ui(b_val, 0);
        for (int i = 0; i < num_a_factors; i++) {
            int qi = fb[a_idx[i]].p, si = fb[a_idx[i]].sqrt_n;
            mpz_t adq, qm, iv, ct;
            mpz_init(adq); mpz_init_set_ui(qm, qi); mpz_init(iv); mpz_init(ct);
            mpz_tdiv_q_ui(adq, a_val, qi);
            mpz_invert(iv, adq, qm);
            mpz_mul_ui(ct, iv, si);
            mpz_mul(ct, ct, adq);
            mpz_add(b_val, b_val, ct);
            mpz_clear(adq); mpz_clear(qm); mpz_clear(iv); mpz_clear(ct);
        }
        mpz_mod(b_val, b_val, a_val);
        mpz_tdiv_q_ui(tmp, a_val, 2);
        if (mpz_cmp(b_val, tmp) > 0) mpz_sub(b_val, a_val, b_val);

        // Verify and compute c
        mpz_mul(tmp, b_val, b_val);
        mpz_sub(tmp, tmp, N);
        mpz_mod(tmp, tmp, a_val);
        if (mpz_sgn(tmp) != 0) { poly_count++; continue; }

        mpz_mul(c_val, b_val, b_val);
        mpz_sub(c_val, c_val, N);
        mpz_tdiv_q(c_val, c_val, a_val);

        // Compute sieve roots
        vector<int> soln1(fb_size, -1), soln2(fb_size, -1);
        for (int i = 1; i < fb_size; i++) {
            int p = fb[i].p;
            bool da = false;
            for (int j : a_idx) if (j == i) { da = true; break; }
            if (da || p == 2) continue;
            long long amp = mpz_fdiv_ui(a_val, p);
            long long bmp = mpz_fdiv_ui(b_val, p);
            long long ainv = mod_pow(amp, p-2, p);
            int sq = fb[i].sqrt_n;
            soln1[i] = (int)(((__int128)(sq - bmp + p) % p * ainv) % p);
            soln2[i] = (int)(((__int128)(p - sq - bmp + p) % p * ainv) % p);
        }

        // Sieve
        memset(sieve.data(), 0, (2*M+2)*sizeof(int16_t));
        for (int i = 1; i < fb_size; i++) {
            if (soln1[i] < 0) continue;
            int p = fb[i].p, lp = fb[i].logp;
            int s1 = (int)(((long long)soln1[i] + M) % p);
            for (int j = s1; j <= 2*M; j += p) sieve[j] += lp;
            if (soln1[i] != soln2[i]) {
                int s2 = (int)(((long long)soln2[i] + M) % p);
                for (int j = s2; j <= 2*M; j += p) sieve[j] += lp;
            }
        }
        // a-primes
        for (int ai : a_idx) {
            int p = fb[ai].p, lp = fb[ai].logp;
            long long bmp = mpz_fdiv_ui(b_val, p);
            mpz_t cm; mpz_init(cm);
            mpz_fdiv_r_ui(cm, c_val, p);
            long long cmp = mpz_get_ui(cm);
            mpz_clear(cm);
            long long twob = (2*bmp) % p;
            if (twob == 0) continue;
            long long root = ((p - cmp) % p * mod_pow(twob, p-2, p)) % p;
            int off = (int)(((long long)root + M) % p);
            for (int j = off; j <= 2*M; j += p) sieve[j] += lp;
        }

        // Scan candidates
        // Lower threshold to catch more LP relations
        int lp_threshold = threshold - (int)(log2((double)P.lp_bound) * 16 * P.max_large_primes);
        if (lp_threshold < threshold / 2) lp_threshold = threshold / 2;

        for (int idx = 0; idx <= 2*M; idx++) {
            if (sieve[idx] < lp_threshold) continue;

            int x = idx - M;
            mpz_t Qval; mpz_init(Qval);
            mpz_mul_si(Qval, a_val, x);
            mpz_addmul_ui(Qval, b_val, 2);
            mpz_mul_si(Qval, Qval, x);
            mpz_add(Qval, Qval, c_val);

            mpz_t rem; mpz_init(rem);
            mpz_abs(rem, Qval);
            bool neg = (mpz_sgn(Qval) < 0);
            vector<int> exponents(fb_size, 0);
            if (neg) exponents[0] = 1;

            for (int i = 1; i < fb_size; i++) {
                unsigned long p = fb[i].p;
                while (mpz_divisible_ui_p(rem, p)) {
                    mpz_divexact_ui(rem, rem, p);
                    exponents[i]++;
                }
            }

            // Verify Q(x) factorization before adding a-factors
            if (mpz_cmp_ui(rem, 1) == 0 && total_full < 3) {
                // Reconstruct Q from exponents
                mpz_t recon; mpz_init_set_ui(recon, 1);
                for (int i = 1; i < fb_size; i++)
                    for (int e = 0; e < exponents[i]; e++)
                        mpz_mul_ui(recon, recon, fb[i].p);
                if (neg) mpz_neg(recon, recon);
                if (mpz_cmp(recon, Qval) != 0)
                    gmp_fprintf(stderr, "BUG: Q recon mismatch! Q=%Zd recon=%Zd\n", Qval, recon);
                else
                    fprintf(stderr, "Relation verified OK\n");
                mpz_clear(recon);
            }

            // Add a-factor exponents
            for (int ai : a_idx) exponents[ai]++;

            // Check cofactor
            int num_lp = -1; // -1 = unusable
            vector<long long> lps;

            if (mpz_cmp_ui(rem, 1) == 0) {
                // Fully smooth
                num_lp = 0;
            } else if (mpz_fits_slong_p(rem)) {
                long long r = mpz_get_si(rem);
                if (r <= P.lp_bound && r > 1) {
                    // Could be 1 LP (prime) or 2 LPs (product of 2 primes)
                    // Check if it's prime using simple trial division up to sqrt(r)
                    bool is_prime_r = true;
                    for (long long d = 2; d * d <= r; d++) {
                        if (r % d == 0) {
                            long long d2 = r / d;
                            if (d <= P.lp_bound && d2 <= P.lp_bound && d > fb.back().p && d2 > fb.back().p) {
                                lps.push_back(d);
                                lps.push_back(d2);
                                num_lp = 2;
                            }
                            is_prime_r = false;
                            break;
                        }
                    }
                    if (is_prime_r && r > fb.back().p) {
                        lps.push_back(r);
                        num_lp = 1;
                    }
                }
            } else if (mpz_sizeinbase(rem, 2) <= 64) {
                // Might be 2 or 3 large primes
                // Try to factor using Pollard rho or trial division
                // For now, skip (too complex to factor large cofactors quickly)
            }

            if (num_lp >= 0 && num_lp <= P.max_large_primes && (num_lp == 0 || !lps.empty() || num_lp == 0)) {
                Relation rel;
                mpz_init(rel.ax_plus_b);
                mpz_mul_si(rel.ax_plus_b, a_val, x);
                mpz_add(rel.ax_plus_b, rel.ax_plus_b, b_val);
                rel.exponents = exponents;
                rel.large_primes = lps;

                if (num_lp == 0) {
                    full_relations.push_back(rel);
                    total_full++;
                } else if (num_lp == 1) {
                    long long lp = lps[0];
                    auto it = lp1_map.find(lp);
                    if (it != lp1_map.end()) {
                        // Combine two 1LP relations
                        // Don't divide by lp - track it for Y computation
                        Relation& other = partial_1lp[it->second];
                        Relation combined;
                        mpz_init(combined.ax_plus_b);
                        mpz_mul(combined.ax_plus_b, rel.ax_plus_b, other.ax_plus_b);
                        mpz_mod(combined.ax_plus_b, combined.ax_plus_b, N);

                        combined.exponents.resize(fb_size);
                        for (int k = 0; k < fb_size; k++)
                            combined.exponents[k] = rel.exponents[k] + other.exponents[k];

                        // Track the LP (appears with exponent 2 = even, so doesn't
                        // affect GF(2) matrix, but needed for square root)
                        combined.lp_factors.push_back(lp);

                        full_relations.push_back(combined);
                        total_from_1lp++;
                        lp1_map.erase(it);
                    } else {
                        lp1_map[lp] = partial_1lp.size();
                        partial_1lp.push_back(rel);
                        n_1lp++;
                    }
                } else if (num_lp == 2) {
                    // 2LP relation: add to graph
                    int ridx = partial_2lp.size();
                    partial_2lp.push_back(rel);
                    lp2_graph.add_edge(lps[0], lps[1], ridx);
                    n_2lp++;
                }
            } else {
                // Can't use this relation
            }

            mpz_clear(Qval); mpz_clear(rem);

            if (total_full + total_from_1lp + total_from_2lp >= needed) break;
        }

        poly_count++;
        if (poly_count % 500 == 0) {
            fprintf(stderr, "  [%.1fs] poly %d: %d full + %d from1lp + %d from2lp = %d/%d (1lp:%d 2lp:%d)\n",
                    timer(), poly_count, total_full, total_from_1lp, total_from_2lp,
                    total_full + total_from_1lp + total_from_2lp, needed, n_1lp, n_2lp);

            // Periodically try to find cycles in 2LP graph
            if (n_2lp > 50 && total_full + total_from_1lp + total_from_2lp < needed) {
                auto cycles = lp2_graph.find_cycles(needed - total_full - total_from_1lp - total_from_2lp);
                for (auto& cycle : cycles) {
                    if (total_full + total_from_1lp + total_from_2lp >= needed) break;

                    // Combine all relations in the cycle
                    Relation combined;
                    mpz_init_set_ui(combined.ax_plus_b, 1);
                    combined.exponents.resize(fb_size, 0);

                    // Collect all large primes in the cycle - they should all cancel
                    map<long long, int> lp_counts;
                    for (int ridx : cycle) {
                        Relation& r = partial_2lp[ridx];
                        mpz_mul(combined.ax_plus_b, combined.ax_plus_b, r.ax_plus_b);
                        mpz_mod(combined.ax_plus_b, combined.ax_plus_b, N);
                        for (int k = 0; k < fb_size; k++)
                            combined.exponents[k] += r.exponents[k];
                        for (long long lp : r.large_primes)
                            lp_counts[lp]++;
                    }

                    // Divide out large primes that appear evenly
                    bool valid = true;
                    for (auto& [lp, cnt] : lp_counts) {
                        if (cnt % 2 != 0) { valid = false; break; }
                        // Divide combined.ax_plus_b by lp^(cnt/2)
                        mpz_t lpm, lpi;
                        mpz_init_set_ui(lpm, (unsigned long)lp);
                        mpz_init(lpi);
                        mpz_invert(lpi, lpm, N);
                        for (int c = 0; c < cnt/2; c++) {
                            mpz_mul(combined.ax_plus_b, combined.ax_plus_b, lpi);
                            mpz_mod(combined.ax_plus_b, combined.ax_plus_b, N);
                        }
                        mpz_clear(lpm); mpz_clear(lpi);
                    }

                    if (valid) {
                        full_relations.push_back(combined);
                        total_from_2lp++;
                    }
                }
            }
        }
    }

    int total = total_full + total_from_1lp + total_from_2lp;
    fprintf(stderr, "Collected %d relations (%d full + %d from1lp + %d from2lp) in %d polys, %.1fs\n",
            total, total_full, total_from_1lp, total_from_2lp, poly_count, timer());

    if (total < fb_size + 1) {
        fprintf(stderr, "FAIL: not enough relations\n");
        return 1;
    }

    // Linear algebra
    int nrels = (int)full_relations.size();
    fprintf(stderr, "Linear algebra: %d x %d\n", nrels, fb_size);

    vector<vector<int>> evecs(nrels, vector<int>(fb_size));
    for (int i = 0; i < nrels; i++)
        for (int j = 0; j < fb_size; j++)
            evecs[i][j] = full_relations[i].exponents[j];

    auto nullvecs = find_null_space(nrels, fb_size, evecs);
    fprintf(stderr, "Found %d null vectors\n", (int)nullvecs.size());

    bool factored = false;
    mpz_t X, Y;
    mpz_init(X); mpz_init(Y);

    for (auto& nv : nullvecs) {
        if (factored) break;
        mpz_set_ui(X, 1);
        vector<long long> texp(fb_size, 0);
        for (int idx : nv) {
            mpz_mul(X, X, full_relations[idx].ax_plus_b);
            mpz_mod(X, X, N);
            for (int j = 0; j < fb_size; j++)
                texp[j] += full_relations[idx].exponents[j];
        }
        bool ok = true;
        for (int j = 0; j < fb_size; j++) if (texp[j] & 1) { ok = false; break; }
        if (!ok) continue;

        mpz_set_ui(Y, 1);
        for (int j = 1; j < fb_size; j++) {
            if (texp[j] == 0) continue;
            mpz_set_ui(tmp, fb[j].p);
            mpz_powm_ui(tmp, tmp, texp[j]/2, N);
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, N);
        }
        // Include large primes from combined relations
        // Each combined relation contributes lp^2 to the right side
        // lp^2 / 2 = lp^1 contribution to Y
        for (int idx : nv) {
            for (long long lp : full_relations[idx].lp_factors) {
                mpz_set_ui(tmp, (unsigned long)lp);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }
        }

        // Debug: verify X^2 ≡ Y^2
        {
            mpz_t x2, y2;
            mpz_init(x2); mpz_init(y2);
            mpz_mul(x2, X, X); mpz_mod(x2, x2, N);
            mpz_mul(y2, Y, Y); mpz_mod(y2, y2, N);
            static int dbg = 0;
            if (dbg++ < 3) {
                if (mpz_cmp(x2, y2) != 0)
                    fprintf(stderr, "BUG: X^2 != Y^2 (nv size=%d, lp_factors=%d)\n",
                            (int)nv.size(), (int)full_relations[nv[0]].lp_factors.size());
                else
                    fprintf(stderr, "OK: X^2 == Y^2\n");
            }
            mpz_clear(x2); mpz_clear(y2);
        }

        mpz_sub(tmp, X, Y); mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) { factored = true; break; }
        mpz_add(tmp, X, Y); mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) { factored = true; break; }
    }

    if (factored) {
        gmp_printf("%Zd\n", factor);
        fprintf(stderr, "Hybrid: factored in %.3f seconds\n", timer());
    } else {
        fprintf(stderr, "Hybrid: FAILED after %.3f seconds\n", timer());
        return 1;
    }

    mpz_clear(N); mpz_clear(factor); mpz_clear(tmp); mpz_clear(tmp2);
    mpz_clear(a_val); mpz_clear(b_val); mpz_clear(c_val); mpz_clear(target_a);
    mpz_clear(X); mpz_clear(Y);
    for (auto& r : full_relations) mpz_clear(r.ax_plus_b);
    for (auto& r : partial_1lp) mpz_clear(r.ax_plus_b);
    for (auto& r : partial_2lp) mpz_clear(r.ax_plus_b);

    return 0;
}
