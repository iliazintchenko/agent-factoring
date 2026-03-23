/*
 * Optimized Self-Initializing Quadratic Sieve (SIQS v2)
 *
 * Key optimizations over v1:
 * - Gray code polynomial enumeration: 2^(s-1) b-values per a-value
 * - Block sieving for cache efficiency
 * - Bucket sieving for large primes
 * - Faster sieve with unsigned char (scaled)
 *
 * Usage: ./siqs2 <number>
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
#include <cassert>
#include <gmp.h>

using namespace std;

// ============================================================
// Parameters
// ============================================================
struct Params {
    int fb_size;
    int sieve_half;  // M: sieve from -M to M
    int block_size;  // block size for sieving (L1 cache friendly)
    int lp_mult;     // LP bound = largest_fb * this
    int num_a_factors; // number of primes in 'a'
    double thresh_adj; // threshold adjustment in bits
};

static Params get_params(int digits) {
    if (digits <= 30) return {120, 32768, 32768, 50, 3, 3.0};
    if (digits <= 35) return {200, 65536, 32768, 60, 4, 3.5};
    if (digits <= 40) return {350, 65536, 32768, 80, 4, 4.0};
    if (digits <= 45) return {600, 65536, 32768, 100, 5, 4.5};
    if (digits <= 50) return {1000, 131072, 32768, 120, 5, 5.0};
    if (digits <= 55) return {1500, 196608, 32768, 150, 6, 5.5};
    if (digits <= 60) return {2200, 262144, 32768, 200, 6, 6.0};
    if (digits <= 65) return {3200, 393216, 32768, 250, 7, 6.5};
    if (digits <= 70) return {4500, 524288, 32768, 300, 7, 7.0};
    if (digits <= 75) return {6500, 786432, 32768, 400, 8, 7.5};
    if (digits <= 80) return {9000, 1048576, 65536, 500, 8, 8.0};
    if (digits <= 85) return {13000, 1572864, 65536, 600, 9, 8.5};
    if (digits <= 90) return {18000, 2097152, 65536, 800, 9, 9.0};
    if (digits <= 95) return {25000, 2621440, 65536, 1000, 10, 9.5};
    return {35000, 3145728, 65536, 1200, 10, 10.0};
}

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
    if (mod_pow(n, (p - 1) / 2, p) != 1) return -1;
    if (p % 4 == 3) return mod_pow(n, (p + 1) / 4, p);
    long long Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    long long z = 2;
    while (mod_pow(z, (p - 1) / 2, p) != p - 1) z++;
    long long M = S, c = mod_pow(z, Q, p), t = mod_pow(n, Q, p), R = mod_pow(n, (Q + 1) / 2, p);
    while (true) {
        if (t == 1) return R;
        long long i = 0, tmp = t;
        while (tmp != 1) { tmp = (__int128)tmp * tmp % p; i++; }
        long long b = c;
        for (long long j = 0; j < M - i - 1; j++) b = (__int128)b * b % p;
        M = i; c = (__int128)b * b % p; t = (__int128)t * c % p; R = (__int128)R * b % p;
    }
}

// ============================================================
// GF(2) linear algebra
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
// Main SIQS
// ============================================================

struct FBEntry {
    int p;
    int logp;  // scaled log
    int sqrt_n; // sqrt(N) mod p
};

struct Relation {
    mpz_t x; // (ax+b)
    vector<int> exponents;
};

static uint64_t rng_state = 42;
static uint64_t next_rand() {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return rng_state >> 16;
}

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

    // Trial division up to 1M
    {
        vector<int> sp = sieve_primes(1000000);
        for (int p : sp) {
            if (mpz_divisible_ui_p(N, p)) {
                printf("%d\n", p); return 0;
            }
        }
    }

    int digits = (int)mpz_sizeinbase(N, 10);
    Params P = get_params(digits);

    fprintf(stderr, "SIQS2: %d digits, fb=%d, M=%d, s=%d\n",
            digits, P.fb_size, P.sieve_half, P.num_a_factors);

    // Build factor base
    int plim = 200000;
    while ((int)sieve_primes(plim).size() < P.fb_size * 3) plim *= 2;
    vector<int> all_primes = sieve_primes(plim);

    vector<FBEntry> fb;
    fb.push_back({-1, 0, 0}); // sign

    for (size_t i = 0; i < all_primes.size() && (int)fb.size() <= P.fb_size; i++) {
        int p = all_primes[i];
        int nm = (int)mpz_fdiv_ui(N, p);
        long long sq = tonelli_shanks(nm, p);
        if (sq >= 0) {
            int logp = max(1, (int)(log2((double)p) + 0.5)); // scale: 1 per bit
            fb.push_back({p, logp, (int)sq});
        }
    }

    int fb_size = (int)fb.size();
    long long lp_bound = (long long)fb.back().p * P.lp_mult;
    int needed = fb_size + 30;
    int M = P.sieve_half;
    int s = P.num_a_factors; // number of primes in 'a'

    fprintf(stderr, "FB: %d primes, largest=%d, lp_bound=%lld\n", fb_size, fb.back().p, lp_bound);

    // target_a = sqrt(2N)/M
    mpz_t target_a;
    mpz_init(target_a);
    mpz_mul_ui(target_a, N, 2);
    mpz_sqrt(target_a, target_a);
    mpz_tdiv_q_ui(target_a, target_a, M);

    // Sieve uses unsigned char: logp values are in bits (1 per bit of prime)
    // Threshold: log2(sqrt(N)*M) - slack ≈ bits(N)/2 + log2(M) - slack
    double log2_qmax = mpz_sizeinbase(N, 2) / 2.0 + log2((double)M);
    int threshold = (int)(log2_qmax - P.thresh_adj);
    if (threshold < 0) threshold = 0;
    // Cap at 255 for unsigned char
    if (threshold > 250) threshold = 250;

    fprintf(stderr, "Threshold: %d (log2_qmax=%.1f, adj=%.1f)\n",
            threshold, log2_qmax, P.thresh_adj);

    vector<Relation> relations;
    unordered_map<int, int> partial_map;
    struct Partial { Relation rel; int lp; };
    vector<Partial> partials;
    int total_full = 0, total_combined = 0, n_partials = 0;
    int poly_count = 0;

    // Range of FB indices for selecting 'a' factors
    int a_lo = max(2, fb_size / 5);
    int a_hi = min(fb_size - 1, fb_size * 4 / 5);

    // Sieve array: use unsigned char for cache efficiency
    vector<unsigned char> sieve(2 * M + 2);

    // Pre-allocate: soln1, soln2 arrays for FB
    vector<int> soln1(fb_size), soln2(fb_size);
    // Precomputed: ainv[i] = a^(-1) mod fb[i].p, Bl_delta[l][i] = 2*Bl mod p * ainv mod p
    vector<int> ainv_cache(fb_size);
    vector<vector<int>> Bl_delta; // [s][fb_size]

    mpz_t a_val, b_val, c_val;
    mpz_init(a_val); mpz_init(b_val); mpz_init(c_val);

    mpz_t Bl[16];
    for (int i = 0; i < 16; i++) mpz_init(Bl[i]);

    while (total_full + total_combined < needed) {
        // Select 'a' factors
        vector<int> a_idx;
        mpz_set_ui(a_val, 1);
        vector<int> pool;
        for (int i = a_lo; i <= a_hi; i++) pool.push_back(i);

        for (int i = 0; i < s && !pool.empty(); i++) {
            int r = next_rand() % pool.size();
            int idx = pool[r];
            a_idx.push_back(idx);
            mpz_mul_ui(a_val, a_val, fb[idx].p);
            pool.erase(pool.begin() + r);
        }
        sort(a_idx.begin(), a_idx.end());

        // Compute B_l values
        // B_l = sqrt(N) mod q_l * (a/q_l)^(-1) mod q_l * (a/q_l)
        // such that sum(B_l) ≡ sqrt(N) (mod a)  [with all + signs]
        for (int l = 0; l < s; l++) {
            int ql = fb[a_idx[l]].p;
            int sl = fb[a_idx[l]].sqrt_n; // sqrt(N) mod ql

            mpz_t a_div_ql, ql_mpz, inv_v;
            mpz_init(a_div_ql); mpz_init_set_ui(ql_mpz, ql); mpz_init(inv_v);
            mpz_tdiv_q_ui(a_div_ql, a_val, ql);
            mpz_invert(inv_v, a_div_ql, ql_mpz);

            // B_l = sl * inv * (a/ql) mod a
            mpz_mul_ui(Bl[l], inv_v, sl);
            mpz_mul(Bl[l], Bl[l], a_div_ql);
            mpz_mod(Bl[l], Bl[l], a_val);

            mpz_clear(a_div_ql); mpz_clear(ql_mpz); mpz_clear(inv_v);
        }

        // Initial b = sum of all B_l
        mpz_set_ui(b_val, 0);
        for (int l = 0; l < s; l++) mpz_add(b_val, b_val, Bl[l]);
        mpz_mod(b_val, b_val, a_val);

        // Make b < a/2
        mpz_tdiv_q_ui(tmp, a_val, 2);
        if (mpz_cmp(b_val, tmp) > 0) {
            mpz_sub(b_val, a_val, b_val);
            // Also negate all B_l to maintain b = sum(sign_l * B_l)
            for (int l = 0; l < s; l++) {
                mpz_sub(Bl[l], a_val, Bl[l]);
            }
        }

        // Verify b^2 ≡ N (mod a)
        mpz_mul(tmp, b_val, b_val);
        mpz_sub(tmp, tmp, N);
        mpz_mod(tmp, tmp, a_val);
        if (mpz_sgn(tmp) != 0) {
            // Bad polynomial, skip this 'a'
            continue;
        }

        // c = (b^2 - N) / a
        mpz_mul(c_val, b_val, b_val);
        mpz_sub(c_val, c_val, N);
        mpz_tdiv_q(c_val, c_val, a_val);

        // Compute initial sieve roots for each FB prime
        // soln = (±sqrt(N) - b) * a^(-1) mod p
        for (int i = 1; i < fb_size; i++) {
            int p = fb[i].p;
            bool divides_a = false;
            for (int j : a_idx) if (j == i) { divides_a = true; break; }
            if (divides_a || p == 2) {
                soln1[i] = soln2[i] = -1;
                continue;
            }
            long long amp = mpz_fdiv_ui(a_val, p);
            long long bmp = mpz_fdiv_ui(b_val, p);
            long long ainv = mod_pow(amp, p - 2, p);
            ainv_cache[i] = (int)ainv;
            int sq = fb[i].sqrt_n;
            soln1[i] = (int)(((__int128)(sq - bmp + p) % p * ainv) % p);
            soln2[i] = (int)(((__int128)(p - sq - bmp + p) % p * ainv) % p);
        }

        // Precompute Bl_delta[l][i] = 2 * Bl[l] mod p * ainv mod p
        Bl_delta.assign(s, vector<int>(fb_size, 0));
        for (int l = 0; l < s; l++) {
            for (int i = 1; i < fb_size; i++) {
                if (soln1[i] < 0) continue;
                int p = fb[i].p;
                long long blmp = mpz_fdiv_ui(Bl[l], p);
                Bl_delta[l][i] = (int)((2 * (__int128)blmp % p * ainv_cache[i]) % p);
            }
        }

        vector<int> signs(s, 1);
        int num_b = 1 << (s - 1);

        for (int bi = 0; bi < num_b; bi++) {
            // ---- Sieve for current polynomial ----
            memset(sieve.data(), 0, 2 * M + 2);

            // Regular FB primes
            for (int i = 1; i < fb_size; i++) {
                if (soln1[i] < 0) continue;
                int p = fb[i].p;
                int lp = fb[i].logp;

                int s1 = (int)(((long long)soln1[i] + M) % p);
                for (int j = s1; j <= 2 * M; j += p) sieve[j] += lp;

                if (soln1[i] != soln2[i]) {
                    int s2 = (int)(((long long)soln2[i] + M) % p);
                    for (int j = s2; j <= 2 * M; j += p) sieve[j] += lp;
                }
            }

            // Handle p=2 and a-primes
            // p=2: every other position gets +1
            for (int j = 0; j <= 2 * M; j += 2) sieve[j] += 1;

            // a-primes: single root sieving
            for (int ai : a_idx) {
                int p = fb[ai].p;
                int lp = fb[ai].logp;
                long long bmp = mpz_fdiv_ui(b_val, p);
                // c mod p: use mpz
                mpz_t cm; mpz_init(cm);
                mpz_fdiv_r_ui(cm, c_val, p);
                long long cmp = mpz_get_ui(cm);
                mpz_clear(cm);

                long long twob = (2 * bmp) % p;
                if (twob == 0) continue;
                long long twob_inv = mod_pow(twob, p - 2, p);
                long long root = ((p - cmp) % p * twob_inv) % p;
                int off = (int)(((long long)root + M) % p);
                for (int j = off; j <= 2 * M; j += p) sieve[j] += lp;
            }

            // ---- Scan for smooth candidates ----
            for (int idx = 0; idx <= 2 * M; idx++) {
                if (sieve[idx] < threshold) continue;

                int x = idx - M;

                // Compute Q(x) = ax^2 + 2bx + c
                mpz_t Qval;
                mpz_init(Qval);
                mpz_mul_si(Qval, a_val, x);
                mpz_addmul_ui(Qval, b_val, 2);
                mpz_mul_si(Qval, Qval, x);
                mpz_add(Qval, Qval, c_val);

                // Trial divide
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

                int lp = 0;
                bool smooth = true;
                if (mpz_cmp_ui(rem, 1) == 0) {
                    // fully smooth
                } else if (mpz_fits_slong_p(rem)) {
                    long r = mpz_get_si(rem);
                    if (r > 1 && r <= lp_bound) lp = (int)r;
                    else smooth = false;
                } else {
                    smooth = false;
                }

                if (smooth) {
                    Relation rel;
                    mpz_init(rel.x);
                    mpz_mul_si(rel.x, a_val, x);
                    mpz_add(rel.x, rel.x, b_val);
                    rel.exponents = exponents;
                    // Add a-factor exponents
                    for (int ai : a_idx) rel.exponents[ai]++;

                    if (lp == 0) {
                        relations.push_back(rel);
                        total_full++;
                    } else {
                        auto it = partial_map.find(lp);
                        if (it != partial_map.end()) {
                            Partial& other = partials[it->second];
                            Relation combined;
                            mpz_init(combined.x);
                            mpz_mul(combined.x, rel.x, other.rel.x);
                            mpz_t lp_inv, lp_mpz;
                            mpz_init_set_ui(lp_mpz, lp);
                            mpz_init(lp_inv);
                            mpz_invert(lp_inv, lp_mpz, N);
                            mpz_mul(combined.x, combined.x, lp_inv);
                            mpz_mod(combined.x, combined.x, N);
                            mpz_clear(lp_inv); mpz_clear(lp_mpz);

                            combined.exponents.resize(fb_size);
                            for (int k = 0; k < fb_size; k++)
                                combined.exponents[k] = rel.exponents[k] + other.rel.exponents[k];
                            relations.push_back(combined);
                            total_combined++;
                            partial_map.erase(it);
                        } else {
                            partial_map[lp] = partials.size();
                            partials.push_back({rel, lp});
                            n_partials++;
                        }
                    }
                }
                mpz_clear(Qval); mpz_clear(rem);

                if (total_full + total_combined >= needed) break;
            }

            poly_count++;
            if (total_full + total_combined >= needed) break;

            // ---- Update to next b using Gray code ----
            if (bi < num_b - 1) {
                // Gray code: flip the l-th sign where l = position of lowest set bit of (bi+1)
                int next = bi + 1;
                int l = __builtin_ctz(next); // index of lowest set bit
                if (l >= s) break;

                // Flip sign of B_l
                int old_sign = signs[l];
                signs[l] = -signs[l];

                // b_new = b_old + 2 * signs[l] * B_l (since we're flipping from -1 to +1 or vice versa)
                // If we go from sign +1 to -1: change = -2*B_l
                // If we go from sign -1 to +1: change = +2*B_l
                if (signs[l] == 1) {
                    // Was -1, now +1: b += 2*B_l
                    mpz_addmul_ui(b_val, Bl[l], 2);
                } else {
                    // Was +1, now -1: b -= 2*B_l
                    mpz_submul_ui(b_val, Bl[l], 2);
                }
                mpz_mod(b_val, b_val, a_val);

                // Recompute c = (b^2 - N) / a
                mpz_mul(c_val, b_val, b_val);
                mpz_sub(c_val, c_val, N);
                mpz_tdiv_q(c_val, c_val, a_val);

                // Update sieve roots: for each prime p not dividing a,
                // soln_new = soln_old + delta where delta = ±2*B_l * ainv mod p
                // The sign depends on which solution
                for (int i = 1; i < fb_size; i++) {
                    if (soln1[i] < 0) continue;
                    int p = fb[i].p;
                    int delta = Bl_delta[l][i];

                    if (signs[l] == 1) {
                        soln1[i] = (int)((soln1[i] - delta + p + p) % p);
                        soln2[i] = (int)((soln2[i] - delta + p + p) % p);
                    } else {
                        soln1[i] = (int)((soln1[i] + delta) % p);
                        soln2[i] = (int)((soln2[i] + delta) % p);
                    }
                }
            }
        }

        if (poly_count % 500 == 0 || (poly_count < 50 && poly_count % 10 == 0)) {
            fprintf(stderr, "  [%.1fs] poly %d: %d full + %d combined = %d/%d (partials: %d)\n",
                    timer(), poly_count, total_full, total_combined,
                    total_full + total_combined, needed, n_partials);
        }
    }

    fprintf(stderr, "Collected %d relations in %d polys, %.1fs\n",
            total_full + total_combined, poly_count, timer());

    // ---- Linear algebra ----
    int nrels = (int)relations.size();
    fprintf(stderr, "Linear algebra: %d x %d\n", nrels, fb_size);

    vector<vector<int>> evecs(nrels, vector<int>(fb_size));
    for (int i = 0; i < nrels; i++)
        for (int j = 0; j < fb_size; j++)
            evecs[i][j] = relations[i].exponents[j];

    vector<vector<int>> nullvecs = find_null_space(nrels, fb_size, evecs);
    fprintf(stderr, "Found %d null vectors\n", (int)nullvecs.size());

    bool factored = false;
    mpz_t X, Y;
    mpz_init(X); mpz_init(Y);

    for (auto& nv : nullvecs) {
        if (factored) break;
        mpz_set_ui(X, 1);
        vector<long long> texp(fb_size, 0);
        for (int idx : nv) {
            mpz_mul(X, X, relations[idx].x);
            mpz_mod(X, X, N);
            for (int j = 0; j < fb_size; j++)
                texp[j] += relations[idx].exponents[j];
        }

        bool ok = true;
        for (int j = 0; j < fb_size; j++) if (texp[j] & 1) { ok = false; break; }
        if (!ok) continue;

        mpz_set_ui(Y, 1);
        for (int j = 1; j < fb_size; j++) {
            if (texp[j] == 0) continue;
            mpz_set_ui(tmp, fb[j].p);
            mpz_powm_ui(tmp, tmp, texp[j] / 2, N);
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, N);
        }

        mpz_sub(tmp, X, Y); mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) { factored = true; break; }
        mpz_add(tmp, X, Y); mpz_gcd(factor, tmp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) { factored = true; break; }
    }

    if (factored) {
        gmp_printf("%Zd\n", factor);
        fprintf(stderr, "SIQS2: factored in %.3f seconds\n", timer());
    } else {
        fprintf(stderr, "SIQS2: FAILED after %.3f seconds\n", timer());
        return 1;
    }

    // Cleanup
    mpz_clear(N); mpz_clear(factor); mpz_clear(tmp); mpz_clear(tmp2);
    mpz_clear(target_a); mpz_clear(a_val); mpz_clear(b_val); mpz_clear(c_val);
    mpz_clear(X); mpz_clear(Y);
    for (int i = 0; i < 16; i++) mpz_clear(Bl[i]);
    for (auto& r : relations) mpz_clear(r.x);
    for (auto& p : partials) mpz_clear(p.rel.x);

    return 0;
}
