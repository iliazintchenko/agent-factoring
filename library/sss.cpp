/*
 * Smooth Subsum Search (SSS) for integer factoring.
 *
 * Based on Hittmeir (2023): "Smooth Subsum Search: A Heuristic for Practical
 * Integer Factorization". Int. J. Found. Comp. Sc.
 *
 * Key idea: Instead of sieving over a polynomial, use CRT to construct
 * x values where pol(x) = (x+b)^2 - N is guaranteed divisible by several
 * small factor base primes. The remaining cofactor is smaller, so it's more
 * likely to be smooth over the rest of the factor base.
 *
 * This replaces sieving with structured candidate generation + batch
 * smoothness testing (product/remainder trees).
 *
 * Usage: ./sss <number>
 * Single-threaded, seed=42.
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
#include <unordered_set>
#include <gmp.h>

/* ---- helpers ---- */

static long power_mod_l(long base, long exp, long mod) {
    long r=1; base=((base%mod)+mod)%mod;
    while(exp>0){if(exp&1)r=(__int128)r*base%mod;base=(__int128)base*base%mod;exp>>=1;}
    return r;
}
static int legendre_sym(long a, long p) {
    long r = power_mod_l(a,(p-1)/2,p); return r==p-1?-1:(int)r;
}
static long tonelli(long n, long p) {
    if(p==2) return n&1; n=((n%p)+p)%p; if(!n)return 0;
    long Q=p-1,S=0; while(Q%2==0){Q/=2;S++;}
    if(S==1) return power_mod_l(n,(p+1)/4,p);
    long z=2; while(power_mod_l(z,(p-1)/2,p)!=p-1)z++;
    long M2=S,c=power_mod_l(z,Q,p),t=power_mod_l(n,Q,p),R=power_mod_l(n,(Q+1)/2,p);
    while(t!=1){long i=0,v=t;while(v!=1){v=(__int128)v*v%p;i++;}
    long b2=c;for(long j=0;j<M2-i-1;j++)b2=(__int128)b2*b2%p;
    M2=i;c=(__int128)b2*b2%p;t=(__int128)t*c%p;R=(__int128)R*b2%p;}
    return R;
}
static std::vector<long> sieve_primes(long lim) {
    std::vector<bool> s(lim+1,true); s[0]=s[1]=false;
    for(long i=2;i*i<=lim;i++) if(s[i]) for(long j=i*i;j<=lim;j+=i) s[j]=false;
    std::vector<long> p; for(long i=2;i<=lim;i++) if(s[i]) p.push_back(i); return p;
}

struct BitRow {
    std::vector<uint64_t> bits,hist; int wb,wh;
    void init(int nb,int nh){wb=(nb+63)/64;wh=(nh+63)/64;bits.assign(wb,0);hist.assign(wh,0);}
    void sbit(int i){bits[i/64]|=1ULL<<(i%64);}
    int gbit(int i)const{return(bits[i/64]>>(i%64))&1;}
    void shist(int i){hist[i/64]|=1ULL<<(i%64);}
    int ghist(int i)const{return(hist[i/64]>>(i%64))&1;}
    void xr(const BitRow&o){for(int i=0;i<wb;i++)bits[i]^=o.bits[i];for(int i=0;i<wh;i++)hist[i]^=o.hist[i];}
    bool zb()const{for(int i=0;i<wb;i++)if(bits[i])return false;return true;}
};

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <number>\n", argv[0]); return 1; }

    mpz_t N, tmp, g;
    mpz_init(N); mpz_init(tmp); mpz_init(g);
    if (mpz_set_str(N, argv[1], 10) != 0) { fprintf(stderr, "Bad number\n"); return 1; }

    struct timespec tstart;
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    size_t digits = mpz_sizeinbase(N, 10);

    /* Parameter selection (matching Hittmeir's tuning) */
    int m; /* number of factor base primes */
    if (digits <= 18) m = 120;
    else if (digits <= 25) m = 300;
    else if (digits <= 34) m = 400;
    else if (digits <= 36) m = 600;
    else if (digits <= 38) m = 800;
    else if (digits <= 40) m = 1000;
    else if (digits <= 42) m = 1200;
    else if (digits <= 44) m = 1400;
    else if (digits <= 48) m = 2000;
    else if (digits <= 52) m = 2400;
    else if (digits <= 56) m = 4000;
    else if (digits <= 60) m = 8000;
    else if (digits <= 66) m = 12000;
    else if (digits <= 74) m = 20000;
    else if (digits <= 80) m = 60000;
    else if (digits <= 88) m = 100000;
    else if (digits <= 94) m = 120000;
    else m = 200000;

    /* Polynomial: pol(j) = (j + b)^2 - N, where b = ceil(sqrt(N)) */
    mpz_t b_val;
    mpz_init(b_val);
    mpz_sqrt(b_val, N);
    mpz_mul(tmp, b_val, b_val);
    if (mpz_cmp(tmp, N) < 0) mpz_add_ui(b_val, b_val, 1);

    /* Build factor base: first m primes where N is QR */
    std::vector<long> all_primes = sieve_primes(m * 20); /* enough to get m QR primes */
    std::vector<long> fbase;
    for (long p : all_primes) {
        mpz_mod_ui(tmp, N, p);
        long nm = mpz_get_ui(tmp);
        if (p == 2 || legendre_sym(nm, p) == 1) {
            fbase.push_back(p);
            if ((int)fbase.size() >= m) break;
        }
    }

    int fb_size = fbase.size();
    int plist_size = fb_size / 5; /* small primes for CRT */
    if (plist_size < 10) plist_size = 10;
    if (plist_size > fb_size) plist_size = fb_size;

    long last_prime = fbase.back();
    long partial_bound = 128 * last_prime;
    int target = fb_size + 10;

    fprintf(stderr, "SSS: %zu dig, fb=%d, plist=%d, target=%d\n",
            digits, fb_size, plist_size, target);

    /* Precompute square roots of N mod each fb prime */
    std::vector<long> sqrt_r0(fb_size), sqrt_r1(fb_size);
    for (int i = 0; i < fb_size; i++) {
        long p = fbase[i];
        mpz_mod_ui(tmp, N, p);
        long nm = mpz_get_ui(tmp);
        long r = tonelli(nm, p);
        /* Roots of (j+b)^2 ≡ N mod p → j+b ≡ ±sqrt(N) mod p → j ≡ ±sqrt(N) - b mod p */
        mpz_mod_ui(tmp, b_val, p);
        long bm = mpz_get_ui(tmp);
        sqrt_r0[i] = ((r - bm) % p + p) % p;
        sqrt_r1[i] = ((p - r - bm) % p + p) % p;
        if (sqrt_r0[i] > sqrt_r1[i]) std::swap(sqrt_r0[i], sqrt_r1[i]);
    }

    /* CRT coefficients for plist */
    /* nval = product of plist primes */
    mpz_t nval;
    mpz_init(nval);
    mpz_set_ui(nval, 1);
    for (int i = 0; i < plist_size; i++)
        mpz_mul_ui(nval, nval, fbase[i]);

    /* coeffs[i] = (nval/p_i) * (nval/p_i)^{-1} mod p_i */
    std::vector<mpz_t> coeffs(plist_size);
    for (int i = 0; i < plist_size; i++) {
        mpz_init(coeffs[i]);
        long p = fbase[i];
        mpz_t nval_over_p;
        mpz_init(nval_over_p);
        mpz_divexact_ui(nval_over_p, nval, p);
        mpz_mod_ui(tmp, nval_over_p, p);
        long inv = power_mod_l(mpz_get_ui(tmp), p - 2, p);
        mpz_mul_ui(coeffs[i], nval_over_p, inv);
        mpz_clear(nval_over_p);
    }

    /* difflist[i] = coeffs[i] * (sqrt_r1[i] - sqrt_r0[i]) for i in plist */
    std::vector<mpz_t> difflist(plist_size);
    for (int i = 0; i < plist_size; i++) {
        mpz_init(difflist[i]);
        mpz_mul_ui(difflist[i], coeffs[i], sqrt_r1[i] - sqrt_r0[i]);
    }

    /* Precompute mval = product of fbase^k for batch smoothness */
    mpz_t mval;
    mpz_init(mval);
    mpz_set_ui(mval, 1);
    int nbits = mpz_sizeinbase(N, 2);
    for (int i = 0; i < fb_size; i++) {
        long p = fbase[i];
        int k = (int)(15.0 * log(2) / log((double)p));
        if (k < 1) k = 1;
        mpz_t pk;
        mpz_init(pk);
        mpz_ui_pow_ui(pk, p, k);
        mpz_mul(mval, mval, pk);
        mpz_clear(pk);
    }

    /* RNG seeded with 42 */
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    /* Main SSS loop */
    std::unordered_set<long long> smooth_x_set;
    std::vector<long long> smooth_x_list;
    std::unordered_map<long, std::vector<long long>> lp_map; /* large_prime -> list of x values */
    int usable_count = 0;
    int search_count = 0;
    int length = 6; /* number of CRT primes per search */

    while (usable_count < target) {
        struct timespec tnow;
        clock_gettime(CLOCK_MONOTONIC, &tnow);
        double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
        if (elapsed > 280) {
            fprintf(stderr, "Timeout after %.1fs\n", elapsed);
            break;
        }

        /* Pick random subset of 'length' primes from plist */
        std::vector<int> chosen;
        {
            std::vector<int> pool(plist_size);
            for (int i = 0; i < plist_size; i++) pool[i] = i;
            for (int i = 0; i < length && !pool.empty(); i++) {
                int idx = gmp_urandomm_ui(rstate, pool.size());
                chosen.push_back(pool[idx]);
                pool.erase(pool.begin() + idx);
            }
        }
        if ((int)chosen.size() < length) continue;

        /* M = product of chosen primes */
        mpz_t M;
        mpz_init(M);
        mpz_set_ui(M, 1);
        for (int i : chosen) mpz_mul_ui(M, M, fbase[i]);

        /* xval = sum of coeffs[i] * sqrt_r0[i] for i in chosen, mod M */
        /* This gives x where pol(x) ≡ 0 mod each chosen prime (root 0) */
        mpz_t xval_mpz;
        mpz_init(xval_mpz);
        mpz_set_ui(xval_mpz, 0);
        for (int i : chosen) {
            mpz_mul_ui(tmp, coeffs[i], sqrt_r0[i]);
            mpz_add(xval_mpz, xval_mpz, tmp);
        }
        mpz_mod(xval_mpz, xval_mpz, M);

        /* For each prime i in chosen, try switching to root 1 */
        for (int sel = 0; sel < (int)chosen.size(); sel++) {
            int switched_prime = chosen[sel];

            /* Switch from root 0 to root 1 for this prime */
            mpz_add(xval_mpz, xval_mpz, difflist[switched_prime]);
            mpz_mod(xval_mpz, xval_mpz, M);

            /* For each sub-product M/p_j (j != switched_prime), generate candidates */
            for (int j_idx = 0; j_idx <= (int)chosen.size(); j_idx++) {
                /* j_idx = 0 means we use M itself (all primes)
                 * j_idx > 0 means we drop chosen[j_idx-1] */
                if (j_idx > 0 && chosen[j_idx - 1] == switched_prime) continue;

                long dropped_prime = 1;
                if (j_idx > 0) dropped_prime = fbase[chosen[j_idx - 1]];

                mpz_t m_sub;
                mpz_init(m_sub);
                if (j_idx == 0) mpz_set(m_sub, M);
                else mpz_divexact_ui(m_sub, M, dropped_prime);

                /* For remaining fb primes (not in plist), compute sieve positions
                 * x ≡ xval mod m_sub, and also x ≡ sqrt_r0[v] or sqrt_r1[v] mod fbase[v]
                 * Use CRT to combine: offset = (root - xval) * m_sub^{-1} mod fbase[v]
                 * Then k = xval + offset * m_sub */

                /* Compute inverse of m_sub mod each remaining fb prime */
                /* Count collisions: for each k value, how many fb primes p have a root at that position */
                std::unordered_map<long long, int> k_count;

                for (int v = plist_size; v < fb_size; v++) {
                    long p = fbase[v];
                    mpz_mod_ui(tmp, m_sub, p);
                    long ms_mod_p = mpz_get_ui(tmp);
                    if (ms_mod_p == 0) continue;
                    long inv = power_mod_l(ms_mod_p, p - 2, p);

                    mpz_mod_ui(tmp, xval_mpz, p);
                    long xm = mpz_get_ui(tmp);

                    long r0 = sqrt_r0[v], r1 = sqrt_r1[v];
                    long k0 = (__int128)(r0 - xm + p) % p * inv % p;
                    long k1 = (__int128)(r1 - xm + p) % p * inv % p;

                    /* Both positive and negative offsets */
                    k_count[k0]++;
                    k_count[k0 - p]++;
                    k_count[k1]++;
                    k_count[k1 - p]++;
                }

                /* For positions with multiple hits (>2), test for smoothness */
                for (auto &[k, cnt] : k_count) {
                    if (cnt <= 2) continue;

                    /* Compute x = xval + k * m_sub */
                    mpz_t x_val, pol_val, pol_abs;
                    mpz_init(x_val); mpz_init(pol_val); mpz_init(pol_abs);

                    mpz_set_si(tmp, k);
                    mpz_mul(x_val, tmp, m_sub);
                    mpz_add(x_val, x_val, xval_mpz);

                    /* pol(x) = (x + b)^2 - N */
                    mpz_add(pol_val, x_val, b_val);
                    mpz_mul(pol_val, pol_val, pol_val);
                    mpz_sub(pol_val, pol_val, N);

                    /* pol(x) should be divisible by m_sub (by construction) */
                    if (mpz_sgn(pol_val) == 0) {
                        /* Perfect square - check for direct factor */
                        mpz_add(tmp, x_val, b_val);
                        mpz_gcd(g, tmp, N);
                        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                            struct timespec tnow2;
                            clock_gettime(CLOCK_MONOTONIC, &tnow2);
                            double el = (tnow2.tv_sec - tstart.tv_sec) + (tnow2.tv_nsec - tstart.tv_nsec) / 1e9;
                            gmp_printf("FACTOR: %Zd\n", g);
                            mpz_t cof; mpz_init(cof); mpz_divexact(cof, N, g);
                            gmp_printf("COFACTOR: %Zd\n", cof);
                            fprintf(stderr, "SSS perfect square in %.3fs\n", el);
                            return 0;
                        }
                        mpz_clear(x_val); mpz_clear(pol_val); mpz_clear(pol_abs);
                        continue;
                    }

                    mpz_abs(pol_abs, pol_val);

                    /* Divide by m_sub to get the reduced value */
                    if (mpz_divisible_p(pol_abs, m_sub)) {
                        mpz_divexact(pol_abs, pol_abs, m_sub);
                    }

                    /* Batch smoothness test: gcd with mval */
                    mpz_mod(tmp, mval, pol_abs);
                    mpz_gcd(g, tmp, pol_abs);

                    /* Iteratively extract smooth part */
                    mpz_t remainder;
                    mpz_init_set(remainder, pol_abs);
                    while (mpz_cmp_ui(g, 1) > 0) {
                        while (mpz_divisible_p(remainder, g))
                            mpz_divexact(remainder, remainder, g);
                        if (mpz_cmp_ui(remainder, 1) == 0) break;
                        mpz_mod(tmp, mval, remainder);
                        mpz_gcd(g, tmp, remainder);
                    }

                    bool is_smooth = (mpz_cmp_ui(remainder, 1) == 0);
                    bool lp_ok = false;
                    long lp_val = 0;

                    if (!is_smooth && mpz_fits_ulong_p(remainder)) {
                        lp_val = mpz_get_ui(remainder);
                        if (lp_val > 1 && lp_val <= (unsigned long)partial_bound) {
                            if (mpz_probab_prime_p(remainder, 3))
                                lp_ok = true;
                        }
                    }
                    mpz_clear(remainder);

                    if (is_smooth) {
                        /* Check for duplicate */
                        long long xll = 0;
                        if (mpz_fits_slong_p(x_val)) xll = mpz_get_si(x_val);
                        else { /* hash it */ mpz_mod_ui(tmp, x_val, 1000000007); xll = mpz_get_si(tmp); }

                        if (smooth_x_set.find(xll) == smooth_x_set.end()) {
                            smooth_x_set.insert(xll);
                            smooth_x_list.push_back(xll);
                            usable_count++;
                        }
                    } else if (lp_ok) {
                        long long xll = 0;
                        if (mpz_fits_slong_p(x_val)) xll = mpz_get_si(x_val);
                        else { mpz_mod_ui(tmp, x_val, 1000000007); xll = mpz_get_si(tmp); }

                        auto it = lp_map.find(lp_val);
                        if (it != lp_map.end()) {
                            if (it->second.size() >= 2) {
                                /* Can combine */
                                usable_count++;
                            }
                            it->second.push_back(xll);
                        } else {
                            lp_map[lp_val] = {xll};
                        }
                    }

                    mpz_clear(x_val); mpz_clear(pol_val); mpz_clear(pol_abs);
                }

                mpz_clear(m_sub);
            }
        }

        mpz_clear(M); mpz_clear(xval_mpz);
        search_count++;

        if (search_count % 50 == 0) {
            int lp_count = 0;
            for (auto &[k, v] : lp_map)
                if (v.size() >= 2) lp_count += (v.size() - 1) / 2;
            clock_gettime(CLOCK_MONOTONIC, &tnow);
            elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
            fprintf(stderr, "  search %d: %d smooth + %d LP = %d/%d (%.1fs)\n",
                    search_count, (int)smooth_x_list.size(), lp_count,
                    (int)smooth_x_list.size() + lp_count, target, elapsed);
        }
    }

    fprintf(stderr, "SSS relation finding: %d smooth, %d searches\n",
            (int)smooth_x_list.size(), search_count);

    /* Now we need to do trial division on the smooth values to get exponent vectors,
     * then linear algebra. */

    /* Collect all relations with full exponent vectors */
    struct Relation {
        mpz_t x_plus_b; /* (x + b) */
        std::vector<int> exps;
    };
    std::vector<Relation> rels;

    for (long long xll : smooth_x_list) {
        /* Reconstruct pol(x) = (x+b)^2 - N */
        mpz_t xv, polv;
        mpz_init_set_si(xv, xll);
        mpz_init(polv);
        mpz_add(polv, xv, b_val);
        mpz_mul(polv, polv, polv);
        mpz_sub(polv, polv, N);

        std::vector<int> exps(fb_size, 0);
        int neg = (mpz_sgn(polv) < 0);
        if (neg) mpz_neg(polv, polv);
        /* Note: we don't have -1 in our fb. We'll handle sign separately. */

        for (int j = 0; j < fb_size; j++) {
            long p = fbase[j];
            while (mpz_divisible_ui_p(polv, p)) {
                mpz_divexact_ui(polv, polv, p);
                exps[j]++;
            }
        }

        if (mpz_cmp_ui(polv, 1) == 0) {
            Relation rel;
            mpz_init(rel.x_plus_b);
            mpz_set_si(tmp, xll);
            mpz_add(rel.x_plus_b, tmp, b_val);
            /* Add sign as extra column */
            if (neg) exps.push_back(1);
            else exps.push_back(0);
            rel.exps = exps;
            rels.push_back(rel);
        }

        mpz_clear(xv); mpz_clear(polv);
    }

    /* Also add combined LP relations */
    for (auto &[lp, xvals] : lp_map) {
        if (xvals.size() < 2) continue;
        for (size_t i = 0; i + 1 < xvals.size(); i += 2) {
            mpz_t x1, x2, polv1, polv2;
            mpz_init_set_si(x1, xvals[i]);
            mpz_init_set_si(x2, xvals[i + 1]);
            mpz_init(polv1); mpz_init(polv2);

            mpz_add(polv1, x1, b_val);
            mpz_mul(polv1, polv1, polv1);
            mpz_sub(polv1, polv1, N);

            mpz_add(polv2, x2, b_val);
            mpz_mul(polv2, polv2, polv2);
            mpz_sub(polv2, polv2, N);

            std::vector<int> exps(fb_size, 0);
            int neg = 0;
            if (mpz_sgn(polv1) < 0) { neg++; mpz_neg(polv1, polv1); }
            if (mpz_sgn(polv2) < 0) { neg++; mpz_neg(polv2, polv2); }

            /* Trial divide both */
            for (int j = 0; j < fb_size; j++) {
                long p = fbase[j];
                while (mpz_divisible_ui_p(polv1, p)) { mpz_divexact_ui(polv1, polv1, p); exps[j]++; }
                while (mpz_divisible_ui_p(polv2, p)) { mpz_divexact_ui(polv2, polv2, p); exps[j]++; }
            }

            /* Both remainders should be the LP value, so product is LP^2 (even exponent) */
            exps.push_back(neg & 1); /* sign */

            Relation rel;
            mpz_init(rel.x_plus_b);
            mpz_set_si(tmp, xvals[i]);
            mpz_add(rel.x_plus_b, tmp, b_val);
            mpz_set_si(tmp, xvals[i + 1]);
            mpz_add(tmp, tmp, b_val);
            mpz_mul(rel.x_plus_b, rel.x_plus_b, tmp);
            /* Also need LP^{-1} mod N for the x side */
            mpz_set_ui(tmp, lp);
            mpz_t lp_inv;
            mpz_init(lp_inv);
            if (mpz_invert(lp_inv, tmp, N)) {
                mpz_mul(rel.x_plus_b, rel.x_plus_b, lp_inv);
            }
            mpz_mod(rel.x_plus_b, rel.x_plus_b, N);
            mpz_clear(lp_inv);

            rel.exps = exps;
            rels.push_back(rel);

            mpz_clear(x1); mpz_clear(x2); mpz_clear(polv1); mpz_clear(polv2);
        }
    }

    int nrels = rels.size();
    int ncols = fb_size + 1; /* +1 for sign */
    fprintf(stderr, "SSS: %d relations, %d columns\n", nrels, ncols);

    if (nrels <= ncols) {
        fprintf(stderr, "Not enough relations.\n");
        printf("FAIL\n");
        goto cleanup;
    }

    {
        /* Linear algebra */
        if (nrels > ncols + 200) nrels = ncols + 200;
        fprintf(stderr, "LinAlg: %dx%d\n", nrels, ncols);

        std::vector<BitRow> mat(nrels);
        for (int i = 0; i < nrels; i++) {
            mat[i].init(ncols, nrels); mat[i].shist(i);
            for (int j = 0; j < ncols; j++)
                if (rels[i].exps[j] & 1) mat[i].sbit(j);
        }

        int rank = 0;
        for (int col = 0; col < ncols && rank < nrels; col++) {
            int pivot = -1;
            for (int row = rank; row < nrels; row++)
                if (mat[row].gbit(col)) { pivot = row; break; }
            if (pivot == -1) continue;
            if (pivot != rank) std::swap(mat[pivot], mat[rank]);
            for (int row = 0; row < nrels; row++)
                if (row != rank && mat[row].gbit(col))
                    mat[row].xr(mat[rank]);
            rank++;
        }
        fprintf(stderr, "Rank %d, null %d\n", rank, nrels - rank);

        int found = 0;
        mpz_t px, py;
        mpz_init(px); mpz_init(py);

        for (int row = 0; row < nrels && !found; row++) {
            if (!mat[row].zb()) continue;

            mpz_set_ui(px, 1);
            std::vector<long> se(ncols, 0);

            for (int i = 0; i < nrels; i++) {
                if (!mat[row].ghist(i)) continue;
                mpz_mul(px, px, rels[i].x_plus_b);
                mpz_mod(px, px, N);
                for (int j = 0; j < ncols; j++)
                    se[j] += rels[i].exps[j];
            }

            bool ok = true;
            for (int j = 0; j < ncols; j++) if (se[j] & 1) { ok = false; break; }
            if (!ok) continue;

            mpz_set_ui(py, 1);
            for (int j = 0; j < fb_size; j++) {
                if (se[j] == 0) continue;
                mpz_set_ui(tmp, fbase[j]);
                mpz_powm_ui(tmp, tmp, se[j] / 2, N);
                mpz_mul(py, py, tmp);
                mpz_mod(py, py, N);
            }

            for (int sign = 0; sign < 2 && !found; sign++) {
                if (sign == 0) mpz_sub(tmp, px, py);
                else mpz_add(tmp, px, py);
                mpz_gcd(g, tmp, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    struct timespec tnow;
                    clock_gettime(CLOCK_MONOTONIC, &tnow);
                    double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                    gmp_printf("FACTOR: %Zd\n", g);
                    mpz_t cof; mpz_init(cof); mpz_divexact(cof, N, g);
                    gmp_printf("COFACTOR: %Zd\n", cof);
                    fprintf(stderr, "SSS factored in %.3fs (%d searches)\n", elapsed, search_count);
                    mpz_clear(cof);
                    found = 1;
                }
            }
        }
        mpz_clear(px); mpz_clear(py);
        if (!found) { fprintf(stderr, "No factor from deps.\n"); printf("FAIL\n"); }
    }

cleanup:
    for (auto &r : rels) mpz_clear(r.x_plus_b);
    mpz_clear(N); mpz_clear(tmp); mpz_clear(g); mpz_clear(b_val);
    mpz_clear(nval); mpz_clear(mval);
    for (int i = 0; i < plist_size; i++) { mpz_clear(coeffs[i]); mpz_clear(difflist[i]); }
    gmp_randclear(rstate);
    return 0;
}
