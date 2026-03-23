/*
 * CFRAC - Continued Fraction Factoring Algorithm.
 *
 * Uses the continued fraction expansion of sqrt(N) to generate smooth values.
 * Each convergent p_k/q_k of sqrt(N) satisfies: p_k^2 - N*q_k^2 = (-1)^{k+1} * A_k
 * where |A_k| < 2*sqrt(N). These A_k values are much smaller than N,
 * giving good smoothness probability.
 *
 * The relation is: p_k^2 ≡ (-1)^{k+1} * A_k (mod N)
 * where A_k needs to factor over the factor base.
 *
 * Key advantage over QS: no polynomial generation needed, the CF expansion
 * naturally produces good candidates. Disadvantage: only one "polynomial"
 * (the CF expansion), so no parallelism in relation generation.
 *
 * Usage: ./cfrac <number>
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
#include <gmp.h>

static long power_mod_l(long base, long exp, long mod) {
    long r=1; base=((base%mod)+mod)%mod;
    while(exp>0){if(exp&1)r=(__int128)r*base%mod;base=(__int128)base*base%mod;exp>>=1;}
    return r;
}
static int legendre_sym(long a, long p) {
    long r = power_mod_l(a,(p-1)/2,p); return r==p-1?-1:(int)r;
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

struct Relation {
    mpz_t p_val; /* convergent numerator mod N */
    std::vector<int> exps;
    long large_prime;
};

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <number>\n", argv[0]); return 1; }

    mpz_t N, tmp, g;
    mpz_init(N); mpz_init(tmp); mpz_init(g);
    if (mpz_set_str(N, argv[1], 10) != 0) { fprintf(stderr, "Bad number\n"); return 1; }

    struct timespec tstart;
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    size_t digits = mpz_sizeinbase(N, 10);
    double ln_n = log(2.0) * mpz_sizeinbase(N, 2);
    double ln_ln_n = log(ln_n);
    double opt = exp(0.5 * sqrt(ln_n * ln_ln_n));

    long B = (long)(opt * 0.9);
    if (B < 2000) B = 2000;
    if (B > 8000000) B = 8000000;
    long lp_mult = 80;

    std::vector<long> all_primes = sieve_primes(B);
    std::vector<long> fb; fb.push_back(-1);
    {mpz_t nm; mpz_init(nm);
    for (long p : all_primes) {
        mpz_mod_ui(nm, N, p);
        long n2 = mpz_get_ui(nm);
        if (p == 2 || legendre_sym(n2, p) == 1)
            fb.push_back(p);
    }
    mpz_clear(nm);}

    int fb_size = fb.size();
    int target = fb_size + 30;
    long lp_bound = B * lp_mult;

    fprintf(stderr, "CFRAC: %zu dig, B=%ld, fb=%d, target=%d\n", digits, B, fb_size, target);

    /* Continued fraction expansion of sqrt(N):
     * a_0 = floor(sqrt(N))
     * Then iteratively:
     *   m_{k+1} = a_k * d_k - m_k
     *   d_{k+1} = (N - m_{k+1}^2) / d_k
     *   a_{k+1} = floor((a_0 + m_{k+1}) / d_{k+1})
     *
     * Convergents p_k/q_k:
     *   p_{-1} = 1, p_0 = a_0
     *   q_{-1} = 0, q_0 = 1
     *   p_k = a_k * p_{k-1} + p_{k-2}
     *   q_k = a_k * q_{k-1} + q_{k-2}
     *
     * The key identity: p_k^2 - N * q_k^2 = (-1)^{k+1} * d_{k+1}
     * So p_k^2 ≡ (-1)^{k+1} * d_{k+1} (mod N)
     * We need d_{k+1} (or its absolute value) to factor over our base.
     */

    /* CF state variables */
    mpz_t a0, m_cf, d_cf, a_cf;
    mpz_t p_prev, p_curr, q_prev, q_curr;
    mpz_init(a0); mpz_init(m_cf); mpz_init(d_cf); mpz_init(a_cf);
    mpz_init(p_prev); mpz_init(p_curr); mpz_init(q_prev); mpz_init(q_curr);

    mpz_sqrt(a0, N);

    /* Check if N is a perfect square */
    mpz_mul(tmp, a0, a0);
    if (mpz_cmp(tmp, N) == 0) {
        gmp_printf("FACTOR: %Zd\n", a0);
        fprintf(stderr, "Perfect square!\n");
        return 0;
    }

    /* Initialize CF */
    mpz_set_ui(m_cf, 0);
    mpz_set_ui(d_cf, 1);
    mpz_set(a_cf, a0);

    mpz_set_ui(p_prev, 1); /* p_{-1} */
    mpz_set(p_curr, a0);   /* p_0 */
    mpz_set_ui(q_prev, 0); /* q_{-1} */
    mpz_set_ui(q_curr, 1); /* q_0 */

    std::vector<Relation> rels;
    std::unordered_map<long, int> lp_map;
    int usable_count = 0;
    int step = 0;

    while (usable_count < target) {
        struct timespec tnow;
        clock_gettime(CLOCK_MONOTONIC, &tnow);
        double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
        if (elapsed > 280) {
            fprintf(stderr, "Timeout at step %d, %d usable\n", step, usable_count);
            break;
        }

        /* Advance CF by one step */
        /* m_{k+1} = a_k * d_k - m_k */
        mpz_mul(tmp, a_cf, d_cf);
        mpz_sub(tmp, tmp, m_cf);
        mpz_set(m_cf, tmp);

        /* d_{k+1} = (N - m_{k+1}^2) / d_k */
        mpz_mul(tmp, m_cf, m_cf);
        mpz_sub(tmp, N, tmp);
        mpz_divexact(d_cf, tmp, d_cf);

        /* a_{k+1} = floor((a0 + m_{k+1}) / d_{k+1}) */
        mpz_add(tmp, a0, m_cf);
        mpz_tdiv_q(a_cf, tmp, d_cf);

        /* Update convergents: p_k = a_k * p_{k-1} + p_{k-2} */
        mpz_mul(tmp, a_cf, p_curr);
        mpz_add(tmp, tmp, p_prev);
        mpz_set(p_prev, p_curr);
        mpz_set(p_curr, tmp);

        mpz_mul(tmp, a_cf, q_curr);
        mpz_add(tmp, tmp, q_prev);
        mpz_set(q_prev, q_curr);
        mpz_set(q_curr, tmp);

        step++;

        /* The relation: p_curr^2 ≡ (-1)^(step+1) * d_cf (mod N)
         * More precisely: p_{k-1}^2 - N * q_{k-1}^2 = (-1)^k * d_k
         * After step k: p_curr = p_k, d_cf = d_{k+1}
         * The identity using the NEW d_cf: p_prev^2 - N*q_prev^2 = (-1)^step * d_cf
         * Wait, let me get this right.
         *
         * Actually, the CF identity is:
         * p_k^2 - N * q_k^2 = (-1)^{k+1} * A_{k+1}
         * where A_{k+1} is the denominator d at step k+1.
         *
         * After our step k (step variable = k), d_cf = d_{k+1}, and p_curr = p_k.
         * So: p_k^2 ≡ (-1)^{k+1} * d_{k+1} mod N
         * = p_curr^2 ≡ (-1)^{step+1} * d_cf mod N
         */

        /* We need to factor d_cf (= |A|) over the factor base */
        /* d_cf should be positive and < 2*sqrt(N) */
        if (mpz_sgn(d_cf) <= 0) continue; /* shouldn't happen */

        mpz_t d_copy;
        mpz_init_set(d_copy, d_cf);

        std::vector<int> exps(fb_size, 0);
        /* Identity: p_prev^2 ≡ (-1)^step * d_cf (mod N)
         * We use p_prev (the previous convergent), not p_curr. */
        if (step & 1) exps[0] = 1; /* sign is (-1)^step: odd step = negative */

        for (int j = 1; j < fb_size; j++) {
            long p = fb[j];
            while (mpz_divisible_ui_p(d_copy, p)) {
                mpz_divexact_ui(d_copy, d_copy, p);
                exps[j]++;
            }
        }

        bool smooth = (mpz_cmp_ui(d_copy, 1) == 0);
        long rem = 0;
        bool lp_ok = false;

        if (!smooth && mpz_fits_ulong_p(d_copy)) {
            rem = mpz_get_ui(d_copy);
            if (rem > 1 && rem <= (unsigned long)lp_bound) {
                if (mpz_probab_prime_p(d_copy, 3))
                    lp_ok = true;
            }
        }

        if (smooth) {
            Relation rel;
            mpz_init(rel.p_val);
            mpz_mod(rel.p_val, p_prev, N); /* p_{k-1} mod N (the previous convergent) */
            rel.exps = exps;
            rel.large_prime = 0;
            rels.push_back(rel);
            usable_count++;
        } else if (lp_ok) {
            auto it = lp_map.find(rem);
            if (it != lp_map.end()) {
                int oi = it->second;
                Relation &orel = rels[oi];
                Relation rel;
                mpz_init(rel.p_val);
                mpz_mul(rel.p_val, p_prev, orel.p_val);
                mpz_mod(rel.p_val, rel.p_val, N);
                rel.exps.resize(fb_size);
                for (int j = 0; j < fb_size; j++)
                    rel.exps[j] = exps[j] + orel.exps[j];
                rel.large_prime = -rem;
                rels.push_back(rel);
                usable_count++;
            } else {
                Relation rel;
                mpz_init(rel.p_val);
                mpz_mod(rel.p_val, p_prev, N);
                rel.exps = exps;
                rel.large_prime = rem;
                rels.push_back(rel);
                lp_map[rem] = rels.size() - 1;
            }
        }

        mpz_clear(d_copy);

        if (step % 50000 == 0) {
            fprintf(stderr, "  step %d: %d usable/%zu total rels (%.1fs)\n",
                    step, usable_count, rels.size(), elapsed);
        }
    }

    fprintf(stderr, "%d usable from %d steps\n", usable_count, step);

    mpz_clear(a0); mpz_clear(m_cf); mpz_clear(d_cf); mpz_clear(a_cf);
    mpz_clear(p_prev); mpz_clear(p_curr); mpz_clear(q_prev); mpz_clear(q_curr);

    if (usable_count <= fb_size) {
        fprintf(stderr, "Not enough.\n");
        printf("FAIL\n");
        goto cleanup;
    }

    {
        std::vector<int> usable;
        for (int i = 0; i < (int)rels.size(); i++)
            if (rels[i].large_prime == 0 || rels[i].large_prime < 0)
                usable.push_back(i);
        if ((int)usable.size() > target + 100) usable.resize(target + 100);

        int nrels = usable.size();
        fprintf(stderr, "LinAlg: %dx%d\n", nrels, fb_size);

        std::vector<BitRow> mat(nrels);
        for (int i = 0; i < nrels; i++) {
            mat[i].init(fb_size, nrels); mat[i].shist(i);
            auto &exps = rels[usable[i]].exps;
            for (int j = 0; j < fb_size; j++)
                if (exps[j] & 1) mat[i].sbit(j);
        }

        int rank = 0;
        for (int col = 0; col < fb_size && rank < nrels; col++) {
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
            std::vector<long> se(fb_size, 0);
            std::vector<long> lps;

            for (int i = 0; i < nrels; i++) {
                if (!mat[row].ghist(i)) continue;
                int ri = usable[i];
                mpz_mul(px, px, rels[ri].p_val);
                mpz_mod(px, px, N);
                for (int j = 0; j < fb_size; j++)
                    se[j] += rels[ri].exps[j];
                if (rels[ri].large_prime < 0) lps.push_back(-rels[ri].large_prime);
            }

            bool ok = true;
            for (int j = 0; j < fb_size; j++) if (se[j] & 1) { ok = false; break; }
            if (!ok) continue;

            mpz_set_ui(py, 1);
            for (int j = 1; j < fb_size; j++) {
                if (se[j] == 0) continue;
                mpz_set_ui(tmp, fb[j]);
                mpz_powm_ui(tmp, tmp, se[j] / 2, N);
                mpz_mul(py, py, tmp);
                mpz_mod(py, py, N);
            }
            for (long lp : lps) {
                mpz_set_ui(tmp, lp);
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
                    fprintf(stderr, "CFRAC factored in %.3fs (%d steps)\n", elapsed, step);
                    mpz_clear(cof);
                    found = 1;
                }
            }
        }
        mpz_clear(px); mpz_clear(py);
        if (!found) { fprintf(stderr, "No factor.\n"); printf("FAIL\n"); }
    }

cleanup:
    for (auto &r : rels) mpz_clear(r.p_val);
    mpz_clear(N); mpz_clear(tmp); mpz_clear(g);
    return 0;
}
