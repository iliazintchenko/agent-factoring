/*
 * Dixon's random squares with batch smoothness detection.
 *
 * Instead of sieving, we:
 * 1. Generate many random x values
 * 2. Compute x^2 mod N for each
 * 3. Use product trees + remainder trees to batch-test smoothness
 *    (Bernstein's method: compute gcd of product of all candidates with
 *     the product of primes in the factor base)
 * 4. Collect smooth relations and do linear algebra
 *
 * This is O(B * polylog(B)) per smoothness test batch, vs O(B) for sieving.
 * The main advantage: we can test candidates that don't come from a polynomial,
 * potentially finding more relations per candidate.
 *
 * Usage: ./dixon_batch <number>
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
#include <gmp.h>

static long power_mod_l(long base, long exp, long mod) {
    long r = 1; base = ((base%mod)+mod)%mod;
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
    long b=c;for(long j=0;j<M2-i-1;j++)b=(__int128)b*b%p;
    M2=i;c=(__int128)b*b%p;t=(__int128)t*c%p;R=(__int128)R*b%p;}
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

struct Relation {
    mpz_t x;
    std::vector<int> exps;
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

    long B = (long)(opt * 0.7);
    if (B < 1000) B = 1000;
    if (B > 5000000) B = 5000000;

    /* Build factor base */
    std::vector<long> all_primes = sieve_primes(B);
    std::vector<long> fb;
    fb.push_back(-1);
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

    fprintf(stderr, "Dixon batch: %zu dig, B=%ld, fb=%d, target=%d\n", digits, B, fb_size, target);

    /* Compute product of factor base primes raised to appropriate powers */
    /* P = prod(p^k) where p^k <= N for each fb prime p */
    mpz_t P;
    mpz_init(P);
    mpz_set_ui(P, 1);
    int nbits = mpz_sizeinbase(N, 2);
    for (int i = 1; i < fb_size; i++) {
        long p = fb[i];
        /* p^k where k = floor(log_p(N)) = floor(nbits / log2(p)) */
        int k = (int)(nbits / log2((double)p));
        if (k < 1) k = 1;
        mpz_t pk;
        mpz_init(pk);
        mpz_ui_pow_ui(pk, p, k);
        mpz_mul(P, P, pk);
        mpz_clear(pk);
    }

    fprintf(stderr, "Product of fb prime powers: %zu bits\n", mpz_sizeinbase(P, 2));

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    /* Main loop: generate batches of candidates, test for smoothness */
    std::vector<Relation> rels;
    int batch_size = 10000;
    int total_tested = 0;

    mpz_t sqrtN;
    mpz_init(sqrtN);
    mpz_sqrt(sqrtN, N);

    while ((int)rels.size() < target) {
        struct timespec tnow;
        clock_gettime(CLOCK_MONOTONIC, &tnow);
        double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
        if (elapsed > 280) {
            fprintf(stderr, "Timeout after %.1fs\n", elapsed);
            break;
        }

        /* Generate batch of x values and compute x^2 mod N */
        std::vector<mpz_t> x_vals(batch_size), q_vals(batch_size);
        for (int i = 0; i < batch_size; i++) {
            mpz_init(x_vals[i]);
            mpz_init(q_vals[i]);
            /* x = sqrtN + random offset */
            mpz_urandomm(x_vals[i], rstate, N);
            /* q = x^2 mod N */
            mpz_mul(q_vals[i], x_vals[i], x_vals[i]);
            mpz_mod(q_vals[i], q_vals[i], N);
        }

        /* Batch smoothness test using GCD with product of fb prime powers.
         * For each q_val, compute gcd(q_val, P). If gcd == q_val, it's smooth.
         * But P is huge, so computing gcd directly is expensive.
         *
         * Better approach: compute R = P mod q_val for each candidate,
         * then gcd(R, q_val) reveals the smooth part.
         *
         * Even better: use a product tree + remainder tree.
         * But for simplicity, just use the direct approach with repeated squaring:
         * smooth(q) iff gcd(q, P^k mod q) = q for large enough k */

        for (int i = 0; i < batch_size; i++) {
            if (mpz_cmp_ui(q_vals[i], 0) == 0) continue;

            /* Compute P mod q_val, then gcd */
            mpz_mod(tmp, P, q_vals[i]);
            mpz_gcd(g, tmp, q_vals[i]);

            /* Extract smooth part: repeatedly divide by gcd until stable */
            mpz_t remainder;
            mpz_init_set(remainder, q_vals[i]);

            /* Iteratively remove smooth factors */
            while (mpz_cmp_ui(g, 1) > 0) {
                while (mpz_divisible_p(remainder, g)) {
                    mpz_divexact(remainder, remainder, g);
                }
                if (mpz_cmp_ui(remainder, 1) == 0) break;
                mpz_mod(tmp, P, remainder);
                mpz_gcd(g, tmp, remainder);
            }

            if (mpz_cmp_ui(remainder, 1) == 0) {
                /* q is smooth! Trial divide to get exponents */
                std::vector<int> exps(fb_size, 0);
                mpz_set(tmp, q_vals[i]);

                for (int j = 1; j < fb_size; j++) {
                    long p = fb[j];
                    while (mpz_divisible_ui_p(tmp, p)) {
                        mpz_divexact_ui(tmp, tmp, p);
                        exps[j]++;
                    }
                }

                if (mpz_cmp_ui(tmp, 1) == 0) {
                    Relation rel;
                    mpz_init_set(rel.x, x_vals[i]);
                    rel.exps = exps;
                    rels.push_back(rel);
                }
            }

            mpz_clear(remainder);
        }

        total_tested += batch_size;

        fprintf(stderr, "  tested %d: %zu/%d rels (%.1fs)\n",
                total_tested, rels.size(), target, elapsed);

        for (int i = 0; i < batch_size; i++) {
            mpz_clear(x_vals[i]);
            mpz_clear(q_vals[i]);
        }
    }

    fprintf(stderr, "%zu rels from %d candidates (%.4f hit rate)\n",
            rels.size(), total_tested, (double)rels.size() / total_tested);

    if ((int)rels.size() <= fb_size) {
        fprintf(stderr, "Not enough.\n");
        printf("FAIL\n");
        goto cleanup;
    }

    {
        /* Linear algebra */
        int nrels = std::min((int)rels.size(), target + 100);
        fprintf(stderr, "LinAlg: %dx%d\n", nrels, fb_size);

        std::vector<BitRow> mat(nrels);
        for (int i = 0; i < nrels; i++) {
            mat[i].init(fb_size, nrels);
            mat[i].shist(i);
            for (int j = 0; j < fb_size; j++)
                if (rels[i].exps[j] & 1) mat[i].sbit(j);
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

            for (int i = 0; i < nrels; i++) {
                if (!mat[row].ghist(i)) continue;
                mpz_mul(px, px, rels[i].x);
                mpz_mod(px, px, N);
                for (int j = 0; j < fb_size; j++)
                    se[j] += rels[i].exps[j];
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
                    fprintf(stderr, "Dixon batch factored in %.3fs\n", elapsed);
                    mpz_clear(cof);
                    found = 1;
                }
            }
        }

        mpz_clear(px); mpz_clear(py);
        if (!found) { fprintf(stderr, "No factor from deps.\n"); printf("FAIL\n"); }
    }

cleanup:
    for (auto &r : rels) mpz_clear(r.x);
    mpz_clear(N); mpz_clear(tmp); mpz_clear(g); mpz_clear(P); mpz_clear(sqrtN);
    gmp_randclear(rstate);
    return 0;
}
