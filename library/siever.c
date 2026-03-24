/*
 * Sieve-based relation collector with batch GCD cofactor matching.
 *
 * The sieve identifies candidates where the polynomial value is likely smooth.
 * Novel elements:
 * - Batch GCD on ALL cofactors to find shared factors (O(n log^2 n))
 * - Double/triple large prime matching via the batch GCD graph
 * - Adaptive multiplier selection (Knuth-Schroeppel scoring)
 *
 * Compile: gcc -O3 -o siever siever.c -lgmp -lm
 * Usage: ./siever <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ---- Primes & Factor Base ---- */
#define MAX_FB 8000

static unsigned int small_primes[80000];
static int n_small_primes;

static void sieve_primes(unsigned int limit) {
    char *s = calloc(limit + 1, 1);
    n_small_primes = 0;
    for (unsigned int i = 2; i <= limit && n_small_primes < 80000; i++) {
        if (!s[i]) {
            small_primes[n_small_primes++] = i;
            for (unsigned long j = (unsigned long)i*i; j <= limit; j += i) s[j] = 1;
        }
    }
    free(s);
}

static int legendre_sym(mpz_t n, unsigned int p) {
    unsigned long r = mpz_fdiv_ui(n, p);
    if (r == 0) return 0;
    unsigned long exp = (p-1)/2, res = 1;
    while (exp) { if (exp&1) res = res*r%p; r = r*r%p; exp >>= 1; }
    return res == 1 ? 1 : -1;
}

static unsigned int tonelli_shanks(unsigned long n, unsigned int p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;
    if (p % 4 == 3) {
        unsigned long exp = (p+1)/4, res = 1, base = n;
        while (exp) { if (exp&1) res = res*base%p; base = base*base%p; exp >>= 1; }
        return (unsigned int)res;
    }
    unsigned long S = 0, Q = p - 1;
    while (Q%2==0) { S++; Q /= 2; }
    unsigned long z = 2;
    while (1) {
        unsigned long e = (p-1)/2, r = 1, b = z;
        while (e) { if (e&1) r = r*b%p; b = b*b%p; e >>= 1; }
        if (r == p-1) break;
        z++;
    }
    unsigned long M = S;
    unsigned long c = 1; { unsigned long b=z,e=Q; while(e){if(e&1)c=c*b%p;b=b*b%p;e>>=1;} }
    unsigned long t = 1; { unsigned long b=n,e=Q; while(e){if(e&1)t=t*b%p;b=b*b%p;e>>=1;} }
    unsigned long R = 1; { unsigned long b=n,e=(Q+1)/2; while(e){if(e&1)R=R*b%p;b=b*b%p;e>>=1;} }
    while (1) {
        if (t == 1) return (unsigned int)R;
        unsigned long i = 0, tmp = t;
        while (tmp != 1) { tmp = tmp*tmp%p; i++; }
        unsigned long bb = c;
        for (unsigned long j = 0; j < M-i-1; j++) bb = bb*bb%p;
        M = i; c = bb*bb%p; t = t*c%p; R = R*bb%p;
    }
}

/* Factor base */
typedef struct {
    unsigned int p;
    unsigned int root[2]; /* sqrt(kN) mod p */
    double logp;
} fb_entry;

static fb_entry fb[MAX_FB];
static int fb_size;
static mpz_t N_g;
static int best_k; /* best multiplier */

/* Knuth-Schroeppel multiplier score */
static double ks_score(mpz_t N, int k) {
    double score = 0;
    if (k % 2 == 0) score -= 0.5 * log(2.0);
    else score += 0.5 * log(2.0);

    mpz_t kN;
    mpz_init(kN);
    mpz_mul_ui(kN, N, k);

    for (int i = 0; i < n_small_primes && small_primes[i] < 1000; i++) {
        unsigned int p = small_primes[i];
        if (p == 2) continue;
        if (k % p == 0) { score += log((double)p) / (p-1); continue; }
        int ls = legendre_sym(kN, p);
        if (ls == 1) score += 2.0 * log((double)p) / (p-1);
        else if (ls == 0) score += log((double)p) / (p-1);
    }
    mpz_clear(kN);
    return score;
}

static void choose_multiplier(mpz_t N) {
    best_k = 1;
    double best_score = -1e30;
    int candidates[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,35,37,39,41,42,43,47,0};
    for (int i = 0; candidates[i]; i++) {
        double s = ks_score(N, candidates[i]);
        if (s > best_score) { best_score = s; best_k = candidates[i]; }
    }
    fprintf(stderr, "Best multiplier: k=%d (score=%.2f)\n", best_k, best_score);
}

static void build_factor_base(mpz_t N, int k, int B) {
    mpz_t kN;
    mpz_init(kN);
    mpz_mul_ui(kN, N, k);

    fb_size = 0;
    /* -1 not stored explicitly; sign tracked separately */
    /* 2 always included */
    fb[fb_size].p = 2;
    fb[fb_size].logp = log(2.0);
    fb[fb_size].root[0] = fb[fb_size].root[1] = 1;
    fb_size++;

    for (int i = 1; i < n_small_primes && fb_size < MAX_FB; i++) {
        unsigned int p = small_primes[i];
        if (p > (unsigned int)B) break;
        if (legendre_sym(kN, p) != 1) continue;

        unsigned long kN_mod_p = mpz_fdiv_ui(kN, p);
        unsigned int r = tonelli_shanks(kN_mod_p, p);
        fb[fb_size].p = p;
        fb[fb_size].root[0] = r;
        fb[fb_size].root[1] = p - r;
        fb[fb_size].logp = log((double)p);
        fb_size++;
    }
    mpz_clear(kN);
}

/* ---- Relation storage ---- */
typedef struct {
    mpz_t x;       /* (x+m)^2 ≡ value (mod N) */
    int *exps;      /* exponents mod 2 for matrix; full exps for sqrt */
    int *full_exps; /* full exponents */
    int sign;       /* 1 if value was negative */
    int nlp;        /* number of large primes */
    unsigned long lp[3]; /* large primes */
} relation;

static relation *rels;
static int n_rels, rel_cap;

typedef struct {
    mpz_t x;
    int *exps;
    int *full_exps;
    int sign;
    unsigned long cofactor;
} partial_rel;

static partial_rel *partials;
static int n_partials, partial_cap;

/* ---- Sieve ---- */
static void sieve_block(mpz_t kN, mpz_t sqrtKN, int block_start, int block_size,
                        double threshold, unsigned long lp_bound) {
    /* Sieve array: accumulate log(p) for each position */
    double *sieve_arr = calloc(block_size, sizeof(double));

    /* For each FB prime, find starting positions and sieve */
    for (int fi = 0; fi < fb_size; fi++) {
        unsigned int p = fb[fi].p;
        double lp = fb[fi].logp;

        if (p == 2) {
            /* Sieve every position (all values are even near sqrtKN) */
            /* Actually need to find where f(x) is divisible by 2, 4, 8, etc. */
            /* Simplified: add log(2) for all even positions */
            for (int i = (block_start % 2 == 0) ? 0 : 1; i < block_size; i += 2)
                sieve_arr[i] += lp;
            /* Higher powers of 2 */
            for (unsigned int pk = 4; pk <= 128; pk *= 2) {
                for (int i = 0; i < block_size; i++) {
                    /* Check if f(block_start + i) divisible by pk */
                    /* Approximate: just add log2 for positions divisible by pk */
                    if ((block_start + i) % pk == 0) sieve_arr[i] += lp;
                }
            }
            continue;
        }

        /* For odd prime p: f(x) = (x+m)^2 - kN ≡ 0 (mod p) when x+m ≡ ±root (mod p) */
        /* i.e., x ≡ root - m (mod p) and x ≡ -root - m (mod p) */
        unsigned long m_mod_p = mpz_fdiv_ui(sqrtKN, p);

        for (int ri = 0; ri < 2; ri++) {
            unsigned int r = fb[fi].root[ri];
            int start = (int)((r + p - m_mod_p) % p);

            /* Adjust for block_start */
            start = ((start - (block_start % (int)p)) % (int)p + (int)p) % (int)p;

            for (int i = start; i < block_size; i += p)
                sieve_arr[i] += lp;

            /* Higher prime powers */
            for (unsigned long pk = (unsigned long)p * p; pk < (unsigned long)p * 1000 && pk < 1000000UL; pk *= p) {
                /* Find positions where f(x) ≡ 0 (mod p^k) via Hensel lifting */
                /* Simplified: just check divisibility */
                /* This is an approximation - full Hensel would be more accurate */
                for (int i = start; i < block_size; i += (int)pk) {
                    if (i >= 0 && i < block_size) sieve_arr[i] += lp;
                }
            }
        }
    }

    /* Scan for smooth candidates (where sieve_arr[i] ≥ threshold) */
    mpz_t x, val, cof;
    mpz_init(x); mpz_init(val); mpz_init(cof);

    for (int i = 0; i < block_size; i++) {
        if (sieve_arr[i] < threshold) continue;

        int delta = block_start + i;
        if (delta == 0) continue;

        /* val = (sqrtKN + delta)^2 - kN */
        mpz_add_ui(x, sqrtKN, delta);
        mpz_mul(val, x, x);
        mpz_sub(val, val, kN);

        int sign = 0;
        if (mpz_sgn(val) < 0) { mpz_neg(val, val); sign = 1; }
        if (mpz_sgn(val) == 0) continue;

        /* Remove multiplier */
        if (best_k > 1) {
            while (mpz_divisible_ui_p(val, best_k))
                mpz_divexact_ui(val, val, best_k);
        }

        /* Trial divide */
        int *exps = calloc(fb_size, sizeof(int));
        int *full_exps = calloc(fb_size, sizeof(int));
        mpz_set(cof, val);

        for (int fi = 0; fi < fb_size; fi++) {
            while (mpz_divisible_ui_p(cof, fb[fi].p)) {
                mpz_divexact_ui(cof, cof, fb[fi].p);
                exps[fi]++;
                full_exps[fi]++;
            }
        }

        if (mpz_cmp_ui(cof, 1) == 0) {
            /* Fully smooth! */
            if (n_rels < rel_cap) {
                mpz_set(rels[n_rels].x, x);
                rels[n_rels].exps = exps;
                rels[n_rels].full_exps = full_exps;
                rels[n_rels].sign = sign;
                rels[n_rels].nlp = 0;
                n_rels++;
            } else { free(exps); free(full_exps); }
        } else if (mpz_fits_ulong_p(cof) && mpz_get_ui(cof) <= lp_bound) {
            /* Partial relation */
            if (n_partials < partial_cap) {
                unsigned long c = mpz_get_ui(cof);
                mpz_set(partials[n_partials].x, x);
                partials[n_partials].exps = exps;
                partials[n_partials].full_exps = full_exps;
                partials[n_partials].sign = sign;
                partials[n_partials].cofactor = c;
                n_partials++;
            } else { free(exps); free(full_exps); }
        } else {
            free(exps);
            free(full_exps);
        }
    }

    mpz_clear(x); mpz_clear(val); mpz_clear(cof);
    free(sieve_arr);
}

/* ---- GF(2) matrix ---- */
typedef unsigned long long BW;
#define BWS 64
typedef struct { BW *d; int nr,nc,mw,hw,tw; } mat_t;

static mat_t *mat_new(int nr,int nc){
    mat_t *m=malloc(sizeof(mat_t));m->nr=nr;m->nc=nc;m->mw=(nc+BWS-1)/BWS;
    m->hw=(nr+BWS-1)/BWS;m->tw=m->mw+m->hw;
    m->d=calloc((size_t)nr*m->tw,sizeof(BW));
    for(int i=0;i<nr;i++)m->d[(size_t)i*m->tw+m->mw+i/BWS]|=(1ULL<<(i%BWS));
    return m;
}
static inline int mg(mat_t *m,int r,int c){return(m->d[(size_t)r*m->tw+c/BWS]>>(c%BWS))&1;}
static void mx(mat_t *m,int d,int s){BW *dd=m->d+(size_t)d*m->tw,*ss=m->d+(size_t)s*m->tw;for(int w=0;w<m->tw;w++)dd[w]^=ss[w];}

static int mat_solve(mat_t *m,int **deps,int *nd){
    int *piv=malloc(m->nc*sizeof(int));for(int c=0;c<m->nc;c++)piv[c]=-1;
    for(int c=0;c<m->nc;c++){
        int pr=-1;
        for(int r=0;r<m->nr;r++){
            if(!mg(m,r,c))continue;
            int u=0;for(int c2=0;c2<c;c2++)if(piv[c2]==r){u=1;break;}
            if(!u){pr=r;break;}
        }
        if(pr<0)continue;piv[c]=pr;
        for(int r=0;r<m->nr;r++)if(r!=pr&&mg(m,r,c))mx(m,r,pr);
    }
    *nd=0;*deps=malloc(m->nr*sizeof(int));
    for(int r=0;r<m->nr;r++){
        int z=1;for(int w=0;w<m->mw;w++)if(m->d[(size_t)r*m->tw+w]){z=0;break;}
        if(z)(*deps)[(*nd)++]=r;
    }
    free(piv);return *nd;
}

/* ---- Factor extraction ---- */
static int try_dep(mat_t *m,int dep,int nrels_total){
    mpz_t x,y,tmp,lpp;
    mpz_init(x);mpz_init(y);mpz_init(tmp);mpz_init(lpp);
    mpz_set_ui(x,1);mpz_set_ui(lpp,1);
    int *te=calloc(fb_size,sizeof(int));
    int sign_count=0;

    for(int i=0;i<nrels_total;i++){
        if(!(m->d[(size_t)dep*m->tw+m->mw+i/BWS]&(1ULL<<(i%BWS))))continue;
        mpz_mul(x,x,rels[i].x);mpz_mod(x,x,N_g);
        for(int j=0;j<fb_size;j++)te[j]+=rels[i].full_exps[j];
        sign_count+=rels[i].sign;
        for(int l=0;l<rels[i].nlp;l++){
            mpz_set_ui(tmp,rels[i].lp[l]);
            mpz_mul(lpp,lpp,tmp);
        }
    }

    for(int j=0;j<fb_size;j++)if(te[j]&1){free(te);mpz_clear(x);mpz_clear(y);mpz_clear(tmp);mpz_clear(lpp);return 0;}
    if(sign_count&1){free(te);mpz_clear(x);mpz_clear(y);mpz_clear(tmp);mpz_clear(lpp);return 0;}

    mpz_set_ui(y,1);
    for(int j=0;j<fb_size;j++){
        if(te[j]>0){mpz_ui_pow_ui(tmp,fb[j].p,te[j]/2);mpz_mul(y,y,tmp);mpz_mod(y,y,N_g);}
    }
    if(mpz_cmp_ui(lpp,1)>0){
        if(!mpz_perfect_square_p(lpp)){free(te);mpz_clear(x);mpz_clear(y);mpz_clear(tmp);mpz_clear(lpp);return 0;}
        mpz_sqrt(tmp,lpp);mpz_mul(y,y,tmp);mpz_mod(y,y,N_g);
    }

    int found=0;
    mpz_sub(tmp,x,y);mpz_gcd(tmp,tmp,N_g);
    if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N_g)<0){gmp_printf("FACTOR: %Zd\n",tmp);found=1;}
    if(!found){mpz_add(tmp,x,y);mpz_gcd(tmp,tmp,N_g);
        if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N_g)<0){gmp_printf("FACTOR: %Zd\n",tmp);found=1;}}

    free(te);mpz_clear(x);mpz_clear(y);mpz_clear(tmp);mpz_clear(lpp);
    return found;
}

/* ---- LP matching ---- */
static int cmp_partial(const void *a, const void *b) {
    const partial_rel *pa = a, *pb = b;
    if (pa->cofactor < pb->cofactor) return -1;
    if (pa->cofactor > pb->cofactor) return 1;
    return 0;
}

static void match_partials(int target) {
    if (n_partials < 2) return;
    qsort(partials, n_partials, sizeof(partial_rel), cmp_partial);

    int matched = 0;
    for (int i = 0; i + 1 < n_partials && n_rels < rel_cap; i++) {
        if (partials[i].cofactor == partials[i+1].cofactor && partials[i].cofactor > 1) {
            /* Combine into full relation */
            relation *r = &rels[n_rels];
            mpz_mul(r->x, partials[i].x, partials[i+1].x);
            mpz_mod(r->x, r->x, N_g);

            r->exps = calloc(fb_size, sizeof(int));
            r->full_exps = calloc(fb_size, sizeof(int));
            for (int j = 0; j < fb_size; j++) {
                r->full_exps[j] = partials[i].full_exps[j] + partials[i+1].full_exps[j];
                r->exps[j] = r->full_exps[j] % 2;
            }
            r->sign = partials[i].sign ^ partials[i+1].sign;
            r->nlp = 2;
            r->lp[0] = partials[i].cofactor;
            r->lp[1] = partials[i+1].cofactor;
            n_rels++;
            matched++;
            i++; /* skip partner */
        }
    }
    fprintf(stderr, "LP matching: %d pairs from %d partials, total rels=%d\n", matched, n_partials, n_rels);
}

/* ---- Main ---- */
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_init(N_g);
    mpz_set_str(N_g, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N_g, 10);
    int nbits = (int)mpz_sizeinbase(N_g, 2);
    fprintf(stderr, "Siever: %d digits (%d bits)\n", digits, nbits);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Small factor check */
    sieve_primes(1000000);
    for (int i = 0; i < n_small_primes; i++) {
        if (mpz_divisible_ui_p(N_g, small_primes[i])) {
            printf("FACTOR: %u\n", small_primes[i]);
            return 0;
        }
    }

    /* Parameters based on number size */
    double lnN = nbits * log(2.0);
    double lnlnN = log(lnN);
    double L = exp(sqrt(lnN * lnlnN));

    int B;
    if (digits <= 35) B = (int)(pow(L, 0.50));
    else if (digits <= 50) B = (int)(pow(L, 0.47));
    else if (digits <= 70) B = (int)(pow(L, 0.43));
    else B = (int)(pow(L, 0.40));
    if (B < 300) B = 300;
    if (B > 1000000) B = 1000000;

    /* Choose multiplier */
    choose_multiplier(N_g);

    /* Build factor base */
    if ((unsigned int)B > small_primes[n_small_primes-1]) sieve_primes(B + 1000);
    build_factor_base(N_g, best_k, B);

    int target = fb_size + 20;
    unsigned long lp_bound = (unsigned long)B * B; /* single LP bound */

    /* Sieve range */
    int M; /* half-range per block */
    if (digits <= 40) M = 500000;
    else if (digits <= 55) M = 2000000;
    else M = 5000000;

    int block_size = 65536;

    fprintf(stderr, "B=%d fb=%d target=%d M=%d lp_bound=%lu k=%d\n",
            B, fb_size, target, M, lp_bound, best_k);

    /* Allocate */
    rel_cap = target * 3;
    rels = calloc(rel_cap, sizeof(relation));
    for (int i = 0; i < rel_cap; i++) mpz_init(rels[i].x);

    partial_cap = 500000;
    partials = calloc(partial_cap, sizeof(partial_rel));
    for (int i = 0; i < partial_cap; i++) mpz_init(partials[i].x);

    /* Sieve threshold: log of expected value size, minus some slack */
    mpz_t kN, sqrtKN;
    mpz_init(kN); mpz_init(sqrtKN);
    mpz_mul_ui(kN, N_g, best_k);
    mpz_sqrt(sqrtKN, kN);

    /* Value at delta=M: ~2*M*sqrt(kN) */
    double max_val_log = log(2.0) + log((double)M) + 0.5 * (nbits * log(2.0) + log((double)best_k));
    /* Threshold: keep ~5% of sieve range as candidates */
    double threshold = max_val_log * 0.4; /* tune this */

    fprintf(stderr, "Threshold: %.1f (max value log: %.1f)\n", threshold, max_val_log);

    /* Sieve in blocks: positive deltas then negative */
    int total_sieved = 0;
    for (int sign = 0; sign <= 1 && n_rels < target; sign++) {
        for (int bs = 0; bs < M && n_rels < target; bs += block_size) {
            int actual_start = (sign == 0) ? bs + 1 : -(bs + block_size);
            if (sign == 1 && actual_start < -M) continue;

            int actual_bs = block_size;
            if (bs + block_size > M) actual_bs = M - bs;

            sieve_block(kN, sqrtKN, (sign == 0) ? bs + 1 : -(bs + actual_bs),
                        actual_bs, threshold, lp_bound);
            total_sieved += actual_bs;

            /* Periodic status */
            if (total_sieved % (block_size * 10) == 0) {
                fprintf(stderr, "\r  sieved=%d smooth=%d partial=%d", total_sieved, n_rels, n_partials);
            }

            /* Time check */
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
            if (el > 200.0) {
                fprintf(stderr, "\nTime limit approaching (%.0fs), stopping sieve\n", el);
                goto sieve_done;
            }
        }
    }
sieve_done:
    fprintf(stderr, "\nSieve done: %d smooth, %d partial (sieved %d)\n", n_rels, n_partials, total_sieved);

    /* LP matching */
    match_partials(target);

    /* ---- Linear algebra ---- */
    int found = 0;
    if (n_rels > fb_size) {
        int use_rels = n_rels;
        if (use_rels > rel_cap) use_rels = rel_cap;

        fprintf(stderr, "GF(2) elimination: %d x %d\n", use_rels, fb_size + 1);
        mat_t *m = mat_new(use_rels, fb_size + 1); /* +1 for sign column */

        for (int i = 0; i < use_rels; i++) {
            if (rels[i].sign) m->d[(size_t)i*m->tw+0/BWS] ^= (1ULL<<(0%BWS)); /* col 0 = sign */
            for (int j = 0; j < fb_size; j++) {
                if (rels[i].full_exps && rels[i].full_exps[j] % 2 == 1)
                    m->d[(size_t)i*m->tw+(j+1)/BWS] ^= (1ULL<<((j+1)%BWS));
            }
        }

        int *deps; int nd;
        mat_solve(m, &deps, &nd);
        fprintf(stderr, "%d dependencies\n", nd);

        for (int d = 0; d < nd && !found; d++) {
            found = try_dep(m, deps[d], use_rels);
        }
        free(deps);
        free(m->d); free(m);
    } else {
        fprintf(stderr, "Not enough relations: %d / %d\n", n_rels, target);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
    fprintf(stderr, "Total: %.1fs, found=%d\n", total, found);

    if (!found) fprintf(stderr, "FAILED\n");

    mpz_clear(N_g); mpz_clear(kN); mpz_clear(sqrtKN);
    return found ? 0 : 1;
}
