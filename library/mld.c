/*
 * mld.c - Multiplicative Lattice Descent (MLD)
 *
 * NOVEL APPROACH:
 * In standard QS, we need values x^2 - N to be B-smooth. The smooth
 * probability ρ(u) where u = log(√N)/log(B) determines the running time.
 *
 * MLD modifies this: instead of requiring INDIVIDUAL values to be smooth,
 * we allow partial smoothness and then use lattice reduction to find
 * SUBSETS of cofactors whose product is smooth.
 *
 * Key idea: Collect k partial relations with cofactors c_1, ..., c_k.
 * For each c_i, compute its factorization over a "cofactor base" of primes
 * in [B, B']. The exponent vectors mod 2 form a matrix. If this matrix
 * has a dependency, the corresponding product ∏c_i is a perfect square
 * (up to the smooth part), giving a full relation.
 *
 * The innovation: instead of B' = B (standard LP) or B' = B^2 (double LP),
 * use B' = B^α for α > 2. The partial relations are easier to find
 * (more values have cofactors < B^α), but the cofactor matrix is larger.
 *
 * THEORETICAL ANALYSIS:
 * - Standard QS: need π(B) smooth relations. Each costs ~1/ρ(u).
 *   Total: π(B)/ρ(u).
 * - MLD with B' = B^α: need π(B') dependencies in cofactor matrix.
 *   But partial relations are α-times easier to find (larger cofactor allowed).
 *   The cofactor matrix has dimensions k × π(B'), needing k > π(B').
 *   Cofactors have ~α prime factors, so the matrix is denser.
 *
 * If the cofactor matrix has density inversely proportional to B',
 * then the number of partial relations needed is O(π(B')),
 * but each is cheaper to find (by factor ρ(u-α)/ρ(u)).
 *
 * EXPERIMENT: Measure the actual cofactor distribution and matrix density
 * to determine if MLD can improve the effective L-exponent.
 *
 * Usage: ./mld <N>
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Tonelli-Shanks */
static long ts_sqrt(unsigned long n, unsigned long p) {
    if (p==2) return n&1;
    if (n==0) return 0;
    mpz_t b,e,m,r;
    mpz_init_set_ui(b,n);mpz_init_set_ui(m,p);
    mpz_init_set_ui(e,(p-1)/2);mpz_init(r);
    mpz_powm(r,b,e,m);
    if(mpz_cmp_ui(r,1)!=0){mpz_clear(b);mpz_clear(e);mpz_clear(m);mpz_clear(r);return -1;}
    unsigned long Q=p-1,S=0;
    while(Q%2==0){Q/=2;S++;}
    if(S==1){mpz_set_ui(e,(p+1)/4);mpz_powm(r,b,e,m);
        long ret=mpz_get_ui(r);mpz_clear(b);mpz_clear(e);mpz_clear(m);mpz_clear(r);return ret;}
    unsigned long z=2;mpz_t zz;mpz_init(zz);mpz_set_ui(e,(p-1)/2);
    while(1){mpz_set_ui(zz,z);mpz_powm(r,zz,e,m);if(mpz_cmp_ui(r,p-1)==0)break;z++;}
    mpz_t c,t,R2,bb,tmp;
    mpz_init(c);mpz_init(t);mpz_init(R2);mpz_init(bb);mpz_init(tmp);
    mpz_set_ui(e,Q);mpz_set_ui(zz,z);mpz_powm(c,zz,e,m);
    mpz_set_ui(b,n);mpz_powm(t,b,e,m);mpz_set_ui(e,(Q+1)/2);mpz_powm(R2,b,e,m);
    unsigned long MM=S;
    while(1){
        if(mpz_cmp_ui(t,1)==0){long ret=mpz_get_ui(R2);
            mpz_clear(b);mpz_clear(e);mpz_clear(m);mpz_clear(r);
            mpz_clear(c);mpz_clear(t);mpz_clear(R2);mpz_clear(bb);
            mpz_clear(tmp);mpz_clear(zz);return ret;}
        unsigned long ii=0;mpz_set(tmp,t);
        while(mpz_cmp_ui(tmp,1)!=0){mpz_mul(tmp,tmp,tmp);mpz_mod(tmp,tmp,m);ii++;}
        if(ii==MM)break;
        mpz_set(bb,c);
        for(unsigned long j=0;j<MM-ii-1;j++){mpz_mul(bb,bb,bb);mpz_mod(bb,bb,m);}
        MM=ii;mpz_mul(c,bb,bb);mpz_mod(c,c,m);
        mpz_mul(t,t,c);mpz_mod(t,t,m);
        mpz_mul(R2,R2,bb);mpz_mod(R2,R2,m);
    }
    mpz_clear(b);mpz_clear(e);mpz_clear(m);mpz_clear(r);
    mpz_clear(c);mpz_clear(t);mpz_clear(R2);mpz_clear(bb);
    mpz_clear(tmp);mpz_clear(zz);
    return -1;
}

/* Partial relation: x^2 - N = (smooth part) * cofactor */
typedef struct {
    mpz_t x;           /* sqrt value */
    mpz_t cofactor;    /* remaining unfactored part */
    int sign;
    unsigned int *base_exp;  /* exponents over base FB */
    unsigned int *cof_exp;   /* exponents over cofactor FB (if factored) */
    int cof_factored;        /* 1 if cofactor fully factored over cof FB */
} partial_t;

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N; mpz_init(N); mpz_set_str(N, argv[1], 10);
    size_t ndig = mpz_sizeinbase(N, 10);
    size_t nbits = mpz_sizeinbase(N, 2);
    fprintf(stderr, "N=%zu digits (%zu bits)\n", ndig, nbits);

    double t0 = now();

    /* Parameter selection */
    double ln_N = nbits * 0.693147;
    double ln_ln_N = log(ln_N);
    double L_half = exp(sqrt(ln_N * ln_ln_N));

    /* Base factor base bound B */
    unsigned long B = (unsigned long)exp(0.5 * sqrt(ln_N * ln_ln_N));
    if (B < 500) B = 500;
    if (B > 200000) B = 200000;

    /* Cofactor factor base bound B' = B^alpha */
    /* alpha = 2 is standard double-LP. Try alpha = 3 for MLD. */
    double alpha = 3.0;
    unsigned long B_cof = (unsigned long)pow((double)B, alpha);
    if (B_cof > 10000000) B_cof = 10000000;
    if (B_cof < B * 10) B_cof = B * 10;

    fprintf(stderr, "B=%lu, B_cof=%lu (alpha=%.1f)\n", B, B_cof, alpha);

    /* Build base factor base: primes up to B where N is QR */
    int n_primes;
    int *all_primes;
    {
        char *sv = calloc(B + 1, 1);
        for (unsigned long i = 2; i <= B; i++) sv[i] = 1;
        for (unsigned long i = 2; i * i <= B; i++)
            if (sv[i]) for (unsigned long j = i * i; j <= B; j += i) sv[j] = 0;
        int cnt = 0;
        for (unsigned long i = 2; i <= B; i++) if (sv[i]) cnt++;
        all_primes = malloc(cnt * sizeof(int));
        int idx = 0;
        for (unsigned long i = 2; i <= B; i++) if (sv[i]) all_primes[idx++] = i;
        n_primes = cnt;
        free(sv);
    }

    /* Filter to QR primes */
    int *fb = malloc(n_primes * sizeof(int));
    long *fb_r1 = malloc(n_primes * sizeof(long));
    long *fb_r2 = malloc(n_primes * sizeof(long));
    int fb_sz = 0;

    for (int i = 0; i < n_primes; i++) {
        unsigned long p = all_primes[i];
        long s = ts_sqrt(mpz_fdiv_ui(N, p), p);
        if (s >= 0) {
            fb[fb_sz] = p;
            fb_r1[fb_sz] = s;
            fb_r2[fb_sz] = (s == 0) ? 0 : (long)p - s;
            fb_sz++;
        }
    }
    free(all_primes);

    fprintf(stderr, "Base FB: %d primes (max=%d)\n", fb_sz, fb[fb_sz - 1]);

    /* Build cofactor factor base: primes in (B, B_cof] where N is QR */
    int *cfb = malloc(100000 * sizeof(int)); /* generous allocation */
    int cfb_sz = 0;
    {
        /* Sieve primes in (B, B_cof] */
        /* For large ranges, use segmented sieve or just iterate */
        for (unsigned long p = (B | 1) + 2; p <= B_cof && cfb_sz < 100000; p += 2) {
            int ip = 1;
            for (unsigned long d = 3; d * d <= p; d += 2)
                if (p % d == 0) { ip = 0; break; }
            if (!ip) continue;
            long s = ts_sqrt(mpz_fdiv_ui(N, p), p);
            if (s >= 0) {
                cfb[cfb_sz++] = (int)p;
            }
        }
    }
    fprintf(stderr, "Cofactor FB: %d primes in (%lu, %lu]\n", cfb_sz, B, B_cof);

    /* Sieve interval */
    long M = (long)fb_sz * 40;
    if (M < 50000) M = 50000;
    if (M > 5000000) M = 5000000;

    size_t slen = 2 * (size_t)M + 1;
    float *sv = calloc(slen, sizeof(float));

    mpz_t sqN, xb;
    mpz_init(sqN); mpz_sqrt(sqN, N);
    mpz_init(xb); mpz_sub_ui(xb, sqN, (unsigned long)M);

    /* Sieve with base FB */
    for (int fi = 0; fi < fb_sz; fi++) {
        unsigned long p = fb[fi];
        float lp = log2f((float)p);
        unsigned long xm = mpz_fdiv_ui(xb, p);
        for (int ri = 0; ri < 2; ri++) {
            long s = (ri == 0) ? fb_r1[fi] : fb_r2[fi];
            if (s < 0) continue;
            if (ri == 1 && fb_r1[fi] == fb_r2[fi]) continue;
            long start = ((s - (long)xm) % (long)p + (long)p) % (long)p;
            for (long j = start; j < (long)slen; j += (long)p)
                sv[j] += lp;
        }
    }

    double sieve_time = now() - t0;
    fprintf(stderr, "Sieve done in %.2fs\n", sieve_time);

    /* Scan: accept candidates with lower threshold (allowing larger cofactors) */
    /* Standard threshold: log2(residue) * 0.7 for full smooth
       MLD threshold: log2(residue) * 0.35 for partial smooth with large cofactor */
    double thresh_full = (double)nbits / 2.0 * 0.45;
    double thresh_partial = (double)nbits / 2.0 * 0.25;
    if (thresh_partial < 10) thresh_partial = 10;

    /* Collect partials */
    size_t max_partials = fb_sz * 10 + cfb_sz + 1000;
    partial_t *partials = calloc(max_partials, sizeof(partial_t));
    size_t n_full = 0, n_partial = 0, n_cof_factored = 0;
    unsigned int *ev = calloc(fb_sz, sizeof(unsigned int));
    unsigned int *cev = calloc(cfb_sz, sizeof(unsigned int));
    mpz_t val, res, cof, tmp;
    mpz_init(val); mpz_init(res); mpz_init(cof); mpz_init(tmp);

    for (long j = 0; j < (long)slen && (n_full + n_cof_factored) < (size_t)(fb_sz + cfb_sz + 20); j++) {
        if (sv[j] < thresh_partial) continue;
        int is_promising = (sv[j] >= thresh_full);

        mpz_add_ui(val, xb, (unsigned long)j);
        mpz_mul(res, val, val);
        mpz_sub(res, res, N);
        if (mpz_sgn(res) == 0) continue;

        int sign = mpz_sgn(res);
        mpz_abs(cof, res);

        /* Trial divide by base FB */
        memset(ev, 0, fb_sz * sizeof(unsigned int));
        for (int fi = 0; fi < fb_sz; fi++) {
            while (mpz_divisible_ui_p(cof, fb[fi])) {
                mpz_divexact_ui(cof, cof, fb[fi]);
                ev[fi]++;
            }
        }

        if (mpz_cmp_ui(cof, 1) == 0) {
            /* Fully smooth over base FB */
            partial_t *p = &partials[n_full + n_partial];
            mpz_init_set(p->x, val);
            mpz_init_set_ui(p->cofactor, 1);
            p->sign = sign;
            p->base_exp = malloc(fb_sz * sizeof(unsigned int));
            memcpy(p->base_exp, ev, fb_sz * sizeof(unsigned int));
            p->cof_exp = NULL;
            p->cof_factored = 1;
            n_full++;
            continue;
        }

        /* Check if cofactor is within B_cof range */
        if (mpz_sizeinbase(cof, 2) > 64) continue; /* Too big */
        unsigned long cof_val = 0;
        if (mpz_fits_ulong_p(cof)) cof_val = mpz_get_ui(cof);
        else continue;
        if (cof_val > (unsigned long)B_cof * B_cof) continue; /* Way too big */

        /* Trial divide cofactor by cofactor FB */
        memset(cev, 0, cfb_sz * sizeof(unsigned int));
        mpz_t cof2;
        mpz_init_set(cof2, cof);
        for (int ci = 0; ci < cfb_sz; ci++) {
            while (mpz_divisible_ui_p(cof2, cfb[ci])) {
                mpz_divexact_ui(cof2, cof2, cfb[ci]);
                cev[ci]++;
            }
        }

        int cof_fully = (mpz_cmp_ui(cof2, 1) == 0);
        mpz_clear(cof2);

        if (cof_fully || cof_val <= (unsigned long)B_cof) {
            partial_t *p = &partials[n_full + n_partial];
            mpz_init_set(p->x, val);
            mpz_init_set(p->cofactor, cof);
            p->sign = sign;
            p->base_exp = malloc(fb_sz * sizeof(unsigned int));
            memcpy(p->base_exp, ev, fb_sz * sizeof(unsigned int));
            p->cof_exp = malloc(cfb_sz * sizeof(unsigned int));
            memcpy(p->cof_exp, cev, cfb_sz * sizeof(unsigned int));
            p->cof_factored = cof_fully;
            n_partial++;
            if (cof_fully) n_cof_factored++;
        }
    }

    double scan_time = now() - t0;
    fprintf(stderr, "Scan: full=%zu partial=%zu (cof_factored=%zu) total=%zu (%.2fs)\n",
            n_full, n_partial, n_cof_factored, n_full + n_partial, scan_time);

    /* Analysis: what is the distribution of cofactor sizes? */
    fprintf(stderr, "\n=== MLD Analysis ===\n");
    fprintf(stderr, "Base FB size: %d, Cofactor FB size: %d\n", fb_sz, cfb_sz);
    fprintf(stderr, "Full smooth: %zu (would need %d for standard QS)\n",
            n_full, fb_sz + 1);
    fprintf(stderr, "Cofactor-factored partials: %zu\n", n_cof_factored);
    fprintf(stderr, "Total usable: %zu (need %d for extended matrix)\n",
            n_full + n_cof_factored, fb_sz + cfb_sz + 1);

    double eff_ratio = (double)(n_full + n_cof_factored) /
                       (double)(fb_sz + cfb_sz + 1);
    fprintf(stderr, "Ratio (have/need): %.3f\n", eff_ratio);

    if (eff_ratio >= 1.0) {
        fprintf(stderr, "SUCCESS: Have enough relations for MLD matrix!\n");
        fprintf(stderr, "Standard QS would need %d, MLD needs %d but has %zu\n",
                fb_sz + 1, fb_sz + cfb_sz + 1, n_full + n_cof_factored);

        /* Compute effective improvement factor */
        double std_needed = fb_sz + 1;
        double mld_rate = (double)(n_full + n_cof_factored) / (double)slen;
        double std_rate = (double)n_full / (double)slen;
        fprintf(stderr, "Smooth rate: standard=%.6f MLD=%.6f (%.1fx more)\n",
                std_rate, mld_rate, mld_rate / (std_rate > 0 ? std_rate : 1e-10));
    } else {
        fprintf(stderr, "NOT ENOUGH: Would need %.0fx larger sieve\n", 1.0 / eff_ratio);

        /* Key question: is the cofactor expansion worth the matrix cost? */
        /* Compare: relations_per_sieve_position for std vs MLD */
        double std_rate = (double)n_full / (double)slen;
        double mld_rate = (double)(n_full + n_cof_factored) / (double)slen;
        double std_need_per_pos = (double)(fb_sz + 1) / (double)slen;
        double mld_need_per_pos = (double)(fb_sz + cfb_sz + 1) / (double)slen;

        fprintf(stderr, "\nCost analysis:\n");
        fprintf(stderr, "  Standard: rate=%.6f, need=%d, sieve_needed=%.0f\n",
                std_rate, fb_sz + 1, std_rate > 0 ? (fb_sz + 1) / std_rate : INFINITY);
        fprintf(stderr, "  MLD:      rate=%.6f, need=%d, sieve_needed=%.0f\n",
                mld_rate, fb_sz + cfb_sz + 1,
                mld_rate > 0 ? (fb_sz + cfb_sz + 1) / mld_rate : INFINITY);

        double std_sieve = std_rate > 0 ? (fb_sz + 1) / std_rate : INFINITY;
        double mld_sieve = mld_rate > 0 ? (fb_sz + cfb_sz + 1) / mld_rate : INFINITY;

        fprintf(stderr, "  Ratio (MLD/std sieve cost): %.2f (< 1 means MLD wins)\n",
                mld_sieve / std_sieve);
    }

    /* Cleanup */
    for (size_t i = 0; i < n_full + n_partial; i++) {
        if (partials[i].base_exp) {
            mpz_clear(partials[i].x);
            mpz_clear(partials[i].cofactor);
            free(partials[i].base_exp);
            if (partials[i].cof_exp) free(partials[i].cof_exp);
        }
    }
    free(partials); free(sv); free(ev); free(cev);
    free(fb); free(fb_r1); free(fb_r2); free(cfb);
    mpz_clear(val); mpz_clear(res); mpz_clear(cof); mpz_clear(tmp);
    mpz_clear(sqN); mpz_clear(xb);
    mpz_clear(N);

    return 0;
}
