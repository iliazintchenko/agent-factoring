/*
 * hybrid_qs_nfs.c — Hybrid QS-NFS: QS-style 1D sieving on NFS polynomials
 *
 * Novel idea: Use a degree-d NFS polynomial for value generation but
 * apply standard QS 1D sieving for smoothness detection. Fix b=1 and
 * sieve over a, then try b=2, b=3, etc.
 *
 * Advantages over pure QS:
 * - Algebraic norm ≈ m * a^d ≈ N^{1/(d+1)} * a^d (smaller than QS's √N * a for small a)
 * - Rational norm ≈ |a - m| (very small near a=m)
 *
 * Advantages over pure NFS:
 * - 1D sieve (not 2D lattice sieve) — much simpler implementation
 * - Standard QS infrastructure for relation collection
 *
 * The key question: does the hybrid outperform pure QS for 40-70 digit N?
 *
 * Usage: ./hybrid_qs_nfs <N>
 * Output: FACTOR: <p>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

static double wall(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

#define MAX_FB 4096
#define MAX_RELS 16384

typedef struct { int p; int r; } afb_entry_t; /* algebraic factor base: (prime, root) */
typedef struct { int p; } rfb_entry_t; /* rational factor base: just primes */

static mpz_t gN, co[8]; /* N and polynomial coefficients */
static mpz_t gm;        /* base m */
static int gdeg;

/* Evaluate homogeneous polynomial: sum c_i * a^i * b^(d-i) */
static void poly_hom(mpz_t res, long a, long b) {
    mpz_t term, tmp;
    mpz_init(term); mpz_init(tmp);
    mpz_set_ui(res, 0);
    for (int i = 0; i <= gdeg; i++) {
        mpz_set(term, co[i]);
        for (int j = 0; j < i; j++) mpz_mul_si(term, term, a);
        for (int j = 0; j < gdeg - i; j++) mpz_mul_si(term, term, b);
        mpz_add(res, res, term);
    }
    mpz_clear(term); mpz_clear(tmp);
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    double start = wall();

    mpz_init(gN); mpz_init(gm);
    mpz_set_str(gN, argv[1], 10);
    int digits = strlen(argv[1]);

    /* Trial division */
    for (unsigned long p = 2; p < 1000000; p++)
        if (mpz_divisible_ui_p(gN, p)) {
            printf("FACTOR: %lu\n", p); return 0;
        }

    /* Choose degree */
    gdeg = (digits < 35) ? 4 : 5;

    /* Base-m polynomial selection */
    mpz_root(gm, gN, gdeg + 1);
    for (int i = 0; i <= gdeg; i++) mpz_init(co[i]);
    { mpz_t rem; mpz_init(rem); mpz_set(rem, gN);
      for (int i = 0; i < gdeg; i++) mpz_fdiv_qr(rem, co[i], rem, gm);
      mpz_set(co[gdeg], rem); mpz_clear(rem); }

    /* Verify f(m) = N */
    { mpz_t check; mpz_init(check);
      mpz_set(check, co[gdeg]);
      for (int i = gdeg-1; i >= 0; i--) { mpz_mul(check, check, gm); mpz_add(check, check, co[i]); }
      if (mpz_cmp(check, gN) != 0) { fprintf(stderr, "Poly verification FAILED\n"); return 1; }
      mpz_clear(check); }

    /* Factor base */
    double lnN = digits * 2.302585;
    double lnlnN = log(lnN);
    double B_val = exp(0.45 * sqrt(lnN * lnlnN));
    unsigned long B = (unsigned long)B_val;
    if (B < 500) B = 500;
    if (B > 50000) B = 50000;

    /* Rational FB: all primes up to B */
    rfb_entry_t *rfb = malloc(MAX_FB * sizeof(rfb_entry_t));
    int rfb_n = 0;
    for (unsigned long p = 2; p <= B && rfb_n < MAX_FB; p++) {
        int ip = 1;
        for (unsigned long d = 2; d*d <= p; d++) if (p%d==0) { ip=0; break; }
        if (!ip) continue;
        rfb[rfb_n++].p = (int)p;
    }

    /* Algebraic FB: primes with roots of f mod p */
    afb_entry_t *afb = malloc(MAX_FB * sizeof(afb_entry_t));
    int afb_n = 0;
    for (int i = 0; i < rfb_n && afb_n < MAX_FB; i++) {
        int p = rfb[i].p;
        for (int r = 0; r < p && afb_n < MAX_FB; r++) {
            long v = 0, rp = 1;
            for (int j = 0; j <= gdeg; j++) {
                long cm = ((long)mpz_fdiv_ui(co[j], p) + p) % p;
                v = (v + cm * rp) % p;
                rp = rp * r % p;
            }
            if (v == 0) { afb[afb_n].p = p; afb[afb_n].r = r; afb_n++; }
        }
    }

    int target_rels = rfb_n + afb_n + 64;
    fprintf(stderr, "HQNFS: %d-digit, deg=%d, B=%lu, RFB=%d, AFB=%d, target=%d\n",
            digits, gdeg, B, rfb_n, afb_n, target_rels);

    /* Sieve range */
    int sieve_half = 50000;
    int sieve_len = 2 * sieve_half + 1;

    float *rsieve = calloc(sieve_len, sizeof(float));
    float *asieve = calloc(sieve_len, sizeof(float));

    /* Pre-compute log(p) for sieving */
    float *rlog = malloc(rfb_n * sizeof(float));
    float *alog = malloc(afb_n * sizeof(float));
    for (int i = 0; i < rfb_n; i++) rlog[i] = logf((float)rfb[i].p);
    for (int i = 0; i < afb_n; i++) alog[i] = logf((float)afb[i].p);

    /* Relations storage */
    typedef struct {
        mpz_t x_val; /* a - b*m for congruence of squares */
        int *ev;     /* exponent vector: [rsign, r_exp[rfb_n], asign, a_exp[afb_n]] */
    } rel_t;

    int vlen = 1 + rfb_n + 1 + afb_n;
    rel_t *rels = malloc(MAX_RELS * sizeof(rel_t));
    for (int i = 0; i < MAX_RELS; i++) {
        mpz_init(rels[i].x_val);
        rels[i].ev = calloc(vlen, sizeof(int));
    }
    int nrels = 0;

    /* Main sieve loop: for each b, sieve over a */
    mpz_t rn, an, td, bz, az;
    mpz_inits(rn, an, td, bz, az, NULL);

    for (long b = 1; b <= 5000 && nrels < target_rels; b++) {
        if (wall() - start > 250.0) break;

        memset(rsieve, 0, sieve_len * sizeof(float));
        memset(asieve, 0, sieve_len * sizeof(float));

        /* Rational sieve: p | (a - b*m) when a ≡ b*m (mod p) */
        for (int i = 0; i < rfb_n; i++) {
            int p = rfb[i].p;
            long bm = ((b % p) * mpz_fdiv_ui(gm, p)) % p;
            long s = ((bm + sieve_half) % p + p) % p;
            for (long j = s; j < sieve_len; j += p)
                rsieve[j] += rlog[i];
        }

        /* Algebraic sieve: p | f_hom(a,b) when a ≡ b*r (mod p) */
        for (int i = 0; i < afb_n; i++) {
            int p = afb[i].p;
            long br = ((b % p) * afb[i].r) % p;
            long s = ((br + sieve_half) % p + p) % p;
            for (long j = s; j < sieve_len; j += p)
                asieve[j] += alog[i];
        }

        /* Estimate max norms for threshold */
        double max_a = sieve_half;
        double rne = log(max_a) + mpz_sizeinbase(gm, 2) * log(2.0); /* log(a*m) ≈ log(rat_norm) */
        double ane = gdeg * log(max_a) + mpz_sizeinbase(co[gdeg], 2) * log(2.0); /* log(alg_norm) */
        float rthr = (float)(rne * 0.25); /* aggressive threshold */
        float athr = (float)(ane * 0.25);

        /* Scan for candidates passing both thresholds */
        for (int idx = 0; idx < sieve_len && nrels < target_rels; idx++) {
            if (rsieve[idx] < rthr || asieve[idx] < athr) continue;

            long a = (long)idx - sieve_half;
            if (a == 0) continue;
            /* gcd(a, b) must be 1 */
            { long g = a < 0 ? -a : a, t = b;
              while (t) { long r = g % t; g = t; t = r; }
              if (g != 1) continue; }

            /* Trial divide rational norm */
            mpz_set_si(az, a);
            mpz_set_si(bz, b);
            mpz_mul(rn, bz, gm);
            mpz_sub(rn, az, rn); /* a - b*m */
            mpz_abs(td, rn);

            int *ev = rels[nrels].ev;
            memset(ev, 0, vlen * sizeof(int));
            ev[0] = (mpz_sgn(rn) < 0) ? 1 : 0;

            int r_smooth = 1;
            for (int i = 0; i < rfb_n; i++) {
                while (mpz_divisible_ui_p(td, rfb[i].p)) {
                    mpz_divexact_ui(td, td, rfb[i].p);
                    ev[1 + i]++;
                }
            }
            if (mpz_cmp_ui(td, 1) != 0) continue; /* rational side not smooth */

            /* Trial divide algebraic norm */
            poly_hom(an, a, b);
            mpz_abs(td, an);
            ev[1 + rfb_n] = (mpz_sgn(an) < 0) ? 1 : 0;

            for (int i = 0; i < afb_n; i++) {
                while (mpz_divisible_ui_p(td, afb[i].p)) {
                    mpz_divexact_ui(td, td, afb[i].p);
                    ev[1 + rfb_n + 1 + i]++;
                }
            }
            if (mpz_cmp_ui(td, 1) != 0) continue; /* algebraic side not smooth */

            /* Both smooth! Record relation. */
            mpz_set(rels[nrels].x_val, rn); /* rational value a - b*m */
            nrels++;

            if (nrels % 20 == 0) {
                fprintf(stderr, "HQNFS: [%.1fs] %d/%d rels (b=%ld)\n",
                        wall() - start, nrels, target_rels, b);
            }
        }
    }

    fprintf(stderr, "HQNFS: collected %d relations (need %d)\n", nrels, target_rels);

    /* GF(2) linear algebra */
    if (nrels >= rfb_n + afb_n + 2) {
        int m = vlen, n = nrels;
        int words = (m + n + 63) / 64;
        unsigned long *mat = calloc((size_t)n * words, sizeof(unsigned long));

        #define MG(r,c) ((mat[(size_t)(r)*words+(c)/64]>>((c)%64))&1)
        #define MS(r,c) (mat[(size_t)(r)*words+(c)/64]|=(1UL<<((c)%64)))
        #define MX(d,s) do{for(int _w=0;_w<words;_w++) mat[(size_t)(d)*words+_w]^=mat[(size_t)(s)*words+_w];}while(0)

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++)
                if (rels[i].ev[j] & 1) MS(i, j);
            MS(i, m + i);
        }

        int rank = 0;
        for (int j = 0; j < m && rank < n; j++) {
            int piv = -1;
            for (int i = rank; i < n; i++) if (MG(i,j)) { piv = i; break; }
            if (piv < 0) continue;
            if (piv != rank) {
                for (int w = 0; w < words; w++) {
                    unsigned long t = mat[(size_t)rank*words+w];
                    mat[(size_t)rank*words+w] = mat[(size_t)piv*words+w];
                    mat[(size_t)piv*words+w] = t;
                }
            }
            for (int i = 0; i < n; i++)
                if (i != rank && MG(i,j)) MX(i, rank);
            rank++;
        }

        fprintf(stderr, "HQNFS: rank=%d, %d null vectors\n", rank, n - rank);

        /* Extract factor from null vectors */
        mpz_t X, Y, g, tmp;
        mpz_inits(X, Y, g, tmp, NULL);
        int found = 0;

        for (int row = rank; row < n && !found; row++) {
            mpz_set_ui(X, 1);
            unsigned long *tot = calloc(vlen, sizeof(unsigned long));
            int cnt = 0;

            for (int i = 0; i < n; i++) {
                if (!MG(row, m + i)) continue;
                cnt++;
                mpz_mul(X, X, rels[i].x_val);
                mpz_mod(X, X, gN);
                for (int k = 0; k < vlen; k++)
                    tot[k] += rels[i].ev[k];
            }
            if (cnt == 0) { free(tot); continue; }

            /* Y = product of rational primes^(half_exp) mod N */
            mpz_set_ui(Y, 1);
            for (int k = 0; k < rfb_n; k++) {
                if (tot[1 + k] == 0) continue;
                unsigned long half = tot[1 + k] / 2;
                if (half == 0) continue;
                mpz_set_ui(tmp, rfb[k].p);
                mpz_powm_ui(tmp, tmp, half, gN);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, gN);
            }

            /* Note: for FULL NFS, we'd also need the algebraic square root.
             * This simplified version only uses the rational side.
             * It works when the rational exponents alone form a dependency. */

            /* Check rational exponents are all even */
            int ok = 1;
            if (tot[0] & 1) ok = 0;
            for (int k = 0; k < rfb_n; k++)
                if (tot[1 + k] & 1) { ok = 0; break; }
            if (!ok) { free(tot); continue; }

            mpz_sub(tmp, X, Y);
            mpz_gcd(g, tmp, gN);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
                mpz_t cof; mpz_init(cof);
                mpz_divexact(cof, gN, g);
                if (mpz_cmp(g, cof) > 0) mpz_set(g, cof);
                gmp_printf("FACTOR: %Zd\n", g);
                fprintf(stderr, "HQNFS: success in %.2fs\n", wall() - start);
                found = 1;
                mpz_clear(cof);
            }
            if (!found) {
                mpz_add(tmp, X, Y);
                mpz_gcd(g, tmp, gN);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
                    mpz_t cof; mpz_init(cof);
                    mpz_divexact(cof, gN, g);
                    if (mpz_cmp(g, cof) > 0) mpz_set(g, cof);
                    gmp_printf("FACTOR: %Zd\n", g);
                    fprintf(stderr, "HQNFS: success in %.2fs\n", wall() - start);
                    found = 1;
                    mpz_clear(cof);
                }
            }
            free(tot);
        }

        mpz_clears(X, Y, g, tmp, NULL);
        free(mat);

        if (found) goto cleanup;

        #undef MG
        #undef MS
        #undef MX
    }

    /* ECM fallback */
    fprintf(stderr, "HQNFS: ECM fallback\n");
    {
        ecm_params params; ecm_init(params);
        params->param = ECM_PARAM_SUYAMA;
        mpz_t factor; mpz_init(factor);
        for (unsigned long b1 = 2000; b1 < 300000000; b1 = (unsigned long)(b1 * 2.5)) {
            for (int c = 0; c < 200; c++) {
                if (wall() - start > 290.0) goto ecm_done;
                mpz_set_ui(factor, 0);
                params->B1done = 1.0;
                mpz_set_ui(params->sigma, 42 + (unsigned long)(b1/1000)*1000 + c);
                if (ecm_factor(factor, gN, (double)b1, params) > 0 &&
                    mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, gN) < 0) {
                    mpz_t cof; mpz_init(cof);
                    mpz_divexact(cof, gN, factor);
                    if (mpz_cmp(factor, cof) > 0) mpz_set(factor, cof);
                    gmp_printf("FACTOR: %Zd\n", factor);
                    mpz_clear(cof); ecm_clear(params); mpz_clear(factor);
                    goto cleanup;
                }
            }
        }
        ecm_done:
        ecm_clear(params); mpz_clear(factor);
    }

    fprintf(stderr, "HQNFS: failed\n");

cleanup:
    for (int i = 0; i < MAX_RELS; i++) { mpz_clear(rels[i].x_val); free(rels[i].ev); }
    free(rels); free(rfb); free(afb); free(rsieve); free(asieve);
    free(rlog); free(alog);
    mpz_clears(rn, an, td, bz, az, gN, gm, NULL);
    for (int i = 0; i <= gdeg; i++) mpz_clear(co[i]);
    return 0;
}
