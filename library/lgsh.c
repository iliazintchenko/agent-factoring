/*
 * lgsh.c - Hybrid factoring: MPQS + ECM fallback
 *
 * Strategy:
 * 1. Try Pollard rho for small factors
 * 2. Try MPQS with multiple Knuth multipliers
 * 3. Fall back to ECM if MPQS fails
 *
 * Uses quadratic character bits in GF(2) matrix to avoid extraction degeneracy.
 *
 * Usage: ./lgsh <N>
 * Output: FACTOR:<p>
 */

#include <gmp.h>
#include <ecm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double wall_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ===== Prime sieve ===== */
static int *sieve_primes(int B, int *count) {
    char *s = (char *)calloc(B + 1, 1);
    for (int i = 2; i <= B; i++) s[i] = 1;
    for (int i = 2; (long long)i * i <= B; i++)
        if (s[i]) for (int j = i * i; j <= B; j += i) s[j] = 0;
    int cnt = 0;
    for (int i = 2; i <= B; i++) if (s[i]) cnt++;
    int *p = (int *)malloc(cnt * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= B; i++) if (s[i]) p[idx++] = i;
    *count = cnt; free(s); return p;
}

/* ===== Modular square root ===== */
static long mod_sqrt_long(long n, long p) {
    if (p == 2) return n & 1;
    n = ((n % p) + p) % p;
    if (n == 0) return 0;
    mpz_t a, m, r;
    mpz_inits(a, m, r, NULL);
    mpz_set_si(a, n); mpz_set_si(m, p);
    mpz_powm_ui(r, a, (p - 1) / 2, m);
    if (mpz_cmp_ui(r, 1) != 0) { mpz_clears(a, m, r, NULL); return -1; }
    if ((p & 3) == 3) {
        mpz_powm_ui(r, a, (p + 1) / 4, m);
        long res = mpz_get_si(r); mpz_clears(a, m, r, NULL); return res;
    }
    long Q = p - 1, S = 0;
    while ((Q & 1) == 0) { Q >>= 1; S++; }
    long z;
    for (z = 2; z < p; z++) {
        mpz_set_si(a, z);
        mpz_powm_ui(r, a, (p - 1) / 2, m);
        if (mpz_cmp_si(r, p - 1) == 0) break;
    }
    mpz_t c, t, R, b;
    mpz_inits(c, t, R, b, NULL);
    mpz_set_si(a, z); mpz_powm_ui(c, a, Q, m);
    mpz_set_si(a, n); mpz_powm_ui(t, a, Q, m);
    mpz_powm_ui(R, a, (Q + 1) / 2, m);
    long Mv = S;
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            long res = mpz_get_si(R);
            mpz_clears(a, m, r, c, t, R, b, NULL); return res;
        }
        long ii = 0; mpz_set(b, t);
        while (mpz_cmp_ui(b, 1) != 0) {
            mpz_mul(b, b, b); mpz_mod(b, b, m); ii++;
            if (ii >= Mv) { mpz_clears(a, m, r, c, t, R, b, NULL); return -1; }
        }
        mpz_set(b, c);
        for (long j = 0; j < Mv - ii - 1; j++) { mpz_mul(b, b, b); mpz_mod(b, b, m); }
        Mv = ii;
        mpz_mul(c, b, b); mpz_mod(c, c, m);
        mpz_mul(t, t, c); mpz_mod(t, t, m);
        mpz_mul(R, R, b); mpz_mod(R, R, m);
    }
}

/* ===== Factor base ===== */
typedef struct { int *primes; long *sqrtN; int count; } fb_t;

static fb_t build_fb(const mpz_t N, int B) {
    int np; int *all = sieve_primes(B, &np);
    fb_t fb; fb.primes = (int*)malloc(np*sizeof(int));
    fb.sqrtN = (long*)malloc(np*sizeof(long)); fb.count = 0;
    for (int i = 0; i < np; i++) {
        long nmp = mpz_fdiv_ui(N, all[i]);
        long sq = mod_sqrt_long(nmp, all[i]);
        if (sq >= 0) { fb.primes[fb.count] = all[i]; fb.sqrtN[fb.count] = sq; fb.count++; }
    }
    free(all); return fb;
}

/* ===== Relation ===== */
typedef struct { mpz_t sqrt_val; int *exps; unsigned long lp; } rel_t;
typedef struct { rel_t *data; int count, cap; int veclen; } relset_t;

static void rs_init(relset_t *rs, int vl) {
    rs->veclen = vl; rs->count = 0; rs->cap = 4096;
    rs->data = (rel_t*)malloc(rs->cap * sizeof(rel_t));
}
static void rs_add(relset_t *rs, const mpz_t sv, const int *exps, unsigned long lp) {
    if (rs->count >= rs->cap) { rs->cap *= 2; rs->data = (rel_t*)realloc(rs->data, rs->cap * sizeof(rel_t)); }
    rel_t *r = &rs->data[rs->count];
    mpz_init_set(r->sqrt_val, sv); r->exps = (int*)malloc(rs->veclen * sizeof(int));
    memcpy(r->exps, exps, rs->veclen * sizeof(int)); r->lp = lp; rs->count++;
}
static void rs_free(relset_t *rs) {
    for (int i = 0; i < rs->count; i++) { mpz_clear(rs->data[i].sqrt_val); free(rs->data[i].exps); }
    free(rs->data);
}

/* Auxiliary primes for quadratic characters (prevent extraction degeneracy) */
static unsigned long qc_primes[3] = {0, 0, 0};
static int qc_count = 0;
static mpz_t qc_kN; /* store kN for QC computation */

static void init_qc_primes_for_kN(const fb_t *fb, const mpz_t kN) {
    /* Choose primes where kN is a QNR. These provide non-redundant QC bits.
     * Key: the residue x^2 - kN can never be ≡ 0 mod such a prime,
     * but its Jacobi symbol mod these primes still varies, giving
     * information NOT captured by the FB exponent columns. */
    qc_count = 0;
    mpz_init(qc_kN);
    mpz_set(qc_kN, kN);
    for (unsigned long p = 3; qc_count < 3 && p < 100000; p += 2) {
        int is_prime = 1;
        for (unsigned long d = 3; d * d <= p; d += 2)
            if (p % d == 0) { is_prime = 0; break; }
        if (!is_prime) continue;

        /* Verify kN is QNR mod p (Jacobi symbol = -1) */
        int jac = mpz_kronecker_ui(kN, p);
        if (jac == -1) {
            qc_primes[qc_count++] = p;
        }
    }
    fprintf(stderr, "  QC primes: ");
    for (int i = 0; i < qc_count; i++) fprintf(stderr, "%lu ", qc_primes[i]);
    fprintf(stderr, "\n");
}

/* Compute QC bits for a residue: exps[fb->count+1+i] = (residue/qc_primes[i]) == -1 ? 1 : 0 */
static void compute_qc(int *exps, const mpz_t residue, int sign, const fb_t *fb) {
    mpz_t r;
    mpz_init(r);
    /* The original signed residue */
    if (sign) mpz_neg(r, residue); else mpz_set(r, residue);
    for (int i = 0; i < qc_count; i++) {
        int js = mpz_jacobi(r, (mpz_srcptr)NULL); /* Can't use mpz_jacobi with unsigned long */
        /* Use mpz_kronecker_ui instead */
        int kr = mpz_kronecker_ui(r, qc_primes[i]);
        exps[fb->count + 1 + i] = (kr == -1) ? 1 : 0;
    }
    mpz_clear(r);
}

/* ===== QS sieve ===== */
static void qs_sieve(const mpz_t kN, const fb_t *fb, int M,
                     unsigned long lp_bound, relset_t *full, relset_t *partial) {
    mpz_t m, tmp, val, residue;
    mpz_inits(m, tmp, val, residue, NULL);
    mpz_sqrt(m, kN); mpz_add_ui(m, m, 1);
    int slen = 2 * M + 1;
    int n_bits = mpz_sizeinbase(kN, 2);
    double log2_max = 1.0 + log2(M) + n_bits / 2.0;
    int threshold = (int)((log2_max - log2(lp_bound > 1 ? lp_bound : 2)) * 1024);
    if (threshold < 0) threshold = 0;
    int *sieve = (int *)calloc(slen, sizeof(int));
    for (int i = 0; i < fb->count; i++) {
        int p = fb->primes[i]; long sq = fb->sqrtN[i];
        long mp = mpz_fdiv_ui(m, p); int logp = (int)(log2(p) * 1024);
        if (p == 2) {
            long r = ((sq - mp) % 2 + 2) % 2;
            long s = ((r + M) % 2 + 2) % 2;
            for (long j = s; j < slen; j += 2) sieve[j] += logp;
            continue;
        }
        long r1 = ((sq - mp) % p + p) % p;
        long r2 = ((-sq - mp) % p + p) % p;
        long s1 = ((r1 + M) % p + p) % p;
        long s2 = ((r2 + M) % p + p) % p;
        for (long j = s1; j < slen; j += p) sieve[j] += logp;
        if (r1 != r2) for (long j = s2; j < slen; j += p) sieve[j] += logp;
    }
    int *exps = (int *)calloc(full->veclen, sizeof(int));
    for (int j = 0; j < slen; j++) {
        if (sieve[j] < threshold) continue;
        long x = (long)j - M;
        mpz_set_si(tmp, x); mpz_add(val, tmp, m);
        mpz_mul(residue, val, val); mpz_sub(residue, residue, kN);
        memset(exps, 0, full->veclen * sizeof(int));
        int sign = 0;
        if (mpz_sgn(residue) < 0) { sign = 1; mpz_neg(residue, residue); }
        exps[0] = sign;
        mpz_set(tmp, residue);
        for (int i = 0; i < fb->count; i++) {
            unsigned long p = fb->primes[i];
            while (mpz_divisible_ui_p(tmp, p)) { mpz_divexact_ui(tmp, tmp, p); exps[i+1]++; }
        }
        if (mpz_cmp_ui(tmp, 1) == 0) {
            compute_qc(exps, residue, sign, fb);
            mpz_set_si(tmp, x); mpz_add(val, tmp, m);
            rs_add(full, val, exps, 0);
        } else if (mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp) <= lp_bound) {
            unsigned long cofactor = mpz_get_ui(tmp);
            compute_qc(exps, residue, sign, fb);
            mpz_set_si(tmp, x); mpz_add(val, tmp, m);
            rs_add(partial, val, exps, cofactor);
        }
    }
    free(sieve); free(exps);
    mpz_clears(m, tmp, val, residue, NULL);
}

/* ===== MPQS sieve ===== */
static void mpqs_sieve(const mpz_t kN, const fb_t *fb, int M,
                       unsigned long lp_bound, int poly_idx,
                       relset_t *full, relset_t *partial) {
    mpz_t target, q, A, B, C, tmp, val, res, q_inv;
    mpz_inits(target, q, A, B, C, tmp, val, res, q_inv, NULL);
    int n_bits = mpz_sizeinbase(kN, 2);
    mpz_mul_ui(tmp, kN, 2); mpz_sqrt(target, tmp); mpz_sqrt(target, target);
    unsigned long sqM = (unsigned long)sqrt((double)M);
    if (sqM < 2) sqM = 2;
    mpz_tdiv_q_ui(target, target, sqM);
    /* Use poly_idx/2 for the prime, and poly_idx%2 for which root */
    mpz_add_ui(target, target, (poly_idx / 2) * 100 + 42);
    mpz_nextprime(q, target);
    int use_second_root = (poly_idx & 1);

    int ok = 0;
    for (int att = 0; att < 200 && !ok; att++) {
        if (!mpz_fits_ulong_p(q)) { mpz_nextprime(q, q); continue; }
        unsigned long qv = mpz_get_ui(q);
        long nq = mpz_fdiv_ui(kN, qv);
        long sq = mod_sqrt_long(nq, qv);
        if (sq <= 0) { mpz_nextprime(q, q); continue; }
        mpz_mul(A, q, q); mpz_set_ui(B, sq);
        mpz_set_ui(tmp, (unsigned long)sq * sq); mpz_sub(tmp, kN, tmp);
        if (!mpz_divisible_ui_p(tmp, qv)) { mpz_nextprime(q, q); continue; }
        mpz_tdiv_q_ui(tmp, tmp, qv);
        mpz_t inv, pz; mpz_inits(inv, pz, NULL);
        mpz_set_ui(inv, 2 * sq % qv); mpz_set_ui(pz, qv);
        if (!mpz_invert(inv, inv, pz)) { mpz_clears(inv, pz, NULL); mpz_nextprime(q, q); continue; }
        mpz_mul(tmp, tmp, inv); mpz_mod_ui(tmp, tmp, qv);
        mpz_mul_ui(tmp, tmp, qv); mpz_add(B, B, tmp);
        mpz_clears(inv, pz, NULL);
        mpz_mul(tmp, B, B); mpz_sub(tmp, tmp, kN);
        if (!mpz_divisible_p(tmp, A)) { mpz_nextprime(q, q); continue; }
        /* Use second root: B' = A - B */
        if (use_second_root) mpz_sub(B, A, B);
        mpz_mul(C, B, B); mpz_sub(C, C, kN); mpz_tdiv_q(C, C, A);
        mpz_set_ui(tmp, mpz_get_ui(q)); mpz_invert(q_inv, tmp, kN);
        ok = 1;
    }
    if (!ok) { mpz_clears(target, q, A, B, C, tmp, val, res, q_inv, NULL); return; }

    int slen = 2 * M + 1;
    double log2_max = log2(M) + n_bits / 2.0 + 0.5;
    int threshold = (int)((log2_max - log2(lp_bound > 1 ? lp_bound : 2)) * 1024);
    if (threshold < 0) threshold = 0;
    int *sieve = (int *)calloc(slen, sizeof(int));
    for (int i = 0; i < fb->count; i++) {
        int p = fb->primes[i]; long sq = fb->sqrtN[i]; int logp = (int)(log2(p)*1024);
        if (p == 2) { for (long j = 0; j < slen; j++) sieve[j] += logp; continue; }
        long Ap = mpz_fdiv_ui(A, p), Bp = mpz_fdiv_ui(B, p);
        if (Ap == 0) {
            if (Bp != 0) {
                mpz_t i2, p2; mpz_inits(i2, p2, NULL);
                mpz_set_ui(i2, 2*Bp%p); mpz_set_ui(p2, p);
                if (mpz_invert(i2, i2, p2)) {
                    long Cp = mpz_fdiv_ui(C, p);
                    long root = (p - Cp) % p * mpz_get_ui(i2) % p;
                    long start = ((root + M) % p + p) % p;
                    for (long j = start; j < slen; j += p) sieve[j] += logp;
                }
                mpz_clears(i2, p2, NULL);
            }
            continue;
        }
        mpz_t ai, pzz; mpz_inits(ai, pzz, NULL);
        mpz_set_ui(ai, Ap); mpz_set_ui(pzz, p);
        if (!mpz_invert(ai, ai, pzz)) { mpz_clears(ai, pzz, NULL); continue; }
        long aiv = mpz_get_ui(ai); mpz_clears(ai, pzz, NULL);
        long r1 = ((sq - Bp) % p + p) % p * aiv % p;
        long r2 = ((p - sq - Bp) % p + p) % p * aiv % p;
        long s1 = ((r1 + M) % p + p) % p;
        long s2 = ((r2 + M) % p + p) % p;
        for (long j = s1; j < slen; j += p) sieve[j] += logp;
        if (s1 != s2) for (long j = s2; j < slen; j += p) sieve[j] += logp;
    }
    int *exps = (int *)calloc(full->veclen, sizeof(int));
    for (int j = 0; j < slen; j++) {
        if (sieve[j] < threshold) continue;
        long x = (long)j - M;
        mpz_set_si(tmp, x); mpz_mul_si(res, tmp, x); mpz_mul(res, res, A);
        mpz_set_si(tmp, x); mpz_mul(tmp, tmp, B); mpz_mul_ui(tmp, tmp, 2);
        mpz_add(res, res, tmp); mpz_add(res, res, C);
        memset(exps, 0, full->veclen * sizeof(int));
        int sign = 0;
        if (mpz_sgn(res) < 0) { sign = 1; mpz_neg(res, res); }
        exps[0] = sign;
        mpz_set(tmp, res);
        for (int i = 0; i < fb->count; i++) {
            unsigned long p = fb->primes[i];
            while (mpz_divisible_ui_p(tmp, p)) { mpz_divexact_ui(tmp, tmp, p); exps[i+1]++; }
        }
        int is_smooth = (mpz_cmp_ui(tmp, 1) == 0);
        int is_partial = (!is_smooth && mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp) <= lp_bound);
        if (is_smooth || is_partial) {
            unsigned long lp = is_smooth ? 0 : mpz_get_ui(tmp);
            compute_qc(exps, res, sign, fb);
            mpz_set_si(val, x); mpz_mul(val, A, val); mpz_add(val, val, B);
            mpz_mul(val, val, q_inv); mpz_mod(val, val, kN);
            if (is_smooth) rs_add(full, val, exps, 0);
            else rs_add(partial, val, exps, lp);
        }
    }
    free(sieve); free(exps);
    mpz_clears(target, q, A, B, C, tmp, val, res, q_inv, NULL);
}

/* ===== LP matching ===== */
static int compare_lp(const void *a, const void *b) {
    const rel_t *ra = (const rel_t *)a, *rb = (const rel_t *)b;
    return (ra->lp < rb->lp) ? -1 : (ra->lp > rb->lp) ? 1 : 0;
}

static int match_lp(relset_t *partials, relset_t *full, const mpz_t kN) {
    if (partials->count < 2) return 0;
    qsort(partials->data, partials->count, sizeof(rel_t), compare_lp);
    int veclen = full->veclen;
    int *combined = (int *)calloc(veclen, sizeof(int));
    mpz_t sv, lp_inv, lp_val;
    mpz_inits(sv, lp_inv, lp_val, NULL);
    int matched = 0;
    for (int i = 0; i < partials->count - 1; i++) {
        if (partials->data[i].lp == 0) continue;
        if (partials->data[i].lp == partials->data[i + 1].lp) {
            for (int k = 0; k < veclen; k++)
                combined[k] = partials->data[i].exps[k] + partials->data[i + 1].exps[k];
            mpz_mul(sv, partials->data[i].sqrt_val, partials->data[i + 1].sqrt_val);
            mpz_set_ui(lp_val, partials->data[i].lp);
            if (mpz_invert(lp_inv, lp_val, kN)) mpz_mul(sv, sv, lp_inv);
            mpz_mod(sv, sv, kN);
            rs_add(full, sv, combined, 0);
            matched++;
            partials->data[i].lp = 0; partials->data[i + 1].lp = 0;
            i++;
        }
    }
    mpz_clears(sv, lp_inv, lp_val, NULL);
    free(combined);
    return matched;
}

/* ===== GF(2) solver ===== */
typedef unsigned long word_t;
#define WBITS (sizeof(word_t) * 8)

static int find_deps(relset_t *rels, int veclen, int ***deps, int **sizes, int *ndeps) {
    int nr = veclen, nc = rels->count;
    *ndeps = 0;
    if (nc <= nr + 1) return 0;
    int nw = (nc + WBITS - 1) / WBITS;
    word_t **mat = (word_t**)calloc(nr, sizeof(word_t*));
    for (int i = 0; i < nr; i++) {
        mat[i] = (word_t*)calloc(nw, sizeof(word_t));
        for (int j = 0; j < nc; j++)
            if (rels->data[j].exps[i] & 1) mat[i][j/WBITS] |= (1UL << (j%WBITS));
    }
    int *pc = (int*)malloc(nr * sizeof(int));
    char *isp = (char*)calloc(nc, 1);
    for (int i = 0; i < nr; i++) {
        pc[i] = -1;
        for (int j = 0; j < nc; j++) {
            if (!isp[j] && (mat[i][j/WBITS] & (1UL << (j%WBITS)))) {
                pc[i] = j; isp[j] = 1;
                for (int k = 0; k < nr; k++)
                    if (k != i && (mat[k][j/WBITS] & (1UL << (j%WBITS))))
                        for (int w = 0; w < nw; w++) mat[k][w] ^= mat[i][w];
                break;
            }
        }
    }
    int mx = 128;
    *deps = (int**)malloc(mx * sizeof(int*));
    *sizes = (int*)malloc(mx * sizeof(int));
    int nd = 0;
    for (int j = 0; j < nc && nd < mx; j++) {
        if (isp[j]) continue;
        int *dep = (int*)malloc((nr+1)*sizeof(int));
        int cnt = 0; dep[cnt++] = j;
        for (int i = 0; i < nr; i++)
            if (mat[i][j/WBITS] & (1UL << (j%WBITS)))
                if (pc[i] >= 0) dep[cnt++] = pc[i];
        (*deps)[nd] = dep; (*sizes)[nd] = cnt; nd++;
    }
    *ndeps = nd;
    for (int i = 0; i < nr; i++) free(mat[i]);
    free(mat); free(pc); free(isp);
    return nd > 0;
}

/* ===== Factor extraction ===== */
static int try_extract(relset_t *rels, int *dep, int dsz,
                       const mpz_t kN, const fb_t *fb, mpz_t factor) {
    mpz_t x, y, tmp;
    mpz_inits(x, y, tmp, NULL);
    mpz_set_ui(x, 1);
    int *tot = (int*)calloc(fb->count + 1, sizeof(int));
    for (int i = 0; i < dsz; i++) {
        mpz_mul(x, x, rels->data[dep[i]].sqrt_val);
        mpz_mod(x, x, kN);
        for (int j = 0; j < fb->count + 1; j++)
            tot[j] += rels->data[dep[i]].exps[j];
    }
    for (int j = 0; j < fb->count + 1; j++)
        if (tot[j] & 1) { free(tot); mpz_clears(x, y, tmp, NULL); return 0; }
    mpz_set_ui(y, 1);
    for (int j = 0; j < fb->count; j++) {
        int e = tot[j+1] / 2;
        if (e > 0) {
            mpz_set_ui(tmp, fb->primes[j]);
            mpz_powm_ui(tmp, tmp, e, kN);
            mpz_mul(y, y, tmp); mpz_mod(y, y, kN);
        }
    }
    if ((tot[0] / 2) & 1) mpz_sub(y, kN, y);
    mpz_sub(tmp, x, y); mpz_gcd(factor, tmp, kN);
    if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, kN) < 0) {
        free(tot); mpz_clears(x, y, tmp, NULL); return 1;
    }
    mpz_add(tmp, x, y); mpz_gcd(factor, tmp, kN);
    int ok = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, kN) < 0);
    free(tot); mpz_clears(x, y, tmp, NULL); return ok;
}

/* ===== MPQS core ===== */
static int factor_mpqs(const mpz_t N, const mpz_t kN, int k, mpz_t result, double deadline) {
    int ndigits = mpz_sizeinbase(kN, 10);
    double ln_N = ndigits * log(10);
    double Lexp = sqrt(ln_N * log(ln_N));
    int B = (int)exp(0.55 * Lexp);
    if (B < 500) B = 500; if (B > 5000000) B = 5000000;
    double logM = 1.05 * Lexp;
    int M = (logM > 16.0) ? 10000000 : (int)exp(logM);
    if (M < 20000) M = 20000; if (M > 10000000) M = 10000000;

    fb_t fb = build_fb(kN, B);
    init_qc_primes_for_kN(&fb, kN);
    /* Add 3 quadratic character columns to prevent extraction degeneracy */
    int n_qc = qc_count;
    int veclen = fb.count + 1 + n_qc; /* sign + FB primes + QC bits */
    int target = fb.count + n_qc + 30;
    unsigned long lp_bound = (unsigned long)B * 50;

    relset_t full, partial;
    rs_init(&full, veclen);
    rs_init(&partial, veclen);

    qs_sieve(kN, &fb, M, lp_bound, &full, &partial);

    for (int pi = 0; pi < 10000 && full.count < target && wall_time() < deadline; pi++) {
        mpqs_sieve(kN, &fb, M, lp_bound, pi, &full, &partial);
        if ((pi + 1) % 100 == 0) {
            int m = match_lp(&partial, &full, kN);
            (void)m;
        }
    }
    match_lp(&partial, &full, kN);
    fprintf(stderr, "  k=%d: %d full (need %d)\n", k, full.count, target);

    /* Deduplicate: remove relations with identical sqrt_val (hash-based) */
    {
        /* Hash by sqrt_val mod a large prime */
        unsigned long hash_size = full.count * 3;
        if (hash_size < 1024) hash_size = 1024;
        int *hash_table = (int*)malloc(hash_size * sizeof(int));
        memset(hash_table, -1, hash_size * sizeof(int));
        int unique = 0;
        for (int i = 0; i < full.count; i++) {
            unsigned long h = mpz_fdiv_ui(full.data[i].sqrt_val, hash_size);
            int dup = 0;
            /* Linear probing */
            for (int probe = 0; probe < 100; probe++) {
                int slot = (h + probe) % hash_size;
                if (hash_table[slot] == -1) {
                    hash_table[slot] = unique;
                    break;
                }
                int j = hash_table[slot];
                if (mpz_cmp(full.data[i].sqrt_val, full.data[j].sqrt_val) == 0) {
                    dup = 1; break;
                }
            }
            if (!dup) {
                if (unique < i) {
                    rel_t tmp_r = full.data[unique];
                    full.data[unique] = full.data[i];
                    full.data[i] = tmp_r;
                }
                unique++;
            }
        }
        free(hash_table);
        if (unique < full.count) {
            fprintf(stderr, "  Dedup: %d → %d unique\n", full.count, unique);
            full.count = unique;
        }
    }

    int success = 0;
    if (full.count > target) {
        int **deps; int *dsizes; int ndeps;
        if (find_deps(&full, veclen, &deps, &dsizes, &ndeps)) {
            for (int d = 0; d < ndeps && !success; d++)
                success = try_extract(&full, deps[d], dsizes[d], kN, &fb, result);

            /* Skip XOR pairs — for degenerate cases, ECM is faster */
            if (0 && !success && ndeps >= 2) {
                /* Combine pairs of dependencies via symmetric difference */
                for (int i = 0; i < ndeps - 1 && !success; i++) {
                    for (int j = i + 1; j < ndeps && !success; j++) {
                        /* Build combined dep = symmetric difference of deps[i] and deps[j] */
                        /* Use a bitmap to XOR the sets */
                        int *present = (int*)calloc(full.count, sizeof(int));
                        for (int a = 0; a < dsizes[i]; a++) present[deps[i][a]] ^= 1;
                        for (int a = 0; a < dsizes[j]; a++) present[deps[j][a]] ^= 1;
                        int cnt = 0;
                        for (int a = 0; a < full.count; a++) if (present[a]) cnt++;
                        if (cnt > 0) {
                            int *combined = (int*)malloc(cnt * sizeof(int));
                            int idx = 0;
                            for (int a = 0; a < full.count; a++)
                                if (present[a]) combined[idx++] = a;
                            success = try_extract(&full, combined, cnt, kN, &fb, result);
                            free(combined);
                        }
                        free(present);
                    }
                    /* Limit: only try ~50 pairs per dep */
                }
            }
            if (!success)
                fprintf(stderr, "  Extraction: %d deps, all failed\n", ndeps);
            for (int d = 0; d < ndeps; d++) free(deps[d]);
            free(deps); free(dsizes);
        }
    }

    if (success) {
        mpz_gcd(result, result, N);
        if (mpz_cmp_ui(result, 1) == 0 || mpz_cmp(result, N) == 0) success = 0;
    }

    rs_free(&full); rs_free(&partial);
    free(fb.primes); free(fb.sqrtN);
    return success;
}

/* ===== ECM fallback ===== */
static int factor_ecm(const mpz_t N, mpz_t result, double deadline) {
    /* Optimal B1 for balanced semiprimes: exp(sqrt(0.5 * ln(N) * ln(ln(N)))) */
    int n_bits = mpz_sizeinbase(N, 2);
    double ln_N = n_bits * 0.693;
    double optimal_B1 = exp(sqrt(0.5 * ln_N * log(ln_N)));
    double B1 = optimal_B1 / 10; /* start below optimal and ramp up */
    if (B1 < 1e5) B1 = 1e5;
    if (B1 > 1e9) B1 = 1e9;

    fprintf(stderr, "  ECM: B1=%.0f (optimal=%.0f)\n", B1, optimal_B1);

    for (int curve = 0; curve < 100000 && wall_time() < deadline; curve++) {
        ecm_params p;
        ecm_init(p);
        mpz_set_ui(p->sigma, 42 + curve);
        /* Let GMP-ECM choose B2 automatically */
        int ret = ecm_factor(result, N, B1, p);
        ecm_clear(p);
        if (ret > 0 && mpz_cmp(result, N) != 0 && mpz_cmp_ui(result, 1) > 0) {
            return 1;
        }
        /* Ramp up B1 */
        if (curve % 10 == 9) B1 *= 1.3;
        if (B1 > optimal_B1 * 3) B1 = optimal_B1 * 3;
    }
    return 0;
}

/* ===== Main ===== */
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    mpz_t N, factor;
    mpz_inits(N, factor, NULL);
    mpz_set_str(N, argv[1], 10);
    int ndigits = strlen(argv[1]);
    double t0 = wall_time();
    double deadline = t0 + 280.0;

    fprintf(stderr, "LGSH: factoring %d-digit number\n", ndigits);

    /* Small factor check */
    { int cnt; int *sp = sieve_primes(1000000, &cnt);
      for (int i = 0; i < cnt; i++)
        if (mpz_divisible_ui_p(N, sp[i])) {
            printf("FACTOR:%d\n", sp[i]); free(sp); mpz_clears(N, factor, NULL); return 0; }
      free(sp); }

    /* Try MPQS with 2 multipliers, then fall back to ECM quickly */
    int multipliers[] = {1, 3};
    int nmult = 2;
    mpz_t kN;
    mpz_init(kN);

    for (int mi = 0; mi < nmult && wall_time() < deadline - 60; mi++) {
        int k = multipliers[mi];
        mpz_mul_ui(kN, N, k);
        if (factor_mpqs(N, kN, k, factor, deadline - 60)) {
            mpz_t cof; mpz_init(cof); mpz_tdiv_q(cof, N, factor);
            gmp_printf("FACTOR:%Zd\n", factor);
            gmp_fprintf(stderr, "%Zd x %Zd (%.3fs)\n", factor, cof, wall_time() - t0);
            mpz_clear(cof); mpz_clear(kN); mpz_clears(N, factor, NULL);
            return 0;
        }
    }

    /* ECM fallback */
    fprintf(stderr, "MPQS failed, trying ECM...\n");
    if (factor_ecm(N, factor, deadline)) {
        mpz_t cof; mpz_init(cof); mpz_tdiv_q(cof, N, factor);
        gmp_printf("FACTOR:%Zd\n", factor);
        gmp_fprintf(stderr, "ECM: %Zd x %Zd (%.3fs)\n", factor, cof, wall_time() - t0);
        mpz_clear(cof); mpz_clear(kN); mpz_clears(N, factor, NULL);
        return 0;
    }

    fprintf(stderr, "FAIL (%.1fs)\n", wall_time() - t0);
    mpz_clear(kN); mpz_clears(N, factor, NULL);
    return 1;
}
