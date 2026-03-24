/*
 * hsd.c - Hierarchical Smooth Decomposition
 *
 * Novel factoring approach:
 * 1. Generate QS-style candidates of size ~√N
 * 2. Trial divide by small primes to find "partially smooth" values
 *    where the cofactor C has size N^α for α < 1/2
 * 3. Collect cofactors. Build a relation graph among cofactors by finding
 *    pairs/triples whose products/ratios are smooth.
 * 4. Use ECM + batch GCD to find cofactor relations.
 * 5. Any cofactor relation, combined with the original partial relations,
 *    gives a full relation for the matrix.
 *
 * Hypothesis: by allowing larger cofactors (more partial relations),
 * we get a denser cofactor graph, potentially improving scaling.
 *
 * Usage: ./hsd <N>
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
    *count = cnt;
    free(s);
    return p;
}

/* ===== Modular sqrt ===== */
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
    long Mval = S;
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            long res = mpz_get_si(R);
            mpz_clears(a, m, r, c, t, R, b, NULL); return res;
        }
        long i = 0; mpz_set(b, t);
        while (mpz_cmp_ui(b, 1) != 0) {
            mpz_mul(b, b, b); mpz_mod(b, b, m); i++;
            if (i >= Mval) { mpz_clears(a, m, r, c, t, R, b, NULL); return -1; }
        }
        mpz_set(b, c);
        for (long j = 0; j < Mval - i - 1; j++) { mpz_mul(b, b, b); mpz_mod(b, b, m); }
        Mval = i;
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

/* ===== Partial relation: sqrt_val^2 ≡ sign * smooth_part * cofactor (mod N) ===== */
typedef struct {
    mpz_t sqrt_val;
    int *exps;        /* exponent vector [sign, p0, p1, ...], length fb_count+1 */
    mpz_t cofactor;   /* remaining part after trial division (1 if fully smooth) */
    int cofactor_bits;
} partial_rel_t;

typedef struct {
    partial_rel_t *data;
    int count, cap;
    int veclen;
} prel_set_t;

static void prs_init(prel_set_t *s, int veclen) {
    s->veclen = veclen; s->count = 0; s->cap = 8192;
    s->data = (partial_rel_t*)malloc(s->cap * sizeof(partial_rel_t));
}

static void prs_add(prel_set_t *s, const mpz_t sv, const int *exps, const mpz_t cof) {
    if (s->count >= s->cap) {
        s->cap *= 2;
        s->data = (partial_rel_t*)realloc(s->data, s->cap * sizeof(partial_rel_t));
    }
    partial_rel_t *r = &s->data[s->count];
    mpz_init_set(r->sqrt_val, sv);
    r->exps = (int*)malloc(s->veclen * sizeof(int));
    memcpy(r->exps, exps, s->veclen * sizeof(int));
    mpz_init_set(r->cofactor, cof);
    r->cofactor_bits = mpz_sizeinbase(cof, 2);
    s->count++;
}

static void prs_free(prel_set_t *s) {
    for (int i = 0; i < s->count; i++) {
        mpz_clear(s->data[i].sqrt_val);
        mpz_clear(s->data[i].cofactor);
        free(s->data[i].exps);
    }
    free(s->data);
}

/* ===== Full relation (for matrix) ===== */
typedef struct { mpz_t sqrt_val; int *exps; } frel_t;
typedef struct { frel_t *data; int count, cap; int veclen; } frel_set_t;

static void frs_init(frel_set_t *s, int veclen) {
    s->veclen = veclen; s->count = 0; s->cap = 4096;
    s->data = (frel_t*)malloc(s->cap * sizeof(frel_t));
}

static void frs_add(frel_set_t *s, const mpz_t sv, const int *exps) {
    if (s->count >= s->cap) { s->cap *= 2; s->data = (frel_t*)realloc(s->data, s->cap * sizeof(frel_t)); }
    frel_t *r = &s->data[s->count];
    mpz_init_set(r->sqrt_val, sv);
    r->exps = (int*)malloc(s->veclen * sizeof(int));
    memcpy(r->exps, exps, s->veclen * sizeof(int));
    s->count++;
}

static void frs_free(frel_set_t *s) {
    for (int i = 0; i < s->count; i++) { mpz_clear(s->data[i].sqrt_val); free(s->data[i].exps); }
    free(s->data);
}

/* ===== QS sieve producing partial relations ===== */
static void qs_sieve_partial(const mpz_t N, const fb_t *fb, int M,
                             int max_cofactor_bits,
                             prel_set_t *partials, frel_set_t *fulls) {
    mpz_t m, tmp, val, residue;
    mpz_inits(m, tmp, val, residue, NULL);
    mpz_sqrt(m, N); mpz_add_ui(m, m, 1);

    int sieve_len = 2 * M + 1;
    int n_bits = mpz_sizeinbase(N, 2);

    /* Sieve threshold: allow cofactors up to max_cofactor_bits */
    double log2_max = 1.0 + log2(M) + n_bits / 2.0;
    int threshold = (int)((log2_max - max_cofactor_bits) * 1024);
    if (threshold < 512) threshold = 512;

    int *sieve = (int *)calloc(sieve_len, sizeof(int));

    for (int i = 0; i < fb->count; i++) {
        int p = fb->primes[i];
        long sq = fb->sqrtN[i];
        long m_mod = mpz_fdiv_ui(m, p);
        int logp = (int)(log2(p) * 1024);
        if (p == 2) {
            long r = ((sq - m_mod) % 2 + 2) % 2;
            long start = ((r + M) % 2 + 2) % 2;
            for (long j = start; j < sieve_len; j += 2) sieve[j] += logp;
            continue;
        }
        long r1 = ((sq - m_mod) % p + p) % p;
        long r2 = ((-sq - m_mod) % p + p) % p;
        long s1 = ((r1 + M) % p + p) % p;
        long s2 = ((r2 + M) % p + p) % p;
        for (long j = s1; j < sieve_len; j += p) sieve[j] += logp;
        if (r1 != r2) for (long j = s2; j < sieve_len; j += p) sieve[j] += logp;
    }

    int *exps = (int*)calloc(fb->count + 1, sizeof(int));
    mpz_t cofactor;
    mpz_init(cofactor);

    for (int j = 0; j < sieve_len; j++) {
        if (sieve[j] < threshold) continue;
        long x = (long)j - M;

        mpz_set_si(tmp, x); mpz_add(val, tmp, m);
        mpz_mul(residue, val, val); mpz_sub(residue, residue, N);

        memset(exps, 0, (fb->count + 1) * sizeof(int));
        int sign = 0;
        if (mpz_sgn(residue) < 0) { sign = 1; mpz_neg(residue, residue); }
        exps[0] = sign;

        mpz_set(cofactor, residue);
        for (int i = 0; i < fb->count; i++) {
            unsigned long p = fb->primes[i];
            while (mpz_divisible_ui_p(cofactor, p)) {
                mpz_divexact_ui(cofactor, cofactor, p);
                exps[i + 1]++;
            }
        }

        int cof_bits = mpz_sizeinbase(cofactor, 2);
        if (cof_bits <= 1) {
            /* Fully smooth */
            mpz_set_si(tmp, x); mpz_add(val, tmp, m); mpz_mod(val, val, N);
            frs_add(fulls, val, exps);
        } else if (cof_bits <= max_cofactor_bits) {
            /* Partial: has a cofactor */
            mpz_set_si(tmp, x); mpz_add(val, tmp, m); mpz_mod(val, val, N);
            prs_add(partials, val, exps, cofactor);
        }
    }

    free(sieve); free(exps); mpz_clear(cofactor);
    mpz_clears(m, tmp, val, residue, NULL);
}

/* ===== MPQS sieve producing partial relations ===== */
static void mpqs_sieve_partial(const mpz_t N, const fb_t *fb, int M,
                               int max_cofactor_bits, int poly_idx,
                               prel_set_t *partials, frel_set_t *fulls) {
    mpz_t target, q, A, B, C, tmp, val, residue, cofactor, q_inv;
    mpz_inits(target, q, A, B, C, tmp, val, residue, cofactor, q_inv, NULL);

    int n_bits = mpz_sizeinbase(N, 2);
    mpz_mul_ui(tmp, N, 2); mpz_sqrt(target, tmp); mpz_sqrt(target, target);
    mpz_tdiv_q_ui(target, target, (unsigned long)sqrt(M > 1 ? M : 2));
    mpz_add_ui(target, target, poly_idx * 100);
    mpz_nextprime(q, target);

    int attempts = 0;
    int found_poly = 0;
    while (attempts < 100) {
        if (!mpz_fits_ulong_p(q)) { mpz_nextprime(q, q); attempts++; continue; }
        unsigned long qv = mpz_get_ui(q);
        long nmodq = mpz_fdiv_ui(N, qv);
        long sq = mod_sqrt_long(nmodq, qv);
        if (sq <= 0) { mpz_nextprime(q, q); attempts++; continue; }

        mpz_mul(A, q, q);
        mpz_set_ui(B, sq);
        mpz_set_ui(tmp, (unsigned long)sq * sq);
        mpz_sub(tmp, N, tmp);
        if (!mpz_divisible_ui_p(tmp, qv)) { mpz_nextprime(q, q); attempts++; continue; }
        mpz_tdiv_q_ui(tmp, tmp, qv);

        mpz_t inv, pz;
        mpz_inits(inv, pz, NULL);
        mpz_set_ui(inv, 2 * sq % qv);
        mpz_set_ui(pz, qv);
        if (!mpz_invert(inv, inv, pz)) {
            mpz_clears(inv, pz, NULL); mpz_nextprime(q, q); attempts++; continue;
        }
        mpz_mul(tmp, tmp, inv); mpz_mod_ui(tmp, tmp, qv);
        mpz_mul_ui(tmp, tmp, qv); mpz_add(B, B, tmp);
        mpz_clears(inv, pz, NULL);

        mpz_mul(tmp, B, B); mpz_sub(tmp, tmp, N);
        if (!mpz_divisible_p(tmp, A)) { mpz_nextprime(q, q); attempts++; continue; }

        mpz_mul(C, B, B); mpz_sub(C, C, N); mpz_tdiv_q(C, C, A);

        mpz_set_ui(tmp, mpz_get_ui(q));
        mpz_invert(q_inv, tmp, N);
        found_poly = 1;
        break;
    }

    if (!found_poly) {
        mpz_clears(target, q, A, B, C, tmp, val, residue, cofactor, q_inv, NULL);
        return;
    }

    int sieve_len = 2 * M + 1;
    double log2_max = log2(M) + n_bits / 2.0 + 0.5;
    int threshold = (int)((log2_max - max_cofactor_bits) * 1024);
    if (threshold < 512) threshold = 512;

    int *sieve = (int *)calloc(sieve_len, sizeof(int));

    for (int i = 0; i < fb->count; i++) {
        int p = fb->primes[i];
        long sq = fb->sqrtN[i];
        int logp = (int)(log2(p) * 1024);
        if (p == 2) {
            for (long j = 0; j < sieve_len; j++) sieve[j] += logp;
            continue;
        }
        long Ap = mpz_fdiv_ui(A, p);
        long Bp = mpz_fdiv_ui(B, p);
        if (Ap == 0) {
            if (Bp != 0) {
                mpz_t inv2, pz2;
                mpz_inits(inv2, pz2, NULL);
                mpz_set_ui(inv2, 2 * Bp % p); mpz_set_ui(pz2, p);
                if (mpz_invert(inv2, inv2, pz2)) {
                    long Cp = mpz_fdiv_ui(C, p);
                    long root = (p - Cp) % p * mpz_get_ui(inv2) % p;
                    long start = ((root + M) % p + p) % p;
                    for (long j = start; j < sieve_len; j += p) sieve[j] += logp;
                }
                mpz_clears(inv2, pz2, NULL);
            }
            continue;
        }
        mpz_t Ainv2, pz2;
        mpz_inits(Ainv2, pz2, NULL);
        mpz_set_ui(Ainv2, Ap); mpz_set_ui(pz2, p);
        if (!mpz_invert(Ainv2, Ainv2, pz2)) { mpz_clears(Ainv2, pz2, NULL); continue; }
        long Ainvp = mpz_get_ui(Ainv2);
        mpz_clears(Ainv2, pz2, NULL);

        long r1 = ((sq - Bp) % p + p) % p * Ainvp % p;
        long r2 = ((p - sq - Bp) % p + p) % p * Ainvp % p;
        long s1 = ((r1 + M) % p + p) % p;
        long s2 = ((r2 + M) % p + p) % p;
        for (long j = s1; j < sieve_len; j += p) sieve[j] += logp;
        if (s1 != s2) for (long j = s2; j < sieve_len; j += p) sieve[j] += logp;
    }

    int *exps = (int*)calloc(fb->count + 1, sizeof(int));

    for (int j = 0; j < sieve_len; j++) {
        if (sieve[j] < threshold) continue;
        long x = (long)j - M;

        /* g(x) = Ax^2 + 2Bx + C */
        mpz_set_si(tmp, x); mpz_mul_si(residue, tmp, x);
        mpz_mul(residue, residue, A);
        mpz_set_si(tmp, x); mpz_mul(tmp, tmp, B); mpz_mul_ui(tmp, tmp, 2);
        mpz_add(residue, residue, tmp); mpz_add(residue, residue, C);

        memset(exps, 0, (fb->count + 1) * sizeof(int));
        int sign = 0;
        if (mpz_sgn(residue) < 0) { sign = 1; mpz_neg(residue, residue); }
        exps[0] = sign;

        mpz_set(cofactor, residue);
        for (int i = 0; i < fb->count; i++) {
            unsigned long p = fb->primes[i];
            while (mpz_divisible_ui_p(cofactor, p)) {
                mpz_divexact_ui(cofactor, cofactor, p); exps[i + 1]++;
            }
        }

        int cof_bits = mpz_sizeinbase(cofactor, 2);

        /* sqrt side: (Ax+B) * q^{-1} mod N */
        mpz_set_si(tmp, x); mpz_mul(val, A, tmp); mpz_add(val, val, B);
        mpz_mul(val, val, q_inv); mpz_mod(val, val, N);

        if (cof_bits <= 1) {
            frs_add(fulls, val, exps);
        } else if (cof_bits <= max_cofactor_bits) {
            prs_add(partials, val, exps, cofactor);
        }
    }

    free(sieve); free(exps);
    mpz_clears(target, q, A, B, C, tmp, val, residue, cofactor, q_inv, NULL);
}

/* ===== Cofactor matching with hash table ===== */
/* Sort-based single LP matching + limited double LP */
static int cmp_partial_cofactor(const void *a, const void *b) {
    const partial_rel_t *pa = (const partial_rel_t *)a;
    const partial_rel_t *pb = (const partial_rel_t *)b;
    return mpz_cmp(pa->cofactor, pb->cofactor);
}

static int match_cofactors(prel_set_t *partials, frel_set_t *fulls,
                           const fb_t *fb, unsigned long lp_bound) {
    int matched = 0;
    int veclen = fulls->veclen;
    int n = partials->count;
    if (n < 2) return 0;

    /* Sort by cofactor for O(n log n) single LP matching */
    qsort(partials->data, n, sizeof(partial_rel_t), cmp_partial_cofactor);

    mpz_t combined_sv;
    mpz_init(combined_sv);
    int *combined_exps = (int*)calloc(veclen, sizeof(int));

    /* Single LP: adjacent equal cofactors */
    for (int i = 0; i < n - 1; i++) {
        if (mpz_cmp_ui(partials->data[i].cofactor, 1) <= 0) continue;
        if (mpz_cmp(partials->data[i].cofactor, partials->data[i+1].cofactor) == 0) {
            for (int k = 0; k < veclen; k++)
                combined_exps[k] = partials->data[i].exps[k] + partials->data[i+1].exps[k];
            mpz_mul(combined_sv, partials->data[i].sqrt_val, partials->data[i+1].sqrt_val);
            frs_add(fulls, combined_sv, combined_exps);
            matched++;
            mpz_set_ui(partials->data[i].cofactor, 0);
            mpz_set_ui(partials->data[i+1].cofactor, 0);
            i++; /* skip partner */
        }
    }

    mpz_clear(combined_sv);
    free(combined_exps);
    return matched;
}

/* ===== Cofactor cycle matching (novel: triple LP) ===== */
/* Find triples (i,j,k) where cofactor_i * cofactor_j * cofactor_k is a perfect square
   or where the product is smooth */
static int match_cofactor_triples(prel_set_t *partials, frel_set_t *fulls, const fb_t *fb) {
    int n = partials->count;
    if (n < 3) return 0;
    int matched = 0;
    int veclen = fulls->veclen;

    /* Only try triples among small cofactors */
    /* Collect indices of partials with cofactors < 2^30 */
    int *small_idx = (int*)malloc(n * sizeof(int));
    int nsmall = 0;
    for (int i = 0; i < n; i++) {
        if (mpz_cmp_ui(partials->data[i].cofactor, 0) > 0 &&
            partials->data[i].cofactor_bits <= 30) {
            small_idx[nsmall++] = i;
        }
    }

    if (nsmall < 3) { free(small_idx); return 0; }

    /* Limit the search to avoid O(n^3) blowup */
    int max_triples = 100000;
    int count = 0;

    mpz_t prod, combined_sv;
    mpz_inits(prod, combined_sv, NULL);
    int *combined_exps = (int*)calloc(veclen, sizeof(int));

    for (int a = 0; a < nsmall - 2 && count < max_triples; a++) {
        int i = small_idx[a];
        if (mpz_cmp_ui(partials->data[i].cofactor, 0) <= 0) continue;
        for (int b = a + 1; b < nsmall - 1 && count < max_triples; b++) {
            int j = small_idx[b];
            if (mpz_cmp_ui(partials->data[j].cofactor, 0) <= 0) continue;
            for (int c = b + 1; c < nsmall && count < max_triples; c++) {
                int k = small_idx[c];
                if (mpz_cmp_ui(partials->data[k].cofactor, 0) <= 0) continue;
                count++;

                mpz_mul(prod, partials->data[i].cofactor, partials->data[j].cofactor);
                mpz_mul(prod, prod, partials->data[k].cofactor);

                /* Check if prod is a perfect square times smooth */
                mpz_t rem;
                mpz_init_set(rem, prod);
                for (int f = 0; f < fb->count && mpz_cmp_ui(rem, 1) > 0; f++) {
                    while (mpz_divisible_ui_p(rem, fb->primes[f]))
                        mpz_divexact_ui(rem, rem, fb->primes[f]);
                }
                int is_smooth = (mpz_cmp_ui(rem, 1) == 0);
                mpz_clear(rem);

                if (is_smooth) {
                    /* Triple match! */
                    for (int f = 0; f < veclen; f++)
                        combined_exps[f] = partials->data[i].exps[f] +
                                          partials->data[j].exps[f] +
                                          partials->data[k].exps[f];
                    /* Add prod's factorization */
                    mpz_t p2;
                    mpz_init_set(p2, prod);
                    for (int f = 0; f < fb->count; f++) {
                        while (mpz_divisible_ui_p(p2, fb->primes[f])) {
                            mpz_divexact_ui(p2, p2, fb->primes[f]);
                            combined_exps[f + 1]++;
                        }
                    }
                    mpz_clear(p2);

                    mpz_mul(combined_sv, partials->data[i].sqrt_val, partials->data[j].sqrt_val);
                    mpz_mul(combined_sv, combined_sv, partials->data[k].sqrt_val);
                    frs_add(fulls, combined_sv, combined_exps);
                    matched++;
                    mpz_set_ui(partials->data[i].cofactor, 0);
                    mpz_set_ui(partials->data[j].cofactor, 0);
                    mpz_set_ui(partials->data[k].cofactor, 0);
                    break;
                }
            }
            if (mpz_cmp_ui(partials->data[i].cofactor, 0) <= 0) break;
        }
    }

    free(small_idx);
    free(combined_exps);
    mpz_clears(prod, combined_sv, NULL);
    return matched;
}

/* ===== GF(2) solver ===== */
typedef unsigned long word_t;
#define WBITS (sizeof(word_t) * 8)

static int find_all_deps(frel_set_t *rels, int veclen,
                         int ***all_deps, int **all_sizes, int *ndeps) {
    int nrows = veclen, ncols = rels->count;
    *ndeps = 0;
    if (ncols <= nrows + 1) return 0;

    int nw = (ncols + WBITS - 1) / WBITS;
    word_t **mat = (word_t**)calloc(nrows, sizeof(word_t*));
    for (int i = 0; i < nrows; i++) {
        mat[i] = (word_t*)calloc(nw, sizeof(word_t));
        for (int j = 0; j < ncols; j++)
            if (rels->data[j].exps[i] & 1)
                mat[i][j / WBITS] |= (1UL << (j % WBITS));
    }

    int *pcol = (int*)malloc(nrows * sizeof(int));
    char *is_pc = (char*)calloc(ncols, 1);
    for (int i = 0; i < nrows; i++) {
        pcol[i] = -1;
        for (int j = 0; j < ncols; j++) {
            if (!is_pc[j] && (mat[i][j / WBITS] & (1UL << (j % WBITS)))) {
                pcol[i] = j; is_pc[j] = 1;
                for (int k = 0; k < nrows; k++)
                    if (k != i && (mat[k][j / WBITS] & (1UL << (j % WBITS))))
                        for (int w = 0; w < nw; w++) mat[k][w] ^= mat[i][w];
                break;
            }
        }
    }

    int max_d = 64;
    *all_deps = (int**)malloc(max_d * sizeof(int*));
    *all_sizes = (int*)malloc(max_d * sizeof(int));
    int nd = 0;
    for (int j = 0; j < ncols && nd < max_d; j++) {
        if (is_pc[j]) continue;
        int *dep = (int*)malloc((nrows + 1) * sizeof(int));
        int cnt = 0;
        dep[cnt++] = j;
        for (int i = 0; i < nrows; i++)
            if (mat[i][j / WBITS] & (1UL << (j % WBITS)))
                if (pcol[i] >= 0) dep[cnt++] = pcol[i];
        (*all_deps)[nd] = dep;
        (*all_sizes)[nd] = cnt;
        nd++;
    }
    *ndeps = nd;

    for (int i = 0; i < nrows; i++) free(mat[i]);
    free(mat); free(pcol); free(is_pc);
    return nd > 0;
}

/* ===== Factor extraction ===== */
static int try_extract(frel_set_t *rels, int *dep, int dsz,
                       const mpz_t N, const fb_t *fb, mpz_t factor) {
    mpz_t x, y, tmp;
    mpz_inits(x, y, tmp, NULL);
    mpz_set_ui(x, 1);
    int *tot = (int*)calloc(fb->count + 1, sizeof(int));

    for (int i = 0; i < dsz; i++) {
        mpz_mul(x, x, rels->data[dep[i]].sqrt_val);
        mpz_mod(x, x, N);
        for (int j = 0; j < fb->count + 1; j++)
            tot[j] += rels->data[dep[i]].exps[j];
    }

    for (int j = 0; j < fb->count + 1; j++)
        if (tot[j] & 1) { free(tot); mpz_clears(x, y, tmp, NULL); return 0; }

    mpz_set_ui(y, 1);
    for (int j = 0; j < fb->count; j++) {
        int e = tot[j + 1] / 2;
        if (e > 0) {
            mpz_set_ui(tmp, fb->primes[j]);
            mpz_powm_ui(tmp, tmp, e, N);
            mpz_mul(y, y, tmp); mpz_mod(y, y, N);
        }
    }

    mpz_sub(tmp, x, y); mpz_gcd(factor, tmp, N);
    if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
        free(tot); mpz_clears(x, y, tmp, NULL); return 1;
    }
    mpz_add(tmp, x, y); mpz_gcd(factor, tmp, N);
    int ok = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);
    free(tot); mpz_clears(x, y, tmp, NULL);
    return ok;
}

/* ===== Main ===== */
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N, factor;
    mpz_inits(N, factor, NULL);
    mpz_set_str(N, argv[1], 10);
    int ndigits = strlen(argv[1]);
    double t0 = wall_time();

    fprintf(stderr, "HSD: factoring %d-digit number\n", ndigits);

    /* Small factors */
    { int cnt; int *sp = sieve_primes(1000000, &cnt);
      for (int i = 0; i < cnt; i++)
        if (mpz_divisible_ui_p(N, sp[i])) {
            printf("FACTOR:%d\n", sp[i]); free(sp); mpz_clears(N, factor, NULL); return 0;
        }
      free(sp);
    }

    /* Parameters */
    double n = ndigits, ln_N = n * log(10);
    double L_exp = sqrt(ln_N * log(ln_N));

    int B = (int)exp(0.55 * L_exp);
    if (B < 500) B = 500;
    if (B > 5000000) B = 5000000;

    double log_M = 1.05 * L_exp;
    int M = (log_M > 16.0) ? 10000000 : (int)exp(log_M);
    if (M < 20000) M = 20000;
    if (M > 10000000) M = 10000000;

    /* KEY PARAMETER: allow much larger cofactors than standard QS */
    /* Standard: cofactor < B. HSD: cofactor up to B^2 or N^{1/4} */
    int max_cof_bits = (int)(log2(B) * 2.0);
    if (max_cof_bits > ndigits * 3.32 / 3) /* cap at N^{1/3} */
        max_cof_bits = (int)(ndigits * 3.32 / 3);

    fprintf(stderr, "Params: B=%d, M=%d, max_cof_bits=%d\n", B, M, max_cof_bits);

    fb_t fb = build_fb(N, B);
    fprintf(stderr, "Factor base: %d primes\n", fb.count);
    int target = fb.count + 30;

    prel_set_t partials;
    frel_set_t fulls;
    prs_init(&partials, fb.count + 1);
    frs_init(&fulls, fb.count + 1);

    /* Phase 1: QS sieve */
    fprintf(stderr, "Phase 1: QS sieve...\n");
    qs_sieve_partial(N, &fb, M, max_cof_bits, &partials, &fulls);
    fprintf(stderr, "  QS: %d full, %d partial (%.1fs)\n",
            fulls.count, partials.count, wall_time() - t0);

    /* Phase 2: MPQS polynomials */
    if (fulls.count < target) {
        fprintf(stderr, "Phase 2: MPQS...\n");
        for (int pi = 0; pi < 5000 && fulls.count < target; pi++) {
            mpqs_sieve_partial(N, &fb, M, max_cof_bits, pi, &partials, &fulls);
            if ((pi + 1) % 100 == 0) {
                /* Match cofactors */
                int m1 = match_cofactors(&partials, &fulls, &fb, B);
                fprintf(stderr, "  poly %d: %d full (+%d matched), %d partial (%.1fs)\n",
                        pi + 1, fulls.count, m1, partials.count, wall_time() - t0);
            }
        }
    }

    /* Phase 3: Aggressive cofactor matching (only if needed) */
    if (fulls.count < target) {
        fprintf(stderr, "Phase 3: Cofactor matching (%d partials)...\n", partials.count);
        int m1 = match_cofactors(&partials, &fulls, &fb, B);
        fprintf(stderr, "  Single/double LP: +%d full\n", m1);

        if (fulls.count < target) {
            int m2 = match_cofactor_triples(&partials, &fulls, &fb);
            fprintf(stderr, "  Triple LP: +%d full\n", m2);
        }
    }

    fprintf(stderr, "Total: %d full relations (need %d) (%.1fs)\n",
            fulls.count, target, wall_time() - t0);

    if (fulls.count < target) {
        fprintf(stderr, "FAIL: not enough relations\n");
        prs_free(&partials); frs_free(&fulls);
        free(fb.primes); free(fb.sqrtN);
        mpz_clears(N, factor, NULL);
        return 1;
    }

    /* Phase 4: Matrix + extraction */
    fprintf(stderr, "Phase 4: Matrix (%d x %d)...\n", fb.count + 1, fulls.count);
    int **deps; int *dsizes; int ndeps;
    if (!find_all_deps(&fulls, fb.count + 1, &deps, &dsizes, &ndeps)) {
        fprintf(stderr, "FAIL: no deps\n");
        prs_free(&partials); frs_free(&fulls);
        free(fb.primes); free(fb.sqrtN);
        mpz_clears(N, factor, NULL);
        return 1;
    }

    fprintf(stderr, "%d deps found\n", ndeps);
    int success = 0;
    for (int d = 0; d < ndeps && !success; d++)
        success = try_extract(&fulls, deps[d], dsizes[d], N, &fb, factor);

    if (success) {
        mpz_t cof; mpz_init(cof); mpz_tdiv_q(cof, N, factor);
        gmp_printf("FACTOR:%Zd\n", factor);
        gmp_fprintf(stderr, "%Zd x %Zd (%.3fs)\n", factor, cof, wall_time() - t0);
        mpz_clear(cof);
    } else {
        fprintf(stderr, "FAIL: extraction failed (%.1fs)\n", wall_time() - t0);
    }

    for (int d = 0; d < ndeps; d++) free(deps[d]);
    free(deps); free(dsizes);
    prs_free(&partials); frs_free(&fulls);
    free(fb.primes); free(fb.sqrtN);
    mpz_clears(N, factor, NULL);
    return success ? 0 : 1;
}
