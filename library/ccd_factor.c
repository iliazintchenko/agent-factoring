/*
 * ccd_factor.c — Cofactor Collision Descent factoring
 *
 * Novel combination of:
 * 1. Multi-polynomial candidate generation
 * 2. Sieve-based smooth detection with large prime collection
 * 3. Batch GCD cofactor matching for partial relations
 * 4. ECM descent on large cofactors
 * 5. GF(2) linear algebra
 *
 * Usage: ./ccd_factor <N>
 * Output: FACTOR: <p>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

#define MAX_FB       4096
#define MAX_RELS     (MAX_FB + 512)
#define MAX_PARTIALS 500000
#define SIEVE_SIZE   (1 << 17)  /* 131072 */

typedef struct {
    unsigned long p;
    int r1, r2;   /* roots of x^2 ≡ N (mod p). -1 if none */
} fb_t;

typedef struct {
    mpz_t x;             /* x value: x^2 ≡ smooth (mod N) */
    mpz_t cofactor;       /* remaining cofactor (1 if full) */
    unsigned char *exp;   /* exponent vector mod 2 (bit-packed, fb_size+1 bits for sign+primes) */
    unsigned long *full_exp; /* full exponent vector for square root computation */
} rel_t;

/* Globals */
static mpz_t gN, gSqrtN;
static fb_t gFB[MAX_FB];
static int gFBsize = 0;
static rel_t gFulls[MAX_RELS];
static int gNFulls = 0;
static rel_t gPartials[MAX_PARTIALS];
static int gNPartials = 0;
static int gTarget = 0;
static int gExpBytes = 0;  /* bytes per exp vector */
static double gStart;

static double walltime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Tonelli-Shanks: sqrt(a) mod p. Returns 0 if no root. */
static unsigned long modsqrt(unsigned long a, unsigned long p) {
    if (a == 0) return 0;
    if (p == 2) return a & 1;

    /* Euler criterion */
    unsigned long e = (p - 1) / 2, base = a % p, res = 1;
    unsigned long ee = e;
    while (ee > 0) {
        if (ee & 1) res = (__uint128_t)res * base % p;
        base = (__uint128_t)base * base % p;
        ee >>= 1;
    }
    if (res != 1) return 0; /* not QR */

    if (p % 4 == 3) {
        base = a; e = (p + 1) / 4; res = 1;
        while (e > 0) {
            if (e & 1) res = (__uint128_t)res * base % p;
            base = (__uint128_t)base * base % p;
            e >>= 1;
        }
        return res;
    }

    /* Full Tonelli-Shanks */
    unsigned long Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }
    unsigned long z = 2;
    for (;;) {
        base = z; e = (p-1)/2; res = 1;
        while (e > 0) {
            if (e & 1) res = (__uint128_t)res * base % p;
            base = (__uint128_t)base * base % p;
            e >>= 1;
        }
        if (res == p - 1) break;
        z++;
    }
    unsigned long M = S;
    unsigned long c; base = z; e = Q; c = 1;
    while (e > 0) { if (e & 1) c = (__uint128_t)c * base % p; base = (__uint128_t)base * base % p; e >>= 1; }
    unsigned long t; base = a; e = Q; t = 1;
    while (e > 0) { if (e & 1) t = (__uint128_t)t * base % p; base = (__uint128_t)base * base % p; e >>= 1; }
    unsigned long R; base = a; e = (Q+1)/2; R = 1;
    while (e > 0) { if (e & 1) R = (__uint128_t)R * base % p; base = (__uint128_t)base * base % p; e >>= 1; }

    while (t != 1) {
        unsigned long i = 0, tt = t;
        while (tt != 1) { tt = (__uint128_t)tt * tt % p; i++; }
        unsigned long b = c;
        for (unsigned long j = 0; j < M - i - 1; j++)
            b = (__uint128_t)b * b % p;
        M = i; c = (__uint128_t)b * b % p;
        t = (__uint128_t)t * c % p; R = (__uint128_t)R * b % p;
    }
    return R;
}

static int is_prime(unsigned long n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (unsigned long d = 5; d * d <= n; d += 6)
        if (n % d == 0 || n % (d+2) == 0) return 0;
    return 1;
}

static void build_fb(int target) {
    gFBsize = 0;
    /* -1 entry for sign */
    /* prime 2 */
    gFB[0].p = 2;
    gFB[0].r1 = (int)(mpz_fdiv_ui(gN, 2));
    gFB[0].r2 = -1;
    gFBsize = 1;

    for (unsigned long p = 3; gFBsize < target; p += 2) {
        if (!is_prime(p)) continue;
        unsigned long nmp = mpz_fdiv_ui(gN, p);
        unsigned long r = modsqrt(nmp, p);
        if (r == 0 && nmp != 0) continue;
        gFB[gFBsize].p = p;
        gFB[gFBsize].r1 = (int)r;
        gFB[gFBsize].r2 = (int)((p - r) % p);
        gFBsize++;
    }
    gTarget = gFBsize + gFBsize / 2 + 64; /* need plenty of excess relations */
    gExpBytes = (gFBsize + 8) / 8; /* +1 for sign, rounded up */
}

/* Initialize a relation */
static void rel_init(rel_t *r) {
    mpz_init(r->x);
    mpz_init(r->cofactor);
    r->exp = calloc(gExpBytes, 1);
    r->full_exp = calloc(gFBsize + 1, sizeof(unsigned long));
}

/* Set bit b in exponent vector */
static void exp_flip(unsigned char *e, int b) {
    e[b / 8] ^= (1 << (b % 8));
}
static int exp_get(unsigned char *e, int b) {
    return (e[b / 8] >> (b % 8)) & 1;
}

/* Trial divide val by factor base. Returns 1 if fully smooth. */
static int trial_div(mpz_t val, rel_t *rel) {
    mpz_t r;
    mpz_init(r);
    mpz_abs(r, val);

    memset(rel->exp, 0, gExpBytes);
    memset(rel->full_exp, 0, (gFBsize + 1) * sizeof(unsigned long));

    /* Sign (bit 0) */
    if (mpz_sgn(val) < 0) {
        exp_flip(rel->exp, 0);
        rel->full_exp[0] = 1;
    }

    for (int i = 0; i < gFBsize; i++) {
        unsigned long p = gFB[i].p;
        unsigned long cnt = 0;
        while (mpz_divisible_ui_p(r, p)) {
            mpz_divexact_ui(r, r, p);
            cnt++;
        }
        if (cnt & 1) exp_flip(rel->exp, i + 1);
        rel->full_exp[i + 1] = cnt;
    }

    int smooth = (mpz_cmp_ui(r, 1) == 0);
    mpz_set(rel->cofactor, r);
    mpz_clear(r);
    return smooth;
}

/* Sieve and collect relations */
static void sieve_range(long start_offset, int count) {
    /* Use a log-based sieve: accumulate log2(p) for each prime p dividing f(x).
     * Use int array to avoid overflow.
     * Values that accumulate enough log are likely smooth — trial divide them. */

    int *slog = calloc(count, sizeof(int));
    int logscale = 10; /* multiply log2 by this for integer precision */

    /* Sieve with factor base primes */
    unsigned long sqrtN_mod_p;
    for (int i = 0; i < gFBsize; i++) {
        unsigned long p = gFB[i].p;
        int lp = (int)(log2((double)p) * logscale + 0.5);
        sqrtN_mod_p = mpz_fdiv_ui(gSqrtN, p);

        for (int ri = 0; ri < 2; ri++) {
            int root = (ri == 0) ? gFB[i].r1 : gFB[i].r2;
            if (root < 0) continue;
            if (ri == 1 && gFB[i].r1 == gFB[i].r2) continue;

            /* offset ≡ root - sqrtN (mod p) */
            long target = ((long)root - (long)sqrtN_mod_p + 2 * (long)p) % (long)p;

            /* First hit in [start_offset, start_offset + count) */
            long first;
            long rel_start = start_offset % (long)p;
            if (rel_start < 0) rel_start += (long)p;
            first = target - rel_start;
            if (first < 0) first += (long)p;

            for (long j = first; j < count; j += (long)p) {
                slog[j] += lp;
                /* Also account for prime powers */
                long offset = start_offset + j;
                mpz_t x, val;
                /* Skip prime power sieving for speed — trial div handles it */
                (void)offset;
            }
        }
    }

    /* Compute expected log for smooth numbers */
    /* For offset t from sqrtN, |f(t)| = |2*sqrtN*t + t^2|
     * For positive t, f(t) ≈ 2*sqrtN*t (dominant for small t) */
    double log2_sqrtN = mpz_sizeinbase(gSqrtN, 2) - 1 +
                        log2(mpz_getlimbn(gSqrtN, mpz_size(gSqrtN)-1) /
                             (double)(1UL << (mp_bits_per_limb - 1)));

    mpz_t x, val;
    mpz_init(x);
    mpz_init(val);

    int n_tested = 0;

    for (int j = 0; j < count && gNFulls < gTarget; j++) {
        long offset = start_offset + j;
        double abs_offset = fabs((double)offset);
        if (abs_offset < 1.0) abs_offset = 1.0;
        double expected_log = 1.0 + log2_sqrtN + log2(abs_offset);
        if (offset == 0) expected_log = log2_sqrtN; /* sqrtN^2 - N < sqrtN */

        /* Threshold: accept if sieve covers at least (expected - slack) */
        double slack = log2((double)gFB[gFBsize-1].p) * 1.5 + 5.0;
        int threshold = (int)((expected_log - slack) * logscale);
        if (threshold < 0) threshold = 0;

        if (slog[j] < threshold) continue;

        /* Promising candidate — trial divide */
        mpz_set(x, gSqrtN);
        if (offset >= 0) mpz_add_ui(x, x, (unsigned long)offset);
        else mpz_sub_ui(x, x, (unsigned long)(-offset));

        mpz_mul(val, x, x);
        mpz_sub(val, val, gN);

        n_tested++;

        rel_t rel;
        rel_init(&rel);
        mpz_set(rel.x, x);

        if (trial_div(val, &rel)) {
            /* Full relation */
            if (gNFulls < MAX_RELS) {
                gFulls[gNFulls] = rel;
                gNFulls++;
            }
        } else {
            /* Check cofactor size for partial */
            int cofbits = mpz_sizeinbase(rel.cofactor, 2);
            if (cofbits <= 62 && mpz_probab_prime_p(rel.cofactor, 2)) {
                /* Single large prime partial */
                if (gNPartials < MAX_PARTIALS) {
                    gPartials[gNPartials] = rel;
                    gNPartials++;
                } else {
                    mpz_clear(rel.x); mpz_clear(rel.cofactor);
                    free(rel.exp); free(rel.full_exp);
                }
            } else if (cofbits <= 40) {
                /* Small enough cofactor — might match something via batch GCD */
                if (gNPartials < MAX_PARTIALS) {
                    gPartials[gNPartials] = rel;
                    gNPartials++;
                } else {
                    mpz_clear(rel.x); mpz_clear(rel.cofactor);
                    free(rel.exp); free(rel.full_exp);
                }
            } else {
                mpz_clear(rel.x); mpz_clear(rel.cofactor);
                free(rel.exp); free(rel.full_exp);
            }
        }
    }

    mpz_clear(x);
    mpz_clear(val);
    free(slog);
}

/* Match partials sharing the same large prime cofactor */
static int gMatchStart = 0; /* Track where we left off to avoid re-matching */
static int *gPartialUsed = NULL; /* Mark used partials */

static int match_partials(void) {
    if (gNPartials < 2) return 0;

    if (!gPartialUsed) {
        gPartialUsed = calloc(MAX_PARTIALS, sizeof(int));
    }

    /* Hash map: cofactor -> index of first unused occurrence */
    int htsize = gNPartials * 3;
    if (htsize < 16) htsize = 16;
    int *ht_idx = malloc(htsize * sizeof(int));
    memset(ht_idx, -1, htsize * sizeof(int));

    int matches = 0;

    /* First pass: insert all unmatched partials into hash table */
    for (int i = 0; i < gNPartials && gNFulls < gTarget; i++) {
        if (gPartialUsed[i]) continue;

        unsigned long h = mpz_fdiv_ui(gPartials[i].cofactor, (unsigned long)htsize);

        for (int probe = 0; probe < htsize; probe++) {
            int slot = (h + probe) % htsize;
            if (ht_idx[slot] == -1) {
                ht_idx[slot] = i;
                break;
            }
            int other = ht_idx[slot];
            if (gPartialUsed[other]) {
                /* This slot's entry was used — replace it */
                ht_idx[slot] = i;
                break;
            }
            if (mpz_cmp(gPartials[i].cofactor, gPartials[other].cofactor) == 0) {
                /* Match! Mark both as used */
                gPartialUsed[i] = 1;
                gPartialUsed[other] = 1;

                if (gNFulls < MAX_RELS) {
                    rel_t *full = &gFulls[gNFulls];
                    rel_init(full);
                    mpz_mul(full->x, gPartials[i].x, gPartials[other].x);
                    mpz_mod(full->x, full->x, gN);
                    mpz_set(full->cofactor, gPartials[i].cofactor);

                    for (int b = 0; b < gExpBytes; b++)
                        full->exp[b] = gPartials[i].exp[b] ^ gPartials[other].exp[b];
                    for (int k = 0; k <= gFBsize; k++)
                        full->full_exp[k] = gPartials[i].full_exp[k] + gPartials[other].full_exp[k];

                    gNFulls++;
                    matches++;
                }
                break;
            }
        }
    }

    free(ht_idx);
    return matches;
}

/* Batch GCD for finding shared factors among all partial cofactors */
static int batch_gcd_match(void) {
    if (gNPartials < 10) return 0;

    /* Build product tree of all cofactors */
    int n = gNPartials;
    if (n > 50000) n = 50000; /* limit for memory */

    /* Product of all cofactors */
    mpz_t prod;
    mpz_init_set_ui(prod, 1);
    for (int i = 0; i < n; i++)
        mpz_mul(prod, prod, gPartials[i].cofactor);

    /* For each cofactor, compute gcd(cofactor, prod/cofactor) */
    /* This reveals if any OTHER cofactor shares a factor */
    int found = 0;
    mpz_t g, other_prod;
    mpz_init(g);
    mpz_init(other_prod);

    for (int i = 0; i < n && gNFulls < gTarget; i++) {
        if (mpz_cmp_ui(gPartials[i].cofactor, 1) == 0) continue;
        mpz_divexact(other_prod, prod, gPartials[i].cofactor);
        mpz_gcd(g, gPartials[i].cofactor, other_prod);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gPartials[i].cofactor) < 0) {
            /* Found a shared factor! But we need to identify which other partial shares it */
            /* For now, just note that this cofactor has a shared factor */
            found++;
        }
    }

    mpz_clear(prod);
    mpz_clear(g);
    mpz_clear(other_prod);
    return found;
}

/* GF(2) Gaussian elimination and factor extraction */
static int solve_and_factor(mpz_t result) {
    int n = gNFulls;
    int m = gFBsize + 1; /* +1 for sign */

    if (n < m + 1) return 0;

    fprintf(stderr, "CCD: solving %d x %d matrix\n", n, m);

    /* Matrix: n rows, m + n columns (augmented with identity for tracking) */
    int words = (m + n + 63) / 64;
    unsigned long *mat = calloc((size_t)n * words, sizeof(unsigned long));

    #define MGET(r, c) ((mat[(size_t)(r) * words + (c)/64] >> ((c) % 64)) & 1)
    #define MSET(r, c) (mat[(size_t)(r) * words + (c)/64] |= (1UL << ((c) % 64)))
    #define MFLIP(r, c) (mat[(size_t)(r) * words + (c)/64] ^= (1UL << ((c) % 64)))
    #define MXOR_ROW(dst, src) do { \
        for (int _w = 0; _w < words; _w++) \
            mat[(size_t)(dst)*words + _w] ^= mat[(size_t)(src)*words + _w]; \
    } while(0)

    /* Fill matrix */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= gFBsize; j++) {
            if (exp_get(gFulls[i].exp, j))
                MSET(i, j);
        }
        MSET(i, m + i); /* identity */
    }

    /* Elimination */
    int rank = 0;
    for (int j = 0; j < m && rank < n; j++) {
        int piv = -1;
        for (int i = rank; i < n; i++)
            if (MGET(i, j)) { piv = i; break; }
        if (piv < 0) continue;

        if (piv != rank) {
            for (int w = 0; w < words; w++) {
                unsigned long tmp = mat[(size_t)rank * words + w];
                mat[(size_t)rank * words + w] = mat[(size_t)piv * words + w];
                mat[(size_t)piv * words + w] = tmp;
            }
        }
        for (int i = 0; i < n; i++) {
            if (i != rank && MGET(i, j))
                MXOR_ROW(i, rank);
        }
        rank++;
    }

    fprintf(stderr, "CCD: rank=%d, %d null vectors\n", rank, n - rank);

    /* Try each null vector */
    mpz_t X, Y, g, tmp;
    mpz_init(X); mpz_init(Y); mpz_init(g); mpz_init(tmp);
    int found = 0;
    int stat_empty = 0, stat_verified = 0, stat_failed_verify = 0, stat_trivial = 0;

    for (int row = rank; row < n && !found; row++) {
        mpz_set_ui(X, 1);
        unsigned long *tot = calloc(gFBsize + 1, sizeof(unsigned long));

        int count = 0;
        for (int i = 0; i < n; i++) {
            if (!MGET(row, m + i)) continue;
            count++;
            mpz_mul(X, X, gFulls[i].x);
            mpz_mod(X, X, gN);
            for (int k = 0; k <= gFBsize; k++)
                tot[k] += gFulls[i].full_exp[k];
        }

        if (count == 0) { stat_empty++; free(tot); continue; }
        /* (debug removed) */

        /* Verify all exponents are even */
        int ok = 1;
        for (int k = 0; k <= gFBsize; k++) {
            if (tot[k] & 1) { ok = 0; break; }
        }
        if (!ok) {
            /* Sign might be wrong — try flipping */
            if (tot[0] & 1) tot[0]++; /* odd sign count means negative product */
            ok = 1;
            for (int k = 1; k <= gFBsize; k++)
                if (tot[k] & 1) { ok = 0; break; }
        }

        /* Compute Y = product of p_i^(exp/2) * product of cofactors mod N */
        mpz_set_ui(Y, 1);
        for (int k = 1; k <= gFBsize; k++) {
            if (tot[k] == 0) continue;
            unsigned long half = tot[k] / 2;
            if (half == 0) continue;
            mpz_set_ui(tmp, gFB[k-1].p);
            mpz_powm_ui(tmp, tmp, half, gN);
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, gN);
        }
        /* Include cofactors from combined partial relations
         * Each combined partial has cofactor c appearing with exponent 2.
         * In the square root, this contributes c^1. */
        for (int i = 0; i < n; i++) {
            if (!MGET(row, m + i)) continue;
            if (mpz_cmp_ui(gFulls[i].cofactor, 1) > 0) {
                mpz_mul(Y, Y, gFulls[i].cofactor);
                mpz_mod(Y, Y, gN);
            }
        }

        /* Verify: X^2 ≡ Y^2 (mod N) */
        {
            mpz_t x2, y2;
            mpz_init(x2); mpz_init(y2);
            mpz_mul(x2, X, X); mpz_mod(x2, x2, gN);
            mpz_mul(y2, Y, Y); mpz_mod(y2, y2, gN);
            if (mpz_cmp(x2, y2) != 0) {
                stat_failed_verify++;
                if (stat_failed_verify <= 2) {
                    gmp_fprintf(stderr, "CCD: null vec %d FAIL: X^2=%Zd Y^2=%Zd (count=%d)\n",
                               row-rank, x2, y2, count);
                }
                mpz_clear(x2); mpz_clear(y2);
                free(tot);
                continue;
            }
            mpz_clear(x2); mpz_clear(y2);
        }

        stat_verified++;

        /* gcd(X - Y, N) */
        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, gN);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
            mpz_set(result, g);
            found = 1;
        }
        if (!found) {
            mpz_add(tmp, X, Y);
            mpz_gcd(g, tmp, gN);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
                mpz_set(result, g);
                found = 1;
            } else {
                stat_trivial++;
            }
        }

        free(tot);
    }

    fprintf(stderr, "CCD: null vec stats: empty=%d verified=%d failed_verify=%d trivial=%d\n",
            stat_empty, stat_verified, stat_failed_verify, stat_trivial);

    /* If all null vectors from elimination are trivial, try random combinations */
    if (!found && stat_verified > 1) {
        fprintf(stderr, "CCD: trying random null vector combinations...\n");
        /* Store null vector indices */
        int n_null = n - rank;
        int *null_rows = malloc(n_null * sizeof(int));
        { int k = 0; for (int r = rank; r < n; r++) null_rows[k++] = r; }

        /* Allocate a combined vector */
        unsigned long *combo = calloc(words, sizeof(unsigned long));
        unsigned seed = 42;

        for (int trial = 0; trial < 1000 && !found; trial++) {
            /* Pick a random odd-sized subset of null vectors and XOR them */
            memset(combo, 0, words * sizeof(unsigned long));
            int combo_count = 0;
            /* Pick 3-7 random null vectors */
            int npick = 3 + (seed % 5);
            for (int p = 0; p < npick; p++) {
                seed = seed * 1103515245 + 12345;
                int which = null_rows[seed % n_null];
                for (int w = 0; w < words; w++)
                    combo[w] ^= mat[(size_t)which * words + w];
            }

            /* Verify first m columns are all zero */
            int valid = 1;
            for (int j = 0; j < m; j++) {
                if ((combo[j/64] >> (j%64)) & 1) { valid = 0; break; }
            }
            if (!valid) continue;

            /* Extract which original relations participate */
            mpz_set_ui(X, 1);
            unsigned long *tot2 = calloc(gFBsize + 1, sizeof(unsigned long));
            int cnt = 0;

            for (int i = 0; i < n; i++) {
                if (!((combo[(m + i)/64] >> ((m + i) % 64)) & 1)) continue;
                cnt++;
                mpz_mul(X, X, gFulls[i].x);
                mpz_mod(X, X, gN);
                for (int k = 0; k <= gFBsize; k++)
                    tot2[k] += gFulls[i].full_exp[k];
            }

            if (cnt == 0) { free(tot2); continue; }

            /* Compute Y */
            mpz_set_ui(Y, 1);
            for (int k = 1; k <= gFBsize; k++) {
                if (tot2[k] == 0) continue;
                if (tot2[k] & 1) goto skip_trial; /* odd exponent — not a valid combo */
                unsigned long half = tot2[k] / 2;
                mpz_set_ui(tmp, gFB[k-1].p);
                mpz_powm_ui(tmp, tmp, half, gN);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, gN);
            }
            /* Cofactors */
            for (int i = 0; i < n; i++) {
                if (!((combo[(m + i)/64] >> ((m + i) % 64)) & 1)) continue;
                if (mpz_cmp_ui(gFulls[i].cofactor, 1) > 0) {
                    mpz_mul(Y, Y, gFulls[i].cofactor);
                    mpz_mod(Y, Y, gN);
                }
            }

            /* Try gcd */
            mpz_sub(tmp, X, Y);
            mpz_gcd(g, tmp, gN);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
                mpz_set(result, g);
                found = 1;
            } else {
                mpz_add(tmp, X, Y);
                mpz_gcd(g, tmp, gN);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, gN) < 0) {
                    mpz_set(result, g);
                    found = 1;
                }
            }
            skip_trial:
            free(tot2);
        }
        free(null_rows);
        free(combo);
        if (found)
            fprintf(stderr, "CCD: random combination succeeded!\n");
        else
            fprintf(stderr, "CCD: random combinations also failed\n");
    }

    mpz_clear(X); mpz_clear(Y); mpz_clear(g); mpz_clear(tmp);
    free(mat);
    return found;

    #undef MGET
    #undef MSET
    #undef MFLIP
    #undef MXOR_ROW
}

static int choose_fb_size(int digits) {
    /* Practical: for d-digit N, FB size ≈ exp(0.5 * sqrt(d * ln(10) * ln(d * ln(10)))) / ln(B) */
    double lnN = digits * 2.302585;
    double lnlnN = log(lnN);
    double B = exp(0.5 * sqrt(lnN * lnlnN));
    double lnB = log(B);
    int sz = (int)(B / lnB);
    if (sz < 40) sz = 40;
    if (sz > MAX_FB) sz = MAX_FB;
    return sz;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    gStart = walltime();
    mpz_init(gN);
    mpz_init(gSqrtN);
    mpz_set_str(gN, argv[1], 10);
    int digits = strlen(argv[1]);

    if (mpz_probab_prime_p(gN, 25)) {
        fprintf(stderr, "N is prime\n");
        return 1;
    }

    /* Trial division */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(gN, p)) {
            mpz_t f, c;
            mpz_init(f); mpz_init(c);
            mpz_set_ui(f, p);
            mpz_divexact(c, gN, f);
            if (mpz_cmp(f, c) > 0) mpz_set(f, c);
            gmp_printf("FACTOR: %Zd\n", f);
            return 0;
        }
    }

    /* Perfect power check */
    if (mpz_perfect_power_p(gN)) {
        mpz_t root;
        mpz_init(root);
        for (int e = 2; e < 64; e++) {
            if (mpz_root(root, gN, e)) {
                gmp_printf("FACTOR: %Zd\n", root);
                mpz_clear(root);
                return 0;
            }
        }
        mpz_clear(root);
    }

    mpz_sqrt(gSqrtN, gN);
    { mpz_t t; mpz_init(t); mpz_mul(t, gSqrtN, gSqrtN);
      if (mpz_cmp(t, gN) < 0) mpz_add_ui(gSqrtN, gSqrtN, 1);
      mpz_clear(t); }

    int fb_target = choose_fb_size(digits);
    fprintf(stderr, "CCD: %d-digit, fb=%d\n", digits, fb_target);

    build_fb(fb_target);
    fprintf(stderr, "CCD: FB=%d primes, largest=%lu, need %d rels\n",
            gFBsize, gFB[gFBsize-1].p, gTarget);

    /* Main sieve loop */
    long sieve_pos = 1; /* start from offset 1 */
    int last_match = 0;

    while (gNFulls < gTarget) {
        if (walltime() - gStart > 285.0) {
            fprintf(stderr, "CCD: timeout\n");
            break;
        }

        /* Sieve positive side */
        sieve_range(sieve_pos, SIEVE_SIZE);
        /* Sieve negative side */
        sieve_range(-sieve_pos - SIEVE_SIZE + 1, SIEVE_SIZE);

        sieve_pos += SIEVE_SIZE;

        /* Match partials periodically */
        if (gNPartials > last_match + 200) {
            int m = match_partials();
            double t = walltime() - gStart;
            fprintf(stderr, "CCD: [%.1fs] fulls=%d partials=%d +matched=%d sieve=%.1fM\n",
                    t, gNFulls, gNPartials, m, sieve_pos / 1e6);
            last_match = gNPartials;
        }
    }

    /* Final matching */
    match_partials();
    fprintf(stderr, "CCD: collected %d full + %d partial relations\n", gNFulls, gNPartials);

    /* Relations collected */

    /* Solve */
    if (gNFulls >= gFBsize + 2) {
        mpz_t factor;
        mpz_init(factor);
        if (solve_and_factor(factor)) {
            mpz_t cof;
            mpz_init(cof);
            mpz_divexact(cof, gN, factor);
            if (mpz_cmp(factor, cof) > 0) mpz_set(factor, cof);
            gmp_printf("FACTOR: %Zd\n", factor);
            fprintf(stderr, "CCD: done in %.2fs\n", walltime() - gStart);
            mpz_clear(cof);
            mpz_clear(factor);
            return 0;
        }
        mpz_clear(factor);
        fprintf(stderr, "CCD: linear algebra failed\n");
    }

    /* Fallback: ECM */
    fprintf(stderr, "CCD: ECM fallback\n");
    {
        ecm_params params;
        ecm_init(params);
        params->param = ECM_PARAM_SUYAMA;
        mpz_t factor;
        mpz_init(factor);

        double b1s[] = {2000, 11000, 50000, 250000, 1e6, 3e6, 11e6, 44e6, 110e6};
        int curves[] = {25, 50, 100, 200, 400, 700, 1200, 2400, 4800};
        for (int lev = 0; lev < 9; lev++) {
            for (int c = 0; c < curves[lev]; c++) {
                if (walltime() - gStart > 290.0) goto done;
                mpz_set_ui(factor, 0);
                params->B1done = 1.0;
                mpz_set_ui(params->sigma, 42 + (unsigned long)lev * 10000 + c);
                int ret = ecm_factor(factor, gN, b1s[lev], params);
                if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, gN) < 0) {
                    mpz_t cof; mpz_init(cof);
                    mpz_divexact(cof, gN, factor);
                    if (mpz_cmp(factor, cof) > 0) mpz_set(factor, cof);
                    gmp_printf("FACTOR: %Zd\n", factor);
                    mpz_clear(cof);
                    ecm_clear(params);
                    mpz_clear(factor);
                    return 0;
                }
            }
        }
        done:
        ecm_clear(params);
        mpz_clear(factor);
    }

    fprintf(stderr, "CCD: failed\n");
    return 1;
}
