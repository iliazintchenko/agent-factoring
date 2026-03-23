/*
 * specialq_siqs.c - Special-Q SIQS: novel QS variant using lattice sieving
 *
 * Key idea: Apply NFS-style "special-Q" technique to QS.
 * For each large prime Q not in the factor base:
 *   1. Find positions x where Q | QS_poly(x)
 *   2. Divide out Q from QS_poly(x) at those positions
 *   3. The quotient QS_poly(x)/Q is smaller -> more likely B-smooth
 *   4. Relations from special-Q: QS_poly(x) = Q * (smooth cofactor)
 *
 * This is different from standard large prime variation:
 * - LP variation: sieve normally, accept relations with 1-2 large prime factors
 * - Special-Q: CHOOSE a large prime Q upfront, sieve only positions divisible by Q
 *
 * The advantage: for each special-Q, the effective polynomial values are
 * QS_poly(x)/Q ~ sqrt(N)/Q, which is Q times smaller than standard QS values.
 * This dramatically increases smoothness probability per candidate.
 *
 * The cost: only ~2/Q fraction of sieve positions are divisible by Q,
 * so we need Q different special-Q values to cover the same sieve space.
 *
 * Net effect: smoothness increases by Q^u (where u = ln(sqrt(N))/ln(B)),
 * but we process Q times fewer candidates. For u > 1, the smoothness
 * gain dominates, giving a net speedup.
 *
 * Compile: gcc -O3 -march=native -o specialq_siqs library/specialq_siqs.c -lgmp -lm
 * Usage: ./specialq_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_FB 50000
#define MAX_RELS 200000
#define SIEVE_BLOCK 32768
#define SEED 42

static int primes[MAX_FB];
static int fb_roots[MAX_FB]; /* sqrt(N) mod p */
static int fb_size;
static int fb_logp[MAX_FB];

/* ==================== Relation storage ==================== */

typedef struct {
    mpz_t Y;            /* Y = (x + sqrt_kN) mod N, so Y^2 = QS_val (mod N) */
    int *exponents;     /* exponent vector over factor base (fb_size entries) */
    int special_q;      /* the special-Q prime for this relation */
    int sq_exp;         /* exponent of special-Q (always 1) */
    int large_prime;    /* 0 for full, >0 for SLP on the reduced part */
} rel_t;

static rel_t rels[MAX_RELS];
static int nrels = 0;

/* LP hash for SLP combining */
#define LP_HASH_SZ (1 << 18)
static int lp_hash[LP_HASH_SZ];
static int lp_next[MAX_RELS];
static int ncombined = 0;

static void init_lp(void) { memset(lp_hash, -1, sizeof(lp_hash)); }

/* ==================== Distinct special-Q tracking ==================== */

#define MAX_SQ 100000
static int distinct_sq[MAX_SQ];
static int n_distinct_sq = 0;

static int find_or_add_sq(int q) {
    for (int i = 0; i < n_distinct_sq; i++)
        if (distinct_sq[i] == q) return i;
    if (n_distinct_sq < MAX_SQ)
        distinct_sq[n_distinct_sq++] = q;
    return n_distinct_sq - 1;
}

/* ==================== Primes and math ==================== */

static void gen_primes(int limit) {
    char *s = calloc(limit + 1, 1);
    fb_size = 0;
    for (int i = 2; i <= limit; i++) s[i] = 1;
    for (int i = 2; (long)i * i <= limit; i++)
        if (s[i]) for (int j = i * i; j <= limit; j += i) s[j] = 0;
    for (int i = 2; i <= limit && fb_size < MAX_FB; i++)
        if (s[i]) primes[fb_size++] = i;
    free(s);
}

static uint32_t modsqrt(uint32_t a, uint32_t p) {
    if (a == 0) return 0;
    if (p == 2) return a & 1;
    /* Euler criterion */
    uint64_t test = 1, base = a % p;
    uint32_t exp = (p - 1) / 2;
    while (exp > 0) {
        if (exp & 1) test = test * base % p;
        base = base * base % p;
        exp >>= 1;
    }
    if (test != 1) return UINT32_MAX; /* not QR */

    /* Tonelli-Shanks */
    uint32_t Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }
    if (S == 1) {
        base = a; exp = (p + 1) / 4;
        uint64_t r = 1;
        while (exp > 0) {
            if (exp & 1) r = r * base % p;
            base = base * base % p;
            exp >>= 1;
        }
        return (uint32_t)r;
    }
    uint32_t z = 2;
    while (1) {
        base = z; exp = (p - 1) / 2; test = 1;
        while (exp > 0) {
            if (exp & 1) test = test * base % p;
            base = base * base % p;
            exp >>= 1;
        }
        if (test == p - 1) break;
        z++;
    }
    uint32_t M = S;
    uint64_t c = 1; base = z; exp = Q;
    while (exp > 0) { if (exp & 1) c = c * base % p; base = base * base % p; exp >>= 1; }
    uint64_t t = 1; base = a; exp = Q;
    while (exp > 0) { if (exp & 1) t = t * base % p; base = base * base % p; exp >>= 1; }
    uint64_t R = 1; base = a; exp = (Q + 1) / 2;
    while (exp > 0) { if (exp & 1) R = R * base % p; base = base * base % p; exp >>= 1; }

    while (t != 1) {
        uint32_t i = 1; uint64_t tmp = t * t % p;
        while (tmp != 1) { tmp = tmp * tmp % p; i++; }
        uint64_t b = c;
        for (uint32_t j = 0; j < M - i - 1; j++) b = b * b % p;
        M = i; c = b * b % p; t = t * c % p; R = R * b % p;
    }
    return (uint32_t)R;
}

static uint32_t modinv(uint32_t a, uint32_t m) {
    int64_t g = m, x = 0, y = 1, a1 = a;
    while (a1) { int64_t q = g / a1, t = g - q * a1; g = a1; a1 = t; t = x - q * y; x = y; y = t; }
    return (uint32_t)((x % (int64_t)m + m) % m);
}

/* Build factor base for kN */
static int build_fb(mpz_t kN, int target) {
    gen_primes(target * 15);
    int count = 0;
    int *new_primes = malloc(target * sizeof(int));
    int *new_roots = malloc(target * sizeof(int));
    int *new_logp = malloc(target * sizeof(int));

    for (int i = 0; i < fb_size && count < target; i++) {
        int p = primes[i];
        uint32_t kn_mod = mpz_fdiv_ui(kN, p);
        uint32_t s = modsqrt(kn_mod, p);
        if (s == UINT32_MAX) continue;
        new_primes[count] = p;
        new_roots[count] = s;
        new_logp[count] = (int)(log2(p) * 1.44 + 0.5);
        count++;
    }

    memcpy(primes, new_primes, count * sizeof(int));
    memcpy(fb_roots, new_roots, count * sizeof(int));
    memcpy(fb_logp, new_logp, count * sizeof(int));
    fb_size = count;
    free(new_primes); free(new_roots); free(new_logp);
    return count;
}

/* ==================== Add relation ==================== */

static void add_relation(mpz_t Y, int *exp_vec, int sq, int lp_val) {
    if (nrels >= MAX_RELS) return;
    rel_t *r = &rels[nrels];
    mpz_init_set(r->Y, Y);
    r->exponents = malloc(fb_size * sizeof(int));
    memcpy(r->exponents, exp_vec, fb_size * sizeof(int));
    r->special_q = sq;
    r->sq_exp = 1;
    r->large_prime = lp_val;

    if (lp_val > 1) {
        /* LP hash: track pairs with same large prime for combining */
        int h = lp_val % LP_HASH_SZ;
        int idx = lp_hash[h];
        while (idx >= 0) {
            if (rels[idx].large_prime == lp_val) { ncombined++; break; }
            idx = lp_next[idx];
        }
        lp_next[nrels] = lp_hash[h];
        lp_hash[h] = nrels;
    }
    nrels++;
}

/*
 * Special-Q sieving:
 * For special prime Q with root r (QS_poly(r) = 0 mod Q):
 *   Sieve positions: x = r + k*Q for k = 0, +/-1, +/-2, ...
 *   At each position: evaluate QS_poly(x) / Q, test for B-smoothness
 *
 * QS polynomial: Q(x) = (x + ceil(sqrt(kN)))^2 - kN
 * Q(r) = 0 mod Q means (r + s)^2 = kN mod Q where s = ceil(sqrt(kN))
 * So r = +/-sqrt(kN mod Q) - s mod Q
 *
 * Relation: Y^2 = QS_val (mod N) where Y = (x + sqrt_kN) mod N
 *   QS_val = special_q * product(p_i^e_i) * [large_prime]
 */
static int specialq_sieve(mpz_t N, mpz_t kN, int special_q, int q_root,
                           int sieve_half, mpz_t sqrt_kN) {
    int rels_found = 0;
    uint8_t *sieve = calloc(SIEVE_BLOCK, 1);

    uint32_t Q_inv[MAX_FB]; /* Q^(-1) mod p for each FB prime */
    for (int i = 0; i < fb_size; i++) {
        int p = primes[i];
        if (p == special_q) { Q_inv[i] = 0; continue; }
        Q_inv[i] = modinv(special_q % p, p);
    }

    /* Sieve threshold */
    int bits_kN = mpz_sizeinbase(kN, 2);
    int bits_val = bits_kN / 2 + (int)log2(sieve_half) + (int)log2(special_q);
    int bits_reduced = bits_val - (int)log2(special_q);
    int threshold = (int)(bits_reduced * 0.7);

    memset(sieve, 0, SIEVE_BLOCK);
    int sieve_len = 2 * sieve_half;
    if (sieve_len > SIEVE_BLOCK) sieve_len = SIEVE_BLOCK;

    /* Sieve: for each FB prime, add logp at divisible positions */
    for (int i = 0; i < fb_size; i++) {
        int p = primes[i];
        if (p == special_q || p < 3) continue;
        int logp = fb_logp[i];
        int s = fb_roots[i];

        int sqrt_kN_mod_p = mpz_fdiv_ui(sqrt_kN, p);
        int root1 = ((s - sqrt_kN_mod_p) % p + p) % p;
        int root2 = ((-s - sqrt_kN_mod_p) % p + p) % p;

        int k1 = ((int64_t)(root1 - q_root % p + p) % p * Q_inv[i]) % p;
        int k2 = ((int64_t)(root2 - q_root % p + p) % p * Q_inv[i]) % p;

        int pos1 = (k1 + sieve_half) % p;
        int pos2 = (k2 + sieve_half) % p;

        for (int j = pos1; j < sieve_len; j += p) sieve[j] += logp;
        if (root1 != root2) {
            for (int j = pos2; j < sieve_len; j += p) sieve[j] += logp;
        }
    }

    /* Scan for candidates and build full relations */
    mpz_t qval, reduced, x_val, cofactor, Y_val;
    mpz_inits(qval, reduced, x_val, cofactor, Y_val, NULL);
    int *exp_vec = malloc(fb_size * sizeof(int));

    long lp_bound = (long)primes[fb_size - 1] * 30;

    for (int j = 0; j < sieve_len; j++) {
        if (sieve[j] < threshold) continue;
        if (nrels >= MAX_RELS) break;

        /* Compute x = q_root + (j - sieve_half) * special_q */
        int64_t k = j - sieve_half;
        mpz_set_si(x_val, k);
        mpz_mul_ui(x_val, x_val, special_q);
        mpz_add_ui(x_val, x_val, q_root);

        /* Y = x + sqrt_kN (mod N); QS_val = Y^2 - kN */
        mpz_add(Y_val, x_val, sqrt_kN);
        mpz_mul(qval, Y_val, Y_val);
        mpz_sub(qval, qval, kN);

        if (mpz_sgn(qval) <= 0) continue;

        /* Check Q divides QS_val */
        if (!mpz_divisible_ui_p(qval, special_q)) continue;

        /* Reduced value = QS_val / Q */
        mpz_divexact_ui(reduced, qval, special_q);

        /* Trial divide reduced value by factor base to get exponent vector */
        memset(exp_vec, 0, fb_size * sizeof(int));
        mpz_set(cofactor, reduced);

        for (int i = 0; i < fb_size && mpz_cmp_ui(cofactor, 1) > 0; i++) {
            while (mpz_divisible_ui_p(cofactor, primes[i])) {
                mpz_divexact_ui(cofactor, cofactor, primes[i]);
                exp_vec[i]++;
            }
        }

        int lp_val = 0;
        int is_smooth = 0;

        if (mpz_cmp_ui(cofactor, 1) == 0) {
            is_smooth = 1; /* fully smooth */
        } else if (mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= (unsigned long)lp_bound) {
            unsigned long cof = mpz_get_ui(cofactor);
            if (mpz_probab_prime_p(cofactor, 2)) {
                lp_val = (int)cof;
                is_smooth = 1; /* SLP */
            }
        }

        if (is_smooth) {
            /* Y mod N */
            mpz_mod(Y_val, Y_val, N);
            find_or_add_sq(special_q);
            add_relation(Y_val, exp_vec, special_q, lp_val);
            rels_found++;
        }
    }

    free(exp_vec);
    free(sieve);
    mpz_clears(qval, reduced, x_val, cofactor, Y_val, NULL);
    return rels_found;
}

/* ==================== LA: GF(2) Gaussian Elimination ==================== */

static int do_linalg_and_sqrt(mpz_t N) {
    /*
     * Each relation: Y_i^2 = special_q_i * prod(p_j^e_ij) * [lp_i] (mod N)
     *
     * Special-Q appears with exponent 1 (odd), so like LPs, we need >= 2
     * relations sharing the same special-Q for the column to be zeroed in
     * a dependency. Only include matched SQs and matched LPs as columns.
     *
     * Matrix columns: fb_size (FB) + matched_SQs + matched_LPs
     */

    /* Count special-Q occurrences */
    int *sq_cnt = calloc(n_distinct_sq, sizeof(int));
    for (int i = 0; i < nrels; i++)
        for (int j = 0; j < n_distinct_sq; j++)
            if (distinct_sq[j] == rels[i].special_q) { sq_cnt[j]++; break; }

    int *matched_sqs = malloc(n_distinct_sq * sizeof(int));
    int n_msqs = 0;
    for (int j = 0; j < n_distinct_sq; j++)
        if (sq_cnt[j] >= 2) matched_sqs[n_msqs++] = distinct_sq[j];
    free(sq_cnt);

    /* Collect distinct LPs appearing >= 2 times */
    int *distinct_lps = malloc(nrels * sizeof(int));
    int n_dlps = 0;
    for (int i = 0; i < nrels; i++) {
        if (rels[i].large_prime <= 1) continue;
        int lp = rels[i].large_prime;
        int found = 0;
        for (int j = 0; j < n_dlps; j++)
            if (distinct_lps[j] == lp) { found = 1; break; }
        if (!found) {
            int cnt = 0;
            for (int j = i; j < nrels; j++)
                if (rels[j].large_prime == lp) { cnt++; if (cnt >= 2) break; }
            if (cnt >= 2) distinct_lps[n_dlps++] = lp;
        }
    }

    /* A relation is usable if its SQ is matched AND its LP is matched (or full) */
    int *usable = calloc(nrels, sizeof(int));
    int n_usable = 0;
    for (int i = 0; i < nrels; i++) {
        int sq_ok = 0;
        for (int j = 0; j < n_msqs; j++)
            if (matched_sqs[j] == rels[i].special_q) { sq_ok = 1; break; }
        if (!sq_ok) continue;

        if (rels[i].large_prime <= 1) {
            usable[i] = 1; n_usable++;
        } else {
            for (int j = 0; j < n_dlps; j++)
                if (rels[i].large_prime == distinct_lps[j]) {
                    usable[i] = 1; n_usable++; break;
                }
        }
    }

    int ncols = fb_size + n_msqs + n_dlps;
    fprintf(stderr, "LA: %d usable rels, %d cols (fb=%d + %d matched SQs + %d LPs)\n",
            n_usable, ncols, fb_size, n_msqs, n_dlps);

    if (n_usable < ncols + 1) {
        fprintf(stderr, "Not enough relations for LA (%d < %d+1)\n", n_usable, ncols);
        free(matched_sqs); free(distinct_lps); free(usable);
        return 0;
    }

    /* Build GF(2) matrix with augmented identity */
    int rw = (ncols + 63) / 64;
    int aw = (n_usable + 63) / 64;
    uint64_t **mat = malloc(n_usable * sizeof(uint64_t *));
    uint64_t **aug = malloc(n_usable * sizeof(uint64_t *));
    int *rel_idx = malloc(n_usable * sizeof(int));

    int row = 0;
    for (int i = 0; i < nrels; i++) {
        if (!usable[i]) continue;
        mat[row] = calloc(rw, sizeof(uint64_t));
        aug[row] = calloc(aw, sizeof(uint64_t));
        aug[row][row / 64] |= (1ULL << (row % 64));

        /* Set FB exponent parities */
        for (int j = 0; j < fb_size; j++) {
            if (rels[i].exponents[j] & 1)
                mat[row][j / 64] |= (1ULL << (j % 64));
        }

        /* Set special-Q column (exponent 1 = always odd parity) */
        for (int j = 0; j < n_msqs; j++) {
            if (matched_sqs[j] == rels[i].special_q) {
                int col = fb_size + j;
                mat[row][col / 64] |= (1ULL << (col % 64));
                break;
            }
        }

        /* Set LP column */
        if (rels[i].large_prime > 1) {
            for (int j = 0; j < n_dlps; j++) {
                if (distinct_lps[j] == rels[i].large_prime) {
                    int col = fb_size + n_msqs + j;
                    mat[row][col / 64] |= (1ULL << (col % 64));
                    break;
                }
            }
        }

        rel_idx[row] = i;
        row++;
    }

    /* Gaussian elimination */
    int *pcol = malloc(n_usable * sizeof(int));
    memset(pcol, -1, n_usable * sizeof(int));

    for (int col = 0; col < ncols; col++) {
        int piv = -1;
        for (int r = 0; r < n_usable; r++) {
            if (pcol[r] >= 0) continue;
            if ((mat[r][col / 64] >> (col % 64)) & 1) { piv = r; break; }
        }
        if (piv < 0) continue;
        pcol[piv] = col;
        for (int r = 0; r < n_usable; r++) {
            if (r == piv) continue;
            if (!((mat[r][col / 64] >> (col % 64)) & 1)) continue;
            for (int w = 0; w < rw; w++) mat[r][w] ^= mat[piv][w];
            for (int w = 0; w < aw; w++) aug[r][w] ^= aug[piv][w];
        }
    }

    /* Extract and try dependencies (null space vectors) */
    mpz_t X, Z, factor, temp;
    mpz_init(X); mpz_init(Z); mpz_init(factor); mpz_init(temp);

    int factored = 0;
    int deps_tried = 0;

    for (int r = 0; r < n_usable && !factored; r++) {
        /* Check if zero row */
        int zero = 1;
        for (int w = 0; w < rw && zero; w++) if (mat[r][w]) zero = 0;
        if (!zero) continue;

        /* Extract dependency: which relations are in this null vector */
        int *total_exp = calloc(fb_size, sizeof(int));
        int *sq_counts = calloc(n_msqs, sizeof(int));
        int *lps = malloc(n_usable * sizeof(int));
        int nlps = 0;
        int dep_cnt = 0;

        mpz_set_ui(X, 1);

        for (int i = 0; i < n_usable; i++) {
            if (!((aug[r][i / 64] >> (i % 64)) & 1)) continue;
            dep_cnt++;
            int ri = rel_idx[i];

            /* X = product of Y_i mod N */
            mpz_mul(X, X, rels[ri].Y);
            mpz_mod(X, X, N);

            /* Accumulate exponents */
            for (int j = 0; j < fb_size; j++)
                total_exp[j] += rels[ri].exponents[j];

            /* Track special-Q occurrences */
            for (int j = 0; j < n_msqs; j++) {
                if (matched_sqs[j] == rels[ri].special_q) { sq_counts[j]++; break; }
            }

            /* Track LPs */
            if (rels[ri].large_prime > 1)
                lps[nlps++] = rels[ri].large_prime;
        }

        if (dep_cnt < 2) { free(total_exp); free(sq_counts); free(lps); continue; }

        /* Verify all exponents are even */
        int ok = 1;
        for (int j = 0; j < fb_size && ok; j++)
            if (total_exp[j] % 2) ok = 0;
        for (int j = 0; j < n_msqs && ok; j++)
            if (sq_counts[j] % 2) ok = 0;

        /* Check LP exponents even */
        if (ok && nlps > 0) {
            /* Sort LPs */
            for (int i = 0; i < nlps - 1; i++)
                for (int j = i + 1; j < nlps; j++)
                    if (lps[i] > lps[j]) { int t = lps[i]; lps[i] = lps[j]; lps[j] = t; }
            for (int i = 0; i < nlps && ok; ) {
                int cnt = 1;
                while (i + cnt < nlps && lps[i + cnt] == lps[i]) cnt++;
                if (cnt % 2) { ok = 0; break; }
                i += cnt;
            }
        }

        if (!ok) { free(total_exp); free(sq_counts); free(lps); continue; }

        /* Compute Z = sqrt(product of QS_vals) mod N
         * = prod(p_j^(total_exp[j]/2)) * prod(sq_k^(sq_count[k]/2)) * prod(lp^(count/2)) mod N */
        mpz_set_ui(Z, 1);

        /* FB primes contribution */
        for (int j = 0; j < fb_size; j++) {
            if (total_exp[j] == 0) continue;
            mpz_set_ui(temp, primes[j]);
            mpz_powm_ui(temp, temp, total_exp[j] / 2, N);
            mpz_mul(Z, Z, temp);
            mpz_mod(Z, Z, N);
        }

        /* Special-Q contribution */
        for (int j = 0; j < n_msqs; j++) {
            if (sq_counts[j] == 0) continue;
            mpz_set_ui(temp, matched_sqs[j]);
            mpz_powm_ui(temp, temp, sq_counts[j] / 2, N);
            mpz_mul(Z, Z, temp);
            mpz_mod(Z, Z, N);
        }

        /* LP contribution */
        if (nlps > 0) {
            for (int i = 0; i < nlps; ) {
                int cnt = 1;
                while (i + cnt < nlps && lps[i + cnt] == lps[i]) cnt++;
                mpz_set_ui(temp, lps[i]);
                mpz_powm_ui(temp, temp, cnt / 2, N);
                mpz_mul(Z, Z, temp);
                mpz_mod(Z, Z, N);
                i += cnt;
            }
        }

        /* Try gcd(X - Z, N) and gcd(X + Z, N) */
        deps_tried++;
        mpz_sub(temp, X, Z);
        mpz_gcd(factor, temp, N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            factored = 1;
        } else {
            mpz_add(temp, X, Z);
            mpz_gcd(factor, temp, N);
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0)
                factored = 1;
        }

        free(total_exp); free(sq_counts); free(lps);
    }

    fprintf(stderr, "LA: tried %d dependencies\n", deps_tried);

    /* Cleanup */
    for (int i = 0; i < n_usable; i++) { free(mat[i]); free(aug[i]); }
    free(mat); free(aug); free(rel_idx); free(pcol);
    free(matched_sqs); free(distinct_lps); free(usable);

    if (factored) {
        mpz_t cofn; mpz_init(cofn);
        mpz_divexact(cofn, N, factor);
        if (mpz_cmp(factor, cofn) > 0) mpz_swap(factor, cofn);
        gmp_printf("Factor: %Zd\n", factor);
        gmp_printf("Cofactor: %Zd\n", cofn);
        mpz_clear(cofn);
    } else {
        fprintf(stderr, "FAILED to find factor from %d dependencies\n", deps_tried);
    }

    mpz_clear(X); mpz_clear(Z); mpz_clear(factor); mpz_clear(temp);
    return factored;
}

/* ==================== Main ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [fb_size]\n", argv[0]);
        return 1;
    }

    mpz_t N, kN, sqrt_kN;
    mpz_inits(N, kN, sqrt_kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);
    int bits = mpz_sizeinbase(N, 2);

    /* Parameters */
    int fb_target = argc > 2 ? atoi(argv[2]) :
        (digits <= 30 ? 200 : digits <= 40 ? 500 : digits <= 50 ? 1500 :
         digits <= 60 ? 4000 : digits <= 70 ? 10000 : 30000);

    /* Use multiplier k=1 for simplicity */
    mpz_set(kN, N);
    mpz_sqrt(sqrt_kN, kN);
    mpz_add_ui(sqrt_kN, sqrt_kN, 1);

    fprintf(stderr, "N: %d digits, %d bits\n", digits, bits);

    /* Build factor base */
    build_fb(kN, fb_target);
    fprintf(stderr, "Factor base: %d primes, max %d\n", fb_size, primes[fb_size-1]);

    init_lp();

    /* We need fb_size + n_distinct_sq + n_dlps columns worth of relations.
     * Each special-Q adds 1 column, so we need many rels per Q.
     * Use a large sieve interval to maximize rels per Q.
     * Dynamically check if we have enough usable rels > columns. */
    int target_rels = fb_size * 3 + 100;

    /* Special-Q range: use primes from fb_max to wider range */
    int q_min = primes[fb_size - 1] + 1;
    int q_max = q_min * 50;
    /* Larger sieve interval = more rels per Q = fewer Q needed = fewer columns */
    int sieve_half = SIEVE_BLOCK / 2;

    fprintf(stderr, "Special-Q range: %d to %d\n", q_min, q_max);
    fprintf(stderr, "Target: %d relations\n", target_rels);

    /* Generate special-Q primes */
    int sq_limit = q_max + 1000;
    char *is_prime = calloc(sq_limit + 1, 1);
    for (int i = 2; i <= sq_limit; i++) is_prime[i] = 1;
    for (int i = 2; (long)i * i <= sq_limit; i++)
        if (is_prime[i]) for (int j = i * i; j <= sq_limit; j += i) is_prime[j] = 0;

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int q_count = 0;
    int total_rels = 0;
    int have_enough = 0;

    for (int q = q_min; q <= q_max && !have_enough; q++) {
        if (!is_prime[q]) continue;

        /* Find roots of kN mod q */
        uint32_t kn_mod_q = mpz_fdiv_ui(kN, q);
        uint32_t sq_root = modsqrt(kn_mod_q, q);
        if (sq_root == UINT32_MAX) continue;

        /* sqrt_kN mod q */
        uint32_t sqrt_mod = mpz_fdiv_ui(sqrt_kN, q);

        /* Two roots: r1 = sq_root - sqrt_mod, r2 = q - sq_root - sqrt_mod */
        int r1 = ((int64_t)sq_root - sqrt_mod + q) % q;
        int r2 = ((int64_t)q - sq_root - sqrt_mod + q) % q;

        int found1 = specialq_sieve(N, kN, q, r1, sieve_half, sqrt_kN);
        int found2 = (r1 != r2) ? specialq_sieve(N, kN, q, r2, sieve_half, sqrt_kN) : 0;
        total_rels += found1 + found2;
        q_count++;

        /* Periodically check if we have enough usable relations.
         * A rel is usable if its SQ has >= 2 rels AND its LP is full or matched.
         * Estimate: count rels whose SQ appears >= 2 times and LP=0. */
        if (q_count % 50 == 0 || total_rels >= target_rels) {
            /* Count SQ occurrences */
            int n_matched_sq_est = 0;
            int est_usable = 0;
            for (int si = 0; si < n_distinct_sq; si++) {
                int cnt = 0;
                for (int ri = 0; ri < nrels; ri++)
                    if (rels[ri].special_q == distinct_sq[si]) cnt++;
                if (cnt >= 2) {
                    n_matched_sq_est++;
                    /* Count usable rels from this SQ: full rels + LP-matched */
                    for (int ri = 0; ri < nrels; ri++) {
                        if (rels[ri].special_q != distinct_sq[si]) continue;
                        if (rels[ri].large_prime <= 1) est_usable++;
                        /* LP matching adds ~ncombined, approximate */
                    }
                }
            }
            est_usable += ncombined; /* rough: combined LP pairs from matched SQs */
            int est_cols = fb_size + n_matched_sq_est + 10;
            if (est_usable > est_cols + 20)
                have_enough = 1;

            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
            fprintf(stderr, "Q=%d, %d SQs(%d matched), %d rels, est_usable=%d, need>%d, %.1fs, %.1f rels/sec\n",
                    q, q_count, n_matched_sq_est, nrels, est_usable, est_cols,
                    elapsed, nrels / elapsed);
            if (elapsed > 280) break;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &now);
    double sieve_time = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "Sieving: %d rels from %d special-Qs in %.1fs (%.1f rels/sec)\n",
            nrels, q_count, sieve_time, nrels / sieve_time);
    fprintf(stderr, "Distinct special-Qs used: %d\n", n_distinct_sq);

    int success = 0;
    fprintf(stderr, "Starting linear algebra...\n");
    success = do_linalg_and_sqrt(N);

    clock_gettime(CLOCK_MONOTONIC, &now);
    double total_time = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "Total time: %.3fs (sieve %.3fs, LA %.3fs)\n",
            total_time, sieve_time, total_time - sieve_time);

    /* Cleanup */
    for (int i = 0; i < nrels; i++) {
        mpz_clear(rels[i].Y);
        free(rels[i].exponents);
    }
    free(is_prime);
    mpz_clears(N, kN, sqrt_kN, NULL);
    return success ? 0 : 1;
}
