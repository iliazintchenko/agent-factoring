/*
 * cf_products.c — Experimental: Continued fraction product candidates
 *
 * HYPOTHESIS: Products of CF convergents of √N have residues that are
 * PRODUCTS of individual CF residues. If individual residues are bounded
 * by 2√N (Thue), then products of k convergents have residues bounded
 * by (2√N)^k. But the SQUARE of the product is:
 *
 *   (p_i * p_j)^2 - N * (q_i * q_j)^2 = (p_i^2 - N*q_i^2)(p_j^2 + N*q_j^2) + ...
 *
 * Wait, that's not right. Let me think more carefully.
 *
 * For convergent k: p_k^2 - N*q_k^2 = (-1)^k * A_k where A_k ≈ √N.
 *
 * For the PRODUCT p_i * p_j:
 *   (p_i * p_j)^2 = p_i^2 * p_j^2
 *   (q_i * q_j)^2 * N = q_i^2 * q_j^2 * N
 *
 *   (p_i*p_j)^2 - N*(q_i*q_j)^2
 *     = (p_i^2 - N*q_i^2)(p_j^2 - N*q_j^2) + N*(p_i*q_j + p_j*q_i)*(p_i*q_j - p_j*q_i) + ...
 *
 * Hmm, this expansion isn't clean. Let me use a different identity.
 *
 * Brahmagupta-Fibonacci identity:
 *   (a² + nb²)(c² + nd²) = (ac ± nbd)² + n(ad ∓ bc)²
 *
 * For our norms: p² - Nq² = A (a "Pell-like" equation)
 * Product of two norms: A₁ * A₂ = (p₁² - Nq₁²)(p₂² - Nq₂²)
 *   = (p₁p₂ + Nq₁q₂)² - N(p₁q₂ + p₂q₁)²
 *   = (p₁p₂ - Nq₁q₂)² - N(p₁q₂ - p₂q₁)²
 *
 * So if we set:
 *   X = p₁p₂ + Nq₁q₂   (or p₁p₂ - Nq₁q₂)
 *   Y = p₁q₂ + p₂q₁     (or p₁q₂ - p₂q₁)
 *
 * Then X² - N*Y² = A₁ * A₂
 *
 * Since A₁ and A₂ are each ≈ √N:
 *   A₁ * A₂ ≈ N
 *
 * And X ≈ p₁p₂ + Nq₁q₂ ≈ 2N (since p ≈ √N*q and p ≈ √N)
 * This gives X² - N*Y² ≈ N, so the residue is about N. WORSE than QS.
 *
 * But what about the OTHER choice: X = p₁p₂ - Nq₁q₂?
 *   X = p₁p₂ - Nq₁q₂, Y = p₁q₂ - p₂q₁
 *
 * p₁p₂ ≈ N (since p ≈ √N), Nq₁q₂ ≈ N (since q ≈ √N/√N = 1... no, q grows)
 * Actually convergents have p_k/q_k → √N, so p_k ≈ √N * q_k.
 * p₁p₂ ≈ N * q₁q₂, so X = p₁p₂ - Nq₁q₂ = (p₁p₂/q₁q₂ - N) * q₁q₂ ≈ 0
 *
 * More precisely: X = p₁p₂ - Nq₁q₂ and Y = p₁q₂ - p₂q₁.
 * p₁/q₁ ≈ √N + ε₁, p₂/q₂ ≈ √N + ε₂ where ε_k = O(1/q_k).
 * So p₁p₂ = q₁q₂(√N + ε₁)(√N + ε₂) = q₁q₂(N + (ε₁+ε₂)√N + ε₁ε₂)
 * X = q₁q₂(ε₁+ε₂)√N + q₁q₂ε₁ε₂ ≈ q₁q₂(ε₁+ε₂)√N
 * Since ε_k ≈ 1/q_k: X ≈ (q₂/q₁ + q₁/q₂)√N ≈ 2√N (for q₁ ≈ q₂)
 *
 * So X ≈ 2√N and Y ≈ ε₁q₂ - ε₂q₁ ≈ (1/q₁)q₂ - (1/q₂)q₁ ≈ 0 (approximately).
 * Actually Y = p₁q₂ - p₂q₁ = det(matrix of convergents). For consecutive
 * convergents: |p_{k+1}q_k - p_kq_{k+1}| = 1 (by CF theory!).
 *
 * So Y = ±1 for consecutive convergents! And X ≈ 2√N.
 * Then X² - N*Y² = A₁*A₂ means (2√N)² - N*1 = 4N - N = 3N ≈ A₁*A₂ ≈ N. ✓
 *
 * For CONSECUTIVE convergents: A₁ * A₂ = X² - N where X ≈ 2√N. Same as QS.
 *
 * For NON-consecutive convergents (far apart): Y = p₁q₂ - p₂q₁ is large.
 *
 * Conclusion: CF product candidates give residues of the same size as QS.
 * No improvement in residue size.
 *
 * But there's a DIFFERENT advantage: A₁*A₂ is AUTOMATICALLY a product of
 * two smaller numbers (A₁ ≈ √N and A₂ ≈ √N). If both A₁ and A₂ are
 * partially smooth (not fully smooth), their PRODUCT might be fully smooth.
 *
 * This is essentially the CFRAC method with product combination.
 * It's the foundation of CFRAC and QS. Not novel.
 *
 * Let me test it anyway to verify the residue sizes experimentally.
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N, sqrtN;
    mpz_init(N); mpz_init(sqrtN);
    mpz_set_str(N, argv[1], 10);
    mpz_sqrt(sqrtN, N);

    int digits = strlen(argv[1]);

    /* Generate CF convergents of √N */
    int max_conv = 200;
    mpz_t *p_arr = malloc(max_conv * sizeof(mpz_t));
    mpz_t *q_arr = malloc(max_conv * sizeof(mpz_t));
    mpz_t *A_arr = malloc(max_conv * sizeof(mpz_t)); /* p^2 - N*q^2 */

    for (int i = 0; i < max_conv; i++) {
        mpz_init(p_arr[i]); mpz_init(q_arr[i]); mpz_init(A_arr[i]);
    }

    /* CF expansion of √N: a0 = ⌊√N⌋, then regular CF */
    mpz_t a, P, Q, tmp;
    mpz_init(a); mpz_init(P); mpz_init(Q); mpz_init(tmp);

    /* Initial: p_{-1} = 1, p_0 = a0; q_{-1} = 0, q_0 = 1 */
    mpz_set(a, sqrtN); /* a0 = ⌊√N⌋ */
    mpz_set(p_arr[0], a);
    mpz_set_ui(q_arr[0], 1);

    /* CF recurrence: a_k, P_k, Q_k */
    mpz_set_ui(P, 0);  /* P_0 = 0 */
    mpz_set_ui(Q, 1);  /* Q_0 = 1 */

    mpz_t p_prev, q_prev;
    mpz_init_set_ui(p_prev, 1);
    mpz_init_set_ui(q_prev, 0);

    int n_conv = 1;
    /* Compute A_0 = p_0^2 - N*q_0^2 */
    mpz_mul(A_arr[0], p_arr[0], p_arr[0]);
    mpz_submul(A_arr[0], N, q_arr[0]);

    for (int k = 1; k < max_conv; k++) {
        /* P_k = a_{k-1} * Q_{k-1} - P_{k-1} */
        mpz_mul(tmp, a, Q);
        mpz_sub(P, tmp, P);

        /* Q_k = (N - P_k^2) / Q_{k-1} */
        mpz_mul(tmp, P, P);
        mpz_sub(tmp, N, tmp);
        mpz_divexact(Q, tmp, Q);

        if (mpz_sgn(Q) == 0) break;

        /* a_k = ⌊(⌊√N⌋ + P_k) / Q_k⌋ */
        mpz_add(tmp, sqrtN, P);
        mpz_fdiv_q(a, tmp, Q);

        /* p_k = a_k * p_{k-1} + p_{k-2} */
        mpz_mul(p_arr[k], a, p_arr[k-1]);
        mpz_add(p_arr[k], p_arr[k], p_prev);

        /* q_k = a_k * q_{k-1} + q_{k-2} */
        mpz_mul(q_arr[k], a, q_arr[k-1]);
        mpz_add(q_arr[k], q_arr[k], q_prev);

        /* Update prev */
        mpz_set(p_prev, p_arr[k-1]);
        mpz_set(q_prev, q_arr[k-1]);

        /* A_k = p_k^2 - N*q_k^2 = (-1)^k * Q_k (a known identity) */
        mpz_mul(A_arr[k], p_arr[k], p_arr[k]);
        mpz_submul(A_arr[k], N, q_arr[k]);

        n_conv++;
    }

    printf("# %d-digit N, %d CF convergents\n", digits, n_conv);
    printf("# Individual |A_k| sizes (bits): ");
    for (int k = 0; k < n_conv && k < 20; k++) {
        printf("%ld ", (long)mpz_sizeinbase(A_arr[k], 2));
    }
    printf("...\n");

    /* Now compute products of PAIRS of A values and check sizes */
    printf("# Product residues |A_i * A_j| sizes for pairs:\n");

    /* Also compute the Brahmagupta identity version */
    mpz_t X, Y, residue;
    mpz_init(X); mpz_init(Y); mpz_init(residue);

    int n_tested = 0, n_smooth = 0;
    for (int i = 0; i < n_conv && i < 100; i++) {
        for (int j = i + 1; j < n_conv && j < 100; j++) {
            /* Brahmagupta: X = p_i*p_j - N*q_i*q_j */
            mpz_mul(X, p_arr[i], p_arr[j]);
            mpz_mul(tmp, q_arr[i], q_arr[j]);
            mpz_submul(X, N, tmp);

            /* Y = p_i*q_j - p_j*q_i */
            mpz_mul(Y, p_arr[i], q_arr[j]);
            mpz_submul(Y, p_arr[j], q_arr[i]);

            /* residue = X^2 - N*Y^2 = A_i * A_j */
            mpz_mul(residue, A_arr[i], A_arr[j]);

            int res_bits = mpz_sizeinbase(residue, 2);
            int x_bits = mpz_sizeinbase(X, 2);
            int y_bits = mpz_sizeinbase(Y, 2);

            n_tested++;

            /* Check if the Brahmagupta candidate (X, Y) gives smaller values */
            if (n_tested <= 10) {
                printf("#   pair (%d,%d): |X|=%d bits, |Y|=%d bits, |A_i*A_j|=%d bits, "
                       "|A_i|=%ld, |A_j|=%ld\n",
                       i, j, x_bits, y_bits, res_bits,
                       (long)mpz_sizeinbase(A_arr[i], 2),
                       (long)mpz_sizeinbase(A_arr[j], 2));
            }

            /* Test smoothness of A_i * A_j with small primes */
            mpz_abs(tmp, residue);
            for (unsigned long p = 2; p < 10000; p++) {
                while (mpz_divisible_ui_p(tmp, p))
                    mpz_divexact_ui(tmp, tmp, p);
            }
            if (mpz_cmp_ui(tmp, 1) == 0) n_smooth++;
        }
    }

    printf("# Tested %d pairs, %d fully 10000-smooth products\n", n_tested, n_smooth);

    /* Compare with individual convergent smoothness */
    int n_indiv_smooth = 0;
    for (int k = 0; k < n_conv && k < 200; k++) {
        mpz_abs(tmp, A_arr[k]);
        for (unsigned long p = 2; p < 10000; p++) {
            while (mpz_divisible_ui_p(tmp, p))
                mpz_divexact_ui(tmp, tmp, p);
        }
        if (mpz_cmp_ui(tmp, 1) == 0) n_indiv_smooth++;
    }
    printf("# Individual: %d/%d convergents are 10000-smooth\n", n_indiv_smooth, n_conv);

    /* Cleanup */
    for (int i = 0; i < max_conv; i++) {
        mpz_clear(p_arr[i]); mpz_clear(q_arr[i]); mpz_clear(A_arr[i]);
    }
    free(p_arr); free(q_arr); free(A_arr);
    mpz_clear(N); mpz_clear(sqrtN);
    mpz_clear(a); mpz_clear(P); mpz_clear(Q); mpz_clear(tmp);
    mpz_clear(p_prev); mpz_clear(q_prev);
    mpz_clear(X); mpz_clear(Y); mpz_clear(residue);
    return 0;
}
