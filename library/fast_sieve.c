/*
 * fast_sieve.c - AVX512BW sieve kernel benchmark
 *
 * Implements YAFU-style batched scattered byte subtraction using AVX512.
 * Processes 32 primes simultaneously, doing 64 byte subtractions per iteration.
 *
 * This is a standalone benchmark to measure raw sieve throughput.
 * Compare with YAFU's ~5000-6500 rels/sec for 90d.
 *
 * Compile: gcc -O3 -march=native -mavx512bw -o fast_sieve library/fast_sieve.c -lgmp -lm
 * Usage: ./fast_sieve <N> [benchmark_only]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <immintrin.h>
#include <gmp.h>

#define BLOCK_SIZE 32768
#define MAX_FB 200000

/* ==================== Sieve kernel ==================== */

/* Process 32 primes at once, YAFU-style:
 * 1. Load 32 root positions into AVX512 vector
 * 2. Store to temp array
 * 3. Do 64 scattered byte subtractions (32 for root1, 32 for root2)
 * 4. Update roots: root += prime, check if < BLOCK_SIZE
 */
static void sieve_batch32(uint8_t *sieve,
                          uint16_t *root1, uint16_t *root2,
                          const uint16_t *primes, uint8_t logp,
                          int count) {
    /* For each batch of 32 primes */
    __m512i vblock = _mm512_set1_epi16(BLOCK_SIZE);
    __m512i vzero = _mm512_setzero_si512();

    for (int i = 0; i < count; i += 32) {
        int batch = (count - i < 32) ? count - i : 32;

        __m512i vprime = _mm512_loadu_si512(primes + i);
        __m512i vroot1 = _mm512_loadu_si512(root1 + i);
        __m512i vroot2 = _mm512_loadu_si512(root2 + i);

        /* Mask for valid primes (non-zero) */
        __mmask32 valid = _mm512_cmpgt_epu16_mask(vprime, vzero);

        /* Set invalid roots to 0 (won't hurt since sieve[0] is expendable) */
        vroot1 = _mm512_mask_set1_epi16(vroot1, ~valid, 0);
        vroot2 = _mm512_mask_set1_epi16(vroot2, ~valid, 0);

        /* Sieve loop: while any root is in range */
        __mmask32 mask1 = _mm512_cmplt_epu16_mask(vroot1, vblock) & valid;
        __mmask32 mask2 = _mm512_cmplt_epu16_mask(vroot2, vblock) & valid;

        while (mask1 | mask2) {
            /* Store roots to temp arrays and do scattered subtractions */
            uint16_t r1[32] __attribute__((aligned(64)));
            uint16_t r2[32] __attribute__((aligned(64)));
            _mm512_store_si512(r1, vroot1);
            _mm512_store_si512(r2, vroot2);

            /* Subtract logp at each valid root position */
            if (mask1) {
                for (int j = 0; j < 32; j++) {
                    if (mask1 & (1u << j))
                        sieve[r1[j]] -= logp;
                }
            }
            if (mask2) {
                for (int j = 0; j < 32; j++) {
                    if (mask2 & (1u << j))
                        sieve[r2[j]] -= logp;
                }
            }

            /* Advance roots */
            vroot1 = _mm512_mask_add_epi16(vroot1, mask1, vroot1, vprime);
            vroot2 = _mm512_mask_add_epi16(vroot2, mask2, vroot2, vprime);

            /* Update masks */
            mask1 = _mm512_cmplt_epu16_mask(vroot1, vblock) & valid;
            mask2 = _mm512_cmplt_epu16_mask(vroot2, vblock) & valid;
        }

        /* Store updated roots back */
        _mm512_storeu_si512(root1 + i, vroot1);
        _mm512_storeu_si512(root2 + i, vroot2);
    }
}

/* Simple scalar sieve for comparison */
static void sieve_scalar(uint8_t *sieve, int block_size,
                         uint16_t *root1, uint16_t *root2,
                         const uint16_t *primes, const uint8_t *logps,
                         int count) {
    for (int i = 0; i < count; i++) {
        int p = primes[i];
        if (p == 0) continue;
        uint8_t logp = logps[i];
        int r1 = root1[i], r2 = root2[i];

        for (int j = r1; j < block_size; j += p)
            sieve[j] -= logp;
        root1[i] = r1 + ((block_size - r1 + p - 1) / p) * p - block_size;

        if (r2 != r1) {
            for (int j = r2; j < block_size; j += p)
                sieve[j] -= logp;
        }
        root2[i] = r2 + ((block_size - r2 + p - 1) / p) * p - block_size;
    }
}

/* Benchmark: measure sieve kernel throughput */
int main(int argc, char **argv) {
    int fb_size = 10000; /* Factor base primes to sieve */
    int num_blocks = 10; /* Number of blocks to sieve */

    if (argc > 1) fb_size = atoi(argv[1]);
    if (argc > 2) num_blocks = atoi(argv[2]);

    printf("Benchmarking sieve kernel: %d primes, %d blocks of %dB\n",
           fb_size, num_blocks, BLOCK_SIZE);

    /* Generate fake factor base */
    /* Round up to multiple of 32 for AVX512 alignment */
    int alloc_size = ((fb_size + 31) / 32) * 32;
    uint16_t *primes = aligned_alloc(64, alloc_size * sizeof(uint16_t));
    uint16_t *root1 = aligned_alloc(64, alloc_size * sizeof(uint16_t));
    uint16_t *root2 = aligned_alloc(64, alloc_size * sizeof(uint16_t));
    uint8_t *logps = aligned_alloc(64, alloc_size);
    memset(primes, 0, alloc_size * sizeof(uint16_t));
    memset(root1, 0, alloc_size * sizeof(uint16_t));
    memset(root2, 0, alloc_size * sizeof(uint16_t));
    memset(logps, 0, alloc_size);

    srand(42);
    int p = 3;
    for (int i = 0; i < fb_size; i++) {
        while (1) { p += 2; int is_prime = 1; for (int d = 3; d*d <= p; d += 2) if (p%d==0) {is_prime=0;break;} if (is_prime) break; }
        primes[i] = (p < 65536) ? p : 0;
        root1[i] = rand() % (primes[i] ? primes[i] : 1);
        root2[i] = rand() % (primes[i] ? primes[i] : 1);
        logps[i] = (int)(log(p) * 1.44 + 0.5);
        if (logps[i] < 1) logps[i] = 1;
    }

    uint8_t *sieve = aligned_alloc(64, BLOCK_SIZE);

    /* Benchmark scalar sieve */
    {
        memset(sieve, 255, BLOCK_SIZE);
        /* Save roots */
        uint16_t *r1_save = malloc(fb_size * sizeof(uint16_t));
        uint16_t *r2_save = malloc(fb_size * sizeof(uint16_t));
        memcpy(r1_save, root1, fb_size * sizeof(uint16_t));
        memcpy(r2_save, root2, fb_size * sizeof(uint16_t));

        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        for (int b = 0; b < num_blocks; b++) {
            memset(sieve, 255, BLOCK_SIZE);
            sieve_scalar(sieve, BLOCK_SIZE, root1, root2, primes, logps, fb_size);
        }
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        printf("Scalar: %.3f ms/block (%.1f blocks/s)\n",
               elapsed * 1000 / num_blocks, num_blocks / elapsed);

        /* Restore roots */
        memcpy(root1, r1_save, fb_size * sizeof(uint16_t));
        memcpy(root2, r2_save, fb_size * sizeof(uint16_t));
        free(r1_save); free(r2_save);
    }

    /* Benchmark AVX512 sieve */
    {
        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);

        /* Count small primes (p < BLOCK_SIZE) */
        int small_count = 0;
        for (int i = 0; i < fb_size; i++)
            if (primes[i] > 0 && primes[i] < BLOCK_SIZE) small_count++;

        for (int b = 0; b < num_blocks; b++) {
            memset(sieve, 255, BLOCK_SIZE);
            /* Use logp from midpoint of batch */
            for (int i = 0; i < fb_size; i += 32) {
                int batch_end = (i + 32 < fb_size) ? i + 32 : fb_size;
                uint8_t logp = logps[i + (batch_end - i) / 2];
                sieve_batch32(sieve, root1 + i, root2 + i, primes + i, logp, batch_end - i);
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        printf("AVX512: %.3f ms/block (%.1f blocks/s)\n",
               elapsed * 1000 / num_blocks, num_blocks / elapsed);
    }

    /* Verify: compare scalar and AVX512 results on single block */
    {
        uint8_t *sieve_scalar_out = aligned_alloc(64, BLOCK_SIZE);
        uint8_t *sieve_avx_out = aligned_alloc(64, BLOCK_SIZE);
        uint16_t *r1a = aligned_alloc(64, alloc_size * sizeof(uint16_t));
        uint16_t *r2a = aligned_alloc(64, alloc_size * sizeof(uint16_t));
        uint16_t *r1b = aligned_alloc(64, alloc_size * sizeof(uint16_t));
        uint16_t *r2b = aligned_alloc(64, alloc_size * sizeof(uint16_t));

        /* Reset roots */
        for (int i = 0; i < fb_size; i++) {
            r1a[i] = r1b[i] = rand() % (primes[i] ? primes[i] : 1);
            r2a[i] = r2b[i] = rand() % (primes[i] ? primes[i] : 1);
        }

        memset(sieve_scalar_out, 255, BLOCK_SIZE);
        sieve_scalar(sieve_scalar_out, BLOCK_SIZE, r1a, r2a, primes, logps, fb_size);

        memset(sieve_avx_out, 255, BLOCK_SIZE);
        for (int i = 0; i < fb_size; i += 32) {
            int batch_end = (i + 32 < fb_size) ? i + 32 : fb_size;
            uint8_t logp = logps[i + (batch_end - i) / 2];
            sieve_batch32(sieve_avx_out, r1b + i, r2b + i, primes + i, logp, batch_end - i);
        }

        int mismatches = 0;
        for (int i = 0; i < BLOCK_SIZE; i++) {
            if (sieve_scalar_out[i] != sieve_avx_out[i]) {
                mismatches++;
                if (mismatches <= 5)
                    printf("MISMATCH at %d: scalar=%d avx=%d\n", i, sieve_scalar_out[i], sieve_avx_out[i]);
            }
        }
        printf("Verification: %d mismatches out of %d positions\n", mismatches, BLOCK_SIZE);

        free(sieve_scalar_out); free(sieve_avx_out);
        free(r1a); free(r2a); free(r1b); free(r2b);
    }

    free(primes); free(root1); free(root2); free(logps); free(sieve);
    return 0;
}
