#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

static inline __m256i avx_mm256_div_epi16_division(
    __m256i a_epi16,
    __m256i b_epi16)

{
    const __m256i a_hi_epi32       = _mm256_srai_epi32(a_epi16, 16);
    const __m256i a_lo_epi32_shift = _mm256_slli_epi32(a_epi16, 16);
    const __m256i a_lo_epi32       = _mm256_srai_epi32(a_lo_epi32_shift, 16);

    const __m256i b_hi_epi32       = _mm256_srai_epi32(b_epi16, 16);
    const __m256i b_lo_epi32_shift = _mm256_slli_epi32(b_epi16, 16);
    const __m256i b_lo_epi32       = _mm256_srai_epi32(b_lo_epi32_shift, 16);

    const __m256 a_hi = _mm256_cvtepi32_ps(a_hi_epi32);
    const __m256 a_lo = _mm256_cvtepi32_ps(a_lo_epi32);
    const __m256 b_hi = _mm256_cvtepi32_ps(b_hi_epi32);
    const __m256 b_lo = _mm256_cvtepi32_ps(b_lo_epi32);

    const __m256 hi = _mm256_div_ps(a_hi, b_hi);
    const __m256 lo = _mm256_div_ps(a_lo, b_lo);

    const __m256i hi_epi32 = _mm256_cvttps_epi32(hi);
    const __m256i lo_epi32 = _mm256_cvttps_epi32(lo);

    const __m256i hi_epi32_shift = _mm256_slli_epi32(hi_epi32, 16);

    return _mm256_blend_epi16(lo_epi32, hi_epi32_shift, 0xAA);
}


static inline uint64_t now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
}


static uint64_t rng_state = 0x123456789abcdefULL;

static uint32_t rng32(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (uint32_t)rng_state;
}


int main(void)
{
    const size_t N = 10000000;

    int16_t *a = aligned_alloc(32, N * sizeof(int16_t));
    int16_t *b = aligned_alloc(32, N * sizeof(int16_t));
    int16_t *out = aligned_alloc(32, N * sizeof(int16_t));

    if (!a || !b || !out) {
        perror("allocation");
        return 1;
    }

    /*
     * Generate random values.
     *
     * Avoid zero divisors.
     */
    for (size_t i = 0; i < N; i++) {
        a[i] = (int16_t)rng32();

        do {
            b[i] = (int16_t)rng32();
        } while (b[i] == 0);
    }

    /*
     * Correctness check.
     */
    for (size_t i = 0; i < N; i += 16) {
        __m256i va = _mm256_load_si256((const __m256i *)&a[i]);
        __m256i vb = _mm256_load_si256((const __m256i *)&b[i]);

        __m256i vr = avx_mm256_div_epi16_division(va, vb);

        _mm256_store_si256((__m256i *)&out[i], vr);
    }

    for (size_t i = 0; i < N; i++) {
        int16_t expected = a[i] / b[i];

        if (out[i] != expected) {
            printf("ERROR at %zu: %d / %d = %d, expected %d\n",
                   i, a[i], b[i], out[i], expected);
            return 1;
        }
    }

    printf("Correctness: PASS\n");

    /*
     * Warm up.
     */
    volatile int16_t sink = 0;

    for (size_t i = 0; i < N; i += 16) {
        __m256i va = _mm256_load_si256((const __m256i *)&a[i]);
        __m256i vb = _mm256_load_si256((const __m256i *)&b[i]);

        __m256i vr = avx_mm256_div_epi16_division(va, vb);

        _mm256_store_si256((__m256i *)&out[i], vr);
    }

    /*
     * SIMD benchmark.
     */
    uint64_t start = now_ns();

    for (size_t i = 0; i < N; i += 16) {
        __m256i va = _mm256_load_si256((const __m256i *)&a[i]);
        __m256i vb = _mm256_load_si256((const __m256i *)&b[i]);

        __m256i vr = avx_mm256_div_epi16_division(va, vb);

        _mm256_store_si256((__m256i *)&out[i], vr);
    }

    uint64_t end = now_ns();

    sink ^= out[N - 1];

    double simd_seconds = (end - start) / 1e9;

    /*
     * Scalar benchmark.
     */
    start = now_ns();

    for (size_t i = 0; i < N; i++) {
        out[i] = a[i] / b[i];
    }

    end = now_ns();

    sink ^= out[N - 1];

    double scalar_seconds = (end - start) / 1e9;

    printf("\n");
    printf("N:              %zu\n", N);
    printf("SIMD:           %.6f s\n", simd_seconds);
    printf("Scalar:         %.6f s\n", scalar_seconds);
    printf("SIMD / scalar:  %.3fx\n",
           simd_seconds / scalar_seconds);
    printf("Scalar / SIMD:  %.3fx\n",
           scalar_seconds / simd_seconds);
    printf("checksum:       %d\n", sink);

    free(a);
    free(b);
    free(out);

    return 0;
}
//clear && gcc -O3 -march=native main.c -o m.o && ./m.o
