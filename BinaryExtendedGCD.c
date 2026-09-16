//clear && gcc BinaryExtendedGCD.c -lm -o m.o && ./m.o
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>

typedef uint64_t ulong;
typedef int64_t slong;


ulong n_gcdinv(ulong *s, ulong x, ulong y)
{
    slong v1, v2, t2;
    ulong d, r, quot, rem;

    assert(y > x);

    v1 = 0;
    v2 = 1;
    r = x;
    x = y;

    if ((slong)(x & r) < 0) {
        d = x - r;
        t2 = v2;
        x = r;
        v2 = v1 - v2;
        v1 = t2;
        r = d;
    }

    while ((slong)(r << 1) < 0) {
        d = x - r;

        if (d < r) {
            t2 = v2;
            x = r;
            v2 = v1 - v2;
            v1 = t2;
            r = d;
        } else if (d < (r << 1)) {
            x = r;
            t2 = v2;
            v2 = v1 - ((ulong)v2 << 1);
            v1 = t2;
            r = d - x;
        } else {
            x = r;
            t2 = v2;
            v2 = v1 - 3 * v2;
            v1 = t2;
            r = d - (x << 1);
        }
    }

    while (r) {
        if (x < (r << 2)) {
            d = x - r;

            if (d < r) {
                t2 = v2;
                x = r;
                v2 = v1 - v2;
                v1 = t2;
                r = d;
            } else if (d < (r << 1)) {
                x = r;
                t2 = v2;
                v2 = v1 - ((ulong)v2 << 1);
                v1 = t2;
                r = d - x;
            } else {
                x = r;
                t2 = v2;
                v2 = v1 - 3 * (ulong)v2;
                v1 = t2;
                r = d - (x << 1);
            }
        } else {
            quot = x / r;
            rem = x - r * quot;

            x = r;
            t2 = v2;
            v2 = v1 - quot * v2;
            v1 = t2;
            r = rem;
        }
    }

    if (v1 < 0)
        v1 += (slong)y;

    *s = (ulong)v1;

    return x;
}

int main(void)
{
    long a = 123456789012345L;
    long b = 987654321098765L;

    const unsigned long N = 10000000UL;
    volatile long sink1 = 0;
    volatile ulong sink2 = 0;

    struct timespec start, end;

    /* Binary Extended GCD */
    clock_gettime(CLOCK_MONOTONIC, &start);

    for (unsigned long i = 0; i < N; i++) {
        egcd_result r = extendBinaryGCD(a, b);
        sink1 += r.s + r.t + r.d;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    double t1 = (end.tv_sec - start.tv_sec)
             + (end.tv_nsec - start.tv_nsec) / 1e9;

    /* n_gcdinv */
    clock_gettime(CLOCK_MONOTONIC, &start);

    for (unsigned long i = 0; i < N; i++) {
        ulong inv;
        ulong g = n_gcdinv(&inv, a, b);
        sink2 += inv + g;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);

    double t2 = (end.tv_sec - start.tv_sec)
             + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("Binary GCD : %.3f ns/call\n", t1 * 1e9 / N);
    printf("n_gcdinv   : %.3f ns/call\n", t2 * 1e9 / N);

    printf("Ratio      : %.2fx\n", t1 > t2 ? t1 / t2 : t2 / t1);

    return 0;
}

