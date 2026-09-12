/*
    mv-bench.c -- microbenchmark for the machine_vectors.h backends.

    It times the modular primitives that fft_small spends its time in, and
    flint_dgemm_fallback/flint_sgemm_fallback, with whichever backend the
    compile flags select. gemm.c is compiled into this file (with its
    exported names redefined, so nothing clashes with the library), which is
    what makes the gemm numbers reflect the backend chosen here rather than
    the one FLINT was built with.

    Build it from the top of a FLINT source tree that has been configured and
    built (so that src/flint-config.h exists and libflint is present), once
    per backend:

        SRC="-I src -I . -std=c11 -O3"
        LIB="-L. -lflint -lgmp -lmpfr -lpthread -lm"

        gcc $SRC -march=native            mv-bench.c -o bench-native $LIB
        gcc $SRC -march=native -DFLINT_MACHINE_VECTORS_FORCE_GENERIC \
                                          mv-bench.c -o bench-gnu-native $LIB
        gcc $SRC -march=native -DFLINT_MACHINE_VECTORS_STRICT_C \
                                          mv-bench.c -o bench-strict-native $LIB
        gcc $SRC                          mv-bench.c -o bench-gnu-base $LIB
        gcc $SRC -DFLINT_MACHINE_VECTORS_STRICT_C \
                                          mv-bench.c -o bench-strict-base $LIB

    and run each (LD_LIBRARY_PATH=. ./bench-native, ...). The first three
    compare the three backends at the machine's full instruction set; the last
    two are the portable tiers restricted to the x86-64 baseline, which is what
    a distribution build without -march would get. Each binary prints which
    backend it ended up with, so a mix-up is visible.

    Do not build this with -ffast-math: the modular arithmetic depends on
    exact IEEE rounding, and the timings would not mean anything.

    Optional: -DBENCH_SECONDS=n (default 0.3) for the time given to each
    measurement, -DBENCH_MAXGEMM=n (default 512) for the largest gemm.
*/

#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
    gemm.c compiled into this translation unit under different names, so that
    it uses the backend selected above; the library's own copy is untouched.
*/
#define flint_sgemm            bench_sgemm
#define flint_dgemm            bench_dgemm
#define flint_sgemm_blas       bench_sgemm_blas
#define flint_dgemm_blas       bench_dgemm_blas
#define flint_sgemm_fallback   bench_sgemm_fallback
#define flint_dgemm_fallback   bench_dgemm_fallback
#define flint_gemm_use_blas    bench_gemm_use_blas

#include "machine_vectors/gemm.c"

#include "ulong_extras.h"

#ifndef BENCH_SECONDS
#define BENCH_SECONDS 0.3
#endif

#ifndef BENCH_MAXGEMM
#define BENCH_MAXGEMM 512
#endif

#define LEN 4096

/* a cheap deterministic sequence, so that runs are comparable */
static ulong
bench_rand(ulong i)
{
    ulong x = i * UWORD(0x9e3779b97f4a7c15) + UWORD(12345);
    x ^= x >> 30; x *= UWORD(0xbf58476d1ce4e5b9);
    x ^= x >> 27; x *= UWORD(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static double
wall(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
}

/*
    Run body repeatedly, doubling the repeat count until it has run for
    BENCH_SECONDS, and return the time of one repetition. The bodies below all
    update their input array in place, so nothing can be hoisted out of the
    loop, and each array entry carries its own dependency chain: with LEN
    entries there is enough of them to measure throughput rather than latency.
*/
#define TIME(body, per_rep) \
    do { \
        slong _r, _reps = 1; \
        double _t0, _t1; \
        while (1) \
        { \
            _t0 = wall(); \
            for (_r = 0; _r < _reps; _r++) { body; } \
            _t1 = wall(); \
            if (_t1 - _t0 >= BENCH_SECONDS) \
                break; \
            _reps *= 2; \
        } \
        (per_rep) = (_t1 - _t0) / (double) _reps; \
    } while (0)

/* ns per lane, given the time for one pass over LEN lanes */
#define PER_LANE(t) (1e9 * (t) / (double) LEN)

static void
bench_primitives(void)
{
    double * x = flint_aligned_alloc(64, LEN * sizeof(double));
    double * y = flint_aligned_alloc(64, LEN * sizeof(double));
    ulong p = (UWORD(1) << 49) + 2 * 1234567 + 1;
    vec4d n = vec4d_set_d((double) p);
    vec4d ninv = vec4d_set_d(1.0 / (double) p);
    vec8d n8 = vec8d_set_d((double) p);
    vec8d ninv8 = vec8d_set_d(1.0 / (double) p);
    double t, checksum = 0;
    slong i;

    for (i = 0; i < LEN; i++)
    {
        x[i] = (double) (bench_rand(i) % p);
        y[i] = (double) (bench_rand(i + 1) % p);
    }

    flint_printf("  %-18s %10s %10s\n", "", "ns/lane", "ns/vector");

    TIME(for (i = 0; i < LEN; i += 4)
            vec4d_store(x + i, vec4d_mulmod(vec4d_load(x + i),
                        vec4d_load(y + i), n, ninv)), t);
    flint_printf("  %-18s %10.3f %10.3f\n", "vec4d_mulmod",
            PER_LANE(t), 4 * PER_LANE(t));

    TIME(for (i = 0; i < LEN; i += 4)
            vec4d_store(x + i, vec4d_reduce_to_pm1n(vec4d_load(x + i),
                        n, ninv)), t);
    flint_printf("  %-18s %10.3f %10.3f\n", "vec4d_reduce_to_pm1n",
            PER_LANE(t), 4 * PER_LANE(t));

    TIME(for (i = 0; i < LEN; i += 4)
            vec4d_store(x + i, vec4d_reduce_to_0n(vec4d_load(x + i),
                        n, ninv)), t);
    flint_printf("  %-18s %10.3f %10.3f\n", "vec4d_reduce_to_0n",
            PER_LANE(t), 4 * PER_LANE(t));

    TIME(for (i = 0; i < LEN; i += 4)
            vec4d_store(x + i, vec4d_round(vec4d_mul(vec4d_load(x + i),
                        ninv))), t);
    flint_printf("  %-18s %10.3f %10.3f\n", "vec4d_round(mul)",
            PER_LANE(t), 4 * PER_LANE(t));

    for (i = 0; i < LEN; i++)
        x[i] = (double) (bench_rand(i) % p);

    TIME(for (i = 0; i < LEN; i += 4)
            vec4d_store(x + i, vec4d_fmadd(vec4d_load(x + i),
                        vec4d_set_d(1.0000001), vec4d_load(y + i))), t);
    flint_printf("  %-18s %10.3f %10.3f\n", "vec4d_fmadd",
            PER_LANE(t), 4 * PER_LANE(t));

    for (i = 0; i < LEN; i++)
        x[i] = (double) (bench_rand(i) % p);

    /* the pair type, which is how fft_small gets its instruction level
       parallelism on the transform inner loops */
    TIME(for (i = 0; i < LEN; i += 8)
            vec8d_store(x + i, vec8d_mulmod(vec8d_load(x + i),
                        vec8d_load(y + i), n8, ninv8)), t);
    flint_printf("  %-18s %10.3f %10.3f\n", "vec8d_mulmod",
            PER_LANE(t), 8 * PER_LANE(t));

    TIME(for (i = 0; i < LEN; i += 8)
            vec8d_store(x + i, vec8d_reduce_to_pm1n(vec8d_load(x + i),
                        n8, ninv8)), t);
    flint_printf("  %-18s %10.3f %10.3f\n", "vec8d_reduce_to_pm1n",
            PER_LANE(t), 8 * PER_LANE(t));

    for (i = 0; i < LEN; i++)
        checksum += x[i];
    flint_printf("  (checksum %g)\n", checksum);

    flint_free(x);
    flint_free(y);
}

static void
bench_gemm(void)
{
    slong sz;

    flint_printf("  %-8s %12s %12s %12s %12s\n", "n",
            "dgemm ms", "dgemm GF/s", "sgemm ms", "sgemm GF/s");

    for (sz = 64; sz <= BENCH_MAXGEMM; sz *= 2)
    {
        double * A = flint_malloc(sz * sz * sizeof(double));
        double * B = flint_malloc(sz * sz * sizeof(double));
        double * C = flint_malloc(sz * sz * sizeof(double));
        float * As = flint_malloc(sz * sz * sizeof(float));
        float * Bs = flint_malloc(sz * sz * sizeof(float));
        float * Cs = flint_malloc(sz * sz * sizeof(float));
        double td, ts, flops = 2.0 * (double) sz * (double) sz * (double) sz;
        slong i;

        for (i = 0; i < sz * sz; i++)
        {
            A[i] = (double) (bench_rand(i) % 1000) / 1000.0;
            B[i] = (double) (bench_rand(i + 7) % 1000) / 1000.0;
            As[i] = (float) A[i];
            Bs[i] = (float) B[i];
        }

        TIME(bench_dgemm_fallback(sz, sz, sz, A, sz, B, sz, C, sz), td);
        TIME(bench_sgemm_fallback(sz, sz, sz, As, sz, Bs, sz, Cs, sz), ts);

        flint_printf("  %-8wd %12.4f %12.2f %12.4f %12.2f\n", sz,
                1e3 * td, 1e-9 * flops / td, 1e3 * ts, 1e-9 * flops / ts);
        fflush(stdout);

        flint_free(A); flint_free(B); flint_free(C);
        flint_free(As); flint_free(Bs); flint_free(Cs);
    }
}

int
main(void)
{
    flint_set_num_threads(1);

#if defined(FLINT_MACHINE_VECTORS_AVX2)
    flint_printf("backend: AVX2");
#elif defined(FLINT_MACHINE_VECTORS_NEON)
    flint_printf("backend: NEON");
#elif defined(FLINT_MACHINE_VECTORS_GNU_VECTOR_EXTENSIONS)
    flint_printf("backend: generic, GNU vector extensions");
#else
    flint_printf("backend: generic, strict ISO C");
#endif
    flint_printf("   (sizeof(vec4d) = %d, one thread)\n", (int) sizeof(vec4d));

    flint_printf("\nmodular primitives, 50 bit modulus:\n");
    bench_primitives();

    flint_printf("\nflint_dgemm_fallback / flint_sgemm_fallback, square:\n");
    bench_gemm();

    flint_cleanup_master();
    return 0;
}
