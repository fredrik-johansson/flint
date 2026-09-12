/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "test_helpers.h"
#include "machine_vectors.h"
#include "ulong_extras.h"

/*
    Tests of the non-modular machine_vectors operations against per lane
    references, including the exact semantics that fft_small and nmod_vec
    rely on: bit patterns of comparison masks, the sign bit convention of
    blendv, round to nearest with ties to even, the 128-bit lane pair
    convention of the vec4d unpacks, the low 32 x low 32 semantics of
    vecKn_mul, wraparound of integer add/sub/horizontal_sum, and the full
    and limited input ranges of addmod. The operations exercised are
    gated on the backend since the backends implement slightly different
    supersets of the common interface; see doc/source/machine_vectors.rst.

    The generic backend implements a superset of the AVX2 and NEON
    interfaces, so on FLINT_MACHINE_VECTORS_GENERIC everything is tested.
*/

#define CHECK(cond, name, extra) \
    do { \
        if (!(cond)) \
        { \
            flint_printf("FAIL: %s (%s)\n", name, extra); \
            fflush(stdout); \
            flint_abort(); \
        } \
    } while (0)

/* random double: integer valued, |a| <= b, b < 2^62 */
static double ops_rand_int(flint_rand_t state, ulong b)
{
    ulong u = n_randint(state, 2*b + 1);
    return (double) ((slong) u - (slong) b);
}

/* random double in [-b, b], with quarter unit granularity */
static double ops_rand_q(flint_rand_t state, ulong b)
{
    return 0.25 * ops_rand_int(state, 4*b);
}

static int ops_bit_same(double a, double b)
{
    return memcmp(&a, &b, sizeof(double)) == 0;
}

/* the mask entries produced by cmp are all-zeros or all-ones */
static int ops_is_mask(double a, int expected)
{
    ulong u;
    memcpy(&u, &a, sizeof(ulong));
    return u == (expected ? ~UWORD(0) : UWORD(0));
}

static double ops_ref_blendv(double a, double b, double c)
{
    ulong u;
    memcpy(&u, &c, sizeof(ulong));
    return (u >> 63) ? b : a;
}

TEST_FUNCTION_START(machine_vectors_ops, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        double a[8], b[8], c[8], z[8], w[8];
        slong i;

        for (i = 0; i < 8; i++)
        {
            a[i] = ops_rand_q(state, 1000);
            b[i] = ops_rand_q(state, 1000);
            c[i] = ops_rand_q(state, 1000);
        }

        /* make some interesting special values likely */
        if (n_randint(state, 4) == 0)
        {
            i = n_randint(state, 8);
            a[i] = n_randint(state, 2) ? -0.0 : 0.0;
            i = n_randint(state, 8);
            b[i] = a[n_randint(state, 8)];
            i = n_randint(state, 8);
            c[i] = n_randint(state, 2) ? -0.0 : 0.0;
        }

        /* vec1d ------------------------------------------------------- */
        {
            double x = a[0], y = b[0], t = c[0];

            CHECK(vec1d_same(vec1d_load(a), a[0]), "vec1d_load", "");
            vec1d_store(z, vec1d_add(x, y));
            CHECK(z[0] == x + y, "vec1d_add", "");
            CHECK(vec1d_same(vec1d_sub(x, y), x - y), "vec1d_sub", "");
            CHECK(vec1d_same(vec1d_mul(x, y), x * y), "vec1d_mul", "");
            CHECK(vec1d_same(vec1d_half(x), 0.5 * x), "vec1d_half", "");
            CHECK(vec1d_same(vec1d_neg(x), -x) && ops_bit_same(vec1d_neg(-0.0), 0.0),
                    "vec1d_neg", "");
            CHECK(vec1d_same(vec1d_abs(x), fabs(x)), "vec1d_abs", "");
            CHECK(y == 0.0 || vec1d_same(vec1d_div(x, y), x / y), "vec1d_div", "");
            CHECK(vec1d_same(vec1d_min(x, y), fmin(x, y)), "vec1d_min", "");
            CHECK(vec1d_same(vec1d_max(x, y), fmax(x, y)), "vec1d_max", "");
            CHECK(vec1d_same(vec1d_zero(), 0.0), "vec1d_zero", "");
            CHECK(vec1d_same(vec1d_one(), 1.0), "vec1d_one", "");
            CHECK(vec1d_same(vec1d_fmadd(x, y, t), fma(x, y, t)), "vec1d_fmadd", "");
            CHECK(vec1d_same(vec1d_fmsub(x, y, t), fma(x, y, -t)), "vec1d_fmsub", "");
            CHECK(vec1d_same(vec1d_fnmadd(x, y, t), fma(-x, y, t)), "vec1d_fnmadd", "");
            CHECK(vec1d_same(vec1d_fnmsub(x, y, t), fma(-x, y, -t)), "vec1d_fnmsub", "");
            CHECK(vec1d_same(vec1d_floor(x), floor(x)), "vec1d_floor", "");
            CHECK(vec1d_same(vec1d_round(x), rint(x)), "vec1d_round", "");
            /* ties to even and sign preservation */
            CHECK(vec1d_round(0.5) == 0.0 && vec1d_round(1.5) == 2.0
                    && vec1d_round(2.5) == 2.0 && vec1d_round(-0.5) == 0.0
                    && vec1d_round(-1.5) == -2.0
                    && vec1d_round(0x1.0p52 + 0.0) == 0x1.0p52
                    && vec1d_round(0x1.23p60) == 0x1.23p60,
                    "vec1d_round", "ties");
            /* vec1d_blendv is documented to select on c >= 0, unlike the
               sign bit convention of the vector widths */
            CHECK(vec1d_same(vec1d_blendv(x, y, t), t >= 0 ? x : y),
                    "vec1d_blendv", "");
            /* addsub on vec1d is the n == 1 column of vecKd_addsub: sub */
            CHECK(vec1d_same(vec1d_addsub(x, y), x - y), "vec1d_addsub", "");
        }

        /* vec4d arithmetic -------------------------------------------- */
        {
            vec4d x = vec4d_load_unaligned(a);
            vec4d y = vec4d_load_unaligned(b);
            vec4d t = vec4d_load_unaligned(c);

            vec4d_store_unaligned(z, vec4d_add(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i] + b[i], "vec4d_add", "");

            vec4d_store_unaligned(z, vec4d_sub(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i] - b[i], "vec4d_sub", "");

            vec4d_store_unaligned(z, vec4d_mul(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i] * b[i], "vec4d_mul", "");

            vec4d_store_unaligned(z, vec4d_min(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == fmin(a[i], b[i]), "vec4d_min", "");

            vec4d_store_unaligned(z, vec4d_max(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == fmax(a[i], b[i]), "vec4d_max", "");

            vec4d_store_unaligned(z, vec4d_div(x, y));
            for (i = 0; i < 4; i++)
                CHECK(b[i] == 0.0 || z[i] == a[i] / b[i], "vec4d_div", "");

            vec4d_store_unaligned(z, vec4d_neg(x));
            for (i = 0; i < 4; i++)
                CHECK(ops_bit_same(z[i], -a[i]), "vec4d_neg", "");

            /* fused on AVX2 and NEON, unfused on the generic backend;
               with exactly representable products they must agree */
#if defined(FLINT_MACHINE_VECTORS_GENERIC)
# define FMA_REF(u, v, w) ((u) * (v) + (w))
#else
# define FMA_REF(u, v, w) fma(u, v, w)
#endif
            vec4d_store_unaligned(z, vec4d_fmadd(x, y, t));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == FMA_REF(a[i], b[i], c[i]), "vec4d_fmadd", "");

            vec4d_store_unaligned(z, vec4d_fmsub(x, y, t));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == FMA_REF(a[i], b[i], -c[i]), "vec4d_fmsub", "");

            vec4d_store_unaligned(z, vec4d_fnmadd(x, y, t));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == FMA_REF(-a[i], b[i], c[i]), "vec4d_fnmadd", "");

            vec4d_store_unaligned(z, vec4d_fnmsub(x, y, t));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == FMA_REF(-a[i], b[i], -c[i]), "vec4d_fnmsub", "");
#undef FMA_REF

            vec4d_store_unaligned(z, vec4d_round(x));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == rint(a[i]), "vec4d_round", "");

            vec4d_store_unaligned(z, vec4d_round(
                    vec4d_set_d4(0.5, 2.5, -1.5, 0x1.0p52)));
            CHECK(z[0] == 0.0 && z[1] == 2.0 && z[2] == -2.0 && z[3] == 0x1.0p52,
                    "vec4d_round", "ties");

            vec4d_store_unaligned(z, vec4d_blendv(x, y, t));
            for (i = 0; i < 4; i++)
                CHECK(ops_bit_same(z[i], ops_ref_blendv(a[i], b[i], c[i])),
                        "vec4d_blendv", "");

            vec4d_store_unaligned(z, vec4d_set_d(a[0]));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[0], "vec4d_set_d", "");

            vec4d_store_unaligned(z, vec4d_set_d4(a[0], a[1], a[2], a[3]));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i], "vec4d_set_d4", "");

            vec4d_store_unaligned(z, vec4d_zero());
            vec4d_store_unaligned(w, vec4d_one());
            for (i = 0; i < 4; i++)
                CHECK(z[i] == 0.0 && w[i] == 1.0, "vec4d_zero/one", "");

            CHECK(vec4d_get_index(x, 0) == a[0] && vec4d_get_index(x, 1) == a[1]
                    && vec4d_get_index(x, 2) == a[2] && vec4d_get_index(x, 3) == a[3],
                    "vec4d_get_index", "");

#if !defined(FLINT_MACHINE_VECTORS_NEON)
            /* not provided by the NEON backend */
            vec4d_store_unaligned(z, vec4d_abs(x));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == fabs(a[i]), "vec4d_abs", "");

            vec4d_store_unaligned(z, vec4d_half(x));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == 0.5 * a[i], "vec4d_half", "");

            vec4d_store_unaligned(z, vec4d_addsub(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == ((i % 2) ? a[i] + b[i] : a[i] - b[i]),
                        "vec4d_addsub", "");

            vec4d_store_unaligned(z, vec4d_floor(x));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == floor(a[i]), "vec4d_floor", "");

            vec4d_store_unaligned(z, vec4d_cmp_ge(x, y));
            for (i = 0; i < 4; i++)
                CHECK(ops_is_mask(z[i], a[i] >= b[i]), "vec4d_cmp_ge", "");

            vec4d_store_unaligned(z, vec4d_cmp_gt(x, y));
            for (i = 0; i < 4; i++)
                CHECK(ops_is_mask(z[i], a[i] > b[i]), "vec4d_cmp_gt", "");

            CHECK(vec4d_same(x, x) && !vec4d_same(x, vec4d_add(x, vec4d_one())),
                    "vec4d_same", "");
#endif
#if defined(FLINT_MACHINE_VECTORS_GENERIC)
            vec4d_store_unaligned(z, vec4d_cmp_lt(x, y));
            for (i = 0; i < 4; i++)
                CHECK(ops_is_mask(z[i], a[i] < b[i]), "vec4d_cmp_lt", "");
#endif
        }

        /* vec4d permutations ------------------------------------------ */
        {
            vec4d x = vec4d_load_unaligned(a);
            vec4d y = vec4d_load_unaligned(b);

            vec4d_store_unaligned(z, vec4d_unpacklo(x, y));
            CHECK(z[0] == a[0] && z[1] == b[0] && z[2] == a[2] && z[3] == b[2],
                    "vec4d_unpacklo", "");

            vec4d_store_unaligned(z, vec4d_unpackhi(x, y));
            CHECK(z[0] == a[1] && z[1] == b[1] && z[2] == a[3] && z[3] == b[3],
                    "vec4d_unpackhi", "");

            vec4d_store_unaligned(z, vec4d_permute_0_2_1_3(x));
            CHECK(z[0] == a[0] && z[1] == a[2] && z[2] == a[1] && z[3] == a[3],
                    "vec4d_permute_0_2_1_3", "");

            vec4d_store_unaligned(z, vec4d_permute_3_1_2_0(x));
            CHECK(z[0] == a[3] && z[1] == a[1] && z[2] == a[2] && z[3] == a[0],
                    "vec4d_permute_3_1_2_0", "");

            vec4d_store_unaligned(z, vec4d_permute_3_2_1_0(x));
            CHECK(z[0] == a[3] && z[1] == a[2] && z[2] == a[1] && z[3] == a[0],
                    "vec4d_permute_3_2_1_0", "");

            vec4d_store_unaligned(z, vec4d_permute2_0_2(x, y));
            CHECK(z[0] == a[0] && z[1] == a[1] && z[2] == b[0] && z[3] == b[1],
                    "vec4d_permute2_0_2", "");

            vec4d_store_unaligned(z, vec4d_permute2_1_3(x, y));
            CHECK(z[0] == a[2] && z[1] == a[3] && z[2] == b[2] && z[3] == b[3],
                    "vec4d_permute2_1_3", "");

            vec4d_store_unaligned(z, vec4d_unpack_lo_permute_0_2_1_3(x, y));
            CHECK(z[0] == a[0] && z[1] == a[2] && z[2] == b[0] && z[3] == b[2],
                    "vec4d_unpack_lo_permute_0_2_1_3", "");

            vec4d_store_unaligned(z, vec4d_unpack_hi_permute_0_2_1_3(x, y));
            CHECK(z[0] == a[1] && z[1] == a[3] && z[2] == b[1] && z[3] == b[3],
                    "vec4d_unpack_hi_permute_0_2_1_3", "");

            vec4d_store_unaligned(z, vec4d_unpacklo_permute_3_1_2_0(x, y));
            CHECK(z[0] == b[2] && z[1] == b[0] && z[2] == a[2] && z[3] == a[0],
                    "vec4d_unpacklo_permute_3_1_2_0", "");

            vec4d_store_unaligned(z, vec4d_unpackhi_permute_3_1_2_0(x, y));
            CHECK(z[0] == b[3] && z[1] == b[1] && z[2] == a[3] && z[3] == a[1],
                    "vec4d_unpackhi_permute_3_1_2_0", "");

            {
                vec4d r0, r1, r2, r3;
                vec4d s0 = vec4d_load_unaligned(a);
                vec4d s1 = vec4d_load_unaligned(a + 4);
                vec4d s2 = vec4d_load_unaligned(b);
                vec4d s3 = vec4d_load_unaligned(b + 4);
                double m[16];
                VEC4D_TRANSPOSE(r0, r1, r2, r3, s0, s1, s2, s3);
                vec4d_store_unaligned(m + 0, r0);
                vec4d_store_unaligned(m + 4, r1);
                vec4d_store_unaligned(m + 8, r2);
                vec4d_store_unaligned(m + 12, r3);
                for (i = 0; i < 4; i++)
                    CHECK(m[4*i + 0] == a[i] && m[4*i + 1] == a[4 + i]
                            && m[4*i + 2] == b[i] && m[4*i + 3] == b[4 + i],
                            "VEC4D_TRANSPOSE", "");
            }

#if !defined(FLINT_MACHINE_VECTORS_AVX2)
            {
                vec2d h1 = vec2d_load_unaligned(a);
                vec2d h2 = vec2d_load_unaligned(b);
                vec4d_store_unaligned(z, vec4d_set_vec2d2(h1, h2));
                CHECK(z[0] == a[0] && z[1] == a[1] && z[2] == b[0] && z[3] == b[1],
                        "vec4d_set_vec2d2", "");
            }
#endif
        }

        /* vec8d ------------------------------------------------------- */
        {
            vec8d x = vec8d_load_unaligned(a);
            vec8d y = vec8d_load_unaligned(b);
            vec8d t = vec8d_load_unaligned(c);

            vec8d_store_unaligned(z, vec8d_add(x, y));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == a[i] + b[i], "vec8d_add", "");

            vec8d_store_unaligned(z, vec8d_sub(x, y));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == a[i] - b[i], "vec8d_sub", "");

            vec8d_store_unaligned(z, vec8d_mul(x, y));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == a[i] * b[i], "vec8d_mul", "");

            vec8d_store_unaligned(z, vec8d_neg(x));
            for (i = 0; i < 8; i++)
                CHECK(ops_bit_same(z[i], -a[i]), "vec8d_neg", "");

            vec8d_store_unaligned(z, vec8d_round(x));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == rint(a[i]), "vec8d_round", "");

            vec8d_store_unaligned(z, vec8d_blendv(x, y, t));
            for (i = 0; i < 8; i++)
                CHECK(ops_bit_same(z[i], ops_ref_blendv(a[i], b[i], c[i])),
                        "vec8d_blendv", "");

            vec8d_store_unaligned(z,
                vec8d_set_d8(b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == b[i], "vec8d_set_d8", "");

            vec8d_store_unaligned(z, vec8d_set_d(a[0]));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == a[0], "vec8d_set_d", "");

            for (i = 0; i < 8; i++)
                w[i] = vec8d_get_index(x, i);
            CHECK(memcmp(w, a, sizeof(w)) == 0, "vec8d_get_index", "");

#if !defined(FLINT_MACHINE_VECTORS_AVX2)
            vec8d_store_unaligned(z, vec8d_one());
            for (i = 0; i < 8; i++)
                CHECK(z[i] == 1.0, "vec8d_one", "");

            vec8d_store_unaligned(z, vec8d_unpacklo(x, y));
            for (i = 0; i < 8; i += 2)
                CHECK(z[i] == a[i - i % 2] && z[i + 1] == b[i - i % 2],
                        "vec8d_unpacklo", "");
#endif
#if defined(FLINT_MACHINE_VECTORS_GENERIC)
            vec8d_store_unaligned(z, vec8d_abs(x));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == fabs(a[i]), "vec8d_abs", "");

            vec8d_store_unaligned(z, vec8d_half(x));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == 0.5 * a[i], "vec8d_half", "");

            vec8d_store_unaligned(z, vec8d_floor(x));
            for (i = 0; i < 8; i++)
                CHECK(z[i] == floor(a[i]), "vec8d_floor", "");
#endif
        }
    }

#if FLINT_BITS == 64
    /* integer vectors */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        ulong a[8], b[8], z[8], n[8];
        slong i;

        for (i = 0; i < 8; i++)
        {
            a[i] = n_randtest(state);
            b[i] = n_randtest(state);
        }

        /* vec4n ------------------------------------------------------- */
        {
            vec4n x = vec4n_load_unaligned(a);
            vec4n y = vec4n_load_unaligned(b);

            vec4n_store_unaligned(z, x);
            CHECK(memcmp(z, a, 4*sizeof(ulong)) == 0, "vec4n_load/store", "");

            vec4n_store_unaligned(z, vec4n_add(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i] + b[i], "vec4n_add", "");

            vec4n_store_unaligned(z, vec4n_sub(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i] - b[i], "vec4n_sub", "");

            vec4n_store_unaligned(z, vec4n_bit_and(x, y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == (a[i] & b[i]), "vec4n_bit_and", "");

            vec4n_store_unaligned(z, vec4n_bit_shift_right_32(x));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i] >> 32, "vec4n_bit_shift_right_32", "");

            vec4n_store_unaligned(z, vec4n_set_n(a[0]));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[0], "vec4n_set_n", "");

#if !defined(FLINT_MACHINE_VECTORS_NEON)
            {
                ulong sh = n_randint(state, 64);
                ulong s = 0;

                vec4n_store_unaligned(z, vec4n_bit_shift_right(x, sh));
                for (i = 0; i < 4; i++)
                    CHECK(z[i] == a[i] >> sh, "vec4n_bit_shift_right", "");

                /* low 32 x low 32 -> full 64 product, high halves ignored */
                vec4n_store_unaligned(z, vec4n_mul(x, y));
                for (i = 0; i < 4; i++)
                    CHECK(z[i] == (a[i] & UWORD(0xffffffff)) * (b[i] & UWORD(0xffffffff)),
                            "vec4n_mul", "");

                for (i = 0; i < 4; i++)
                    s += a[i];
                CHECK(vec4n_horizontal_sum(x) == s, "vec4n_horizontal_sum", "");

                vec4n_store_unaligned(z, vec4n_zero());
                for (i = 0; i < 4; i++)
                    CHECK(z[i] == 0, "vec4n_zero", "");

                vec4n_store_unaligned(z, vec4n_set_n4(a[0], a[1], a[2], a[3]));
                CHECK(memcmp(z, a, 4*sizeof(ulong)) == 0, "vec4n_set_n4", "");

                vec4n_store_unaligned(z, vec4n_permute_3_2_1_0(x));
                for (i = 0; i < 4; i++)
                    CHECK(z[i] == a[3 - i], "vec4n_permute_3_2_1_0", "");
            }
#endif
        }

        /* addmod: works for any n != 0 with a, b in [0, n), including
           n > 2^63; addmod_limited requires n < 2^63 */
        {
            vec4n x, y, m;

            for (i = 0; i < 8; i++)
            {
                n[i] = n_randtest_not_zero(state);
                a[i] = n_randlimb(state) % n[i];
                b[i] = n_randlimb(state) % n[i];
            }
            /* the same modulus is broadcast in practice; use lanewise
               moduli here for stronger coverage */
            x = vec4n_load_unaligned(a);
            y = vec4n_load_unaligned(b);
            m = vec4n_load_unaligned(n);

            vec4n_store_unaligned(z, vec4n_addmod(x, y, m));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == n_addmod(a[i], b[i], n[i]), "vec4n_addmod",
                        "full range");

            for (i = 0; i < 8; i++)
            {
                n[i] >>= 1;
                n[i] += (n[i] == 0);
                a[i] %= n[i];
                b[i] %= n[i];
            }
            x = vec4n_load_unaligned(a);
            y = vec4n_load_unaligned(b);
            m = vec4n_load_unaligned(n);

            vec4n_store_unaligned(z, vec4n_addmod_limited(x, y, m));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == n_addmod(a[i], b[i], n[i]),
                        "vec4n_addmod_limited", "n < 2^63");

            /* vec8n */
            {
                vec8n x8, y8, m8, t8;

                /* the vec8n type may not support lanewise moduli, so
                   broadcast n[0] and bring the operands into [0, n[0]) */
                for (i = 0; i < 8; i++)
                {
                    n[i] = n[0];
                    a[i] %= n[0];
                    b[i] %= n[0];
                }
                x8 = vec8n_load_unaligned(a);
                y8 = vec8n_load_unaligned(b);
                m8 = vec8n_set_n(n[0]);

                t8 = vec8n_addmod_limited(x8, y8, m8);
                vec8d_store_unaligned((double *) z, vec8n_convert_limited_vec8d(
                            vec8n_bit_and(t8, vec8n_set_n(UWORD(0xffffffff)))));
                for (i = 0; i < 8; i++)
                    CHECK(((double *) z)[i] ==
                            (double) (n_addmod(a[i], b[i], n[0]) & UWORD(0xffffffff)),
                            "vec8n_addmod_limited", "");

                t8 = vec8n_addmod(x8, y8, m8);
                vec8d_store_unaligned((double *) z, vec8n_convert_limited_vec8d(
                            vec8n_bit_shift_right_32(t8)));
                for (i = 0; i < 8; i++)
                    CHECK(((double *) z)[i] ==
                            (double) (n_addmod(a[i], b[i], n[0]) >> 32),
                            "vec8n_addmod", "");
            }

#if !defined(FLINT_MACHINE_VECTORS_AVX2)
            {
                vec2n x2 = vec2n_load_unaligned(a);
                vec2n y2 = vec2n_load_unaligned(b);
                vec2n m2 = vec2n_load_unaligned(n);

                vec2n_store_unaligned(z, vec2n_addmod(x2, y2, m2));
                for (i = 0; i < 2; i++)
                    CHECK(z[i] == n_addmod(a[i], b[i], n[i]), "vec2n_addmod", "");

                vec2n_store_unaligned(z, vec2n_addmod_limited(x2, y2, m2));
                for (i = 0; i < 2; i++)
                    CHECK(z[i] == n_addmod(a[i], b[i], n[i]),
                            "vec2n_addmod_limited", "");

                vec2n_store_unaligned(z, vec2n_add(x2, vec2n_sub(y2,
                                vec2n_bit_and(x2, vec2n_set_n(0)))));
                for (i = 0; i < 2; i++)
                    CHECK(z[i] == a[i] + b[i], "vec2n_add/sub/bit_and/set_n", "");

                CHECK(vec1n_addmod(a[0], b[0], n[0]) == n_addmod(a[0], b[0], n[0]),
                        "vec1n_addmod", "");
            }
#endif
        }

        /* conversions: exact integers in [0, 2^52) -------------------- */
        {
            double d[8], e[8];
            vec4n x;
            vec4d y;

            for (i = 0; i < 8; i++)
            {
                a[i] = n_randint(state, UWORD(1) << 52);
                if (n_randint(state, 8) == 0)
                    a[i] = n_randint(state, 2) ? 0 : (UWORD(1) << 52) - 1;
                d[i] = (double) a[i];
            }

            x = vec4n_load_unaligned(a);
            vec4d_store_unaligned(e, vec4n_convert_limited_vec4d(x));
            for (i = 0; i < 4; i++)
                CHECK(e[i] == d[i], "vec4n_convert_limited_vec4d", "");

            y = vec4d_load_unaligned(d);
            vec4n_store_unaligned(z, vec4d_convert_limited_vec4n(y));
            for (i = 0; i < 4; i++)
                CHECK(z[i] == a[i], "vec4d_convert_limited_vec4n", "");

            vec8d_store_unaligned(e, vec8n_convert_limited_vec8d(
                        vec8n_load_unaligned(a)));
            for (i = 0; i < 8; i++)
                CHECK(e[i] == d[i], "vec8n_convert_limited_vec8d", "");

#if defined(FLINT_MACHINE_VECTORS_GENERIC)
            {
                vec8n z8 = vec8d_convert_limited_vec8n(vec8d_load_unaligned(d));
                vec8n_store_unaligned(z, z8);
                for (i = 0; i < 8; i++)
                    CHECK(z[i] == a[i], "vec8d_convert_limited_vec8n", "");
            }
            CHECK(vec1n_convert_limited_vec1d(a[0]) == d[0],
                    "vec1n_convert_limited_vec1d", "");
#endif
#if !defined(FLINT_MACHINE_VECTORS_AVX2)
            CHECK(vec1d_convert_limited_vec1n(d[0]) == a[0],
                    "vec1d_convert_limited_vec1n", "");
            {
                vec2n x2 = vec2n_load_unaligned(a);
                vec2d y2;
                vec2d_store_unaligned(e, vec2n_convert_limited_vec2d(x2));
                for (i = 0; i < 2; i++)
                    CHECK(e[i] == d[i], "vec2n_convert_limited_vec2d", "");
                y2 = vec2d_load_unaligned(d);
                vec2n_store_unaligned(z, vec2d_convert_limited_vec2n(y2));
                for (i = 0; i < 2; i++)
                    CHECK(z[i] == a[i], "vec2d_convert_limited_vec2n", "");
            }
#endif
        }
    }

#endif /* FLINT_BITS == 64 */

    /* vec2d, on the backends that provide it -------------------------- */
#if !defined(FLINT_MACHINE_VECTORS_AVX2)
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        double a[2], b[2], c[2], z[2];
        slong i;

        for (i = 0; i < 2; i++)
        {
            a[i] = ops_rand_q(state, 1000);
            b[i] = ops_rand_q(state, 1000);
            c[i] = ops_rand_q(state, 1000);
        }

        {
            vec2d x = vec2d_load_unaligned(a);
            vec2d y = vec2d_load_unaligned(b);
            vec2d t = vec2d_load_unaligned(c);

            vec2d_store_unaligned(z, vec2d_add(x, y));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == a[i] + b[i], "vec2d_add", "");

            vec2d_store_unaligned(z, vec2d_sub(x, y));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == a[i] - b[i], "vec2d_sub", "");

            vec2d_store_unaligned(z, vec2d_mul(x, y));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == a[i] * b[i], "vec2d_mul", "");

            vec2d_store_unaligned(z, vec2d_min(x, y));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == fmin(a[i], b[i]), "vec2d_min", "");

            vec2d_store_unaligned(z, vec2d_max(x, y));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == fmax(a[i], b[i]), "vec2d_max", "");

            vec2d_store_unaligned(z, vec2d_div(x, y));
            for (i = 0; i < 2; i++)
                CHECK(b[i] == 0.0 || z[i] == a[i] / b[i], "vec2d_div", "");

            vec2d_store_unaligned(z, vec2d_neg(x));
            for (i = 0; i < 2; i++)
                CHECK(ops_bit_same(z[i], -a[i]), "vec2d_neg", "");

            vec2d_store_unaligned(z, vec2d_abs(x));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == fabs(a[i]), "vec2d_abs", "");

            vec2d_store_unaligned(z, vec2d_half(x));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == 0.5 * a[i], "vec2d_half", "");

            vec2d_store_unaligned(z, vec2d_round(x));
            for (i = 0; i < 2; i++)
                CHECK(z[i] == rint(a[i]), "vec2d_round", "");

            vec2d_store_unaligned(z, vec2d_blendv(x, y, t));
            for (i = 0; i < 2; i++)
                CHECK(ops_bit_same(z[i], ops_ref_blendv(a[i], b[i], c[i])),
                        "vec2d_blendv", "");

            vec2d_store_unaligned(z, vec2d_unpacklo(x, y));
            CHECK(z[0] == a[0] && z[1] == b[0], "vec2d_unpacklo", "");

            vec2d_store_unaligned(z, vec2d_unpackhi(x, y));
            CHECK(z[0] == a[1] && z[1] == b[1], "vec2d_unpackhi", "");

            vec2d_store_unaligned(z, vec2d_set_d(a[0]));
            CHECK(z[0] == a[0] && z[1] == a[0], "vec2d_set_d", "");

            CHECK(vec2d_get_index(x, 0) == a[0] && vec2d_get_index(x, 1) == a[1],
                    "vec2d_get_index", "");

            vec2d_store_unaligned(z, vec2d_zero());
            CHECK(z[0] == 0.0 && z[1] == 0.0, "vec2d_zero", "");

            vec2d_store_unaligned(z, vec2d_one());
            CHECK(z[0] == 1.0 && z[1] == 1.0, "vec2d_one", "");

            vec2d_store_unaligned(z, vec2d_cmp_gt(x, y));
            for (i = 0; i < 2; i++)
                CHECK(ops_is_mask(z[i], a[i] > b[i]), "vec2d_cmp_gt", "");

            vec2d_store_unaligned(z, vec2d_cmp_lt(x, y));
            for (i = 0; i < 2; i++)
                CHECK(ops_is_mask(z[i], a[i] < b[i]), "vec2d_cmp_lt", "");

#if defined(FLINT_MACHINE_VECTORS_GENERIC)
            vec2d_store_unaligned(z, vec2d_cmp_ge(x, y));
            for (i = 0; i < 2; i++)
                CHECK(ops_is_mask(z[i], a[i] >= b[i]), "vec2d_cmp_ge", "");

            CHECK(vec2d_same(x, x) && !vec2d_same(x, vec2d_add(x, vec2d_one())),
                    "vec2d_same", "");
#endif
        }
    }
#endif

    TEST_FUNCTION_END(state);
}

#undef CHECK
