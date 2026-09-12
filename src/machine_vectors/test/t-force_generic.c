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
    The generic backends of machine_vectors.h can be selected on any
    target with FLINT_MACHINE_VECTORS_FORCE_GENERIC (GNU vector
    extensions tier) or FLINT_MACHINE_VECTORS_STRICT_C (plain C tier).
    This test exercises both tiers even when the library itself was
    built with the AVX2 or NEON backend, by instantiating
    machine_vectors_generic.h two more times in this translation unit
    under renamed identifiers: once as mvg_* (GNU tier) and once as
    mvs_* (strict tier). The renaming is done by the generated
    object like macros in force_generic_renames.h; to regenerate that
    file after changing the generic backend, run from the top level
    directory (the sed pipeline extracts every identifier the header
    defines from preprocessed output, since most names are pasted
    together by macros and never appear literally):

    for tier in FORCE_GENERIC STRICT_C; do
        echo '#include "machine_vectors.h"' |
        gcc -E -std=c11 -DFLINT_MACHINE_VECTORS_$tier -I src -I . -xc -
    done | ... (see the header for the exact recipe)

    The checks:

    1. The two generic tiers must agree bit for bit on every operation,
       for all inputs in the documented domains, including signed zeros
       and range boundaries.

    2. The native backend must agree bit for bit with the generic
       tiers on the operations for which this is guaranteed:

       - add, sub, mul (IEEE operations),
       - round (|a| < 2^52; round to nearest, ties to even),
       - the conversions (exact on [0, 2^52)),
       - mulmod, nmulmod (identical quotient q = round((a b) ninv)
         computed without fused operations, and a remainder a b - q n
         that is exact on all backends),
       - reduce_to_pm1n, reduce_to_pm1no, reduce_to_0n (|a| <= 8n,
         where a - round(a ninv) n is exact),
       - the integer operations and addmod.

       fmadd and friends are intentionally not in this list: they are
       fused on AVX2 and NEON but unfused on the generic backend.
*/

/*
    The native VEC4D_TRANSPOSE and vecKn_bit_shift_right_32 macros are
    replaced by the generic definitions below (a #undef of a name that
    is an inline function on some backend is a no-op). This is why this
    file must be included last in main.c.
*/
#undef VEC4D_TRANSPOSE
#undef vec2n_bit_shift_right_32
#undef vec4n_bit_shift_right_32
#undef vec8n_bit_shift_right_32

/* instantiate the GNU vector extensions tier as mvg_* (if the build
   forces the strict tier globally, this instantiates the strict tier
   again, which keeps the comparisons below valid, just redundant) */
#define FLINT_MV_RENAME_MVG
#include "force_generic_renames.h"
#undef MACHINE_VECTORS_GENERIC_H
#include "machine_vectors_generic.h"
#define FLINT_MV_RENAME_UNDEF
#include "force_generic_renames.h"

/* instantiate the strict C tier as mvs_* */
#ifndef FLINT_MACHINE_VECTORS_STRICT_C
# define FLINT_MACHINE_VECTORS_STRICT_C 1
# define _MV_TEST_UNDO_STRICT
#endif
#undef VEC4D_TRANSPOSE
#undef vec2n_bit_shift_right_32
#undef vec4n_bit_shift_right_32
#undef vec8n_bit_shift_right_32
#define FLINT_MV_RENAME_MVS
#include "force_generic_renames.h"
#undef MACHINE_VECTORS_GENERIC_H
#include "machine_vectors_generic.h"
#define FLINT_MV_RENAME_UNDEF
#include "force_generic_renames.h"
#ifdef _MV_TEST_UNDO_STRICT
# undef FLINT_MACHINE_VECTORS_STRICT_C
# undef _MV_TEST_UNDO_STRICT
#endif

#define CHECK(cond, name, extra) \
    do { \
        if (!(cond)) \
        { \
            flint_printf("FAIL: %s (%s)\n", name, extra); \
            fflush(stdout); \
            flint_abort(); \
        } \
    } while (0)

/* random integer valued double in [-b, b] */
static double fg_rand_int(flint_rand_t state, ulong b)
{
    ulong u = n_randint(state, 2*b + 1);
    return (double) ((slong) u - (slong) b);
}

/*
    mvg vs mvs, vec4d, by arity; operands in av, bv, cv, moduli in
    nv/iv; results compared bit for bit
*/
#define TIERS1(op) \
    do { \
        double zg[4], zs[4]; \
        mvg_vec4d_store_unaligned(zg, \
                mvg_vec4d_##op(mvg_vec4d_load_unaligned(av))); \
        mvs_vec4d_store_unaligned(zs, \
                mvs_vec4d_##op(mvs_vec4d_load_unaligned(av))); \
        CHECK(memcmp(zg, zs, sizeof(zg)) == 0, "vec4d_" #op, "mvg vs mvs"); \
    } while (0)

#define TIERS2(op, x, y) \
    do { \
        double zg[4], zs[4]; \
        mvg_vec4d_store_unaligned(zg, \
                mvg_vec4d_##op(mvg_vec4d_load_unaligned(x), \
                               mvg_vec4d_load_unaligned(y))); \
        mvs_vec4d_store_unaligned(zs, \
                mvs_vec4d_##op(mvs_vec4d_load_unaligned(x), \
                               mvs_vec4d_load_unaligned(y))); \
        CHECK(memcmp(zg, zs, sizeof(zg)) == 0, "vec4d_" #op, "mvg vs mvs"); \
    } while (0)

#define TIERS3(op, x, y, w) \
    do { \
        double zg[4], zs[4]; \
        mvg_vec4d_store_unaligned(zg, \
                mvg_vec4d_##op(mvg_vec4d_load_unaligned(x), \
                               mvg_vec4d_load_unaligned(y), \
                               mvg_vec4d_load_unaligned(w))); \
        mvs_vec4d_store_unaligned(zs, \
                mvs_vec4d_##op(mvs_vec4d_load_unaligned(x), \
                               mvs_vec4d_load_unaligned(y), \
                               mvs_vec4d_load_unaligned(w))); \
        CHECK(memcmp(zg, zs, sizeof(zg)) == 0, "vec4d_" #op, "mvg vs mvs"); \
    } while (0)

#define TIERS4(op, x, y, w, v) \
    do { \
        double zg[4], zs[4]; \
        mvg_vec4d_store_unaligned(zg, \
                mvg_vec4d_##op(mvg_vec4d_load_unaligned(x), \
                               mvg_vec4d_load_unaligned(y), \
                               mvg_vec4d_load_unaligned(w), \
                               mvg_vec4d_load_unaligned(v))); \
        mvs_vec4d_store_unaligned(zs, \
                mvs_vec4d_##op(mvs_vec4d_load_unaligned(x), \
                               mvs_vec4d_load_unaligned(y), \
                               mvs_vec4d_load_unaligned(w), \
                               mvs_vec4d_load_unaligned(v))); \
        CHECK(memcmp(zg, zs, sizeof(zg)) == 0, "vec4d_" #op, "mvg vs mvs"); \
    } while (0)

/* native vec4d vs mvg, same shapes */
#define NAT2(op, x, y) \
    do { \
        double zn[4], zg[4]; \
        vec4d_store_unaligned(zn, \
                vec4d_##op(vec4d_load_unaligned(x), \
                           vec4d_load_unaligned(y))); \
        mvg_vec4d_store_unaligned(zg, \
                mvg_vec4d_##op(mvg_vec4d_load_unaligned(x), \
                               mvg_vec4d_load_unaligned(y))); \
        CHECK(memcmp(zn, zg, sizeof(zn)) == 0, "vec4d_" #op, \
                "native vs generic"); \
    } while (0)

#define NAT3(op, x, y, w) \
    do { \
        double zn[4], zg[4]; \
        vec4d_store_unaligned(zn, \
                vec4d_##op(vec4d_load_unaligned(x), \
                           vec4d_load_unaligned(y), \
                           vec4d_load_unaligned(w))); \
        mvg_vec4d_store_unaligned(zg, \
                mvg_vec4d_##op(mvg_vec4d_load_unaligned(x), \
                               mvg_vec4d_load_unaligned(y), \
                               mvg_vec4d_load_unaligned(w))); \
        CHECK(memcmp(zn, zg, sizeof(zn)) == 0, "vec4d_" #op, \
                "native vs generic"); \
    } while (0)

#define NAT4(op, x, y, w, v) \
    do { \
        double zn[4], zg[4]; \
        vec4d_store_unaligned(zn, \
                vec4d_##op(vec4d_load_unaligned(x), \
                           vec4d_load_unaligned(y), \
                           vec4d_load_unaligned(w), \
                           vec4d_load_unaligned(v))); \
        mvg_vec4d_store_unaligned(zg, \
                mvg_vec4d_##op(mvg_vec4d_load_unaligned(x), \
                               mvg_vec4d_load_unaligned(y), \
                               mvg_vec4d_load_unaligned(w), \
                               mvg_vec4d_load_unaligned(v))); \
        CHECK(memcmp(zn, zg, sizeof(zn)) == 0, "vec4d_" #op, \
                "native vs generic"); \
    } while (0)

TEST_FUNCTION_START(machine_vectors_force_generic, state)
{
    slong iter;

    /* real arithmetic, on general values including signed zeros */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        double av[4], bv[4], cv[4];
        slong i;

        for (i = 0; i < 4; i++)
        {
            av[i] = 0.25 * fg_rand_int(state, 4000);
            bv[i] = 0.25 * fg_rand_int(state, 4000);
            cv[i] = 0.25 * fg_rand_int(state, 4000);
        }
        if (n_randint(state, 4) == 0)
        {
            av[n_randint(state, 4)] = n_randint(state, 2) ? -0.0 : 0.0;
            bv[n_randint(state, 4)] = n_randint(state, 2) ? -0.0 : 0.0;
            cv[n_randint(state, 4)] = n_randint(state, 2) ? -0.0 : 0.0;
        }
        if (n_randint(state, 4) == 0)
        {
            /* ties and large values for round */
            av[n_randint(state, 4)] = 0.5 + (double) n_randint(state, 100);
            av[n_randint(state, 4)] = ldexp(
                    (double) (1 + n_randint(state, 100)), 45);
        }

        TIERS1(neg); TIERS1(abs); TIERS1(half); TIERS1(floor);
        TIERS1(round);
        TIERS2(add, av, bv); TIERS2(sub, av, bv); TIERS2(mul, av, bv);
        TIERS2(min, av, bv); TIERS2(max, av, bv); TIERS2(addsub, av, bv);
        TIERS2(unpacklo, av, bv); TIERS2(unpackhi, av, bv);
        TIERS2(permute2_0_2, av, bv); TIERS2(permute2_1_3, av, bv);
        TIERS1(permute_0_2_1_3); TIERS1(permute_3_1_2_0);
        TIERS1(permute_3_2_1_0);
        TIERS3(blendv, av, bv, cv);
        TIERS3(fmadd, av, bv, cv); TIERS3(fmsub, av, bv, cv);
        TIERS3(fnmadd, av, bv, cv); TIERS3(fnmsub, av, bv, cv);
        TIERS2(cmp_ge, av, bv); TIERS2(cmp_gt, av, bv);
        TIERS2(cmp_lt, av, bv);
        if (bv[0] != 0.0 && bv[1] != 0.0 && bv[2] != 0.0 && bv[3] != 0.0)
            TIERS2(div, av, bv);

        /* native agreement for the IEEE operations and round */
        NAT2(add, av, bv); NAT2(sub, av, bv); NAT2(mul, av, bv);
        {
            double zn[4], zg[4];
            vec4d_store_unaligned(zn, vec4d_round(vec4d_load_unaligned(av)));
            mvg_vec4d_store_unaligned(zg,
                    mvg_vec4d_round(mvg_vec4d_load_unaligned(av)));
            CHECK(memcmp(zn, zg, sizeof(zn)) == 0, "vec4d_round",
                    "native vs generic");
        }
    }

#if FLINT_BITS == 64
    /* the modular contract domains: bit for bit across all three */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        double av[4], bv[4], nv[4], iv[4];
        ulong n;
        slong i;

        /* every modulus size the operations are used with, not just the
           50 bit primes of mpn_ctx: nmod_poly transforms over a prime
           input modulus directly from 20 bits up */
        {
            ulong bits = 4 + n_randint(state, 47);
            n = (UWORD(1) << (bits - 1))
                + 2*n_randint(state, UWORD(1) << (bits - 2)) + 1;
            if (n < 3)
                n = 3;
        }

        for (i = 0; i < 4; i++)
        {
            nv[i] = (double) n;
            iv[i] = 1.0 / (double) n;
            if (n_randint(state, 2))
            {
                av[i] = fg_rand_int(state, 2*n);
                bv[i] = fg_rand_int(state, 2*n - 1);
            }
            else
            {
                av[i] = fg_rand_int(state, 4*n - 1);
                bv[i] = fg_rand_int(state, n);
            }
        }

        TIERS4(mulmod, av, bv, nv, iv);
        TIERS4(nmulmod, av, bv, nv, iv);
        NAT4(mulmod, av, bv, nv, iv);
        NAT4(nmulmod, av, bv, nv, iv);

        /* both the range plain arithmetic handles and the range where
           only an exact method agrees with the fma backends */
        for (i = 0; i < 4; i++)
            av[i] = fg_rand_int(state, n_randint(state, 2) ? 8*n
                    : (UWORD(1) << 53) - n);
        if (n_randint(state, 8) == 0)
            av[n_randint(state, 4)] = n_randint(state, 2) ? -0.0 : 0.0;

        TIERS3(reduce_to_pm1n, av, nv, iv);
        TIERS3(reduce_to_pm1no, av, nv, iv);
        TIERS3(reduce_to_0n, av, nv, iv);
        NAT3(reduce_to_pm1n, av, nv, iv);
        NAT3(reduce_to_pm1no, av, nv, iv);
        NAT3(reduce_to_0n, av, nv, iv);

        /* the pure range maps, mvg vs mvs, including the boundaries */
        for (i = 0; i < 4; i++)
            av[i] = 0.5 * fg_rand_int(state, 2*n);
        if (n_randint(state, 4) == 0)
            av[n_randint(state, 4)] = n_randint(state, 2) ? -0.0 : 0.0;

        TIERS2(reduce_pm1no_to_0n, av, nv);
        TIERS2(reduce_pm1n_to_pmhn, av, nv);
        for (i = 0; i < 4; i++)
            av[i] = fabs(av[i]);
        TIERS2(reduce_0n_to_pmhn, av, nv);
        TIERS2(reduce_2n_to_n, av, nv);
    }

    /* integer operations and the conversions */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        ulong a[4], b[4], nn[4], zn[4], zg[4], zs[4];
        double d[4], en[4], eg[4], es[4];
        slong i;

        for (i = 0; i < 4; i++)
        {
            nn[i] = n_randtest_not_zero(state);
            a[i] = n_randlimb(state) % nn[i];
            b[i] = n_randlimb(state) % nn[i];
        }

        vec4n_store_unaligned(zn, vec4n_addmod(vec4n_load_unaligned(a),
                    vec4n_load_unaligned(b), vec4n_load_unaligned(nn)));
        mvg_vec4n_store_unaligned(zg, mvg_vec4n_addmod(
                    mvg_vec4n_load_unaligned(a),
                    mvg_vec4n_load_unaligned(b),
                    mvg_vec4n_load_unaligned(nn)));
        mvs_vec4n_store_unaligned(zs, mvs_vec4n_addmod(
                    mvs_vec4n_load_unaligned(a),
                    mvs_vec4n_load_unaligned(b),
                    mvs_vec4n_load_unaligned(nn)));
        CHECK(memcmp(zn, zg, sizeof(zn)) == 0
                && memcmp(zg, zs, sizeof(zg)) == 0, "vec4n_addmod",
                "native vs generic tiers");

        for (i = 0; i < 4; i++)
        {
            nn[i] >>= 1;
            nn[i] += (nn[i] == 0);
            a[i] %= nn[i];
            b[i] %= nn[i];
        }
        vec4n_store_unaligned(zn, vec4n_addmod_limited(
                    vec4n_load_unaligned(a),
                    vec4n_load_unaligned(b), vec4n_load_unaligned(nn)));
        mvg_vec4n_store_unaligned(zg, mvg_vec4n_addmod_limited(
                    mvg_vec4n_load_unaligned(a),
                    mvg_vec4n_load_unaligned(b),
                    mvg_vec4n_load_unaligned(nn)));
        mvs_vec4n_store_unaligned(zs, mvs_vec4n_addmod_limited(
                    mvs_vec4n_load_unaligned(a),
                    mvs_vec4n_load_unaligned(b),
                    mvs_vec4n_load_unaligned(nn)));
        CHECK(memcmp(zn, zg, sizeof(zn)) == 0
                && memcmp(zg, zs, sizeof(zg)) == 0, "vec4n_addmod_limited",
                "native vs generic tiers");

        /* conversions on [0, 2^52) */
        for (i = 0; i < 4; i++)
        {
            a[i] = n_randint(state, UWORD(1) << 52);
            d[i] = (double) a[i];
        }

        vec4d_store_unaligned(en, vec4n_convert_limited_vec4d(
                    vec4n_load_unaligned(a)));
        mvg_vec4d_store_unaligned(eg, mvg_vec4n_convert_limited_vec4d(
                    mvg_vec4n_load_unaligned(a)));
        mvs_vec4d_store_unaligned(es, mvs_vec4n_convert_limited_vec4d(
                    mvs_vec4n_load_unaligned(a)));
        CHECK(memcmp(en, eg, sizeof(en)) == 0
                && memcmp(eg, es, sizeof(eg)) == 0,
                "vec4n_convert_limited_vec4d", "native vs generic tiers");

        vec4n_store_unaligned(zn, vec4d_convert_limited_vec4n(
                    vec4d_load_unaligned(d)));
        mvg_vec4n_store_unaligned(zg, mvg_vec4d_convert_limited_vec4n(
                    mvg_vec4d_load_unaligned(d)));
        mvs_vec4n_store_unaligned(zs, mvs_vec4d_convert_limited_vec4n(
                    mvs_vec4d_load_unaligned(d)));
        CHECK(memcmp(zn, zg, sizeof(zn)) == 0
                && memcmp(zg, zs, sizeof(zg)) == 0,
                "vec4d_convert_limited_vec4n", "native vs generic tiers");
    }
#endif

    TEST_FUNCTION_END(state);
}

#undef CHECK
