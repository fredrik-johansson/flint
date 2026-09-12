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
    Tests of the documented contracts of the floating point modular
    arithmetic in machine_vectors.h, for moduli 2^49 < n < 2^50 with
    ninv = 1.0/n, as used by fft_small:

    mulmod/nmulmod, for integer valued operands:
      - the result is an integer congruent to +-(a*b) mod n (checked
        exactly against integer arithmetic);
      - for any such n: |a*b| < 2n^2 gives a result in (-9/8 n, 9/8 n)
        and |a*b| < 4n^2 gives a result in (-7/4 n, 7/4 n);
      - if fft_small_mulmod_satisfies_bounds(n): |a*b| < 2n^2 gives
        (-n, n) and |a*b| < 4n^2 gives (-3/2 n, 3/2 n);
      - nmulmod(a, b) is exactly -mulmod(a, b), with +0.0 for both when
        the result is zero.

    reduce_to_pm1n/reduce_to_pm1no, for every integer valued a that a
    double holds exactly and that leaves the quotient times the modulus
    representable, that is |a| <= 2^53 - n:
      the result is a mod n in (-n, n), never -0.0. Note the range. A
      naive floating point a - round(a*ninv)*n is exact only while the
      quotient satisfies |q| <= 8, and fft_small calls these with larger
      quotients, so the tests below go up to |a| = 2^53, where |q| can
      reach 16.

    reduce_to_0n likewise: the result is a mod n in [0, n).

    The pure range->range maps have exact per lane semantics that
    fft_small relies on:
      reduce_pm1no_to_0n(a, n)   adds n when the sign bit of a is set
                                 (so -0.0 maps to +n);
      reduce_2n_to_n(a, n)       subtracts n when a - n >= +0.0;
      reduce_0n_to_pmhn(a, n)    subtracts n when a > n/2;
      reduce_pm1n_to_pmhn(a, n)  adds -+n when |a| > n/2.

    A vector operation must behave as its scalar counterpart in every
    lane, so all widths get the same checks; which widths provide which
    operations depends on the backend (the generic backend provides
    everything, see doc/source/machine_vectors.rst).
*/

#if FLINT_BITS == 64

#include "fft_small.h"

#define CHECK(cond, name, extra) \
    do { \
        if (!(cond)) \
        { \
            flint_printf("FAIL: %s (%s), n = %wu, lane %wd\n", \
                                                    name, extra, n, i); \
            fflush(stdout); \
            flint_abort(); \
        } \
    } while (0)

/* random integer valued double in [-b, b] */
static double mod_rand_int(flint_rand_t state, ulong b)
{
    ulong u = n_randint(state, 2*b + 1);
    return (double) ((slong) u - (slong) b);
}

/* canonical residue of an integer valued |a| < 2^62 */
static ulong mod_res(double a, ulong n)
{
    slong ia = (slong) a;
    slong r = ia % (slong) n;
    return (ulong) (r < 0 ? r + (slong) n : r);
}

static int mod_bit_same(double a, double b)
{
    return memcmp(&a, &b, sizeof(double)) == 0;
}

/*
    mulmod, nmulmod, reduce_to_pm1n, reduce_to_pm1no, reduce_to_0n,
    reduce_pm1no_to_0n, reduce_pm1n_to_pmhn: every backend provides
    these for every width it provides at all.
*/
#define DEFINE_CHECK_CORE(V, K) \
static void V##_check_core(flint_rand_t state, ulong n, int satisfies) \
{ \
    double av[K], bv[K], zv[K], wv[K]; \
    slong i; \
    int which = (int) n_randint(state, 2); \
    V a, b; \
    V vn = V##_set_d((double) n); \
    V vninv = V##_set_d(1.0 / (double) n); \
 \
    /* which = 0: |a*b| < 2n^2;  which = 1: |a*b| < 4n^2 */ \
    for (i = 0; i < K; i++) \
    { \
        if (which == 0) \
        { \
            av[i] = mod_rand_int(state, 2*n); \
            bv[i] = mod_rand_int(state, n - 1); \
        } \
        else if (n_randint(state, 2)) \
        { \
            av[i] = mod_rand_int(state, 2*n); \
            bv[i] = mod_rand_int(state, 2*n - 1); \
        } \
        else \
        { \
            av[i] = mod_rand_int(state, 4*n - 1); \
            bv[i] = mod_rand_int(state, n); \
        } \
        if (n_randint(state, 16) == 0) \
            av[i] = 0.0; \
    } \
 \
    a = V##_load_unaligned(av); \
    b = V##_load_unaligned(bv); \
 \
    V##_store_unaligned(zv, V##_mulmod(a, b, vn, vninv)); \
    V##_store_unaligned(wv, V##_nmulmod(a, b, vn, vninv)); \
    for (i = 0; i < K; i++) \
    { \
        slong r = (slong) zv[i]; \
        ulong m = n_mulmod2_preinv(mod_res(av[i], n), mod_res(bv[i], n), \
                                                n, n_preinvert_limb(n)); \
        CHECK(zv[i] == (double) r, #V "_mulmod", "not integer valued"); \
        CHECK(mod_res(zv[i], n) == m, #V "_mulmod", "wrong residue"); \
        if (which == 0) \
            CHECK(FLINT_ABS(r) * 8 < 9 * (slong) n, #V "_mulmod", \
                    "out of (-9/8 n, 9/8 n)"); \
        else \
            CHECK(FLINT_ABS(r) * 4 < 7 * (slong) n, #V "_mulmod", \
                    "out of (-7/4 n, 7/4 n)"); \
        if (satisfies) \
        { \
            if (which == 0) \
                CHECK(FLINT_ABS(r) < (slong) n, #V "_mulmod", \
                        "out of (-n, n)"); \
            else \
                CHECK(FLINT_ABS(r) * 2 < 3 * (slong) n, #V "_mulmod", \
                        "out of (-3/2 n, 3/2 n)"); \
        } \
        CHECK(mod_bit_same(wv[i], zv[i] == 0.0 ? 0.0 : -zv[i]), \
                #V "_nmulmod", "not -mulmod"); \
    } \
 \
    /* reductions of integer valued a, both in the range that plain \
       arithmetic handles and in the range where it does not */ \
    for (i = 0; i < K; i++) \
        av[i] = mod_rand_int(state, n_randint(state, 2) ? 8*n \
                : (UWORD(1) << 53) - n); \
    if (n_randint(state, 8) == 0) \
        av[n_randint(state, K)] = n_randint(state, 2) ? -0.0 : 0.0; \
    a = V##_load_unaligned(av); \
 \
    V##_store_unaligned(zv, V##_reduce_to_pm1n(a, vn, vninv)); \
    for (i = 0; i < K; i++) \
    { \
        slong r = (slong) zv[i]; \
        CHECK(zv[i] == (double) r && FLINT_ABS(r) < (slong) n \
                && mod_res(zv[i], n) == mod_res(av[i], n) \
                && !mod_bit_same(zv[i], -0.0), \
                #V "_reduce_to_pm1n", ""); \
    } \
 \
    V##_store_unaligned(zv, V##_reduce_to_pm1no(a, vn, vninv)); \
    for (i = 0; i < K; i++) \
    { \
        slong r = (slong) zv[i]; \
        CHECK(zv[i] == (double) r && FLINT_ABS(r) < (slong) n \
                && mod_res(zv[i], n) == mod_res(av[i], n) \
                && !mod_bit_same(zv[i], -0.0), \
                #V "_reduce_to_pm1no", ""); \
    } \
 \
    V##_store_unaligned(zv, V##_reduce_to_0n(a, vn, vninv)); \
    for (i = 0; i < K; i++) \
    { \
        slong r = (slong) zv[i]; \
        CHECK(zv[i] == (double) r && 0 <= r && r < (slong) n \
                && (ulong) r == mod_res(av[i], n), \
                #V "_reduce_to_0n", ""); \
    } \
 \
    /* the range->range maps need not receive integers; feed them \
       half integers and signed zeros */ \
    for (i = 0; i < K; i++) \
        av[i] = 0.5 * mod_rand_int(state, 2*n); \
    if (n_randint(state, 8) == 0) \
        av[n_randint(state, K)] = n_randint(state, 2) ? -0.0 : 0.0; \
    a = V##_load_unaligned(av); \
 \
    V##_store_unaligned(zv, V##_reduce_pm1no_to_0n(a, vn)); \
    for (i = 0; i < K; i++) \
    { \
        /* the vector blendv keys on the sign bit, so -0.0 becomes n, \
           while the scalar one keys on a >= 0 and leaves -0.0 alone */ \
        int neg = (K == 1) ? (av[i] < 0.0) : (signbit(av[i]) != 0); \
        CHECK(mod_bit_same(zv[i], neg ? av[i] + (double) n : av[i]), \
                #V "_reduce_pm1no_to_0n", ""); \
    } \
 \
    V##_store_unaligned(zv, V##_reduce_pm1n_to_pmhn(a, vn)); \
    for (i = 0; i < K; i++) \
    { \
        double e = av[i]; \
        if (fabs(e) > 0.5 * (double) n) \
            e = signbit(e) ? e + (double) n : e - (double) n; \
        CHECK(mod_bit_same(zv[i], e), #V "_reduce_pm1n_to_pmhn", ""); \
    } \
}

/* reduce_2n_to_n */
#define DEFINE_CHECK_2N(V, K) \
static void V##_check_2n(flint_rand_t state, ulong n) \
{ \
    double av[K], zv[K]; \
    slong i; \
    V a, vn = V##_set_d((double) n); \
 \
    for (i = 0; i < K; i++) \
        av[i] = 0.5 * (double) n_randint(state, 4*n); \
    if (n_randint(state, 8) == 0) \
        av[n_randint(state, K)] = (double) n; \
    a = V##_load_unaligned(av); \
    V##_store_unaligned(zv, V##_reduce_2n_to_n(a, vn)); \
    for (i = 0; i < K; i++) \
        CHECK(mod_bit_same(zv[i], av[i] - (double) n >= 0.0 ? \
                                av[i] - (double) n : av[i]), \
                #V "_reduce_2n_to_n", ""); \
}

/* reduce_0n_to_pmhn */
#define DEFINE_CHECK_0N_PMHN(V, K) \
static void V##_check_0n_pmhn(flint_rand_t state, ulong n) \
{ \
    double av[K], zv[K]; \
    slong i; \
    V a, vn = V##_set_d((double) n); \
 \
    for (i = 0; i < K; i++) \
        av[i] = 0.5 * (double) n_randint(state, 2*n + 1); \
    a = V##_load_unaligned(av); \
    V##_store_unaligned(zv, V##_reduce_0n_to_pmhn(a, vn)); \
    for (i = 0; i < K; i++) \
        CHECK(mod_bit_same(zv[i], av[i] > 0.5 * (double) n ? \
                                av[i] - (double) n : av[i]), \
                #V "_reduce_0n_to_pmhn", ""); \
}

/* same_mod on integer valued inputs in (-n, n) */
#define DEFINE_CHECK_SAME_MOD(V, K) \
static void V##_check_same_mod(flint_rand_t state, ulong n) \
{ \
    double av[K], bv[K]; \
    slong i; \
    V a, b; \
    V vn = V##_set_d((double) n); \
    V vninv = V##_set_d(1.0 / (double) n); \
    int s, sref = 1; \
 \
    for (i = 0; i < K; i++) \
    { \
        av[i] = mod_rand_int(state, n - 1); \
        bv[i] = n_randint(state, 2) ? av[i] : mod_rand_int(state, n - 1); \
        if (n_randint(state, 2) && av[i] > 0.0) \
            bv[i] = av[i] - (double) n; \
    } \
    a = V##_load_unaligned(av); \
    b = V##_load_unaligned(bv); \
    s = V##_same_mod(a, b, vn, vninv); \
    for (i = 0; i < K; i++) \
        sref = sref && (mod_res(av[i], n) == mod_res(bv[i], n)); \
    i = -1; \
    CHECK((s != 0) == (sref != 0), #V "_same_mod", ""); \
}

/*
    Which widths exist, and which of them have the three optional
    operations, depends on the backend:

                        AVX2            NEON            generic
    widths              1, 4, 8         1, 2, 4, 8      1, 2, 4, 8
    reduce_2n_to_n      1, 4, 8         -               1, 2, 4, 8
    reduce_0n_to_pmhn   1, 4            1, 2            1, 2, 4, 8
    same_mod            1, 4            1, 2            1, 2, 4, 8
*/

DEFINE_CHECK_CORE(vec1d, 1)
DEFINE_CHECK_CORE(vec4d, 4)
DEFINE_CHECK_CORE(vec8d, 8)
DEFINE_CHECK_0N_PMHN(vec1d, 1)

#if !defined(FLINT_MACHINE_VECTORS_AVX2)
DEFINE_CHECK_CORE(vec2d, 2)
DEFINE_CHECK_0N_PMHN(vec2d, 2)
DEFINE_CHECK_SAME_MOD(vec2d, 2)
#endif

#if !defined(FLINT_MACHINE_VECTORS_NEON)
DEFINE_CHECK_2N(vec1d, 1)
DEFINE_CHECK_2N(vec4d, 4)
DEFINE_CHECK_2N(vec8d, 8)
DEFINE_CHECK_0N_PMHN(vec4d, 4)
DEFINE_CHECK_SAME_MOD(vec4d, 4)
#endif

DEFINE_CHECK_SAME_MOD(vec1d, 1)

#if defined(FLINT_MACHINE_VECTORS_GENERIC)
DEFINE_CHECK_2N(vec2d, 2)
DEFINE_CHECK_0N_PMHN(vec8d, 8)
DEFINE_CHECK_SAME_MOD(vec8d, 8)
#endif

#endif /* FLINT_BITS == 64 */

TEST_FUNCTION_START(machine_vectors_mod, state)
{
#if FLINT_BITS == 64
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        ulong n;
        ulong bits;
        int satisfies;

        /*
            Odd moduli of every size up to the 50 bit limit, prime or
            not. 50 bit primes are what mpn_ctx uses, but a modulus the
            caller chose is also used directly: nmod_poly transforms
            over the input modulus itself whenever it is prime, below
            2^50 and 2-adically deep enough (see
            _nmod_poly_should_directly_fft), which starts at 20 bits.
            Smaller moduli than that are covered too, since nothing in
            these operations knows where the modulus came from.
        */
        bits = 4 + n_randint(state, 47);

        switch (n_randint(state, 6))
        {
            case 0:
                n = UWORD(0x0003f00000000001);   /* the default fft prime */
                break;
            case 1:
                n = n_randprime(state, 50, 0);
                break;
            case 2:
                n = n_randprime(state, FLINT_MAX(bits, 3), 0);
                break;
            case 3:
                /* the extremes of each size */
                n = n_randint(state, 2) ? (UWORD(1) << bits) - 1
                                        : (UWORD(1) << (bits - 1)) + 1;
                break;
            default:
                n = (UWORD(1) << (bits - 1))
                    + 2*n_randint(state, UWORD(1) << (bits - 2)) + 1;
                break;
        }

        if (n < 3 || n >= (UWORD(1) << 50) || n % 2 == 0)
            continue;

        satisfies = fft_small_mulmod_satisfies_bounds(n);

        vec1d_check_core(state, n, satisfies);
        vec4d_check_core(state, n, satisfies);
        vec8d_check_core(state, n, satisfies);
        vec1d_check_0n_pmhn(state, n);
        vec1d_check_same_mod(state, n);

#if !defined(FLINT_MACHINE_VECTORS_AVX2)
        vec2d_check_core(state, n, satisfies);
        vec2d_check_0n_pmhn(state, n);
        vec2d_check_same_mod(state, n);
#endif

#if !defined(FLINT_MACHINE_VECTORS_NEON)
        vec1d_check_2n(state, n);
        vec4d_check_2n(state, n);
        vec8d_check_2n(state, n);
        vec4d_check_0n_pmhn(state, n);
        vec4d_check_same_mod(state, n);
#endif

#if defined(FLINT_MACHINE_VECTORS_GENERIC)
        vec2d_check_2n(state, n);
        vec8d_check_0n_pmhn(state, n);
        vec8d_check_same_mod(state, n);
#endif
    }
#endif

    TEST_FUNCTION_END(state);
}

#undef CHECK
