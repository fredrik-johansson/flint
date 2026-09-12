/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023 Mathieu Gouttenoire

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MACHINE_VECTORS_GENERIC_H
#define MACHINE_VECTORS_GENERIC_H

#ifndef MACHINE_VECTORS_H
# error machine_vectors_generic.h should be included only via machine_vectors.h
#endif

/*
    Generic backends, so that this header is usable on every platform rather
    than only AVX2 and NEON. Two tiers:

    (a) GNU vector extensions (GCC/Clang on any architecture: RISC-V, POWER,
        s390x, ARM without NEON, ...). These are compiler extensions rather
        than target intrinsics, so one source maps onto whatever SIMD the
        target actually has, and -- importantly -- the operations are already
        vector operations, so the compiler does not try to auto-vectorize
        loops over them in unprofitable ways.

    (b) Strict ISO C structs, for any other compiler. Correct everywhere;
        quality of codegen is up to the compiler's own auto-vectorizer.
        Define FLINT_MACHINE_VECTORS_STRICT_C to force this tier.

    On a 64-bit word size these tiers implement the full interface, including
    the integer vectors and the modular operations used by fft_small; on a
    32-bit word size only the floating point subset needed by
    flint_sgemm/flint_dgemm is provided.

    The modular operations must produce exactly the same values as the AVX2
    and NEON backends, whose correctness analysis (see
    fft_small/mulmod_satisfies_bounds.c) depends on the quotient
    q = round(mul(h, ninv)) with h = mul(a, b) being formed with one rounding
    per operation and on the remainder a*b - q*n then being computed exactly.
    AVX2 and NEON get the exact remainder from correctly rounded FMA. Portable
    C has no such primitive (fma() may be an emulation running dozens of times
    slower), so here the same quotient q is formed in pure floating point
    arithmetic -- bit for bit the value the FMA backends compute, since it
    involves no fused operations -- and the remainder a*b - q*n is computed
    exactly in 64-bit integer arithmetic with wraparound, in the style of
    n_mulmod_precomp. The remainder fits far below 2^63, so the wrapped
    difference of the wrapped products recovers it exactly, and the final
    conversion back to double is exact. The double <-> int64 conversions do
    not assume any conversion instruction: they use the usual 2^52 shifting
    tricks so that the GNU vector tier stays vectorizable on targets, such as
    plain x86-64 SSE2, without packed 64-bit conversions.

    Everything here assumes round-to-nearest mode and no
    contraction-with-reassociation (-ffast-math would break the shifting
    tricks; FLINT is not built with it).
*/

typedef double vec1d;
typedef float vec1f;
#if FLINT_BITS == 64
typedef ulong vec1n;
#endif

/* scalar helpers ***********************************************************/

#if FLINT_BITS == 64

/* round to nearest even, exact conversion to slong; assumes |a| < 2^52 */
FLINT_FORCE_INLINE slong _vec_generic_round_i64_52(double a)
{
    double t = fabs(a) + 0x1.0p52;
    slong i;
    memcpy(&i, &t, sizeof(i));
    i -= (slong) WORD(0x4330000000000000);
    return a < 0 ? -i : i;
}

/* exact conversion of a signed integer with |r| < 2^51 to double */
FLINT_FORCE_INLINE double _vec_generic_i64_51_get_d(ulong r)
{
    double d;
    r += UWORD(0x4338000000000000);
    memcpy(&d, &r, sizeof(d));
    return d - 0x1.8p52;
}

/*
    mulmod. See the comment above and doc/source/machine_vectors.rst for the
    contract. The quotient below is identical to the one the FMA backends
    compute: h = fl(a*b) is the same one rounding, x = fl(h*ninv) the same
    second rounding, and rounding x to an integer (nearest, ties to even) is
    the same third. The remainder a*b - q*n is then computed exactly modulo
    2^64; since it lies in (-2n, 2n) (or slightly beyond for a modulus not
    satisfying fft_small_mulmod_satisfies_bounds), it is recovered exactly.
    The casts of a, b and n to slong are exact because the operands are
    integer valued and bounded well below 2^62.
*/
/*
    On a target whose fma is a single instruction, the exact remainder is
    much cheaper to get the way the AVX2 and NEON backends get it: l = a*b - h
    is the part of the product h lost, and (h - q*n) + l is the exact integer.
    __FP_FAST_FMA promises the hardware instruction; where it is absent, fma()
    is a software routine and far slower than the integer detour below, so the
    two paths are not interchangeable. Both compute the same exact integer, so
    the results agree with each other and with the other backends either way.
*/
/* FP_FAST_FMA is the standard spelling, from math.h, which
   machine_vectors.h has already included; __FP_FAST_FMA is what GCC and
   clang predefine and what some libcs derive the former from, so accept
   either. configure tests the same pair. */
#if defined(FP_FAST_FMA) || defined(__FP_FAST_FMA)
# define _VEC_GENERIC_USE_FMA 1
#else
# define _VEC_GENERIC_USE_FMA 0
#endif

FLINT_FORCE_INLINE double _vec_generic_mulmod1(double a, double b, double n, double ninv)
{
    double h = a * b;
#if _VEC_GENERIC_USE_FMA
    double q = rint(h * ninv);
    return fma(-q, n, h) + fma(a, b, -h);
#else
    slong iq = _vec_generic_round_i64_52(h * ninv);
    ulong r = (ulong)(slong) a * (ulong)(slong) b - (ulong) iq * (ulong)(slong) n;
    return _vec_generic_i64_51_get_d(r);
#endif
}

FLINT_FORCE_INLINE double _vec_generic_nmulmod1(double a, double b, double n, double ninv)
{
    double h = a * b;
#if _VEC_GENERIC_USE_FMA
    double q = rint(h * ninv);
    return fma(q, n, -h) - fma(a, b, -h);
#else
    slong iq = _vec_generic_round_i64_52(h * ninv);
    ulong r = (ulong) iq * (ulong)(slong) n - (ulong)(slong) a * (ulong)(slong) b;
    return _vec_generic_i64_51_get_d(r);
#endif
}

/*
    a - round(a*ninv)*n. The plain expression is exact only while q*n stays
    below 2^53, that is for |q| <= 8, and fft_small does call this with
    somewhat larger quotients. The FMA backends are exact for every q, since
    the result is small and the fma rounds only once, so the same has to hold
    here: n is split as nh + nl with nh a multiple of 2^26 and |nl| <= 2^25,
    which makes q*nh, q*nl and both subtractions exact for every |q| < 2^27,
    hence for every integer valued a a double can hold. The split depends
    only on n and is hoisted out of the caller's loops, where n is invariant,
    so the marginal cost is one multiply and one subtract.
*/
FLINT_FORCE_INLINE double _vec_generic_reduce_to_pm1n1(double a, double n, double ninv)
{
    double x = a * ninv;
    double t = x + 0x1.8p52;
    double q = t - 0x1.8p52;
    double nh = ((n * 0x1p-26 + 0x1.8p52) - 0x1.8p52) * 0x1p26;
    double nl = n - nh;
    /* the shift trick returns +0.0 where rint returns -0.0; restoring the
       sign is what makes a = -0.0 give +0.0 here, as on the FMA backends,
       rather than -0.0, which the [0,n) reductions would turn into n */
    q = copysign(q, x);
    return (a - q * nh) - q * nl;
}

#endif

#if defined(__GNUC__) && !defined(FLINT_MACHINE_VECTORS_STRICT_C)

# define FLINT_MACHINE_VECTORS_GNU_VECTOR_EXTENSIONS 1

/*
    These vectors may be wider than the target's native SIMD width, which
    is intentional (extra instruction level parallelism, and the compiler
    splits them). All functions here are force-inlined, so the psabi
    warning about passing such types across function boundaries does not
    apply to any code we actually generate.

    The warning has to be turned off for the rest of the translation unit
    rather than only for this header: GCC reports it while expanding the
    functions of the includer, and therefore at a location in the
    includer, which a push/pop around the header would not cover.
*/
# if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wpsabi"
# endif

typedef double vec2d __attribute__((vector_size(16)));
typedef double vec4d __attribute__((vector_size(32)));
typedef float vec4f __attribute__((vector_size(16)));
typedef float vec8f __attribute__((vector_size(32)));

/* internal: same-sized integer views for sign tricks and shuffle masks */
typedef long long _vec2i __attribute__((vector_size(16)));
typedef long long _vec4i __attribute__((vector_size(32)));
typedef unsigned long long _vec2u __attribute__((vector_size(16)));
typedef unsigned long long _vec4u __attribute__((vector_size(32)));

#if FLINT_BITS == 64
typedef ulong vec2n __attribute__((vector_size(16)));
typedef ulong vec4n __attribute__((vector_size(32)));
#endif

# if defined(__clang__)
#  define _MV_SHUF2(a, b, i0, i1) __builtin_shufflevector(a, b, i0, i1)
#  define _MV_SHUF4(a, b, i0, i1, i2, i3) \
                               __builtin_shufflevector(a, b, i0, i1, i2, i3)
# else
#  define _MV_SHUF2(a, b, i0, i1) \
                    __builtin_shuffle(a, b, (_vec2i){i0, i1})
#  define _MV_SHUF4(a, b, i0, i1, i2, i3) \
                    __builtin_shuffle(a, b, (_vec4i){i0, i1, i2, i3})
# endif

# define VEC_GENERIC_DEF(V, S, N, SUF, FLOOR1) \
FLINT_FORCE_INLINE V V##_load_unaligned(const S* a) { \
    V z; memcpy(&z, a, sizeof(z)); return z; \
} \
FLINT_FORCE_INLINE V V##_load_aligned(const S* a) { \
    V z; memcpy(&z, a, sizeof(z)); return z; \
} \
FLINT_FORCE_INLINE V V##_load(const S* a) { \
    return V##_load_aligned(a); \
} \
FLINT_FORCE_INLINE void V##_store_unaligned(S* z, V a) { \
    memcpy(z, &a, sizeof(a)); \
} \
FLINT_FORCE_INLINE void V##_store_aligned(S* z, V a) { \
    memcpy(z, &a, sizeof(a)); \
} \
FLINT_FORCE_INLINE void V##_store(S* z, V a) { \
    V##_store_aligned(z, a); \
} \
FLINT_FORCE_INLINE V V##_zero(void) { V z = {0}; return z; } \
/* a - z rather than z + a: the latter turns a = -0.0 into +0.0, and \
   -0.0 has to survive, both because callers may broadcast it and \
   because a sign mask is built this way */ \
FLINT_FORCE_INLINE V V##_set_##SUF(S a) { V z = {0}; return a - z; } \
FLINT_FORCE_INLINE V V##_add(V a, V b) { return a + b; } \
FLINT_FORCE_INLINE V V##_sub(V a, V b) { return a - b; } \
FLINT_FORCE_INLINE V V##_mul(V a, V b) { return a * b; } \
FLINT_FORCE_INLINE V V##_fmadd(V a, V b, V c) { return a * b + c; } \
FLINT_FORCE_INLINE V V##_fnmadd(V a, V b, V c) { return c - a * b; } \
FLINT_FORCE_INLINE V V##_floor(V a) { \
    for (int i = 0; i < N; i++) a[i] = FLOOR1(a[i]); \
    return a; \
}

/*
    The double widths get the full interface. IV/UV are the signed/unsigned
    integer views of V used for bit manipulation; a cast between GNU vector
    types of equal size reinterprets the representation.
*/
# define VEC_GENERIC_DEF_D(V, IV, UV, N) \
FLINT_FORCE_INLINE double V##_get_index(V a, const int i) { return a[i]; } \
FLINT_FORCE_INLINE int V##_same(V a, V b) { \
    for (int i = 0; i < N; i++) \
        if (a[i] != b[i]) \
            return 0; \
    return 1; \
} \
FLINT_FORCE_INLINE V V##_one(void) { return V##_set_d(1.0); } \
/* the sign mask is built in the integer view: going through \
   set_d(-0.0) would leave the compiler to fold a floating point \
   broadcast, which it does not always manage */ \
FLINT_FORCE_INLINE UV _##V##_signmask(void) { \
    UV z = {0}; \
    return z + UWORD(0x8000000000000000); \
} \
FLINT_FORCE_INLINE V V##_neg(V a) { \
    return (V) ((UV) a ^ _##V##_signmask()); \
} \
FLINT_FORCE_INLINE V V##_abs(V a) { \
    return (V) ((UV) a & ~_##V##_signmask()); \
} \
FLINT_FORCE_INLINE V V##_half(V a) { return a * 0.5; } \
FLINT_FORCE_INLINE V V##_div(V a, V b) { return a / b; } \
FLINT_FORCE_INLINE V V##_fmsub(V a, V b, V c) { return a * b - c; } \
FLINT_FORCE_INLINE V V##_fnmsub(V a, V b, V c) { return -(a * b) - c; } \
FLINT_FORCE_INLINE V V##_cmp_ge(V a, V b) { return (V) (a >= b); } \
FLINT_FORCE_INLINE V V##_cmp_gt(V a, V b) { return (V) (a > b); } \
FLINT_FORCE_INLINE V V##_cmp_lt(V a, V b) { return (V) (a < b); } \
FLINT_FORCE_INLINE V V##_min(V a, V b) { \
    IV m = (IV) (a < b); \
    return (V) (((UV) a & (UV) m) | ((UV) b & ~(UV) m)); \
} \
FLINT_FORCE_INLINE V V##_max(V a, V b) { \
    IV m = (IV) (b < a); \
    return (V) (((UV) a & (UV) m) | ((UV) b & ~(UV) m)); \
} \
FLINT_FORCE_INLINE V V##_blendv(V a, V b, V c) { \
    IV m = (IV) c >> 63; \
    return (V) (((UV) a & ~(UV) m) | ((UV) b & (UV) m)); \
} \
FLINT_FORCE_INLINE V V##_round(V a) { \
    V big = V##_set_d(0x1.0p52); \
    V ax = V##_abs(a); \
    V r = (ax + big) - big; \
    r = (V) (((UV) r) | ((UV) a & _##V##_signmask())); \
    return V##_blendv(a, r, V##_cmp_lt(ax, big)); \
} \
FLINT_FORCE_INLINE V V##_reduce_0n_to_pmhn(V a, V n) { \
    return V##_blendv(a, a - n, V##_cmp_gt(a, V##_half(n))); \
} \
FLINT_FORCE_INLINE V V##_reduce_pm1n_to_pmhn(V a, V n) { \
    V t = V##_blendv(n, V##_neg(n), a); \
    return V##_blendv(a, a - t, V##_cmp_gt(V##_abs(a), V##_half(n))); \
} \
FLINT_FORCE_INLINE V V##_reduce_2n_to_n(V a, V n) { \
    V s = a - n; \
    return V##_blendv(s, a, s); \
} \
FLINT_FORCE_INLINE V V##_reduce_pm1no_to_0n(V a, V n) { \
    return V##_blendv(a, a + n, a); \
}

VEC_GENERIC_DEF(vec2d, double, 2, d, floor)
VEC_GENERIC_DEF(vec4d, double, 4, d, floor)
VEC_GENERIC_DEF(vec4f, float, 4, f, floorf)
VEC_GENERIC_DEF(vec8f, float, 8, f, floorf)

#undef VEC_GENERIC_DEF

VEC_GENERIC_DEF_D(vec2d, _vec2i, _vec2u, 2)
VEC_GENERIC_DEF_D(vec4d, _vec4i, _vec4u, 4)

#undef VEC_GENERIC_DEF_D

FLINT_FORCE_INLINE vec2d vec2d_set_d2(double a0, double a1) {
    vec2d z = {a0, a1};
    return z;
}

FLINT_FORCE_INLINE vec4d vec4d_set_d4(double a0, double a1, double a2, double a3) {
    vec4d z = {a0, a1, a2, a3};
    return z;
}

FLINT_FORCE_INLINE vec4d vec4d_set_vec2d2(vec2d a, vec2d b) {
    vec4d z = {a[0], a[1], b[0], b[1]};
    return z;
}

FLINT_FORCE_INLINE void vec4d_print(vec4d a) {
    flint_printf("{%f, %f, %f, %f}", a[0], a[1], a[2], a[3]);
}

FLINT_FORCE_INLINE vec2d vec2d_addsub(vec2d a, vec2d b) {
    vec2d s = {-0.0, 0.0};
    return a + (vec2d) ((_vec2u) b ^ (_vec2u) s);
}

FLINT_FORCE_INLINE vec4d vec4d_addsub(vec4d a, vec4d b) {
    vec4d s = {-0.0, 0.0, -0.0, 0.0};
    return a + (vec4d) ((_vec4u) b ^ (_vec4u) s);
}

/* permutations, with the same 128-bit-lane semantics as AVX2 and NEON */

FLINT_FORCE_INLINE vec2d vec2d_unpacklo(vec2d a, vec2d b) {
    return _MV_SHUF2(a, b, 0, 2);
}

FLINT_FORCE_INLINE vec2d vec2d_unpackhi(vec2d a, vec2d b) {
    return _MV_SHUF2(a, b, 1, 3);
}

FLINT_FORCE_INLINE vec4d vec4d_unpacklo(vec4d a, vec4d b) {
    return _MV_SHUF4(a, b, 0, 4, 2, 6);
}

FLINT_FORCE_INLINE vec4d vec4d_unpackhi(vec4d a, vec4d b) {
    return _MV_SHUF4(a, b, 1, 5, 3, 7);
}

#define DEFINE_IT(i0, i1, i2, i3) \
FLINT_FORCE_INLINE vec4d CAT6(vec4d, permute, i0, i1, i2, i3)(vec4d a) { \
    return _MV_SHUF4(a, a, i0, i1, i2, i3); \
}
DEFINE_IT(0,2,1,3)
DEFINE_IT(3,1,2,0)
DEFINE_IT(3,2,1,0)
#undef DEFINE_IT

#define DEFINE_IT(i0, i1) \
FLINT_FORCE_INLINE vec4d CAT4(vec4d, permute2, i0, i1)(vec4d a, vec4d b) { \
    return _MV_SHUF4(a, b, 2*(i0), 2*(i0) + 1, 2*(i1), 2*(i1) + 1); \
}
DEFINE_IT(0,2)
DEFINE_IT(1,3)
#undef DEFINE_IT

#if FLINT_BITS == 64

/* integer vectors -- GNU vector tier **************************************/

# define VEC_GENERIC_DEF_N(V, IV, N) \
FLINT_FORCE_INLINE V V##_load_unaligned(const ulong* a) { \
    V z; memcpy(&z, a, sizeof(z)); return z; \
} \
FLINT_FORCE_INLINE void V##_store_unaligned(ulong* z, V a) { \
    memcpy(z, &a, sizeof(a)); \
} \
FLINT_FORCE_INLINE V V##_zero(void) { V z = {0}; return z; } \
FLINT_FORCE_INLINE V V##_set_n(ulong a) { V z = {0}; return z + a; } \
FLINT_FORCE_INLINE V V##_add(V a, V b) { return a + b; } \
FLINT_FORCE_INLINE V V##_sub(V a, V b) { return a - b; } \
FLINT_FORCE_INLINE V V##_bit_and(V a, V b) { return a & b; } \
FLINT_FORCE_INLINE V V##_bit_shift_right(V a, ulong b) { return a >> b; } \
FLINT_FORCE_INLINE V V##_mul(V a, V b) { \
    return (a & UWORD(0xffffffff)) * (b & UWORD(0xffffffff)); \
} \
FLINT_FORCE_INLINE ulong V##_horizontal_sum(V a) { \
    ulong s = 0; \
    for (int i = 0; i < N; i++) \
        s += a[i]; \
    return s; \
} \
FLINT_FORCE_INLINE V V##_addmod(V a, V b, V n) { \
    V nmb = n - b; \
    V m = (V) (a >= nmb); \
    return a + b - (n & m); \
} \
FLINT_FORCE_INLINE V V##_addmod_limited(V a, V b, V n) { \
    V s = a + b; \
    V m = (V) (s >= n); \
    return s - (n & m); \
}

VEC_GENERIC_DEF_N(vec2n, _vec2i, 2)
VEC_GENERIC_DEF_N(vec4n, _vec4i, 4)

#undef VEC_GENERIC_DEF_N

#define vec2n_bit_shift_right_32(a) vec2n_bit_shift_right((a), 32)
#define vec4n_bit_shift_right_32(a) vec4n_bit_shift_right((a), 32)

FLINT_FORCE_INLINE vec4n vec4n_set_n4(ulong a0, ulong a1, ulong a2, ulong a3) {
    vec4n z = {a0, a1, a2, a3};
    return z;
}

FLINT_FORCE_INLINE vec4n vec4n_permute_3_2_1_0(vec4n a) {
    return _MV_SHUF4(a, a, 3, 2, 1, 0);
}

FLINT_FORCE_INLINE void vec4n_print(vec4n a) {
    flint_printf("[hi %016wx_%016wx_%016wx_%016wx lo]",
                 a[3], a[2], a[1], a[0]);
}

/* conversions and modular arithmetic -- GNU vector tier *******************/

# define VEC_GENERIC_DEF_MOD(V, N, IV) \
FLINT_FORCE_INLINE N V##_convert_limited_##N(V a) { \
    return (N) (a + 0x1.0p52) - UWORD(0x4330000000000000); \
} \
FLINT_FORCE_INLINE V N##_convert_limited_##V(N a) { \
    return (V) (a | UWORD(0x4330000000000000)) - 0x1.0p52; \
} \
/* round to nearest even and convert exactly, per-entry |a| < 2^52 */ \
FLINT_FORCE_INLINE N _##V##_round_i64(V a) { \
    V ax = V##_abs(a); \
    N t = (N) (ax + 0x1.0p52) - UWORD(0x4330000000000000); \
    IV s = (IV) a >> 63; \
    return (N) ((IV) t ^ s) - (N) s; \
} \
/* exact value modulo 2^64 of an integer-valued vector, per-entry |a| < 2^62 */ \
FLINT_FORCE_INLINE N _##V##_i64(V a) { \
    V th = a * 0x1p-32 + 0x1.8p52; \
    N ih = (N) th - UWORD(0x4338000000000000); \
    V lo = a - (th - 0x1.8p52) * 0x1p32; \
    N il = (N) (lo + 0x1.8p52) - UWORD(0x4338000000000000); \
    return (ih << 32) + il; \
} \
/* exact conversion of signed per-entry residues |r| < 2^51 to double */ \
FLINT_FORCE_INLINE V _##N##_small_i64_d(N r) { \
    return (V) (r + UWORD(0x4338000000000000)) - 0x1.8p52; \
} \
/* lane by lane, but the compiler turns these back into vector fmas \
   whenever the target has the instruction, which is the only case in \
   which this path is taken at all */ \
FLINT_FORCE_INLINE V _##V##_fma(V a, V b, V c) { \
    V z = a * b + c; \
    for (int i = 0; i < (int) (sizeof(V)/sizeof(double)); i++) \
        z[i] = fma(a[i], b[i], c[i]); \
    return z; \
} \
/* the integer detour and the fma form compute the same exact integer; \
   see _vec_generic_mulmod1 for why the choice is made on the target */ \
FLINT_FORCE_INLINE V V##_mulmod(V a, V b, V n, V ninv) { \
    V h = a * b; \
    if (_VEC_GENERIC_USE_FMA) \
    { \
        V q = V##_round(h * ninv); \
        return _##V##_fma(-q, n, h) + _##V##_fma(a, b, -h); \
    } \
    else \
    { \
        N iq = _##V##_round_i64(h * ninv); \
        N in = V##_convert_limited_##N(n); \
        N r = _##V##_i64(a) * _##V##_i64(b) - iq * in; \
        return _##N##_small_i64_d(r); \
    } \
} \
FLINT_FORCE_INLINE V V##_nmulmod(V a, V b, V n, V ninv) { \
    V h = a * b; \
    if (_VEC_GENERIC_USE_FMA) \
    { \
        V q = V##_round(h * ninv); \
        return _##V##_fma(q, n, -h) - _##V##_fma(a, b, -h); \
    } \
    else \
    { \
        N iq = _##V##_round_i64(h * ninv); \
        N in = V##_convert_limited_##N(n); \
        N r = iq * in - _##V##_i64(a) * _##V##_i64(b); \
        return _##N##_small_i64_d(r); \
    } \
} \
/* see _vec_generic_reduce_to_pm1n1 for why n is split in two */ \
FLINT_FORCE_INLINE V V##_reduce_to_pm1n(V a, V n, V ninv) { \
    V x = a * ninv; \
    V t = x + 0x1.8p52; \
    V q = t - 0x1.8p52; \
    V nh = ((n * 0x1p-26 + 0x1.8p52) - 0x1.8p52) * 0x1p26; \
    V nl = n - nh; \
    /* the shift trick returns +0.0 where rint returns -0.0; restoring \
       the sign is what makes a = -0.0 give +0.0 here, as on the FMA \
       backends, rather than -0.0, which the [0,n) reductions would \
       turn into n */ \
    q = (V) ((N) q | ((N) x & (N) _##V##_signmask())); \
    return (a - q * nh) - q * nl; \
} \
FLINT_FORCE_INLINE V V##_reduce_to_pm1no(V a, V n, V ninv) { \
    return V##_reduce_to_pm1n(a, n, ninv); \
} \
FLINT_FORCE_INLINE V V##_reduce_to_0n(V a, V n, V ninv) { \
    return V##_reduce_pm1no_to_0n(V##_reduce_to_pm1no(a, n, ninv), n); \
} \
FLINT_FORCE_INLINE int V##_same_mod(V a, V b, V n, V ninv) { \
    return V##_same(V##_reduce_to_0n(a, n, ninv), V##_reduce_to_0n(b, n, ninv)); \
}

VEC_GENERIC_DEF_MOD(vec2d, vec2n, _vec2i)
VEC_GENERIC_DEF_MOD(vec4d, vec4n, _vec4i)

#undef VEC_GENERIC_DEF_MOD

#endif /* FLINT_BITS == 64 */

#else

/*
    The test code instantiates this file more than once in a translation
    unit (under renamed identifiers), possibly with different tiers, so
    the tier indicator must be reset here rather than assumed undefined.
*/
#undef FLINT_MACHINE_VECTORS_GNU_VECTOR_EXTENSIONS

typedef struct {double v[2];} vec2d;
typedef struct {double v[4];} vec4d;
typedef struct {float v[4];} vec4f;
typedef struct {float v[8];} vec8f;
#if FLINT_BITS == 64
typedef struct {ulong v[2];} vec2n;
typedef struct {ulong v[4];} vec4n;
#endif

/*
    The elementwise operations are written as straight-line code rather
    than loops: whether they become SIMD instructions is up to the
    compiler's SLP vectorizer, and an unrolled basic block is
    considerably more reliable across compilers and versions than a loop
    the vectorizer must first unroll. This tier is only reached on
    compilers without GNU vector extensions, so the generated code is
    otherwise outside our control.
*/

# define VEC_SC_EACH2(F) F(0) F(1)
# define VEC_SC_EACH4(F) F(0) F(1) F(2) F(3)
# define VEC_SC_EACH8(F) F(0) F(1) F(2) F(3) F(4) F(5) F(6) F(7)

# define VEC_SC_LOAD(i) z.v[i] = a[i];
# define VEC_SC_STORE(i) z[i] = a.v[i];
# define VEC_SC_ZERO(i) z.v[i] = 0;
# define VEC_SC_SET1(i) z.v[i] = a;
# define VEC_SC_ADD(i) z.v[i] = a.v[i] + b.v[i];
# define VEC_SC_SUB(i) z.v[i] = a.v[i] - b.v[i];
# define VEC_SC_MUL(i) z.v[i] = a.v[i] * b.v[i];
# define VEC_SC_FMADD(i) z.v[i] = a.v[i] * b.v[i] + c.v[i];
# define VEC_SC_FNMADD(i) z.v[i] = c.v[i] - a.v[i] * b.v[i];

# define VEC_GENERIC_DEF(V, S, N, SUF, FLOOR1) \
FLINT_FORCE_INLINE V V##_load_unaligned(const S* a) { \
    V z; VEC_SC_EACH##N(VEC_SC_LOAD) return z; \
} \
FLINT_FORCE_INLINE V V##_load_aligned(const S* a) { \
    return V##_load_unaligned(a); \
} \
FLINT_FORCE_INLINE V V##_load(const S* a) { \
    return V##_load_unaligned(a); \
} \
FLINT_FORCE_INLINE void V##_store_unaligned(S* z, V a) { \
    VEC_SC_EACH##N(VEC_SC_STORE) \
} \
FLINT_FORCE_INLINE void V##_store_aligned(S* z, V a) { \
    V##_store_unaligned(z, a); \
} \
FLINT_FORCE_INLINE void V##_store(S* z, V a) { \
    V##_store_unaligned(z, a); \
} \
FLINT_FORCE_INLINE V V##_zero(void) { \
    V z; VEC_SC_EACH##N(VEC_SC_ZERO) return z; \
} \
FLINT_FORCE_INLINE V V##_set_##SUF(S a) { \
    V z; VEC_SC_EACH##N(VEC_SC_SET1) return z; \
} \
FLINT_FORCE_INLINE V V##_add(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_ADD) return z; \
} \
FLINT_FORCE_INLINE V V##_sub(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_SUB) return z; \
} \
FLINT_FORCE_INLINE V V##_mul(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_MUL) return z; \
} \
FLINT_FORCE_INLINE V V##_fmadd(V a, V b, V c) { \
    V z; VEC_SC_EACH##N(VEC_SC_FMADD) return z; \
} \
FLINT_FORCE_INLINE V V##_fnmadd(V a, V b, V c) { \
    V z; VEC_SC_EACH##N(VEC_SC_FNMADD) return z; \
} \
FLINT_FORCE_INLINE V V##_floor(V a) { \
    int i; \
    for (i = 0; i < N; i++) \
        a.v[i] = FLOOR1(a.v[i]); \
    return a; \
}

VEC_GENERIC_DEF(vec2d, double, 2, d, floor)
VEC_GENERIC_DEF(vec4d, double, 4, d, floor)
VEC_GENERIC_DEF(vec4f, float, 4, f, floorf)
VEC_GENERIC_DEF(vec8f, float, 8, f, floorf)

#undef VEC_GENERIC_DEF

/* the sign-bit select underlying blendv, on the scalar representation */
FLINT_FORCE_INLINE double _vec_generic_blendv1(double a, double b, double c)
{
    slong m;
    memcpy(&m, &c, sizeof(m));
    return m < 0 ? b : a;
}

# define VEC_SC_MSK(i) z.v[i] = a.v[i] OP b.v[i] ? -1.0 : 0.0;
# define VEC_SC_GET(i) if (i == n) return a.v[i];
# define VEC_SC_SAME(i) r = r && (a.v[i] == b.v[i]);
# define VEC_SC_NEG(i) a.v[i] = -a.v[i];
# define VEC_SC_ABS(i) a.v[i] = fabs(a.v[i]);
# define VEC_SC_HALF(i) a.v[i] = a.v[i] * 0.5;
# define VEC_SC_DIV(i) z.v[i] = a.v[i] / b.v[i];
# define VEC_SC_MIN(i) z.v[i] = a.v[i] < b.v[i] ? a.v[i] : b.v[i];
# define VEC_SC_MAX(i) z.v[i] = b.v[i] < a.v[i] ? a.v[i] : b.v[i];
# define VEC_SC_FMSUB(i) z.v[i] = a.v[i] * b.v[i] - c.v[i];
# define VEC_SC_FNMSUB(i) z.v[i] = -(a.v[i] * b.v[i]) - c.v[i];
# define VEC_SC_ADDSUB(i) z.v[i] = (i % 2) ? a.v[i] + b.v[i] : a.v[i] - b.v[i];
# define VEC_SC_BLENDV(i) z.v[i] = _vec_generic_blendv1(a.v[i], b.v[i], c.v[i]);
# define VEC_SC_ROUND(i) a.v[i] = rint(a.v[i]);
# define VEC_SC_R0N(i) z.v[i] = a.v[i] > 0.5 * n.v[i] ? a.v[i] - n.v[i] : a.v[i];
# define VEC_SC_R2N(i) z.v[i] = a.v[i] - n.v[i] >= 0 ? a.v[i] - n.v[i] : a.v[i];
# define VEC_SC_RPM1NO0N(i) z.v[i] = _vec_generic_blendv1(a.v[i], a.v[i] + n.v[i], a.v[i]);

FLINT_FORCE_INLINE double _vec_generic_reduce_pm1n_to_pmhn1(double a, double n)
{
    double t = a + n;
    double halfn = 0.5 * n;
    if (a > halfn)
        return a - n;
    else if (t < halfn)
        return t;
    else
        return a;
}

# define VEC_SC_RPMHN(i) z.v[i] = _vec_generic_reduce_pm1n_to_pmhn1(a.v[i], n.v[i]);

# define VEC_GENERIC_DEF_D(V, N) \
FLINT_FORCE_INLINE double V##_get_index(V a, const int n) { \
    VEC_SC_EACH##N(VEC_SC_GET) \
    return a.v[0]; \
} \
FLINT_FORCE_INLINE int V##_same(V a, V b) { \
    int r = 1; VEC_SC_EACH##N(VEC_SC_SAME) return r; \
} \
FLINT_FORCE_INLINE V V##_one(void) { return V##_set_d(1.0); } \
FLINT_FORCE_INLINE V V##_neg(V a) { VEC_SC_EACH##N(VEC_SC_NEG) return a; } \
FLINT_FORCE_INLINE V V##_abs(V a) { VEC_SC_EACH##N(VEC_SC_ABS) return a; } \
FLINT_FORCE_INLINE V V##_half(V a) { VEC_SC_EACH##N(VEC_SC_HALF) return a; } \
FLINT_FORCE_INLINE V V##_div(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_DIV) return z; \
} \
FLINT_FORCE_INLINE V V##_min(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_MIN) return z; \
} \
FLINT_FORCE_INLINE V V##_max(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_MAX) return z; \
} \
FLINT_FORCE_INLINE V V##_fmsub(V a, V b, V c) { \
    V z; VEC_SC_EACH##N(VEC_SC_FMSUB) return z; \
} \
FLINT_FORCE_INLINE V V##_fnmsub(V a, V b, V c) { \
    V z; VEC_SC_EACH##N(VEC_SC_FNMSUB) return z; \
} \
FLINT_FORCE_INLINE V V##_addsub(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_ADDSUB) return z; \
} \
FLINT_FORCE_INLINE V V##_blendv(V a, V b, V c) { \
    V z; VEC_SC_EACH##N(VEC_SC_BLENDV) return z; \
} \
FLINT_FORCE_INLINE V V##_round(V a) { VEC_SC_EACH##N(VEC_SC_ROUND) return a; } \
FLINT_FORCE_INLINE V V##_reduce_0n_to_pmhn(V a, V n) { \
    V z; VEC_SC_EACH##N(VEC_SC_R0N) return z; \
} \
FLINT_FORCE_INLINE V V##_reduce_pm1n_to_pmhn(V a, V n) { \
    V z; VEC_SC_EACH##N(VEC_SC_RPMHN) return z; \
} \
FLINT_FORCE_INLINE V V##_reduce_2n_to_n(V a, V n) { \
    V z; VEC_SC_EACH##N(VEC_SC_R2N) return z; \
} \
FLINT_FORCE_INLINE V V##_reduce_pm1no_to_0n(V a, V n) { \
    V z; VEC_SC_EACH##N(VEC_SC_RPM1NO0N) return z; \
}

# define VEC_SC_CMP_DEF(V, N, name, OP2) \
FLINT_FORCE_INLINE V V##_cmp_##name(V a, V b) { \
    V z; \
    for (int i = 0; i < N; i++) { \
        slong m = (a.v[i] OP2 b.v[i]) ? -1 : 0; \
        double d; \
        memcpy(&d, &m, sizeof(d)); \
        z.v[i] = d; \
    } \
    return z; \
}

VEC_GENERIC_DEF_D(vec2d, 2)
VEC_GENERIC_DEF_D(vec4d, 4)
VEC_SC_CMP_DEF(vec2d, 2, ge, >=)
VEC_SC_CMP_DEF(vec2d, 2, gt, >)
VEC_SC_CMP_DEF(vec2d, 2, lt, <)
VEC_SC_CMP_DEF(vec4d, 4, ge, >=)
VEC_SC_CMP_DEF(vec4d, 4, gt, >)
VEC_SC_CMP_DEF(vec4d, 4, lt, <)

#undef VEC_GENERIC_DEF_D
#undef VEC_SC_CMP_DEF

FLINT_FORCE_INLINE vec2d vec2d_set_d2(double a0, double a1) {
    vec2d z; z.v[0] = a0; z.v[1] = a1;
    return z;
}

FLINT_FORCE_INLINE vec4d vec4d_set_d4(double a0, double a1, double a2, double a3) {
    vec4d z; z.v[0] = a0; z.v[1] = a1; z.v[2] = a2; z.v[3] = a3;
    return z;
}

FLINT_FORCE_INLINE vec4d vec4d_set_vec2d2(vec2d a, vec2d b) {
    return vec4d_set_d4(a.v[0], a.v[1], b.v[0], b.v[1]);
}

FLINT_FORCE_INLINE void vec4d_print(vec4d a) {
    flint_printf("{%f, %f, %f, %f}", a.v[0], a.v[1], a.v[2], a.v[3]);
}

FLINT_FORCE_INLINE vec2d vec2d_unpacklo(vec2d a, vec2d b) {
    return vec2d_set_d2(a.v[0], b.v[0]);
}

FLINT_FORCE_INLINE vec2d vec2d_unpackhi(vec2d a, vec2d b) {
    return vec2d_set_d2(a.v[1], b.v[1]);
}

FLINT_FORCE_INLINE vec4d vec4d_unpacklo(vec4d a, vec4d b) {
    return vec4d_set_d4(a.v[0], b.v[0], a.v[2], b.v[2]);
}

FLINT_FORCE_INLINE vec4d vec4d_unpackhi(vec4d a, vec4d b) {
    return vec4d_set_d4(a.v[1], b.v[1], a.v[3], b.v[3]);
}

#define DEFINE_IT(i0, i1, i2, i3) \
FLINT_FORCE_INLINE vec4d CAT6(vec4d, permute, i0, i1, i2, i3)(vec4d a) { \
    return vec4d_set_d4(a.v[i0], a.v[i1], a.v[i2], a.v[i3]); \
}
DEFINE_IT(0,2,1,3)
DEFINE_IT(3,1,2,0)
DEFINE_IT(3,2,1,0)
#undef DEFINE_IT

FLINT_FORCE_INLINE vec4d vec4d_permute2_0_2(vec4d a, vec4d b) {
    return vec4d_set_d4(a.v[0], a.v[1], b.v[0], b.v[1]);
}

FLINT_FORCE_INLINE vec4d vec4d_permute2_1_3(vec4d a, vec4d b) {
    return vec4d_set_d4(a.v[2], a.v[3], b.v[2], b.v[3]);
}

#if FLINT_BITS == 64

/* integer vectors and modular arithmetic -- strict tier *******************/

# define VEC_SC_AND(i) z.v[i] = a.v[i] & b.v[i];
# define VEC_SC_SHR(i) a.v[i] = a.v[i] >> b;
# define VEC_SC_MUL32(i) z.v[i] = (a.v[i] & UWORD(0xffffffff)) * (b.v[i] & UWORD(0xffffffff));
# define VEC_SC_HSUM(i) s += a.v[i];
# define VEC_SC_ADDMOD(i) z.v[i] = vec1n_addmod(a.v[i], b.v[i], n.v[i]);
# define VEC_SC_ADDMODL(i) z.v[i] = a.v[i] + b.v[i] >= n.v[i] ? a.v[i] + b.v[i] - n.v[i] : a.v[i] + b.v[i];
# define VEC_SC_D2N(i) z.v[i] = (ulong)(slong) a.v[i];
# define VEC_SC_N2D(i) z.v[i] = (double) a.v[i];
# define VEC_SC_MULMOD(i) z.v[i] = _vec_generic_mulmod1(a.v[i], b.v[i], n.v[i], ninv.v[i]);
# define VEC_SC_NMULMOD(i) z.v[i] = _vec_generic_nmulmod1(a.v[i], b.v[i], n.v[i], ninv.v[i]);
# define VEC_SC_RPM1N(i) z.v[i] = _vec_generic_reduce_to_pm1n1(a.v[i], n.v[i], ninv.v[i]);

FLINT_FORCE_INLINE ulong vec1n_addmod(ulong a, ulong b, ulong n)
{
    ulong nmb = n - b;
    return nmb > a ? a + b : a - nmb;
}

# define VEC_GENERIC_DEF_N(V, D, N) \
FLINT_FORCE_INLINE V V##_load_unaligned(const ulong* a) { \
    V z; VEC_SC_EACH##N(VEC_SC_LOAD) return z; \
} \
FLINT_FORCE_INLINE void V##_store_unaligned(ulong* z, V a) { \
    VEC_SC_EACH##N(VEC_SC_STORE) \
} \
FLINT_FORCE_INLINE V V##_zero(void) { \
    V z; VEC_SC_EACH##N(VEC_SC_ZERO) return z; \
} \
FLINT_FORCE_INLINE V V##_set_n(ulong a) { \
    V z; VEC_SC_EACH##N(VEC_SC_SET1) return z; \
} \
FLINT_FORCE_INLINE V V##_add(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_ADD) return z; \
} \
FLINT_FORCE_INLINE V V##_sub(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_SUB) return z; \
} \
FLINT_FORCE_INLINE V V##_bit_and(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_AND) return z; \
} \
FLINT_FORCE_INLINE V V##_bit_shift_right(V a, ulong b) { \
    VEC_SC_EACH##N(VEC_SC_SHR) return a; \
} \
FLINT_FORCE_INLINE V V##_mul(V a, V b) { \
    V z; VEC_SC_EACH##N(VEC_SC_MUL32) return z; \
} \
FLINT_FORCE_INLINE ulong V##_horizontal_sum(V a) { \
    ulong s = 0; VEC_SC_EACH##N(VEC_SC_HSUM) return s; \
} \
FLINT_FORCE_INLINE V V##_addmod(V a, V b, V n) { \
    V z; VEC_SC_EACH##N(VEC_SC_ADDMOD) return z; \
} \
FLINT_FORCE_INLINE V V##_addmod_limited(V a, V b, V n) { \
    V z; VEC_SC_EACH##N(VEC_SC_ADDMODL) return z; \
} \
FLINT_FORCE_INLINE V D##_convert_limited_##V(D a) { \
    V z; VEC_SC_EACH##N(VEC_SC_D2N) return z; \
} \
FLINT_FORCE_INLINE D V##_convert_limited_##D(V a) { \
    D z; VEC_SC_EACH##N(VEC_SC_N2D) return z; \
} \
FLINT_FORCE_INLINE D D##_mulmod(D a, D b, D n, D ninv) { \
    D z; VEC_SC_EACH##N(VEC_SC_MULMOD) return z; \
} \
FLINT_FORCE_INLINE D D##_nmulmod(D a, D b, D n, D ninv) { \
    D z; VEC_SC_EACH##N(VEC_SC_NMULMOD) return z; \
} \
FLINT_FORCE_INLINE D D##_reduce_to_pm1n(D a, D n, D ninv) { \
    D z; VEC_SC_EACH##N(VEC_SC_RPM1N) return z; \
} \
FLINT_FORCE_INLINE D D##_reduce_to_pm1no(D a, D n, D ninv) { \
    return D##_reduce_to_pm1n(a, n, ninv); \
} \
FLINT_FORCE_INLINE D D##_reduce_to_0n(D a, D n, D ninv) { \
    return D##_reduce_pm1no_to_0n(D##_reduce_to_pm1no(a, n, ninv), n); \
} \
FLINT_FORCE_INLINE int D##_same_mod(D a, D b, D n, D ninv) { \
    return D##_same(D##_reduce_to_0n(a, n, ninv), D##_reduce_to_0n(b, n, ninv)); \
}

VEC_GENERIC_DEF_N(vec2n, vec2d, 2)
VEC_GENERIC_DEF_N(vec4n, vec4d, 4)

#undef VEC_GENERIC_DEF_N

#define vec2n_bit_shift_right_32(a) vec2n_bit_shift_right((a), 32)
#define vec4n_bit_shift_right_32(a) vec4n_bit_shift_right((a), 32)

FLINT_FORCE_INLINE vec4n vec4n_set_n4(ulong a0, ulong a1, ulong a2, ulong a3) {
    vec4n z; z.v[0] = a0; z.v[1] = a1; z.v[2] = a2; z.v[3] = a3;
    return z;
}

FLINT_FORCE_INLINE vec4n vec4n_permute_3_2_1_0(vec4n a) {
    return vec4n_set_n4(a.v[3], a.v[2], a.v[1], a.v[0]);
}

FLINT_FORCE_INLINE void vec4n_print(vec4n a) {
    flint_printf("[hi %016wx_%016wx_%016wx_%016wx lo]",
                 a.v[3], a.v[2], a.v[1], a.v[0]);
}

#endif /* FLINT_BITS == 64 */

#endif

/* combined operations shared by both tiers *********************************/

FLINT_FORCE_INLINE vec4d vec4d_unpack_lo_permute_0_2_1_3(vec4d u, vec4d v) {
    return vec4d_permute_0_2_1_3(vec4d_unpacklo(u, v));
}

FLINT_FORCE_INLINE vec4d vec4d_unpack_hi_permute_0_2_1_3(vec4d u, vec4d v) {
    return vec4d_permute_0_2_1_3(vec4d_unpackhi(u, v));
}

FLINT_FORCE_INLINE vec4d vec4d_unpackhi_permute_3_1_2_0(vec4d u, vec4d v) {
    return vec4d_permute_3_1_2_0(vec4d_unpackhi(u, v));
}

FLINT_FORCE_INLINE vec4d vec4d_unpacklo_permute_3_1_2_0(vec4d u, vec4d v) {
    return vec4d_permute_3_1_2_0(vec4d_unpacklo(u, v));
}

/* view the 4 vectors as the rows of a 4x4 matrix */
#define VEC4D_TRANSPOSE(z0, z1, z2, z3, a0, a1, a2, a3) \
{ \
    vec4d _t0, _t1, _t2, _t3; \
    _t0 = vec4d_unpacklo(a0, a1); \
    _t1 = vec4d_unpackhi(a0, a1); \
    _t2 = vec4d_unpacklo(a2, a3); \
    _t3 = vec4d_unpackhi(a2, a3); \
    z0 = vec4d_permute2_0_2(_t0, _t2); \
    z1 = vec4d_permute2_0_2(_t1, _t3); \
    z2 = vec4d_permute2_1_3(_t0, _t2); \
    z3 = vec4d_permute2_1_3(_t1, _t3); \
}

/* vec1d, vec1f -- generic *************************************************/

FLINT_FORCE_INLINE vec1d vec1d_load(const double* a) { return a[0]; }
FLINT_FORCE_INLINE vec1d vec1d_load_aligned(const double* a) { return a[0]; }
FLINT_FORCE_INLINE vec1d vec1d_load_unaligned(const double* a) { return a[0]; }
FLINT_FORCE_INLINE void vec1d_store(double* z, vec1d a) { z[0] = a; }
FLINT_FORCE_INLINE void vec1d_store_aligned(double* z, vec1d a) { z[0] = a; }
FLINT_FORCE_INLINE void vec1d_store_unaligned(double* z, vec1d a) { z[0] = a; }
FLINT_FORCE_INLINE int vec1d_same(double a, double b) { return a == b; }
FLINT_FORCE_INLINE vec1d vec1d_zero(void) { return 0.0; }
FLINT_FORCE_INLINE vec1d vec1d_one(void) { return 1.0; }
FLINT_FORCE_INLINE vec1d vec1d_set_d(double a) { return a; }
FLINT_FORCE_INLINE vec1d vec1d_round(vec1d a) { return rint(a); }
FLINT_FORCE_INLINE vec1d vec1d_add(vec1d a, vec1d b) { return a + b; }
FLINT_FORCE_INLINE vec1d vec1d_sub(vec1d a, vec1d b) { return a - b; }
FLINT_FORCE_INLINE vec1d vec1d_addsub(vec1d a, vec1d b) { return a - b; }
FLINT_FORCE_INLINE vec1d vec1d_neg(vec1d a) { return -a; }
FLINT_FORCE_INLINE vec1d vec1d_abs(vec1d a) { return fabs(a); }
FLINT_FORCE_INLINE vec1d vec1d_max(vec1d a, vec1d b) { return fmax(a, b); }
FLINT_FORCE_INLINE vec1d vec1d_min(vec1d a, vec1d b) { return fmin(a, b); }
FLINT_FORCE_INLINE vec1d vec1d_mul(vec1d a, vec1d b) { return a * b; }
FLINT_FORCE_INLINE vec1d vec1d_half(vec1d a) { return a * 0.5; }
FLINT_FORCE_INLINE vec1d vec1d_div(vec1d a, vec1d b) { return a / b; }
FLINT_FORCE_INLINE vec1d vec1d_fmadd(vec1d a, vec1d b, vec1d c) {
    return fma(a, b, c);
}
FLINT_FORCE_INLINE vec1d vec1d_fmsub(vec1d a, vec1d b, vec1d c) {
    return fma(a, b, -c);
}
FLINT_FORCE_INLINE vec1d vec1d_fnmadd(vec1d a, vec1d b, vec1d c) {
    return fma(-a, b, c);
}
FLINT_FORCE_INLINE vec1d vec1d_fnmsub(vec1d a, vec1d b, vec1d c) {
    return fma(-a, b, -c);
}
FLINT_FORCE_INLINE vec1d vec1d_floor(vec1d a) { return floor(a); }
FLINT_FORCE_INLINE vec1d vec1d_blendv(vec1d a, vec1d b, vec1d c) {
    return c >= 0 ? a : b;
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_0n_to_pmhn(vec1d a, vec1d n) {
    vec1d halfn = 0.5*n;
    return a > halfn ? a - n : a;
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_pm1n_to_pmhn(vec1d a, vec1d n) {
    vec1d t = a + n;
    vec1d halfn = 0.5*n;
    if (a > halfn)
        return a - n;
    else if (t < halfn)
        return t;
    else
        return a;
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_2n_to_n(vec1d a, vec1d n) {
    return a - n >= 0 ? a - n : a;
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_pm1no_to_0n(vec1d a, vec1d n) {
    return vec1d_blendv(a, vec1d_add(a, n), a);
}

#if FLINT_BITS == 64

FLINT_FORCE_INLINE vec1d vec1d_reduce_to_pm1n(vec1d a, vec1d n, vec1d ninv) {
    return _vec_generic_reduce_to_pm1n1(a, n, ninv);
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_to_pm1no(vec1d a, vec1d n, vec1d ninv) {
    return _vec_generic_reduce_to_pm1n1(a, n, ninv);
}

FLINT_FORCE_INLINE vec1d vec1d_reduce_to_0n(vec1d a, vec1d n, vec1d ninv) {
    return vec1d_reduce_pm1no_to_0n(vec1d_reduce_to_pm1no(a, n, ninv), n);
}

FLINT_FORCE_INLINE int vec1d_same_mod(vec1d a, vec1d b, vec1d n, vec1d ninv) {
    return vec1d_same(vec1d_reduce_to_0n(a, n, ninv), vec1d_reduce_to_0n(b, n, ninv));
}

FLINT_FORCE_INLINE vec1d vec1d_mulmod(vec1d a, vec1d b, vec1d n, vec1d ninv) {
    return _vec_generic_mulmod1(a, b, n, ninv);
}

FLINT_FORCE_INLINE vec1d vec1d_nmulmod(vec1d a, vec1d b, vec1d n, vec1d ninv) {
    return _vec_generic_nmulmod1(a, b, n, ninv);
}

FLINT_FORCE_INLINE vec1n vec1d_convert_limited_vec1n(vec1d a) {
    return (ulong)(slong) a;
}

FLINT_FORCE_INLINE vec1d vec1n_convert_limited_vec1d(vec1n a) {
    return (double) a;
}

FLINT_FORCE_INLINE void vec1n_store_unaligned(ulong* z, vec1n a) {
    z[0] = a;
}

#if defined(FLINT_MACHINE_VECTORS_GNU_VECTOR_EXTENSIONS)
FLINT_FORCE_INLINE ulong vec1n_addmod(ulong a, ulong b, ulong n)
{
    ulong nmb = n - b;
    return nmb > a ? a + b : a - nmb;
}
#endif

#endif /* FLINT_BITS == 64 */

FLINT_FORCE_INLINE vec1f vec1f_load(const float* a) { return a[0]; }
FLINT_FORCE_INLINE vec1f vec1f_load_aligned(const float* a) { return a[0]; }
FLINT_FORCE_INLINE vec1f vec1f_load_unaligned(const float* a) { return a[0]; }
FLINT_FORCE_INLINE void vec1f_store(float* z, vec1f a) { z[0] = a; }
FLINT_FORCE_INLINE void vec1f_store_aligned(float* z, vec1f a) { z[0] = a; }
FLINT_FORCE_INLINE void vec1f_store_unaligned(float* z, vec1f a) { z[0] = a; }
FLINT_FORCE_INLINE vec1f vec1f_zero(void) { return 0.0f; }
FLINT_FORCE_INLINE vec1f vec1f_set_f(float a) { return a; }
FLINT_FORCE_INLINE vec1f vec1f_add(vec1f a, vec1f b) { return a + b; }
FLINT_FORCE_INLINE vec1f vec1f_sub(vec1f a, vec1f b) { return a - b; }
FLINT_FORCE_INLINE vec1f vec1f_mul(vec1f a, vec1f b) { return a * b; }
FLINT_FORCE_INLINE vec1f vec1f_fmadd(vec1f a, vec1f b, vec1f c) {
    return fmaf(a, b, c);
}
FLINT_FORCE_INLINE vec1f vec1f_fnmadd(vec1f a, vec1f b, vec1f c) {
    return fmaf(-a, b, c);
}
FLINT_FORCE_INLINE vec1f vec1f_floor(vec1f a) { return floorf(a); }

/* vec8d, vec16f, vec8n -- generic (pairs) *********************************/

typedef struct {vec4d e1, e2;} vec8d;
typedef struct {vec8f e1, e2;} vec16f;
#if FLINT_BITS == 64
typedef struct {vec4n e1, e2;} vec8n;
#endif

#define VEC_GENERIC_PAIR_DEF(U, V, S, SUF) \
FLINT_FORCE_INLINE V V##_load_unaligned(const S* a) { \
    V z = {U##_load_unaligned(a), U##_load_unaligned(a + sizeof(U)/sizeof(S))}; \
    return z; \
} \
FLINT_FORCE_INLINE V V##_load_aligned(const S* a) { \
    V z = {U##_load_aligned(a), U##_load_aligned(a + sizeof(U)/sizeof(S))}; \
    return z; \
} \
FLINT_FORCE_INLINE V V##_load(const S* a) { return V##_load_aligned(a); } \
FLINT_FORCE_INLINE void V##_store_unaligned(S* z, V a) { \
    U##_store_unaligned(z, a.e1); \
    U##_store_unaligned(z + sizeof(U)/sizeof(S), a.e2); \
} \
FLINT_FORCE_INLINE void V##_store_aligned(S* z, V a) { \
    U##_store_aligned(z, a.e1); \
    U##_store_aligned(z + sizeof(U)/sizeof(S), a.e2); \
} \
FLINT_FORCE_INLINE void V##_store(S* z, V a) { V##_store_aligned(z, a); } \
FLINT_FORCE_INLINE V V##_zero(void) { \
    V z = {U##_zero(), U##_zero()}; return z; \
} \
FLINT_FORCE_INLINE V V##_set_##SUF(S a) { \
    V z = {U##_set_##SUF(a), U##_set_##SUF(a)}; return z; \
} \
FLINT_FORCE_INLINE V V##_add(V a, V b) { \
    V z = {U##_add(a.e1, b.e1), U##_add(a.e2, b.e2)}; return z; \
} \
FLINT_FORCE_INLINE V V##_sub(V a, V b) { \
    V z = {U##_sub(a.e1, b.e1), U##_sub(a.e2, b.e2)}; return z; \
} \
FLINT_FORCE_INLINE V V##_mul(V a, V b) { \
    V z = {U##_mul(a.e1, b.e1), U##_mul(a.e2, b.e2)}; return z; \
} \
FLINT_FORCE_INLINE V V##_fmadd(V a, V b, V c) { \
    V z = {U##_fmadd(a.e1, b.e1, c.e1), U##_fmadd(a.e2, b.e2, c.e2)}; \
    return z; \
} \
FLINT_FORCE_INLINE V V##_fnmadd(V a, V b, V c) { \
    V z = {U##_fnmadd(a.e1, b.e1, c.e1), U##_fnmadd(a.e2, b.e2, c.e2)}; \
    return z; \
} \
FLINT_FORCE_INLINE V V##_floor(V a) { \
    V z = {U##_floor(a.e1), U##_floor(a.e2)}; return z; \
}

VEC_GENERIC_PAIR_DEF(vec4d, vec8d, double, d)
VEC_GENERIC_PAIR_DEF(vec8f, vec16f, float, f)

#undef VEC_GENERIC_PAIR_DEF

#define EXTEND_VEC_DEF1(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a) { \
    U z1 = U##f(a.e1); \
    U z2 = U##f(a.e2); \
    V z = {z1, z2}; return z; \
}

#define EXTEND_VEC_DEF2(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b) { \
    U z1 = U##f(a.e1, b.e1); \
    U z2 = U##f(a.e2, b.e2); \
    V z = {z1, z2}; return z; \
}

#define EXTEND_VEC_DEF3(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b, V c) { \
    U z1 = U##f(a.e1, b.e1, c.e1); \
    U z2 = U##f(a.e2, b.e2, c.e2); \
    V z = {z1, z2}; return z; \
}

#define EXTEND_VEC_DEF4(U, V, f) \
FLINT_FORCE_INLINE V V##f(V a, V b, V c, V d) { \
    U z1 = U##f(a.e1, b.e1, c.e1, d.e1); \
    U z2 = U##f(a.e2, b.e2, c.e2, d.e2); \
    V z = {z1, z2}; return z; \
}

EXTEND_VEC_DEF1(vec4d, vec8d, _neg)
EXTEND_VEC_DEF1(vec4d, vec8d, _abs)
EXTEND_VEC_DEF1(vec4d, vec8d, _half)
EXTEND_VEC_DEF1(vec4d, vec8d, _round)
EXTEND_VEC_DEF2(vec4d, vec8d, _min)
EXTEND_VEC_DEF2(vec4d, vec8d, _max)
EXTEND_VEC_DEF2(vec4d, vec8d, _div)
EXTEND_VEC_DEF2(vec4d, vec8d, _reduce_0n_to_pmhn)
EXTEND_VEC_DEF2(vec4d, vec8d, _reduce_pm1n_to_pmhn)
EXTEND_VEC_DEF2(vec4d, vec8d, _reduce_pm1no_to_0n)
EXTEND_VEC_DEF2(vec4d, vec8d, _reduce_2n_to_n)
EXTEND_VEC_DEF3(vec4d, vec8d, _fmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _fnmsub)
EXTEND_VEC_DEF3(vec4d, vec8d, _blendv)

FLINT_FORCE_INLINE vec8d vec8d_set_vec4d2(vec4d a, vec4d b) {
    vec8d z = {a, b};
    return z;
}

EXTEND_VEC_DEF2(vec4d, vec8d, _unpacklo)
EXTEND_VEC_DEF2(vec4d, vec8d, _unpackhi)

FLINT_FORCE_INLINE vec8d vec8d_one(void) {
    vec8d z = {vec4d_one(), vec4d_one()};
    return z;
}

FLINT_FORCE_INLINE vec8d vec8d_set_d8(double a0, double a1, double a2,
                double a3, double a4, double a5, double a6, double a7) {
    vec8d z = {vec4d_set_d4(a0, a1, a2, a3), vec4d_set_d4(a4, a5, a6, a7)};
    return z;
}

FLINT_FORCE_INLINE double vec8d_get_index(vec8d a, int i) {
    return i < 4 ? vec4d_get_index(a.e1, i) : vec4d_get_index(a.e2, i - 4);
}

FLINT_FORCE_INLINE int vec8d_same(vec8d a, vec8d b) {
    return vec4d_same(a.e1, b.e1) && vec4d_same(a.e2, b.e2);
}

#if FLINT_BITS == 64

EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1n)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_pm1no)
EXTEND_VEC_DEF3(vec4d, vec8d, _reduce_to_0n)
EXTEND_VEC_DEF4(vec4d, vec8d, _mulmod)
EXTEND_VEC_DEF4(vec4d, vec8d, _nmulmod)

FLINT_FORCE_INLINE int vec8d_same_mod(vec8d a, vec8d b, vec8d n, vec8d ninv) {
    return vec4d_same_mod(a.e1, b.e1, n.e1, ninv.e1)
        && vec4d_same_mod(a.e2, b.e2, n.e2, ninv.e2);
}

FLINT_FORCE_INLINE vec8n vec8n_load_unaligned(const ulong* a) {
    vec8n z = {vec4n_load_unaligned(a + 0), vec4n_load_unaligned(a + 4)};
    return z;
}

FLINT_FORCE_INLINE void vec8n_store_unaligned(ulong* z, vec8n a) {
    vec4n_store_unaligned(z + 0, a.e1);
    vec4n_store_unaligned(z + 4, a.e2);
}

FLINT_FORCE_INLINE vec8n vec8n_zero(void) {
    vec8n z = {vec4n_zero(), vec4n_zero()};
    return z;
}

FLINT_FORCE_INLINE vec8n vec8n_set_n(ulong a) {
    vec4n x = vec4n_set_n(a);
    vec8n z = {x, x};
    return z;
}

FLINT_FORCE_INLINE vec8n vec8n_bit_shift_right(vec8n a, ulong b) {
    vec8n z = {vec4n_bit_shift_right(a.e1, b), vec4n_bit_shift_right(a.e2, b)};
    return z;
}

#define vec8n_bit_shift_right_32(a) vec8n_bit_shift_right((a), 32)

EXTEND_VEC_DEF2(vec4n, vec8n, _add)
EXTEND_VEC_DEF2(vec4n, vec8n, _sub)
EXTEND_VEC_DEF2(vec4n, vec8n, _bit_and)
EXTEND_VEC_DEF2(vec4n, vec8n, _mul)
EXTEND_VEC_DEF3(vec4n, vec8n, _addmod)
EXTEND_VEC_DEF3(vec4n, vec8n, _addmod_limited)

FLINT_FORCE_INLINE ulong vec8n_horizontal_sum(vec8n a) {
    return vec4n_horizontal_sum(a.e1) + vec4n_horizontal_sum(a.e2);
}

FLINT_FORCE_INLINE vec8d vec8n_convert_limited_vec8d(vec8n a) {
    vec8d z = {vec4n_convert_limited_vec4d(a.e1),
               vec4n_convert_limited_vec4d(a.e2)};
    return z;
}

FLINT_FORCE_INLINE vec8n vec8d_convert_limited_vec8n(vec8d a) {
    vec8n z = {vec4d_convert_limited_vec4n(a.e1),
               vec4d_convert_limited_vec4n(a.e2)};
    return z;
}

#endif /* FLINT_BITS == 64 */

#undef EXTEND_VEC_DEF4
#undef EXTEND_VEC_DEF3
#undef EXTEND_VEC_DEF2
#undef EXTEND_VEC_DEF1

#endif /* MACHINE_VECTORS_GENERIC_H */
