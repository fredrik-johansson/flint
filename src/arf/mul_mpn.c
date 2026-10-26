/*
    Copyright (C) 2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arf.h"

/*
    Multiplication of finite, nonzero arf values using truncated products.

    Let X and Y be the mantissas (xn and yn limbs, top limbs normalised).
    Rather than forming the full (xn + yn)-limb product, we compute an
    approximation A of the top an limbs of X Y, with

        0 <= P - A < 2^e            (P = X Y / B^(xn + yn - an), B = 2^64)

    for a known e, and then decide whether rounding A to prec bits gives the
    same result as rounding P (Ahlback & Johansson, "Fast basecases for
    arbitrary-size multiplication", Corollary 4). Writing k for the number of
    bits of A below the rounding point, the rounding of P (in any strict
    rounding mode, including the tie and inexact decisions) is determined by
    A as soon as the bits e, ..., k-2 of A are neither all zero nor all one:
    the first condition gives P mod 2^k >= A mod 2^k >= 2^e > 0 (P is
    inexact, and not a tie) and the second gives A mod 2^(k-1) < 2^(k-1) - 2^e,
    so that P cannot cross the rounding boundary 2^(k-1) (round to nearest)
    nor the truncation boundary 2^k. It suffices to inspect any subwindow of
    these bits, so we look at (at most) one limb's worth; the probability of
    an inconclusive check is then about 2^(1-w) for a window of w bits, in
    which case we fall back to the exact product.

    The approximation is the general high product flint_mpn_mulhigh_normalised
    (flint_mpn_sqrhigh_normalised for squares): the top outn = ceil(prec / 64)
    limbs of the product of the operands of any lengths, plus a guard limb,
    normalised so that the top bit is set, with a deficit of at most
    FLINT_MPN_MULHIGH_NORMALISED_GUARD_ERR_BOUND(outn) = 4 outn + 9 ulp of the
    guard limb (it truncates long operands itself, accounting for the
    dropped limbs). With at least 63 bits of A below the rounding point and
    e = bitcount(4 outn + 9) the check fails with probability about
    2^(e - 63).

    Full products are used instead where the high product would compute
    one anyway (operands no longer than the output and either small or far
    from balanced) or where nearly the whole product is wanted.

    With the relaxed rounding modes ARF_RND_FAST and ARF_RND_ACCURATE the
    check is skipped and A is rounded to nearest: since P - A < 2^e is at
    most 2^(e-k) ulp of the result with e - k <= -33, the result is within
    0.5 + 2^-33 ulp of P (when A rounds up to the next power of two, P lies
    within 2^(e-k) ulp above it as well).
*/

/* Bits of a guard window below which the check is not worth attempting. */
#define ARF_MUL_MIN_CHECK_BITS 4

/* Is rounding {a, an} (a lower approximation with 0 <= P - A < 2^e, in
   units of a[0]) to prec bits guaranteed to give the same result as
   rounding P, in every strict rounding mode? */
int
_arf_mul_rounding_certain(nn_srcptr a, slong an, slong e, slong prec)
{
    slong k, w, q, r;
    ulong v, mask;

    /* number of bits of A below the rounding point */
    k = an * FLINT_BITS - flint_clz(a[an - 1]) - prec;
    /* bits [e, k-1) are available; inspect at most one limb of them */
    w = k - 1 - e;

    if (w < ARF_MUL_MIN_CHECK_BITS)
        return 0;

    w = FLINT_MIN(w, FLINT_BITS);

    q = (ulong) e / FLINT_BITS;
    r = (ulong) e % FLINT_BITS;

    v = a[q] >> r;
    if (r != 0 && r + w > FLINT_BITS)
        v |= a[q + 1] << (FLINT_BITS - r);

    mask = (w == FLINT_BITS) ? ~UWORD(0) : ((UWORD(1) << w) - 1);
    v &= mask;

    return (v != 0) && (v != mask);
}


/*
    Lower approximation A of the top limbs of the product of {xptr, xn} and
    {yptr, yn} (top limbs normalised), for use where only the top outn
    limbs of the product matter (outn >= 0). A is written to dst, which must
    have room for ARF_MUL_APPROX_ALLOC(outn, xn, yn) limbs; the number of
    limbs an of A is returned, with 0 <= P - A < 2^e in units of dst[0] for
    the exact top-aligned product P, and the value is (A / B^an) 2^(xexp +
    yexp + expfix). Returns 0 when the full product should be computed
    instead (it is exact and no more expensive). Same method selection
    as _arf_mul_mpn.
*/
slong
_arf_mul_approx(nn_ptr dst, slong * e, slong * expfix,
    nn_srcptr xptr, slong xn, nn_srcptr yptr, slong yn, slong outn, int squaring)
{
    slong n = FLINT_MAX(outn, 1);
    mp_limb_pair_t hi;

    if (_flint_mpn_mulhigh_prefer_full(xn, yn, n))
        return 0;

    if (squaring)
        hi = flint_mpn_sqrhigh_normalised(dst + 1, xptr, xn, n);
    else
        hi = flint_mpn_mulhigh_normalised(dst + 1, xptr, xn, yptr, yn, n);

    dst[0] = hi.m1;
    *expfix = -(slong) hi.m2;
    *e = FLINT_BIT_COUNT(FLINT_MPN_MULHIGH_NORMALISED_GUARD_ERR_BOUND(n));

    return n + 1;
}


/*
    Specialised multiplication for prec <= n * FLINT_BITS with n <= 16: the
    normalised high product of n limbs with guard limb (the operands of any
    lengths, truncated internally), the rounding certified by inspecting
    the guard limb, and the rounding and normalisation done inline. Returns
    -1 (without touching z) if the check is inconclusive, else the inexact
    flag (always 1: a certified result is never exact).
*/
int
_arf_mul_high_small(arf_ptr z, arf_srcptr x, arf_srcptr y, slong n,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    ulong res[ARF_MUL_HIGH_SMALL_MAX], g, v;
    nn_srcptr xptr, yptr;
    slong xn, yn, fix;
    mp_limb_pair_t hi;

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    /* the balanced kernels directly (the general functions reduce to them
       after a dozen checks, which matter at this size) */
    if (x == y)
        hi = (xn == n) ? flint_mpn_sqrhigh_n_normalised(res, xptr, n)
                       : flint_mpn_sqrhigh_normalised(res, xptr, xn, n);
    else if (xn == n && yn == n)
        hi = flint_mpn_mulhigh_n_normalised(res, xptr, yptr, n);
    else
        hi = flint_mpn_mulhigh_normalised(res, xptr, xn, yptr, yn, n);

    g = hi.m1;

    /* Certification (see the top of this file): the deficit is below
       2^ARF_MUL_HIGH_SMALL_ERR_BITS units of g, and at least FLINT_BITS
       bits of the value lie below the rounding point, so the bits
       ARF_MUL_HIGH_SMALL_ERR_BITS .. FLINT_BITS - 2 of g form a valid
       window. */
    v = (g >> ARF_MUL_HIGH_SMALL_ERR_BITS) & ((UWORD(1) << (FLINT_BITS - 1 - ARF_MUL_HIGH_SMALL_ERR_BITS)) - 1);
    if (v == 0 || v == ((UWORD(1) << (FLINT_BITS - 1 - ARF_MUL_HIGH_SMALL_ERR_BITS)) - 1))
        return -1;

    fix = -(slong) hi.m2 + _arf_set_round_small_certified(z, res, n, g, sgnbit, prec, rnd);
    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix);

    return 1;
}


/*
    Rounds the normalised n-limb mantissa res (top bit set, n <= 16) to
    prec <= n * FLINT_BITS bits, where g is the limb below res, the value
    being known to be inexact and not a tie (certified by the caller), and
    writes the mantissa of z. Returns the exponent adjustment (1 if the
    result overflowed to a power of two, else 0). The caller sets the
    exponent.
*/
int
_arf_set_round_small_certified(arf_ptr z, nn_ptr res, slong n, ulong g,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    ulong mask, cy;
    nn_ptr zptr;
    slong zn, sb, i, fix = 0;
    int increment;

    /* Rounding of the n-limb mantissa res to prec bits: keep zn limbs,
       clearing sb bits at the bottom of res[n - zn]. */
    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    sb = zn * FLINT_BITS - prec;

    if (rnd == ARF_RND_DOWN)
        increment = 0;
    else if (rnd == ARF_RND_UP)
        increment = 1;
    else if (rnd == ARF_RND_FLOOR)
        increment = sgnbit;
    else if (rnd == ARF_RND_CEIL)
        increment = !sgnbit;
    else
    {
        /* nearest: the highest discarded bit (no tie is possible) */
        if (sb != 0)
            increment = (res[n - zn] >> (sb - 1)) & 1;
        else if (zn < n)
            increment = res[n - zn - 1] >> (FLINT_BITS - 1);
        else
            increment = g >> (FLINT_BITS - 1);
    }

    mask = (sb == 0) ? ~UWORD(0) : (~UWORD(0) << sb);
    res[n - zn] &= mask;

    if (increment)
    {
        cy = (sb == 0) ? 1 : (UWORD(1) << sb);
        for (i = n - zn; i < n; i++)
        {
            res[i] += cy;
            cy = (res[i] < cy);
            if (!cy)
                break;
        }
        if (cy)
        {
            /* overflow to a power of two */
            res[n - 1] = UWORD(1) << (FLINT_BITS - 1);
            _flint_mpn_zero_small(res + n - zn, zn - 1);
            fix++;
        }
    }

    /* strip trailing zero limbs */
    while (res[n - zn] == 0)
        zn--;

    ARF_GET_MPN_WRITE(zptr, zn, z);
    _flint_mpn_copy_small(zptr, res + n - zn, zn);
    ARF_XSIZE(z) = ARF_MAKE_XSIZE(zn, sgnbit);

    return fix;
}

/* Exact multiplication: full product, then rounding. */
static int
_arf_mul_exact(arf_ptr z, arf_srcptr x, arf_srcptr y, slong xn, slong yn,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xptr, yptr;
    nn_ptr tmp;
    slong zn, alloc, fix;
    int inexact;
    ARF_MUL_TMP_DECL

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    alloc = zn = xn + yn;
    ARF_MUL_TMP_ALLOC(tmp, alloc)

    FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);

    inexact = _arf_set_round_mpn(z, &fix, tmp, zn, sgnbit, prec, rnd);
    _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix);
    ARF_MUL_TMP_FREE(tmp, alloc)

    return inexact;
}

int
_arf_mul_mpn(arf_ptr z, arf_srcptr x, arf_srcptr y, slong xn, slong yn,
    int sgnbit, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xptr, yptr;
    nn_ptr dst, tmp = NULL;
    slong zn, n, an, alloc, fix, e;
    int inexact, relaxed, in_place;
    mp_limb_pair_t hi;
    arf_rnd_t srnd;
    ARF_MUL_TMP_DECL

    FLINT_ASSERT(xn >= yn && yn >= 1);

    /* strict mode used for exact products with a relaxed mode */
    relaxed = (rnd >= ARF_RND_FAST);
    srnd = relaxed ? arf_rnd_relaxed_to_strict(rnd) : rnd;

    zn = xn + yn;

    if (prec >= zn * FLINT_BITS)
        return _arf_mul_exact(z, x, y, xn, yn, sgnbit, prec, srnd);

    n = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (_flint_mpn_mulhigh_prefer_full(xn, yn, n))
        return _arf_mul_exact(z, x, y, xn, yn, sgnbit, prec, srnd);

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    /* The normalised high product (n limbs plus a guard limb) is written
       directly to the output when it is heap allocated and does not hold
       an operand (it is then rounded in place with a plain copy, the top
       bit being set). Besides z == x or z == y, the buffer of z can be
       shared with an operand that is a shallow copy of z (as in
       arf_neg_mul), so we compare the data pointers rather than the
       structs. */
    an = n + 1;
    in_place = (an > ARF_NOPTR_LIMBS) && (z != x) && (z != y) &&
        !(ARF_HAS_PTR(z) && (ARF_PTR_D(z) == xptr || ARF_PTR_D(z) == yptr));

    if (in_place)
    {
        ARF_GET_MPN_WRITE(dst, an, z);
        alloc = 0;
    }
    else
    {
        alloc = an;
        ARF_MUL_TMP_ALLOC(tmp, alloc)
        dst = tmp;
    }

    if (x == y)
        hi = flint_mpn_sqrhigh_normalised(dst + 1, xptr, xn, n);
    else
        hi = flint_mpn_mulhigh_normalised(dst + 1, xptr, xn, yptr, yn, n);

    dst[0] = hi.m1;
    e = FLINT_BIT_COUNT(FLINT_MPN_MULHIGH_NORMALISED_GUARD_ERR_BOUND(n));

    if (relaxed || _arf_mul_rounding_certain(dst, an, e, prec))
    {
        inexact = _arf_set_round_mpn(z, &fix, dst, an, sgnbit, prec,
            relaxed ? ARF_RND_NEAR : rnd);
        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix - (slong) hi.m2);
        ARF_MUL_TMP_FREE(tmp, alloc)
        /* the approximation cannot certify exactness */
        return relaxed ? 1 : inexact;
    }

    ARF_MUL_TMP_FREE(tmp, alloc)

    /* Inconclusive rounding check: compute the exact product. */
    return _arf_mul_exact(z, x, y, xn, yn, sgnbit, prec, srnd);
}
