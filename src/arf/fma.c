/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

/* second scratch buffer for the sum in the fma cores (the product uses
   ARF_MUL_TMP_ALLOC) */
#define SUM_TMP_ALLOC(ptr, alloc) \
    if ((alloc) <= ARF_MUL_STACK_ALLOC) \
        ptr = tmp_stack2; \
    else \
        ptr = flint_malloc(sizeof(ulong) * (alloc));
#define SUM_TMP_FREE(ptr, alloc) \
    if ((alloc) > ARF_MUL_STACK_ALLOC) \
        flint_free(ptr);

#include "mpn_extras.h"


/* Exact product, then exact addition and rounding (the original algorithm). */
int
_arf_fma_exact(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z,
    int tsgnbit, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, zn, tn, alloc;
    nn_srcptr xptr, yptr, zptr;
    nn_ptr tptr, tptr2;
    fmpz_t texp;
    slong shift;
    int inexact;
    ARF_MUL_TMP_DECL

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);

    _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), 0);
    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);

    alloc = tn = xn + yn;
    ARF_MUL_TMP_ALLOC(tptr2, alloc)
    tptr = tptr2;

    FLINT_MPN_MUL_WITH_SPECIAL_CASES(tptr, xptr, xn, yptr, yn);

    tn -= (tptr[0] == 0);
    tptr += (tptr[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(res, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(res, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr2, alloc)
    fmpz_clear(texp);

    return inexact;
}

/* Bits of the guard window needed to skip the rounding check with a
   relaxed rounding mode: the result is then within 0.5 + 2^-30 ulp. */
#define ARF_FMA_RELAXED_MARGIN 30

/*
    Small-size fused multiply-add: prec <= n * FLINT_BITS with n <= 16
    limbs of output, x and y at least n limbs (or one of them longer), the
    product not small (below ARF_FMA_HIGH_SMALL_MIN_PRODUCT limb products
    the full product and exact addition are as fast), LAGOM exponents. The
    normalised high product with guard limb of the (internally truncated)
    operands is formed on the stack, z is added exactly on the stack, the
    rounding is certified from the guard window (as in _arf_fma_mpn) and
    done in registers. Returns -1 (z untouched) when the shape is
    unsuitable or the check is inconclusive.
*/
int
_arf_fma_high_small(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z,
    int tsgnbit, slong prec, arf_rnd_t rnd)
{
    ulong A[ARF_MUL_HIGH_SMALL_MAX + 1], S[3 * ARF_MUL_HIGH_SMALL_MAX + 8], g;
    nn_srcptr xptr, yptr, zptr;
    slong xn, yn, zn, n, an, shift, sn, sn_original, p, e, k, E, lead, fix, nn;
    int ssgnbit, relaxed, zsgnbit;
    mp_limb_pair_t hi;

    if (prec == ARF_PREC_EXACT || !ARF_IS_LAGOM(x) || !ARF_IS_LAGOM(y) || !ARF_IS_LAGOM(z))
        return -1;

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);
    n = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (n < 1 || n > ARF_MUL_HIGH_SMALL_MAX || zn > 2 * n ||
        xn * yn < ARF_FMA_HIGH_SMALL_MIN_PRODUCT || !(FLINT_MIN(xn, yn) >= n || FLINT_MAX(xn, yn) > n))
        return -1;

    relaxed = (rnd >= ARF_RND_FAST);
    zsgnbit = ARF_SGNBIT(z);
    an = n + 1;

    /* exponent of z minus that of the product (up to the normalising
       shift of the high product, at most one bit): both operands must be
       within reach of each other, else the general code handles it */
    shift = ARF_EXP(z) - (ARF_EXP(x) + ARF_EXP(y));
    if (shift > (zn + 2) * FLINT_BITS - 1 || -shift > (an + 2) * FLINT_BITS - 1)
        return -1;


    if (x == y)
        hi = flint_mpn_sqrhigh_normalised(A + 1, xptr, xn, n);
    else
        hi = flint_mpn_mulhigh_normalised(A + 1, xptr, xn, yptr, yn, n);
    A[0] = hi.m1;
    shift += (slong) hi.m2;

    if (shift >= 0)
    {
        slong shift_limbs = (ulong) shift / FLINT_BITS, shift_bits = (ulong) shift % FLINT_BITS;
        slong wn = an + (shift_bits != 0);
        sn = _arf_add_mpn_raw(S, &sn_original, &ssgnbit, zptr, zn, zsgnbit, A, an, tsgnbit, shift);
        p = FLINT_BITS * (sn_original - shift_limbs - wn) + ((shift_bits != 0) ? FLINT_BITS - shift_bits : 0);
        E = ARF_EXP(z);
    }
    else
    {
        sn = _arf_add_mpn_raw(S, &sn_original, &ssgnbit, A, an, tsgnbit, zptr, zn, zsgnbit, -shift);
        p = FLINT_BITS * (sn_original - an);
        E = ARF_EXP(x) + ARF_EXP(y) - (slong) hi.m2;
    }

    if (sn == 0)
        return -1;

    /* the exact result lies within 2^ARF_MUL_HIGH_SMALL_ERR_BITS units of
       A's ulp on one side */
    e = ARF_MUL_HIGH_SMALL_ERR_BITS + p;
    lead = flint_clz(S[sn - 1]);
    k = sn * FLINT_BITS - lead - prec;

    if (relaxed)
    {
        if (k - 1 - e < ARF_FMA_RELAXED_MARGIN)
            return -1;
        rnd = ARF_RND_NEAR;
    }
    else if (!_arf_mul_rounding_certain(S, sn, e, prec))
        return -1;

    /* value = S / B^sn * 2^E' with E' = E + (sn - sn_original) FLINT_BITS;
       normalise and round the top nn limbs with the next limb as guard */
    E += (sn - sn_original) * FLINT_BITS - lead;
    if (lead != 0)
        mpn_lshift(S, S, sn, lead);

    nn = FLINT_MIN(n, sn);
    g = (sn > nn) ? S[sn - nn - 1] : 0;

    fix = _arf_set_round_small_certified(res, S + sn - nn, nn, g, ssgnbit, prec, rnd);
    fmpz_set_si(ARF_EXPREF(res), E + fix);

    return 1;
}

/*
    res = z + (-1)^tsgnbit |x y|, rounded, for finite nonzero x, y, z.

    When the full product would have many more limbs than can affect the
    result, we use a lower approximation A of its top limbs (see
    _arf_mul_approx), add z exactly to A, and certify the rounding of that
    sum S as in _arf_mul_mpn: the exact sum lies within 2^e units of some
    limb position of S on a known side, and rounding S is correct in all
    strict modes if the bits of S from that position up to the rounding
    point are neither all zero nor all one. Cancellation shrinks the
    guard window and the check then fails, in which case (as for any
    inconclusive check, or non-small exponents) we fall back to the exact
    computation. The number of product limbs needed is reduced when z has
    larger magnitude, as in the dot product.
*/
int
_arf_fma_mpn(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z,
    int tsgnbit, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, zn, an, sn, sn_original, outn, eprec, e, expfix, fix, p;
    slong shift, ashift, k;
    nn_srcptr xptr, yptr, zptr;
    nn_ptr aptr, sptr;
    int inexact, ssgnbit, relaxed, zsgnbit, squaring;
    arf_rnd_t srnd;
    fmpz_t texp;
    ARF_MUL_TMP_DECL
    ulong tmp_stack2[ARF_MUL_STACK_ALLOC];

    relaxed = (rnd >= ARF_RND_FAST);
    srnd = relaxed ? arf_rnd_relaxed_to_strict(rnd) : rnd;

    if (prec < 1)
        prec = 1;

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);
    zsgnbit = ARF_SGNBIT(z);

    /* Only bother when the product is long and the exponents are small. */
    if (prec == ARF_PREC_EXACT || xn + yn <= 4 ||
        !ARF_IS_LAGOM(x) || !ARF_IS_LAGOM(y) || !ARF_IS_LAGOM(z))
        return _arf_fma_exact(res, x, y, z, tsgnbit, prec, srnd);


    /* shift = exponent of z minus exponent of the product */
    shift = ARF_EXP(z) - (ARF_EXP(x) + ARF_EXP(y));

    /* Product limbs above the rounding point, assuming no cancellation:
       prec bits below the larger of the two exponents. */
    eprec = prec - FLINT_MAX(shift, 0);

    if (eprec < -FLINT_BITS)
    {
        /* The product is entirely below the guard limb of the rounding
           point, or two limbs below the bottom of z: only its sign matters
           (the "x +/- eps" case of _arf_add_mpn), except when it interacts
           with the low limbs of z through a borrow, in which case we
           still need the top limbs below. */
        if (shift > (zn + 1) * FLINT_BITS && shift > prec + FLINT_BITS + 1)
        {
            ulong one = 1;
            return _arf_add_mpn(res, zptr, zn, zsgnbit, ARF_EXPREF(z),
                &one, 1, tsgnbit, shift, prec, srnd);
        }
        eprec = 0;
    }

    outn = (eprec + FLINT_BITS - 1) / FLINT_BITS;
    if (outn < 0)
        outn = 0;

    /* Quick cancellation check: opposite signs and comparable magnitudes
       with the top limbs nearly equal. */
    if (tsgnbit != zsgnbit && shift >= -2 && shift <= 2)
    {
        ulong hi, lo, zt;
        umul_ppmm(hi, lo, xptr[xn - 1], yptr[yn - 1]);
        /* hi is the top limb of the product (exponent ex + ey), possibly with
           one leading zero; align z's top limb to the same exponent. */
        zt = zptr[zn - 1];
        if (shift > 0) hi >>= shift; else if (shift < 0) zt >>= -shift;
        hi = (hi > zt) ? hi - zt : zt - hi;
        if (hi < (UWORD(1) << (FLINT_BITS - FLINT_MIN(30, FLINT_BITS - 2))))
            return _arf_fma_exact(res, x, y, z, tsgnbit, prec, srnd);
    }

    squaring = (x == y);

    {
        slong alloc, salloc;
        nn_ptr tmp;

        alloc = ARF_MUL_APPROX_ALLOC(outn, xn, yn);
        ARF_MUL_TMP_ALLOC(tmp, alloc)
        aptr = tmp;

        an = _arf_mul_approx(aptr, &e, &expfix, xptr, xn, yptr, yn, outn, squaring);

        if (an == 0)
        {
            ARF_MUL_TMP_FREE(tmp, alloc)
            return _arf_fma_exact(res, x, y, z, tsgnbit, prec, srnd);
        }

        /* Add z exactly to the approximation. */
        ashift = shift - expfix;   /* exponent of z minus exponent of A's frame */

        if (ashift >= 0)
        {
            /* z is the larger operand: z +/- A 2^-ashift */
            if (ashift > (zn + 2) * FLINT_BITS && ashift > prec + FLINT_BITS + 2)
            {
                /* The product is below the two guard limbs appended to z
                   (it is smaller than one unit of the lowest limb), so
                   z +/- eps has the same rounding as the exact sum. */
                sn_original = zn + 2;
                salloc = sn_original + 1;
                SUM_TMP_ALLOC(sptr, salloc)
                sptr[0] = sptr[1] = 0;
                flint_mpn_copyi(sptr + 2, zptr, zn);
                if (tsgnbit == zsgnbit)
                    sptr[0] = 1;
                else
                    mpn_sub_1(sptr, sptr, sn_original, 1);
                sn = sn_original;
                while (sptr[sn - 1] == 0) sn--;
                ssgnbit = zsgnbit;
                /* the exact sum is within one unit of sptr[0] of S */
                p = 0;
                e = 0;
            }
            else
            {
                slong shift_limbs = ashift / FLINT_BITS, shift_bits = ashift % FLINT_BITS;
                slong wn = an + (shift_bits != 0);
                salloc = _ARF_ADD_MPN_SIZE(zn, an, ashift);
                SUM_TMP_ALLOC(sptr, salloc)
                sn = _arf_add_mpn_raw(sptr, &sn_original, &ssgnbit, zptr, zn, zsgnbit, aptr, an, tsgnbit, ashift);
                /* position of A's ulp within the frame */
                p = FLINT_BITS * (sn_original - shift_limbs - wn) + ((shift_bits != 0) ? FLINT_BITS - shift_bits : 0);
            }
            fmpz_init(texp);
            fmpz_set(texp, ARF_EXPREF(z));
        }
        else
        {
            /* the product is the larger operand: A +/- z 2^-(-ashift) */
            ashift = -ashift;
            if (ashift > (an + 2) * FLINT_BITS)
            {
                /* z is smaller than one unit of the lowest of two guard
                   limbs appended to A: the exact sum is within
                   2^e units of A's ulp plus one unit of S. */
                sn_original = an + 2;
                salloc = sn_original + 1;
                SUM_TMP_ALLOC(sptr, salloc)
                sptr[0] = sptr[1] = 0;
                flint_mpn_copyi(sptr + 2, aptr, an);
                if (tsgnbit == zsgnbit)
                    sptr[0] = 1;
                else
                    mpn_sub_1(sptr, sptr, sn_original, 1);
                sn = sn_original;
                while (sptr[sn - 1] == 0) sn--;
                ssgnbit = tsgnbit;
                p = 2 * FLINT_BITS + 1;
            }
            else
            {
                salloc = _ARF_ADD_MPN_SIZE(an, zn, ashift);
                SUM_TMP_ALLOC(sptr, salloc)
                sn = _arf_add_mpn_raw(sptr, &sn_original, &ssgnbit, aptr, an, tsgnbit, zptr, zn, zsgnbit, ashift);
                p = FLINT_BITS * (sn_original - an);
            }
            fmpz_init(texp);
            _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), expfix);
        }

        ARF_MUL_TMP_FREE(tmp, alloc)

        /* The exact sum lies within 2^(e + p) units of sptr[0] of S, on one
           side. Certify the rounding, or with a relaxed mode just check
           that the guard window is wide enough. */
        inexact = -1;

        if (sn > 0)
        {
            e += p;
            k = sn * FLINT_BITS - flint_clz(sptr[sn - 1]) - prec;

            if (relaxed)
            {
                if (k - 1 - e >= ARF_FMA_RELAXED_MARGIN)
                {
                    _arf_set_round_mpn(res, &fix, sptr, sn, ssgnbit, prec, ARF_RND_NEAR);
                    inexact = 1;
                }
            }
            else if (_arf_mul_rounding_certain(sptr, sn, e, prec))
            {
                inexact = _arf_set_round_mpn(res, &fix, sptr, sn, ssgnbit, prec, rnd);
            }

            if (inexact != -1)
            {
                fix += (sn - sn_original) * FLINT_BITS;
                _fmpz_add_fast(ARF_EXPREF(res), texp, fix);
            }
        }

        SUM_TMP_FREE(sptr, salloc)
        fmpz_clear(texp);
    }

    if (inexact == -1)
        return _arf_fma_exact(res, x, y, z, tsgnbit, prec, srnd);

    return inexact;
}

/* Exact products, then exact addition and rounding. */
static int
_arf_fmma_exact(arf_ptr res, arf_srcptr a, arf_srcptr b, int sgn1,
    arf_srcptr c, arf_srcptr d, int sgn2, slong prec, arf_rnd_t rnd)
{
    nn_srcptr ap, bp, cp, dp;
    nn_ptr tmp, t1, t2;
    slong an, bn, cn, dn, n1, n2, alloc, shift;
    int inexact;
    fmpz_t e1, e2;
    ARF_MUL_TMP_DECL

    ARF_GET_MPN_READONLY(ap, an, a);
    ARF_GET_MPN_READONLY(bp, bn, b);
    ARF_GET_MPN_READONLY(cp, cn, c);
    ARF_GET_MPN_READONLY(dp, dn, d);

    fmpz_init(e1);
    fmpz_init(e2);
    _fmpz_add2_fast(e1, ARF_EXPREF(a), ARF_EXPREF(b), 0);
    _fmpz_add2_fast(e2, ARF_EXPREF(c), ARF_EXPREF(d), 0);
    shift = _fmpz_sub_small(e1, e2);

    n1 = an + bn;
    n2 = cn + dn;
    alloc = n1 + n2;
    ARF_MUL_TMP_ALLOC(tmp, alloc)
    t1 = tmp;
    t2 = tmp + n1;

    FLINT_MPN_MUL_WITH_SPECIAL_CASES(t1, ap, an, bp, bn);
    n1 -= (t1[0] == 0);
    t1 += (t1[0] == 0);
    FLINT_MPN_MUL_WITH_SPECIAL_CASES(t2, cp, cn, dp, dn);
    n2 -= (t2[0] == 0);
    t2 += (t2[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(res, t1, n1, sgn1, e1, t2, n2, sgn2, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(res, t2, n2, sgn2, e2, t1, n1, sgn1, -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tmp, alloc)
    fmpz_clear(e1);
    fmpz_clear(e2);
    return inexact;
}

/*
    res = (-1)^sgn1 |a b| + (-1)^sgn2 |c d| for finite nonzero a, b, c, d.
    Both products are approximated from above the rounding point only
    (each to the precision it can contribute, given the other's magnitude)
    and added exactly; the exact result then lies within
    2^(e1 + p1) + 2^(e2 + p2) units of the sum, and the rounding is
    certified as in _arf_fma_mpn, with an exact fallback.
*/
int
_arf_fmma_mpn(arf_ptr res, arf_srcptr a, arf_srcptr b, int sgn1,
    arf_srcptr c, arf_srcptr d, int sgn2, slong prec, arf_rnd_t rnd)
{
    nn_srcptr ap, bp, cp, dp;
    slong an, bn, cn, dn, n1, n2, shift, ashift, out1, out2, e1, e2, fix1, fix2;
    slong p1, p2, e, sn, sn_original, fix, k, alloc, salloc;
    nn_ptr tmp, A1, A2, sptr;
    int inexact, ssgnbit, relaxed;
    arf_rnd_t srnd;
    fmpz_t texp;
    ARF_MUL_TMP_DECL
    ulong tmp_stack2[ARF_MUL_STACK_ALLOC];

    relaxed = (rnd >= ARF_RND_FAST);
    srnd = relaxed ? arf_rnd_relaxed_to_strict(rnd) : rnd;

    if (prec < 1)
        prec = 1;

    ARF_GET_MPN_READONLY(ap, an, a);
    ARF_GET_MPN_READONLY(bp, bn, b);
    ARF_GET_MPN_READONLY(cp, cn, c);
    ARF_GET_MPN_READONLY(dp, dn, d);

    if (prec == ARF_PREC_EXACT || an + bn + cn + dn <= 8 ||
        !ARF_IS_LAGOM(a) || !ARF_IS_LAGOM(b) || !ARF_IS_LAGOM(c) || !ARF_IS_LAGOM(d))
        return _arf_fmma_exact(res, a, b, sgn1, c, d, sgn2, prec, srnd);

    /* small products not much longer than the output: compute in full */
    if (an * bn <= ARF_FMA_SMALL_PRODUCT && cn * dn <= ARF_FMA_SMALL_PRODUCT &&
        FLINT_MAX(an + bn, cn + dn) <= 2 + 2 * ((prec + FLINT_BITS - 1) / FLINT_BITS))
        return _arf_fmma_exact(res, a, b, sgn1, c, d, sgn2, prec, srnd);

    /* exponent of ab minus exponent of cd */
    shift = (ARF_EXP(a) + ARF_EXP(b)) - (ARF_EXP(c) + ARF_EXP(d));

    /* limbs of each product above the rounding point (no cancellation) */
    out1 = prec - FLINT_MAX(-shift, 0);
    out2 = prec - FLINT_MAX(shift, 0);
    out1 = (out1 <= 0) ? 0 : (out1 + FLINT_BITS - 1) / FLINT_BITS;
    out2 = (out2 <= 0) ? 0 : (out2 + FLINT_BITS - 1) / FLINT_BITS;

    /* quick cancellation check on the top limbs */
    if (sgn1 != sgn2 && shift >= -2 && shift <= 2)
    {
        ulong h1, h2, lo;
        umul_ppmm(h1, lo, ap[an - 1], bp[bn - 1]);
        umul_ppmm(h2, lo, cp[cn - 1], dp[dn - 1]);
        if (shift > 0) h2 >>= shift; else if (shift < 0) h1 >>= -shift;
        h1 = (h1 > h2) ? h1 - h2 : h2 - h1;
        if (h1 < (UWORD(1) << (FLINT_BITS - FLINT_MIN(30, FLINT_BITS - 2))))
            return _arf_fmma_exact(res, a, b, sgn1, c, d, sgn2, prec, srnd);
    }

    alloc = ARF_MUL_APPROX_ALLOC(out1, an, bn) + ARF_MUL_APPROX_ALLOC(out2, cn, dn);
    ARF_MUL_TMP_ALLOC(tmp, alloc)
    A1 = tmp;
    A2 = tmp + ARF_MUL_APPROX_ALLOC(out1, an, bn);

    n1 = _arf_mul_approx(A1, &e1, &fix1, ap, an, bp, bn, out1, a == b);
    if (n1 == 0)
    {
        /* exact product */
        n1 = an + bn;
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(A1, ap, an, bp, bn);
        e1 = -1; fix1 = 0;   /* no error */
    }
    n2 = _arf_mul_approx(A2, &e2, &fix2, cp, cn, dp, dn, out2, c == d);
    if (n2 == 0)
    {
        n2 = cn + dn;
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(A2, cp, cn, dp, dn);
        e2 = -1; fix2 = 0;
    }

    /* exponent of A1's frame minus exponent of A2's frame */
    ashift = shift + fix1 - fix2;

    fmpz_init(texp);

    /* Symmetric setup: X is the product with the larger frame exponent. */
    {
        nn_ptr X, Y;
        slong xn, yn, xe, xfix, px, py;
        int xs, ys;
        arf_srcptr xa, xb;

        if (ashift >= 0)
        {
            X = A1; xn = n1; xs = sgn1; xe = e1; xfix = fix1; xa = a; xb = b;
            Y = A2; yn = n2; ys = sgn2;
        }
        else
        {
            X = A2; xn = n2; xs = sgn2; xe = e2; xfix = fix2; xa = c; xb = d;
            Y = A1; yn = n1; ys = sgn1;
            ashift = -ashift;
        }

        if (ashift > (xn + 2) * FLINT_BITS)
        {
            /* Y is smaller than one unit of the lowest of two guard limbs
               appended to X: X +/- eps. */
            sn_original = xn + 2;
            salloc = sn_original + 1;
            SUM_TMP_ALLOC(sptr, salloc)
            sptr[0] = sptr[1] = 0;
            flint_mpn_copyi(sptr + 2, X, xn);
            if (xs == ys)
                sptr[0] = 1;
            else
                mpn_sub_1(sptr, sptr, sn_original, 1);
            sn = sn_original;
            while (sptr[sn - 1] == 0)
                sn--;
            ssgnbit = xs;
            px = 2 * FLINT_BITS + 1;
            py = -1;   /* accounted for by the +1 in px */
        }
        else
        {
            slong shift_limbs = ashift / FLINT_BITS, shift_bits = ashift % FLINT_BITS;
            slong wn = yn + (shift_bits != 0);
            salloc = _ARF_ADD_MPN_SIZE(xn, yn, ashift);
            SUM_TMP_ALLOC(sptr, salloc)
            sn = _arf_add_mpn_raw(sptr, &sn_original, &ssgnbit, X, xn, xs, Y, yn, ys, ashift);
            px = FLINT_BITS * (sn_original - xn);
            py = FLINT_BITS * (sn_original - shift_limbs - wn) + ((shift_bits != 0) ? FLINT_BITS - shift_bits : 0);
        }

        if (X == A1) { p1 = px; p2 = py; }
        else { p2 = px; p1 = py; }

        /* an operand dropped as eps contributes no error of its own */
        if (py < 0) { if (X == A1) e2 = -1; else e1 = -1; }

        _fmpz_add2_fast(texp, ARF_EXPREF(xa), ARF_EXPREF(xb), xfix);
        (void) xe;
    }

    ARF_MUL_TMP_FREE(tmp, alloc)

    /* combined two-sided error bound: 2^(e1+p1) + 2^(e2+p2) <= 2^e */
    e = FLINT_MAX((e1 < 0) ? -1 : e1 + p1, (e2 < 0) ? -1 : e2 + p2) + 1;

    inexact = -1;

    if (sn > 0)
    {
        k = sn * FLINT_BITS - flint_clz(sptr[sn - 1]) - prec;

        if (relaxed)
        {
            if (k - 1 - e >= ARF_FMA_RELAXED_MARGIN)
            {
                _arf_set_round_mpn(res, &fix, sptr, sn, ssgnbit, prec, ARF_RND_NEAR);
                inexact = 1;
            }
        }
        else if (_arf_mul_rounding_certain(sptr, sn, e, prec))
        {
            inexact = _arf_set_round_mpn(res, &fix, sptr, sn, ssgnbit, prec, rnd);
        }

        if (inexact != -1)
        {
            fix += (sn - sn_original) * FLINT_BITS;
            _fmpz_add_fast(ARF_EXPREF(res), texp, fix);
        }
    }

    SUM_TMP_FREE(sptr, salloc)
    fmpz_clear(texp);

    if (inexact == -1)
        return _arf_fmma_exact(res, a, b, sgn1, c, d, sgn2, prec, srnd);

    return inexact;
}

int
arf_fmma(arf_ptr res, arf_srcptr a, arf_srcptr b, arf_srcptr c, arf_srcptr d, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(a) || arf_is_special(b) || arf_is_special(c) || arf_is_special(d))
    {
        arf_t t;
        int inexact;
        if (rnd >= ARF_RND_FAST)
            rnd = arf_rnd_relaxed_to_strict(rnd);
        arf_init(t);
        arf_mul(t, c, d, ARF_PREC_EXACT, ARF_RND_DOWN);
        inexact = arf_fma(res, a, b, t, prec, rnd);
        arf_clear(t);
        return inexact;
    }

    return _arf_fmma_mpn(res, a, b, ARF_SGNBIT(a) ^ ARF_SGNBIT(b),
        c, d, ARF_SGNBIT(c) ^ ARF_SGNBIT(d), prec, rnd);
}

int
arf_fma(arf_ptr res, arf_srcptr x, arf_srcptr y, arf_srcptr z, slong prec, arf_rnd_t rnd)
{
    slong xn, yn, zn, tn, alloc;
    nn_srcptr xptr, yptr, zptr;
    nn_ptr tptr, tptr2;
    fmpz_t texp;
    slong shift;
    int tsgnbit, inexact;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
    {
        if (rnd >= ARF_RND_FAST)
            rnd = arf_rnd_relaxed_to_strict(rnd);

        if (arf_is_zero(z))
        {
            return arf_mul(res, x, y, prec, rnd);
        }
        else if (arf_is_finite(x) && arf_is_finite(y))
        {
            return arf_set_round(res, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_t t;
            arf_init(t);
            arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_add(res, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    /* outputs of up to 16 limbs where the product is not small: high
       product of the (truncated) operands with the addition and certified
       rounding done on the stack */
    if (ARF_SIZE(x) * ARF_SIZE(y) >= ARF_FMA_HIGH_SMALL_MIN_PRODUCT &&
        prec <= ARF_MUL_HIGH_SMALL_MAX * FLINT_BITS)
    {
        int r = _arf_fma_high_small(res, x, y, z, ARF_SGNBIT(x) ^ ARF_SGNBIT(y), prec, rnd);
        if (r != -1)
            return r;
    }

    /* long products: truncated product with certification after the
       addition (see _arf_fma_mpn); small ones are computed in full below */
    if (ARF_SIZE(x) * ARF_SIZE(y) > ARF_FMA_SMALL_PRODUCT || (prec != ARF_PREC_EXACT &&
        ARF_SIZE(x) + ARF_SIZE(y) > 2 + 2 * ((prec + FLINT_BITS - 1) / FLINT_BITS)))
        return _arf_fma_mpn(res, x, y, z, ARF_SGNBIT(x) ^ ARF_SGNBIT(y), prec, rnd);

    if (rnd >= ARF_RND_FAST)
        rnd = arf_rnd_relaxed_to_strict(rnd);

    tsgnbit = ARF_SGNBIT(x) ^ ARF_SGNBIT(y);
    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);

    _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), 0);
    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);

    alloc = tn = xn + yn;
    ARF_MUL_TMP_ALLOC(tptr2, alloc)
    tptr = tptr2;

    FLINT_MPN_MUL_WITH_SPECIAL_CASES(tptr, xptr, xn, yptr, yn);

    tn -= (tptr[0] == 0);
    tptr += (tptr[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(res, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(res, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr2, alloc)
    fmpz_clear(texp);

    return inexact;
}
