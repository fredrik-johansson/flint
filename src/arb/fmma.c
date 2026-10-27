/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb/impl.h"
#include "arb/terms.h"

/* err += |xm| yr + xr |ym| + xr yr */
static void
_arb_fmma_mul_error(mag_t err, const arb_t x, const arb_t y)
{
    mag_t xm, ym, t;
    mag_init(xm);
    mag_init(ym);
    mag_init(t);
    arf_get_mag(xm, arb_midref(x));
    arf_get_mag(ym, arb_midref(y));
    mag_addmul(err, xm, arb_radref(y));
    mag_add(t, ym, arb_radref(y));
    mag_addmul(err, t, arb_radref(x));
    mag_clear(xm);
    mag_clear(ym);
    mag_clear(t);
}

/* the propagated error |xm| yr + |ym| xr + xr yr of a product with a zero
   midpoint in x or y, for LAGOM data */
static void
_arb_fmma_zero_product_error(mag_t e, const arb_t x, const arb_t y)
{
    mag_struct t;
    MAG_MAN(&t) = 0;
    MAG_EXP(&t) = 0;

    if (arf_is_zero(arb_midref(x)))
    {
        mag_fast_init_set_arf(&t, arb_midref(y));
        mag_fast_add(&t, &t, arb_radref(y));
        mag_fast_mul(e, &t, arb_radref(x));
    }
    else
    {
        mag_fast_init_set_arf(&t, arb_midref(x));
        mag_fast_add(&t, &t, arb_radref(x));
        mag_fast_mul(e, &t, arb_radref(y));
    }
}

/* scratch for the two products on the stack */
#define ARB_FMMA_TMP_STACK 128

/*
    res = a b + c d in one pass, as one part of the complex multiplication
    of acb_mul (arb/terms.h): the propagated error |a| br + |b| ar + ar br +
    |c| dr + |d| cr + cr dr as one aligned sum, the two products exact (in
    registers up to two limbs at prec <= 2 FLINT_BITS) or truncated with a
    guard limb, exact when they nearly cancel, added exactly and rounded
    once to nearest with the measured rounding error and deficits. A zero
    midpoint reduces to a product with an extra radius. Returns 0 for
    other inputs (specials, large exponents, prec out of range).
*/
static int
_arb_fmma_onepass(arb_t res, const arb_t x1, const arb_t y1, const arb_t x2, const arb_t y2, slong prec)
{
    arf_srcptr a = arb_midref(x1), b = arb_midref(y1), c = arb_midref(x2), d = arb_midref(y2);
    mag_srcptr ar = arb_radref(x1), br = arb_radref(y1), cr = arb_radref(x2), dr = arb_radref(y2);
    mag_t am, bm, cm, dm, er;
    mag_srcptr erp;
    _arb_term_struct T[2];
    nn_srcptr ap, bp, cp, dp;
    slong an, bn, cn, dn, outn, maxn, slot, alloc, e1, e2;
    ulong tmp_stack[ARB_FMMA_TMP_STACK];
    nn_ptr tmp;
    int s1, s2;

    if (prec < 2 || prec == ARF_PREC_EXACT)
        return 0;

    {
        const ulong M = UWORD(1) << (FLINT_BITS - 6);
        ulong range = (ulong) (ARF_EXP(a) + M) | (ulong) (ARF_EXP(b) + M)
                    | (ulong) (ARF_EXP(c) + M) | (ulong) (ARF_EXP(d) + M)
                    | (ulong) (MAG_EXP(ar) + M) | (ulong) (MAG_EXP(br) + M)
                    | (ulong) (MAG_EXP(cr) + M) | (ulong) (MAG_EXP(dr) + M);

        /* zero midpoints have exponent 0, inside the range */
        if (range >= 2 * M)
            return 0;
    }

    if (ARF_IS_SPECIAL(a) || ARF_IS_SPECIAL(b) || ARF_IS_SPECIAL(c) || ARF_IS_SPECIAL(d))
    {
        int z1, z2;

        /* only zeros (not infinities or NaN) are handled */
        if ((ARF_IS_SPECIAL(a) && !arf_is_zero(a)) || (ARF_IS_SPECIAL(b) && !arf_is_zero(b)) ||
            (ARF_IS_SPECIAL(c) && !arf_is_zero(c)) || (ARF_IS_SPECIAL(d) && !arf_is_zero(d)))
            return 0;

        z1 = arf_is_zero(a) || arf_is_zero(b);
        z2 = arf_is_zero(c) || arf_is_zero(d);

        if (z1 && z2)
            return 0;   /* a zero midpoint: no rounding */

        MAG_MAN(er) = 0;
        MAG_EXP(er) = 0;

        if (z1)
        {
            _arb_fmma_zero_product_error(er, x1, y1);
            return _arb_mul_add_rad(res, x2, y2, 0, prec, er);
        }
        else
        {
            _arb_fmma_zero_product_error(er, x2, y2);
            return _arb_mul_add_rad(res, x1, y1, 0, prec, er);
        }
    }

    /* the propagated error, before the output (which may alias the inputs)
       is written */
    if ((MAG_MAN(ar) | MAG_MAN(br) | MAG_MAN(cr) | MAG_MAN(dr)) != 0)
    {
        mag_fast_init_set_arf(am, a);
        mag_fast_init_set_arf(bm, b);
        mag_fast_init_set_arf(cm, c);
        mag_fast_init_set_arf(dm, d);
        _arb_mag_sum6(er, am, br, bm, ar, ar, br, cm, dr, dm, cr, cr, dr);
        erp = er;
    }
    else
        erp = NULL;

    ARF_GET_MPN_READONLY(ap, an, a);
    ARF_GET_MPN_READONLY(bp, bn, b);
    ARF_GET_MPN_READONLY(cp, cn, c);
    ARF_GET_MPN_READONLY(dp, dn, d);
    e1 = ARF_EXP(a) + ARF_EXP(b);
    e2 = ARF_EXP(c) + ARF_EXP(d);
    s1 = ARF_SGNBIT(a) ^ ARF_SGNBIT(b);
    s2 = ARF_SGNBIT(c) ^ ARF_SGNBIT(d);

#if FLINT_BITS == 64
    if (prec <= 2 * FLINT_BITS && ((an - 1) | (bn - 1) | (cn - 1) | (dn - 1)) == 0)
    {
        ulong h0, l0, h1, l1;
        umul_ppmm(h0, l0, ap[0], bp[0]);
        umul_ppmm(h1, l1, cp[0], dp[0]);
        _arb_round_sum_1limb(res, h0, l0, e1, s1, h1, l1, e2, s2, prec, erp);
        return 1;
    }

    if (prec <= 2 * FLINT_BITS && ((an - 1) | (bn - 1) | (cn - 1) | (dn - 1)) <= 1)
    {
        ulong pt[2][4];
        _arb_term_exact(T + 0, pt[0], ap, an, bp, bn, e1, s1);
        _arb_term_exact(T + 1, pt[1], cp, cn, dp, dn, e2, s2);
        _arb_round_acc4(res, T, 2, prec, erp, 0);
        return 1;
    }
#endif

    outn = (prec + FLINT_BITS - 1) / FLINT_BITS + 1;
    maxn = FLINT_MAX(FLINT_MAX(an, bn), FLINT_MAX(cn, dn));
    slot = (outn >= 2 * maxn) ? 2 * maxn : 2 * maxn + 2;
    alloc = 2 * slot;
    tmp = (alloc <= ARB_FMMA_TMP_STACK) ? tmp_stack
        : (alloc <= ARF_MUL_TLS_ALLOC) ? _arf_mul_tls_tmp(alloc) : flint_malloc(sizeof(ulong) * alloc);

    /* exact products when both are cheap in full, or when they nearly
       cancel (so that their deficits would dominate) */
    if ((_arb_product_full(an, bn, outn) && _arb_product_full(cn, dn, outn)) ||
        (s1 != s2 && _arb_products_cancel(ap + an - 1, bp + bn - 1, e1, cp + cn - 1, dp + dn - 1, e2)))
    {
        _arb_term_exact(T + 0, tmp, ap, an, bp, bn, e1, s1);
        _arb_term_exact(T + 1, tmp + slot, cp, cn, dp, dn, e2, s2);
    }
    else
    {
        _arb_term_high(T + 0, tmp, ap, an, bp, bn, outn, e1, s1);
        _arb_term_high(T + 1, tmp + slot, cp, cn, dp, dn, outn, e2, s2);
    }

    _arb_round_sum2(res, T, T + 1, prec, erp);

    if (alloc > ARF_MUL_TLS_ALLOC)
        flint_free(tmp);

    return 1;
}

void
arb_fmma(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t d, slong prec)
{
    mag_t err;
    int inexact;

    if (_arb_fmma_onepass(res, a, b, c, d, prec))
        return;

    if (!arb_is_finite(a) || !arb_is_finite(b) || !arb_is_finite(c) || !arb_is_finite(d))
    {
        arb_t t;
        arb_init(t);
        arb_mul(t, a, b, prec);
        arb_addmul(t, c, d, prec);
        arb_swap(res, t);
        arb_clear(t);
        return;
    }

    mag_init(err);

    if (!arb_is_exact(a) || !arb_is_exact(b))
        _arb_fmma_mul_error(err, a, b);
    if (!arb_is_exact(c) || !arb_is_exact(d))
        _arb_fmma_mul_error(err, c, d);

    inexact = arf_fmma(arb_midref(res), arb_midref(a), arb_midref(b),
        arb_midref(c), arb_midref(d), prec, ARF_RND_FAST);

    if (inexact)
        arf_mag_add_ulp(err, err, arb_midref(res), prec);

    mag_swap(arb_radref(res), err);
    mag_clear(err);
}

void
arb_sosq(arb_t res, const arb_t a, const arb_t b, slong prec)
{
    arb_fmma(res, a, a, b, b, prec);
}
