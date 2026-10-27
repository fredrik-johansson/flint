/*
    Copyright (C) 2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arb.h"
#include "arb/impl.h"

static void
_arb_mul_arf_arf_rounding(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    mag_t zr, ym;
    int inexact;

    if (arb_is_exact(x))
    {
        inexact = arf_mul(arb_midref(z), arb_midref(x), y, prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
    else if (ARB_IS_LAGOM(x) && ARF_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(ym, y);
        mag_init(zr);

        mag_fast_mul(zr, ym, arb_radref(x));

        inexact = arf_mul(arb_midref(z), arb_midref(x), y, prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_fast_add_ulp(zr, zr, arb_midref(z), prec);

        *arb_radref(z) = *zr;
    }
    else if (arf_is_inf(y) && arb_is_nonzero(x))
    {
        mag_zero(arb_radref(z));

        if (arf_sgn(arb_midref(x)) * arf_sgn(y) > 0)
            arf_pos_inf(arb_midref(z));
        else
            arf_neg_inf(arb_midref(z));
    }
    else
    {
        mag_init_set_arf(ym, y);
        mag_init(zr);

        mag_mul(zr, ym, arb_radref(x));

        inexact = arf_mul(arb_midref(z), arb_midref(x), y, prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(ym);
        mag_clear(zr);

    }
}

void
arb_mul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    /* view y as an exact ball (shallow copy, never cleared) */
    arb_struct yb;
    *arb_midref(&yb) = *y;
    MAG_MAN(arb_radref(&yb)) = 0;
    MAG_EXP(arb_radref(&yb)) = 0;
    arb_mul(z, x, &yb, prec);
}

/* the previous implementation, used as the fallback of arb_mul */
void
_arb_mul_arf_rounding(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    mag_t zr, xm, ym;
    int inexact;

    if (arb_is_exact(x))
    {
        _arb_mul_arf_arf_rounding(z, y, arb_midref(x), prec);
    }
    else if (arb_is_exact(y))
    {
        _arb_mul_arf_arf_rounding(z, x, arb_midref(y), prec);
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(xm, arb_midref(x));
        mag_fast_init_set_arf(ym, arb_midref(y));

        mag_init(zr);
        mag_fast_mul(zr, xm, arb_radref(y));
        mag_fast_addmul(zr, ym, arb_radref(x));
        mag_fast_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_mul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_fast_add_ulp(zr, zr, arb_midref(z), prec);

        *arb_radref(z) = *zr;
    }
    else if ((arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)) && arb_is_nonzero(y)) ||
             (arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)) && arb_is_nonzero(x)))
    {
        mag_zero(arb_radref(z));

        if (arf_sgn(arb_midref(x)) * arf_sgn(arb_midref(y)) > 0)
            arf_pos_inf(arb_midref(z));
        else
            arf_neg_inf(arb_midref(z));
    }
    else
    {
        mag_init_set_arf(xm, arb_midref(x));
        mag_init_set_arf(ym, arb_midref(y));

        mag_init(zr);
        mag_mul(zr, xm, arb_radref(y));
        mag_addmul(zr, ym, arb_radref(x));
        mag_addmul(zr, arb_radref(x), arb_radref(y));

        inexact = arf_mul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(xm);
        mag_clear(ym);
        mag_clear(zr);
    }
}


/*
    Ball multiplication without going through arf rounding: the propagated
    error is computed first (which also bounds the useful precision of the
    output), the midpoint product is formed as an exact or truncated mpn
    product, and arb_set_round_mpn_2exp_si rounds it to nearest and sets
    the rounding error, to which the propagated error (and the deficit of a
    truncated product) is added. Fast paths require LAGOM exponents.
*/
/*
    err = |xm| yr + |ym| xr + xr yr (+ zr, for the fused multiply-add),
    rounded up, for LAGOM operands with at least one of xr, yr nonzero,
    in one pass: each term is the exact product of a 30-bit (31-bit for
    the midpoints, rounded up) mantissa pair with its exponent, the terms
    are aligned to the largest exponent (rounded up) and the sum is
    rounded up to a mantissa. Symmetric in x and y, so the result, and
    the precision clamped from it, do not depend on the operand order.
*/
void
_arb_mul_error(mag_t err, arf_srcptr xm, mag_srcptr xr, arf_srcptr ym, mag_srcptr yr, mag_srcptr zr)
{
    /* products of two 30-bit mantissas: 64-bit arithmetic on 32-bit limbs */
    uint64_t XR = MAG_MAN(xr), YR = MAG_MAN(yr), S, X, Y;
    slong emax, bc;
    nn_srcptr xptr, yptr;
    slong xn, yn;

#define _SHR_UP(p, s) (((p) >> FLINT_MIN((s), 63)) + ((s) != 0))

    /* the top limbs of the midpoints rounded up to MAG_BITS bits */
    ARF_GET_MPN_READONLY(xptr, xn, xm);
    ARF_GET_MPN_READONLY(yptr, yn, ym);
    X = (xptr[xn - 1] >> (FLINT_BITS - MAG_BITS)) + 1;
    Y = (yptr[yn - 1] >> (FLINT_BITS - MAG_BITS)) + 1;

    if (XR == 0)
    {
        S = X * YR;
        emax = ARF_EXP(xm) + MAG_EXP(yr);
    }
    else if (YR == 0)
    {
        S = Y * XR;
        emax = ARF_EXP(ym) + MAG_EXP(xr);
    }
    else
    {
        uint64_t pa = X * YR, pb = Y * XR, pc = XR * YR;
        slong ea = ARF_EXP(xm) + MAG_EXP(yr), eb = ARF_EXP(ym) + MAG_EXP(xr), ec = MAG_EXP(xr) + MAG_EXP(yr);

        emax = FLINT_MAX(FLINT_MAX(ea, eb), ec);
        S = _SHR_UP(pa, emax - ea) + _SHR_UP(pb, emax - eb) + _SHR_UP(pc, emax - ec);
    }

    /* the extra term, in the same units 2^(emax - 2 MAG_BITS) */
    if (zr != NULL && MAG_MAN(zr) != 0)
    {
        uint64_t pz = (uint64_t) MAG_MAN(zr) << MAG_BITS;
        slong ez = MAG_EXP(zr);

        if (ez > emax)
        {
            S = _SHR_UP(S, ez - emax) + pz;
            emax = ez;
        }
        else
            S = S + _SHR_UP(pz, emax - ez);
    }

#undef _SHR_UP

    /* value = S 2^(emax - 2 MAG_BITS), S < 2^63 */
    bc = (S >> 32) ? 32 + FLINT_BIT_COUNT((ulong) (S >> 32)) : FLINT_BIT_COUNT((ulong) S);
    if (bc > MAG_BITS)
        S = (S >> (bc - MAG_BITS)) + 1;
    else
        S = S << (MAG_BITS - bc);
    MAG_MAN(err) = (ulong) S;
    MAG_EXP(err) = emax - 2 * MAG_BITS + bc;
    MAG_FAST_ADJUST_ONE_TOO_LARGE(err);
}

/* res = (z, zn) 2^exp exactly (zn <= 4, top limb nonzero), radius err
   (NULL for zero) */
static void
_arb_set_exact_limbs(arb_t res, nn_srcptr z, slong zn, int sgnbit, slong exp, mag_srcptr err)
{
    arf_ptr mid = arb_midref(res);
    ulong t[4];
    slong lead, low, n;
    nn_ptr zptr;

    /* strip trailing zero limbs, normalise the top limb */
    low = 0;
    while (z[low] == 0)
        low++;
    n = zn - low;
    lead = flint_clz(z[zn - 1]);
    exp = exp + zn * FLINT_BITS - lead;

    if (lead == 0)
    {
        z = z + low;
    }
    else
    {
        mpn_lshift(t, z + low, n, lead);
        z = t;
        if (t[0] == 0)
        {
            z++;
            n--;
        }
    }

    ARF_GET_MPN_WRITE(zptr, n, mid);
    flint_mpn_copyi(zptr, z, n);
    ARF_XSIZE(mid) = ARF_MAKE_XSIZE(n, sgnbit);
    _fmpz_demote(ARF_EXPREF(mid));
    *ARF_EXPREF(mid) = exp;
    if (err == NULL)
        mag_zero(arb_radref(res));
    else
        mag_set(arb_radref(res), err);
}

static void _arb_mul_generic(arb_t res, arf_srcptr xm, slong xn, arf_srcptr ym, slong yn,
    int sgnbit, slong exp, slong prec, mag_srcptr err);

/*
    res = (-1)^neg x y + [+/- zr] (zr NULL for none) for LAGOM operands with
    nonzero midpoints and a LAGOM zr (the caller checks); see arb_mul.
*/
FLINT_FORCE_INLINE void
_arb_mul_main(arb_t res, const arb_t x, const arb_t y, int neg, slong prec, mag_srcptr zr)
{
    arf_srcptr xm = arb_midref(x), ym = arb_midref(y);
    mag_srcptr xr = arb_radref(x), yr = arb_radref(y);
    mag_t err;
    mag_srcptr errp;
    slong xn, yn, exp, acc;
    int sgnbit;
    nn_srcptr xptr, yptr;

    /* propagated error: |xm| yr + |ym| xr + xr yr (rounded up), plus zr,
       and the precision clamped to the bits of the product above it (an
       exact product is never clamped) */
    mag_init(err);
    if ((MAG_MAN(xr) | MAG_MAN(yr)) != 0 || (zr != NULL && MAG_MAN(zr) != 0))
    {
        if ((MAG_MAN(xr) | MAG_MAN(yr)) != 0)
            _arb_mul_error(err, xm, xr, ym, yr, zr);
        else
        {
            MAG_MAN(err) = MAG_MAN(zr);
            MAG_EXP(err) = MAG_EXP(zr);
        }
        errp = err;

        if (prec != ARF_PREC_EXACT)
        {
            acc = (ARF_EXP(xm) + ARF_EXP(ym)) - MAG_EXP(err);
            if (acc < prec - MAG_BITS)
                prec = FLINT_MAX(2, acc + MAG_BITS);
        }
    }
    else
        errp = NULL;

    xn = ARF_SIZE(xm);
    yn = ARF_SIZE(ym);
    sgnbit = ARF_SGNBIT(xm) ^ ARF_SGNBIT(ym) ^ neg;
    exp = ARF_EXP(xm) + ARF_EXP(ym) - (xn + yn) * FLINT_BITS;   /* value = product 2^exp */

    ARF_GET_MPN_READONLY(xptr, xn, xm);
    ARF_GET_MPN_READONLY(yptr, yn, ym);

    if (xn <= 2 && yn <= 2 && (prec <= 2 * FLINT_BITS || prec >= (xn + yn) * FLINT_BITS - 1))
    {
        /* products of up to four limbs, in registers: rounded to at most
           two limbs, or exact (the product of normalised operands has at
           least (xn + yn) FLINT_BITS - 1 bits) */
        ulong zz[4];
        slong zn = xn + yn;

        if (xn == 1 && yn == 1)
        {
            umul_ppmm(zz[1], zz[0], xptr[0], yptr[0]);
        }
        else if (xn == 2 && yn == 2)
        {
            FLINT_MPN_MUL_2X2(zz[3], zz[2], zz[1], zz[0], xptr[1], xptr[0], yptr[1], yptr[0]);
        }
        else
        {
            if (xn == 2)
                FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], xptr[1], xptr[0], yptr[0]);
            else
                FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], yptr[1], yptr[0], xptr[0]);
        }

        if (prec <= 2 * FLINT_BITS)
        {
            /* value = (hi, lo, lo2) 2^(exp + (zn - 3) FLINT_BITS) with sticky below */
            if (zn == 2)
                _arb_set_round_3limbs(res, zz[1], zz[0], 0, 0, sgnbit, exp - FLINT_BITS, prec, errp, NULL);
            else if (zn == 3)
                _arb_set_round_3limbs(res, zz[2], zz[1], zz[0], 0, sgnbit, exp, prec, errp, NULL);
            else
                _arb_set_round_3limbs(res, zz[3], zz[2], zz[1], zz[0], sgnbit, exp + FLINT_BITS, prec, errp, NULL);
        }
        else if (prec >= zn * FLINT_BITS - (int) (zz[zn - 1] >> (FLINT_BITS - 1) == 0))
        {
            if (zn == 2)
                _arb_set_exact_2limbs(res, zz[1], zz[0], sgnbit, exp, errp);
            else
                _arb_set_exact_limbs(res, zz, zn, sgnbit, exp, errp);
        }
        else
        {
            _arb_set_round_mpn_2exp_si_core(res, zz, zn, sgnbit, exp, prec, 0, 0, errp, NULL);
        }
        return;
    }

    _arb_mul_generic(res, xm, xn, ym, yn, sgnbit, exp, prec, errp);
}

void
arb_mul(arb_t res, const arb_t x, const arb_t y, slong prec)
{
    arf_srcptr xm = arb_midref(x), ym = arb_midref(y);
    mag_srcptr xr = arb_radref(x), yr = arb_radref(y);

    if (arf_is_special(xm) || arf_is_special(ym) || mag_is_inf(xr) || mag_is_inf(yr) ||
        !ARF_IS_LAGOM(xm) || !ARF_IS_LAGOM(ym) || !MAG_IS_LAGOM(xr) || !MAG_IS_LAGOM(yr) ||
        prec < 2)
    {
        _arb_mul_arf_rounding(res, x, y, prec);
        return;
    }

    _arb_mul_main(res, x, y, 0, prec, NULL);
}

/*
    res = (-1)^neg x y + [+/- zr] in one pass, for the fused multiply-add
    with a zero midpoint in the addend (zr may be the radius of res).
    Returns 0 (doing nothing) if the operands are not handled.
*/
int
_arb_mul_add_rad(arb_t res, const arb_t x, const arb_t y, int neg, slong prec, mag_srcptr zr)
{
    arf_srcptr xm = arb_midref(x), ym = arb_midref(y);
    mag_srcptr xr = arb_radref(x), yr = arb_radref(y);
    mag_struct z;

    if (arf_is_special(xm) || arf_is_special(ym) || mag_is_inf(xr) || mag_is_inf(yr) ||
        !ARF_IS_LAGOM(xm) || !ARF_IS_LAGOM(ym) || !MAG_IS_LAGOM(xr) || !MAG_IS_LAGOM(yr) ||
        !MAG_IS_LAGOM(zr) || prec < 2)
        return 0;

    /* a copy: zr may be the radius of the output */
    MAG_MAN(&z) = MAG_MAN(zr);
    MAG_EXP(&z) = MAG_EXP(zr);
    _arb_mul_main(res, x, y, neg, prec, &z);
    return 1;
}

/* the paths beyond the register products: high products and exact mpn
   products (kept out of line so that the small case stays lean) */
static void
_arb_mul_generic(arb_t res, arf_srcptr xm, slong xn, arf_srcptr ym, slong yn,
    int sgnbit, slong exp, slong prec, mag_srcptr err)
{
    nn_srcptr xptr, yptr;

    ARF_GET_MPN_READONLY(xptr, xn, xm);
    ARF_GET_MPN_READONLY(yptr, yn, ym);

    /* (prec may be ARF_PREC_EXACT here: compare before counting limbs,
       which would overflow) */
    if (prec <= ARF_MUL_HIGH_SMALL_MAX * FLINT_BITS &&
        (FLINT_MIN(xn, yn) >= (prec + FLINT_BITS - 1) / FLINT_BITS || FLINT_MAX(xn, yn) > (prec + FLINT_BITS - 1) / FLINT_BITS))
    {
        /* outputs of up to 16 limbs from operands at least as long (or one
           of them longer): normalised high product with guard limb on the
           stack; its deficit (below 2^e units of the guard limb) goes into
           the radius */
        ulong t[ARF_MUL_HIGH_SMALL_MAX + ARB_GUARD_EXTRA + 1];
        mp_limb_pair_t hi;
        slong n = (prec + FLINT_BITS - 1) / FLINT_BITS + ARB_GUARD_EXTRA, aexp;

        if (xm == ym)
            hi = flint_mpn_sqrhigh_normalised(t + 1, xptr, xn, n);
        else if (xn == n && yn == n)
            hi = flint_mpn_mulhigh_n_normalised(t + 1, xptr, yptr, n);
        else
            hi = flint_mpn_mulhigh_normalised(t + 1, xptr, xn, yptr, yn, n);
        t[0] = hi.m1;

        aexp = ARF_EXP(xm) + ARF_EXP(ym) - (slong) hi.m2 - (n + 1) * FLINT_BITS;

        /* the deficit, below 2^ARF_MUL_HIGH_SMALL_ERR_BITS units of t[0],
           is below one unit of the radius (which has FLINT_BITS - 1 or
           more discarded bits below the kept ones): one extra unit */
        _arb_set_round_mpn_2exp_si_core(res, t, n + 1, sgnbit, aexp, prec, 0, 1, err, NULL);
    }
    else if (xn + yn <= ARF_MUL_SMALL_LIMBS || prec >= (xn + yn) * FLINT_BITS)
    {
        /* exact product */
        nn_ptr tmp;
        slong alloc = xn + yn;
        ARF_MUL_TMP_DECL

        ARF_MUL_TMP_ALLOC(tmp, alloc)
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);
        _arb_set_round_mpn_2exp_si_long(res, tmp, xn + yn, sgnbit, exp, prec, err, NULL);
        ARF_MUL_TMP_FREE(tmp, alloc)
    }
    else
    {
        /* truncated product: top outn + 1 limbs (plus guard) of the
           operands truncated to outn + 1 limbs, deficit below 2^e units
           of the guard limb, which goes into the radius */
        slong outn = (prec + FLINT_BITS - 1) / FLINT_BITS + ARB_GUARD_EXTRA, an, e, expfix, alloc;
        nn_ptr tmp;
        int in_place;
        ARF_MUL_TMP_DECL

        /* room for the approximation or, if the helper declines, the full product */
        alloc = ARF_MUL_APPROX_ALLOC(outn, xn, yn);

        /* Large products go straight into the output's buffer and are
           rounded in place, when it does not hold an operand. */
        in_place = (alloc > ARF_MUL_STACK_ALLOC) && (xm != arb_midref(res)) && (ym != arb_midref(res)) &&
            !(ARF_HAS_PTR(arb_midref(res)) && (ARF_PTR_D(arb_midref(res)) == xptr || ARF_PTR_D(arb_midref(res)) == yptr));

        if (in_place)
        {
            ARF_GET_MPN_WRITE(tmp, alloc, arb_midref(res));
        }
        else
        {
            ARF_MUL_TMP_ALLOC(tmp, alloc)
        }

        an = _arf_mul_approx(tmp, &e, &expfix, xptr, xn, yptr, yn, outn, xm == ym);

        if (an == 0)
        {
            FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);
            if (in_place)
                _arb_set_round_mpn_2exp_si_inplace_err(res, xn + yn, sgnbit, exp, prec, 0, err);
            else
                _arb_set_round_mpn_2exp_si_err(res, tmp, xn + yn, sgnbit, exp, prec, 0, err);
        }
        else
        {
            /* value = (A / B^an) 2^(xexp + yexp + expfix) = A 2^(xexp + yexp + expfix - an FLINT_BITS);
               the deficit (below 2^e units of tmp[0], e <= 34) is below
               one unit of the radius, which has at least FLINT_BITS
               discarded bits (the guard limb) below the kept ones */
            slong aexp = ARF_EXP(xm) + ARF_EXP(ym) + expfix - an * FLINT_BITS;
            FLINT_ASSERT(e + MAG_BITS <= an * FLINT_BITS - prec);
            if (in_place)
                _arb_set_round_mpn_2exp_si_inplace_err(res, an, sgnbit, aexp, prec, 1, err);
            else
                _arb_set_round_mpn_2exp_si_err(res, tmp, an, sgnbit, aexp, prec, 1, err);
        }

        if (!in_place)
            ARF_MUL_TMP_FREE(tmp, alloc)
    }
}

/* Up to this many limbs the product by a single limb is formed on the
   stack in full; beyond, and when the precision is much smaller than the
   operand, the general code truncates. */
#define ARB_MUL_LIMB_STACK 64

/*
    res = x * (-1)^ysgnbit y for a single nonzero limb y: the mantissa
    product is a single mpn_mul_1 (in registers for x of up to two limbs),
    the propagated error is xr y rounded up, both go through the one-pass
    rounding. Falls back on arb_mul_new (viewing y as an exact ball) for
    special or non-LAGOM x, and when the product should be truncated.
*/
static void
_arb_mul_limb(arb_t res, const arb_t x, ulong y, int ysgnbit, slong prec)
{
    arf_srcptr xm = arb_midref(x);
    mag_srcptr xr = arb_radref(x);
    mag_t err;
    mag_srcptr errp;
    slong xn, exp, acc, c, yexp;
    int sgnbit;
    nn_srcptr xptr;
    ulong ym;

    c = flint_clz(y);
    ym = y << c;
    yexp = FLINT_BITS - c;       /* y = (ym / B) 2^yexp */

    /* (the limb count of prec only once prec != ARF_PREC_EXACT: it would
       overflow) */
    if (arf_is_special(xm) || mag_is_inf(xr) || !ARF_IS_LAGOM(xm) || !MAG_IS_LAGOM(xr) ||
        prec < 2 || prec == ARF_PREC_EXACT ||
        ARF_SIZE(xm) > (prec + FLINT_BITS - 1) / FLINT_BITS + 2 || ARF_SIZE(xm) > ARB_MUL_LIMB_STACK)
    {
        arb_struct yb;
        ARF_EXP(arb_midref(&yb)) = yexp;
        ARF_NOPTR_D(arb_midref(&yb))[0] = ym;
        ARF_XSIZE(arb_midref(&yb)) = ARF_MAKE_XSIZE(1, ysgnbit);
        MAG_MAN(arb_radref(&yb)) = 0;
        MAG_EXP(arb_radref(&yb)) = 0;
        arb_mul(res, x, &yb, prec);
        return;
    }

    /* propagated error xr y, and the precision clamped to the bits of the
       product above it */
    mag_init(err);
    if (MAG_MAN(xr) != 0)
    {
        uint64_t S = (uint64_t) MAG_MAN(xr) * ((ym >> (FLINT_BITS - MAG_BITS)) + 1);
        slong bc = (S >> 32) ? 32 + FLINT_BIT_COUNT((ulong) (S >> 32)) : FLINT_BIT_COUNT((ulong) S);

        /* value = S 2^(MAG_EXP(xr) + yexp - 2 MAG_BITS), S < 2^61 */
        if (bc > MAG_BITS)
            S = (S >> (bc - MAG_BITS)) + 1;
        else
            S = S << (MAG_BITS - bc);
        MAG_MAN(err) = (ulong) S;
        MAG_EXP(err) = MAG_EXP(xr) + yexp - 2 * MAG_BITS + bc;
        MAG_FAST_ADJUST_ONE_TOO_LARGE(err);
        errp = err;

        acc = (ARF_EXP(xm) + yexp) - MAG_EXP(err);
        if (acc < prec - MAG_BITS)
            prec = FLINT_MAX(2, acc + MAG_BITS);
    }
    else
        errp = NULL;

    xn = ARF_SIZE(xm);
    sgnbit = ARF_SGNBIT(xm) ^ ysgnbit;
    exp = ARF_EXP(xm) + yexp - (xn + 1) * FLINT_BITS;   /* value = product 2^exp */
    ARF_GET_MPN_READONLY(xptr, xn, xm);

    if (xn <= 2)
    {
        ulong zz[3];
        slong zn = xn + 1;

        if (xn == 1)
            umul_ppmm(zz[1], zz[0], xptr[0], ym);
        else
            FLINT_MPN_MUL_2X1(zz[2], zz[1], zz[0], xptr[1], xptr[0], ym);

        if (prec <= 2 * FLINT_BITS)
        {
            if (zn == 2)
                _arb_set_round_3limbs(res, zz[1], zz[0], 0, 0, sgnbit, exp - FLINT_BITS, prec, errp, NULL);
            else
                _arb_set_round_3limbs(res, zz[2], zz[1], zz[0], 0, sgnbit, exp, prec, errp, NULL);
        }
        else if (prec >= zn * FLINT_BITS - (int) (zz[zn - 1] >> (FLINT_BITS - 1) == 0))
        {
            if (zn == 2)
                _arb_set_exact_2limbs(res, zz[1], zz[0], sgnbit, exp, errp);
            else
                _arb_set_exact_limbs(res, zz, zn, sgnbit, exp, errp);
        }
        else
        {
            _arb_set_round_mpn_2exp_si_core(res, zz, zn, sgnbit, exp, prec, 0, 0, errp, NULL);
        }
    }
    else
    {
        ulong t[ARB_MUL_LIMB_STACK + 1];

        t[xn] = mpn_mul_1(t, xptr, xn, ym);
        _arb_set_round_mpn_2exp_si_core(res, t, xn + 1, sgnbit, exp, prec, 0, 0, errp, NULL);
    }
}

void
arb_mul_ui(arb_t res, const arb_t x, ulong y, slong prec)
{
    if (y == 0)
    {
        arf_t t;
        arf_init(t);
        arb_mul_arf(res, x, t, prec);   /* the special-value semantics */
    }
    else
        _arb_mul_limb(res, x, y, 0, prec);
}

void
arb_mul_si(arb_t res, const arb_t x, slong y, slong prec)
{
    if (y == 0)
    {
        arf_t t;
        arf_init(t);
        arb_mul_arf(res, x, t, prec);
    }
    else
        _arb_mul_limb(res, x, FLINT_UABS(y), y < 0, prec);
}

void
arb_mul_fmpz(arb_t res, const arb_t x, const fmpz_t y, slong prec)
{
    if (!COEFF_IS_MPZ(*y))
    {
        slong v = *y;
        if (v == 0)
            arb_mul_si(res, x, 0, prec);
        else
            _arb_mul_limb(res, x, FLINT_UABS(v), v < 0, prec);
    }
    else
    {
        arf_t t;
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_mul_arf(res, x, t, prec);
        arf_clear(t);
    }
}
