/*
    Copyright (C) 2012-2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arb.h"
#include "arb/impl.h"

/* the previous implementation: midpoint rounded by arf, radius bumped by
   an ulp when inexact */
void
_arb_add_rounding(arb_t z, const arb_t x, const arb_t y, int ysub, slong prec)
{
    int inexact;

    if (ysub)
        inexact = arf_sub(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);
    else
        inexact = arf_add(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARB_RND);

    mag_add(arb_radref(z), arb_radref(x), arb_radref(y));
    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
}

/*
    res = x + (-1)^ysub y as balls, in one pass: the radius is the sum of
    the radii (rounded up) plus the measured rounding error of the exact
    sum of the midpoints, rounded to nearest. The precision is clamped to
    what the propagated error leaves meaningful. Operands of up to two
    limbs are added in registers (or in a frame of five limbs); otherwise
    the exact sum is formed by _arf_add_mpn_raw and rounded by
    arb_set_round_mpn_2exp_si. Specials and non-LAGOM exponents go through
    the previous implementation.
*/
void
_arb_add_sub(arb_t res, const arb_t x, const arb_t y, int ysub, slong prec)
{
    arf_srcptr xm = arb_midref(x), ym = arb_midref(y);
    mag_srcptr xr = arb_radref(x), yr = arb_radref(y);
    mag_t err;
    mag_srcptr errp;
    slong xn, yn, xexp, shift, sn, sn_original, acc;
    int xs, ys, ssgnbit;
    nn_srcptr xptr, yptr;

    if (arf_is_special(xm) || arf_is_special(ym) || mag_is_inf(xr) || mag_is_inf(yr) ||
        !ARF_IS_LAGOM(xm) || !ARF_IS_LAGOM(ym) || !MAG_IS_LAGOM(xr) || !MAG_IS_LAGOM(yr) ||
        prec < 2)
    {
        /* one midpoint zero: the other rounded, with the sum of the radii */
        if (arf_is_zero(xm) != arf_is_zero(ym) && MAG_IS_LAGOM(xr) && MAG_IS_LAGOM(yr))
        {
            mag_struct e;
            MAG_MAN(&e) = 0;
            MAG_EXP(&e) = 0;
            mag_fast_add(&e, xr, yr);

            if (arf_is_zero(xm) ? _arb_set_round_arf_err(res, ym, ysub, prec, &e)
                                : _arb_set_round_arf_err(res, xm, 0, prec, &e))
                return;
        }

        _arb_add_rounding(res, x, y, ysub, prec);
        return;
    }

    xs = ARF_SGNBIT(xm);
    ys = ARF_SGNBIT(ym) ^ ysub;

    /* x first if it has the larger exponent */
    shift = ARF_EXP(xm) - ARF_EXP(ym);
    if (shift < 0)
    {
        FLINT_SWAP(arf_srcptr, xm, ym);
        FLINT_SWAP(int, xs, ys);
        shift = -shift;
    }
    xexp = ARF_EXP(xm);

    /* propagated error, and the precision clamped to the bits of the sum
       above it (the sum is below 2^(xexp + 1)) */
    if ((MAG_MAN(xr) | MAG_MAN(yr)) == 0)
    {
        errp = NULL;
    }
    else
    {
        mag_init(err);
        mag_fast_add(err, xr, yr);
        errp = err;

        if (prec != ARF_PREC_EXACT)
        {
            acc = xexp + 1 - MAG_EXP(err);
            if (acc < prec - MAG_BITS)
                prec = FLINT_MAX(2, acc + MAG_BITS);
        }
    }

    xn = ARF_SIZE(xm);
    yn = ARF_SIZE(ym);
    ARF_GET_MPN_READONLY(xptr, xn, xm);
    ARF_GET_MPN_READONLY(yptr, yn, ym);

    if (xn <= 2 && yn <= 2)
    {
        if (xn == 1 && yn == 1 && shift < FLINT_BITS - 1)
        {
            /* the sum has fewer than 2 FLINT_BITS bits: in two registers,
               as (hi, lo) 2^(xexp - 2 FLINT_BITS + 1) */
            ulong hi, lo, xhi, xlo, yhi, ylo;

            xhi = xptr[0] >> 1;
            xlo = xptr[0] << (FLINT_BITS - 1);
            yhi = yptr[0] >> (shift + 1);
            ylo = yptr[0] << (FLINT_BITS - (shift + 1));

            if (xs == ys)
            {
                add_ssaaaa(hi, lo, xhi, xlo, yhi, ylo);
            }
            else if (xhi > yhi || (xhi == yhi && xlo >= ylo))
            {
                sub_ddmmss(hi, lo, xhi, xlo, yhi, ylo);
                if ((hi | lo) == 0)
                {
                    arf_zero(arb_midref(res));
                    if (errp == NULL)
                        mag_zero(arb_radref(res));
                    else
                        mag_set(arb_radref(res), errp);
                    return;
                }
            }
            else
            {
                sub_ddmmss(hi, lo, yhi, ylo, xhi, xlo);
                xs = ys;
            }

            if (hi == 0)
            {
                hi = lo;
                lo = 0;
                xexp -= FLINT_BITS;
            }

            if (prec <= 2 * FLINT_BITS)
                _arb_set_round_3limbs(res, hi, lo, 0, 0, xs, xexp - 3 * FLINT_BITS + 1, prec, errp, NULL);
            else
                _arb_set_exact_2limbs(res, hi, lo, xs, xexp - 2 * FLINT_BITS + 1, errp);
            return;
        }
        else if (prec <= 2 * FLINT_BITS)
        {
            /* the top three limbs of the sum in registers, with the part
               of y below them (if any) as a sticky bit */
            ulong t[4], sticky;

            sn = _arb_add_window(t, 3, xptr, xn, xs, yptr, yn, ys, shift, &ssgnbit, &sticky);
            if (sn == 0)
            {
                arf_zero(arb_midref(res));
                if (errp == NULL)
                    mag_zero(arb_radref(res));
                else
                    mag_set(arb_radref(res), errp);
                return;
            }
            _arb_set_round_frame(res, t, sn, ssgnbit, xexp - 3 * FLINT_BITS, prec, errp, sticky);
            return;
        }
        else if (shift <= 3 * FLINT_BITS)
        {
            /* the exact sum in a frame of five limbs: x in the top two, y
               entirely within */
            ulong t[6];

            sn = _arb_add_frame(t, 5, xptr, xn, xs, yptr, yn, ys, shift, &ssgnbit);
            if (sn == 0)
            {
                arf_zero(arb_midref(res));
                if (errp == NULL)
                    mag_zero(arb_radref(res));
                else
                    mag_set(arb_radref(res), errp);
                return;
            }
            _arb_set_round_frame(res, t, sn, ssgnbit, xexp - 5 * FLINT_BITS, prec, errp, 0);
            return;
        }
    }

    if (prec != ARF_PREC_EXACT && shift > prec + FLINT_BITS + 1 && shift > (xn + 1) * FLINT_BITS)
    {
        /* x +/- eps: y lies more than FLINT_BITS bits below the precision */
        _arb_set_round_mpn_2exp_si_eps(res, xptr, xn, xs, xexp - xn * FLINT_BITS, xs != ys, prec, errp);
        return;
    }

    {
        nn_ptr tmp;
        slong alloc;
        ARF_ADD_TMP_DECL

        alloc = _ARF_ADD_MPN_SIZE(xn, yn, shift);
        ARF_ADD_TMP_ALLOC(tmp, alloc)

        sn = _arf_add_mpn_raw(tmp, &sn_original, &ssgnbit, xptr, xn, xs, yptr, yn, ys, shift);

        if (sn == 0)
        {
            arf_zero(arb_midref(res));
            if (errp == NULL)
                mag_zero(arb_radref(res));
            else
                mag_set(arb_radref(res), errp);
        }
        else
        {
            /* value = (S / B^sn_original) 2^xexp */
            _arb_set_round_mpn_2exp_si_err(res, tmp, sn, ssgnbit, xexp - sn_original * FLINT_BITS, prec, 0, errp);
        }

        ARF_ADD_TMP_FREE(tmp, alloc)
    }
}

void
arb_add(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    _arb_add_sub(z, x, y, 0, prec);
}

void
arb_add_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    /* view y as an exact ball (shallow copy, never cleared) */
    arb_struct yb;
    *arb_midref(&yb) = *y;
    MAG_MAN(arb_radref(&yb)) = 0;
    MAG_EXP(arb_radref(&yb)) = 0;
    _arb_add_sub(z, x, &yb, 0, prec);
}

void
arb_add_ui(arb_t z, const arb_t x, ulong y, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_add_arf(z, x, t, prec);
}

void
arb_add_si(arb_t z, const arb_t x, slong y, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_add_arf(z, x, t, prec);
}

void
arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_add_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_add_arf(z, x, t, prec);
        arf_clear(t);
    }
}

void
arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t man, const fmpz_t exp, slong prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz_2exp(t, man, exp);
    arb_add_arf(z, x, t, prec);
    arf_clear(t);
}
