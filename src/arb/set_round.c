/*
    Copyright (C) 2012-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "arb/impl.h"

/* limbs of a value copied to the stack when it is also the output */
#define ARB_SET_ROUND_STACK 64

/*
    res = (-1)^neg {x, xn} 2^exp (top limb nonzero) rounded to nearest at
    prec bits, radius err (NULL for none) plus the measured rounding error;
    the value may be the midpoint of res itself (alias != 0). LAGOM
    exponents, 2 <= prec < ARF_PREC_EXACT.
*/
FLINT_FORCE_INLINE void
_arb_set_round_limbs(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, slong prec, mag_srcptr err, int alias)
{
    slong zn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (xn <= 2)
    {
        /* registers (prec < 2 FLINT_BITS here unless exact) */
        ulong hi = x[xn - 1], lo = (xn == 2) ? x[0] : 0;
        slong bits = FLINT_BITS - flint_clz(hi) + (xn - 1) * FLINT_BITS;

        if (bits <= prec)
        {
            if (lo == 0)
                _arb_set_exact_2limbs(res, hi, 0, sgnbit, exp + (xn - 2) * FLINT_BITS, err);
            else
                _arb_set_exact_2limbs(res, hi, lo, sgnbit, exp, err);
        }
        else if (xn == 2)
            _arb_set_round_3limbs(res, hi, lo, 0, 0, sgnbit, exp - FLINT_BITS, prec, err, NULL);
        else
            _arb_set_round_3limbs(res, hi, 0, 0, 0, sgnbit, exp - 2 * FLINT_BITS, prec, err, NULL);
    }
    else if (!alias || (xn > zn && zn <= ARB_ROUND_WINDOW_MAX))
    {
        /* the window path reads all it needs before writing */
        _arb_set_round_mpn_2exp_si_long(res, x, xn, sgnbit, exp, prec, err, NULL);
    }
    else
    {
        ulong stack[ARB_SET_ROUND_STACK];
        nn_ptr t = (xn <= ARB_SET_ROUND_STACK) ? stack : flint_malloc(xn * sizeof(ulong));
        flint_mpn_copyi(t, x, xn);
        _arb_set_round_mpn_2exp_si_long(res, t, xn, sgnbit, exp, prec, err, NULL);
        if (t != stack)
            flint_free(t);
    }
}

/*
    res = (-1)^neg x rounded to nearest at prec bits with the measured
    rounding error plus err (NULL for none) as radius, the precision clamped
    to what err leaves meaningful as in arb_add; x may be the midpoint of
    res and err its radius. Returns 0 (doing nothing) if x or err is special
    or not LAGOM, or prec is out of range.
*/
FLINT_FORCE_INLINE int
_arb_set_round_arf_err_inline(arb_t res, arf_srcptr x, int neg, slong prec, mag_srcptr err)
{
    mag_struct e;
    mag_srcptr ep = NULL;
    nn_srcptr xp;
    slong xn;

    if (arf_is_special(x) || !ARF_IS_LAGOM(x) || prec < 2 || prec == ARF_PREC_EXACT)
        return 0;

    /* (an infinite radius has a zero mantissa: check the exponent first) */
    if (err != NULL && !MAG_IS_LAGOM(err))
        return 0;

    if (err != NULL && MAG_MAN(err) != 0)
    {
        slong acc;

        /* a copy: err may be the radius of res */
        MAG_MAN(&e) = MAG_MAN(err);
        MAG_EXP(&e) = MAG_EXP(err);
        ep = &e;

        acc = ARF_EXP(x) - MAG_EXP(err);
        if (acc < prec - MAG_BITS)
            prec = FLINT_MAX(2, acc + MAG_BITS);
    }

    ARF_GET_MPN_READONLY(xp, xn, x);

    if (xn * FLINT_BITS - (slong) flint_clz(xp[xn - 1]) <= prec)
    {
        /* exact (the usual case of a value already at this precision) */
        arf_set(arb_midref(res), x);
        if (neg)
            arf_neg(arb_midref(res), arb_midref(res));
        if (ep != NULL)
            mag_set(arb_radref(res), ep);
        else
            mag_zero(arb_radref(res));
        return 1;
    }

    _arb_set_round_limbs(res, xp, xn, ARF_SGNBIT(x) ^ neg, ARF_EXP(x) - xn * FLINT_BITS,
        prec, ep, x == arb_midref(res));
    return 1;
}

int
_arb_set_round_arf_err(arb_t res, arf_srcptr x, int neg, slong prec, mag_srcptr err)
{
    return _arb_set_round_arf_err_inline(res, x, neg, prec, err);
}

void
arb_set_round(arb_t z, const arb_t x, slong prec)
{
    int inexact;

    if (_arb_set_round_arf_err_inline(z, arb_midref(x), 0, prec, arb_radref(x)))
        return;

    inexact = arf_set_round(arb_midref(z), arb_midref(x), prec, ARB_RND);

    if (inexact)
        arf_mag_add_ulp(arb_radref(z), arb_radref(x), arb_midref(z), prec);
    else
        mag_set(arb_radref(z), arb_radref(x));
}

void
arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t exp, slong prec)
{
    int inexact;

    if (!fmpz_is_zero(x) && !COEFF_IS_MPZ(*exp) && prec >= 2 && prec != ARF_PREC_EXACT &&
        *exp > ARF_MIN_LAGOM_EXP / 2 && *exp < ARF_MAX_LAGOM_EXP / 2)
    {
        ulong t;
        nn_srcptr xp;
        slong xn;
        int sgnbit;

        if (!COEFF_IS_MPZ(*x))
        {
            t = FLINT_UABS(*x);
            xp = &t;
            xn = 1;
            sgnbit = (*x < 0);
        }
        else
        {
            mpz_srcptr z = COEFF_TO_PTR(*x);
            xp = z->_mp_d;
            xn = FLINT_ABS(z->_mp_size);
            sgnbit = (z->_mp_size < 0);
        }

        if (xn < ARF_MAX_LAGOM_EXP / (4 * FLINT_BITS))
        {
            /* the value is {xp, xn} 2^exp with an integer mantissa */
            _arb_set_round_limbs(y, xp, xn, sgnbit, *exp, prec, NULL, 0);
            return;
        }
    }

    inexact = arf_set_round_fmpz_2exp(arb_midref(y), x, exp, prec, ARB_RND);

    if (inexact)
        arf_mag_set_ulp(arb_radref(y), arb_midref(y), prec);
    else
        mag_zero(arb_radref(y));
}

void
arb_set_round_fmpz(arb_t y, const fmpz_t x, slong prec)
{
    fmpz e = 0;
    arb_set_round_fmpz_2exp(y, x, &e, prec);
}
