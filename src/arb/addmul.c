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

/* The previous implementations: z + (-1)^neg x y with the midpoint
   rounded by arf, one ulp added to the radius when inexact. */
static void
_arb_addmul_arf_rounding(arb_t z, const arb_t x, const arf_t y, int neg, slong prec)
{
    mag_t ym;
    int inexact;

    if (arb_is_exact(x))
    {
        if (neg)
            inexact = arf_submul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
        else
            inexact = arf_addmul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
    }
    else if (ARB_IS_LAGOM(x) && ARF_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(ym, y);
        mag_fast_addmul(arb_radref(z), ym, arb_radref(x));
        if (neg)
            inexact = arf_submul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
        else
            inexact = arf_addmul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);

        if (inexact)
            arf_mag_fast_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
    }
    else if (arf_is_inf(y) && arb_is_nonzero(x))
    {
        if ((arf_sgn(arb_midref(x)) > 0) ^ neg)
            arb_add_arf(z, z, y, prec);
        else
            arb_sub_arf(z, z, y, prec);
    }
    else
    {
        mag_init_set_arf(ym, y);
        mag_addmul(arb_radref(z), ym, arb_radref(x));

        if (neg)
            inexact = arf_submul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
        else
            inexact = arf_addmul(arb_midref(z), arb_midref(x), y, prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);

        mag_clear(ym);
    }
}

void
_arb_addmul_rounding(arb_t z, const arb_t x, const arb_t y, int neg, slong prec)
{
    mag_t zr, xm, ym;
    int inexact;

    if (arb_is_exact(y))
    {
        _arb_addmul_arf_rounding(z, x, arb_midref(y), neg, prec);
    }
    else if (arb_is_exact(x))
    {
        _arb_addmul_arf_rounding(z, y, arb_midref(x), neg, prec);
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y) && ARB_IS_LAGOM(z))
    {
        mag_fast_init_set_arf(xm, arb_midref(x));
        mag_fast_init_set_arf(ym, arb_midref(y));

        mag_fast_init_set(zr, arb_radref(z));
        mag_fast_addmul(zr, xm, arb_radref(y));
        mag_fast_addmul(zr, ym, arb_radref(x));
        mag_fast_addmul(zr, arb_radref(x), arb_radref(y));

        if (neg)
            inexact = arf_submul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
        else
            inexact = arf_addmul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);

        if (inexact)
            arf_mag_fast_add_ulp(zr, zr, arb_midref(z), prec);

        *arb_radref(z) = *zr;
    }
    else if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)) && arb_is_nonzero(y))
    {
        if ((arf_sgn(arb_midref(y)) > 0) ^ neg)
            arb_add_arf(z, z, arb_midref(x), prec);
        else
            arb_sub_arf(z, z, arb_midref(x), prec);
    }
    else if (arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)) && arb_is_nonzero(x))
    {
        if ((arf_sgn(arb_midref(x)) > 0) ^ neg)
            arb_add_arf(z, z, arb_midref(y), prec);
        else
            arb_sub_arf(z, z, arb_midref(y), prec);
    }
    else
    {
        mag_init_set_arf(xm, arb_midref(x));
        mag_init_set_arf(ym, arb_midref(y));

        mag_init_set(zr, arb_radref(z));
        mag_addmul(zr, xm, arb_radref(y));
        mag_addmul(zr, ym, arb_radref(x));
        mag_addmul(zr, arb_radref(x), arb_radref(y));

        if (neg)
            inexact = arf_submul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
        else
            inexact = arf_addmul(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_set(arb_radref(z), zr);

        mag_clear(zr);
        mag_clear(xm);
        mag_clear(ym);
    }
}

/* res = z + (-1)^neg x y by the previous implementation, for any aliasing */
static void
_arb_addmul_fallback(arb_t res, const arb_t z, const arb_t x, const arb_t y, int neg, slong prec)
{
    if (res == z)
    {
        _arb_addmul_rounding(res, x, y, neg, prec);
    }
    else if (res != x && res != y)
    {
        arb_set(res, z);
        _arb_addmul_rounding(res, x, y, neg, prec);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arb_set(t, z);
        _arb_addmul_rounding(t, x, y, neg, prec);
        arb_swap(res, t);
        arb_clear(t);
    }
}

/* Frame for the register path: the product of up to four limbs and a
   term of up to two limbs, up to this far apart (in limbs) */
#define ARB_ADDMUL_FRAME 7

/*
    res = z + (-1)^neg x y as balls, in one pass: the propagated error
    (|xm| yr + |ym| xr + xr yr + zr, rounded up) first, which also clamps
    the precision; then the midpoint product, exact in registers for
    operands of up to two limbs and added to z in a frame, or formed once
    in full or truncated by _arf_mul_approx, added exactly to z on the
    stack and rounded to nearest with the rounding error (and the deficit
    of a truncated product) measured into the radius. Specials, non-LAGOM
    exponents and terms too far apart for the exact addition to be cheap
    go through the previous implementation.
*/
void
_arb_addmul(arb_t res, const arb_t z, const arb_t x, const arb_t y, int neg, slong prec)
{
    arf_srcptr xm = arb_midref(x), ym = arb_midref(y), zm = arb_midref(z);
    mag_srcptr xr = arb_radref(x), yr = arb_radref(y), zr = arb_radref(z);
    mag_t err;
    mag_srcptr errp;
    slong xn, yn, zn, aexp, an, e, expfix, outn, alloc, salloc, shift, sn, sn_original, acc, pexp;
    int sgnbit, ssgnbit;
    nn_srcptr xptr, yptr, zptr;
    nn_ptr tmp, sptr;
    ARF_MUL_TMP_DECL

    if (arf_is_special(xm) || arf_is_special(ym) || arf_is_special(zm) ||
        mag_is_inf(xr) || mag_is_inf(yr) || mag_is_inf(zr) ||
        !ARF_IS_LAGOM(xm) || !ARF_IS_LAGOM(ym) || !ARF_IS_LAGOM(zm) ||
        !MAG_IS_LAGOM(xr) || !MAG_IS_LAGOM(yr) || !MAG_IS_LAGOM(zr) || prec < 2)
    {
        /* zero midpoints with LAGOM data: one rounding, measured */
        if (MAG_IS_LAGOM(xr) && MAG_IS_LAGOM(yr) && MAG_IS_LAGOM(zr) &&
            (arf_is_zero(xm) || ARF_IS_LAGOM(xm)) && (arf_is_zero(ym) || ARF_IS_LAGOM(ym)) &&
            (arf_is_zero(zm) || ARF_IS_LAGOM(zm)))
        {
            if (arf_is_zero(zm))
            {
                /* the product, with the radius of z */
                if (!arf_is_zero(xm) && !arf_is_zero(ym) && _arb_mul_add_rad(res, x, y, neg, prec, zr))
                    return;
            }
            else if (arf_is_zero(xm) || arf_is_zero(ym))
            {
                /* the product has a zero midpoint: z rounded, with
                   zr + |xm| yr + |ym| xr + xr yr */
                mag_struct e, t;
                MAG_MAN(&e) = 0;
                MAG_EXP(&e) = 0;
                MAG_MAN(&t) = 0;
                MAG_EXP(&t) = 0;
                if (arf_is_zero(xm))
                {
                    mag_fast_init_set_arf(&t, ym);
                    mag_fast_add(&t, &t, yr);
                    mag_fast_mul(&e, &t, xr);
                }
                else
                {
                    mag_fast_init_set_arf(&t, xm);
                    mag_fast_add(&t, &t, xr);
                    mag_fast_mul(&e, &t, yr);
                }
                mag_fast_add(&e, &e, zr);
                if (_arb_set_round_arf_err(res, zm, 0, prec, &e))
                    return;
            }
        }

        _arb_addmul_fallback(res, z, x, y, neg, prec);
        return;
    }

    /* propagated error of the product, plus z's radius */
    if ((MAG_MAN(xr) | MAG_MAN(yr)) != 0)
    {
        mag_init(err);
        _arb_mul_error(err, xm, xr, ym, yr, zr);
        errp = err;
    }
    else if (MAG_MAN(zr) != 0)
    {
        MAG_MAN(err) = MAG_MAN(zr);
        MAG_EXP(err) = MAG_EXP(zr);
        errp = err;
    }
    else
        errp = NULL;

    xn = ARF_SIZE(xm);
    yn = ARF_SIZE(ym);
    zn = ARF_SIZE(zm);
    sgnbit = ARF_SGNBIT(xm) ^ ARF_SGNBIT(ym) ^ neg;
    pexp = ARF_EXP(xm) + ARF_EXP(ym);

    /* precision clamp from the error against the larger of the two terms
       (the sum is below 2^(max + 1)) */
    if (errp != NULL && prec != ARF_PREC_EXACT)
    {
        acc = FLINT_MAX(pexp, ARF_EXP(zm)) + 1 - MAG_EXP(err);
        if (acc < prec - MAG_BITS)
            prec = FLINT_MAX(2, acc + MAG_BITS);
    }

    ARF_GET_MPN_READONLY(xptr, xn, xm);
    ARF_GET_MPN_READONLY(yptr, yn, ym);
    ARF_GET_MPN_READONLY(zptr, zn, zm);

    shift = pexp - ARF_EXP(zm);

    if (xn <= 2 && yn <= 2 && zn <= 2)
    {
        /* the exact product in registers (its top limb may have a leading
           zero bit: it is placed by pexp regardless) */
        ulong zz[4], t[ARB_ADDMUL_FRAME + 1], sticky;
        slong pn = xn + yn, fexp;

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
            /* z added in a window of four limbs headed by the larger
               term, the part of the other below it as a sticky bit */
            if (shift >= 0)
            {
                sn = _arb_add_window(t, 4, zz, pn, sgnbit, zptr, zn, ARF_SGNBIT(zm), shift, &ssgnbit, &sticky);
                fexp = pexp - 4 * FLINT_BITS;
            }
            else
            {
                sn = _arb_add_window(t, 4, zptr, zn, ARF_SGNBIT(zm), zz, pn, sgnbit, -shift, &ssgnbit, &sticky);
                fexp = ARF_EXP(zm) - 4 * FLINT_BITS;
            }
        }
        else if ((shift >= 0 && shift <= (ARB_ADDMUL_FRAME - 3) * FLINT_BITS) ||
                 (shift < 0 && -shift <= (ARB_ADDMUL_FRAME - 5) * FLINT_BITS))
        {
            /* the exact sum in a frame */
            sticky = 0;
            if (shift >= 0)
            {
                sn = _arb_add_frame(t, ARB_ADDMUL_FRAME, zz, pn, sgnbit, zptr, zn, ARF_SGNBIT(zm), shift, &ssgnbit);
                fexp = pexp - ARB_ADDMUL_FRAME * FLINT_BITS;
            }
            else
            {
                sn = _arb_add_frame(t, ARB_ADDMUL_FRAME, zptr, zn, ARF_SGNBIT(zm), zz, pn, sgnbit, -shift, &ssgnbit);
                fexp = ARF_EXP(zm) - ARB_ADDMUL_FRAME * FLINT_BITS;
            }
        }
        else
            goto generic;

        if (sn == 0)
        {
            arf_zero(arb_midref(res));
            if (errp == NULL)
                mag_zero(arb_radref(res));
            else
                mag_set(arb_radref(res), errp);
        }
        else
            _arb_set_round_frame(res, t, sn, ssgnbit, fexp, prec, errp, sticky);
        return;
    }

generic:
    /* the product only needs the bits above the rounding point, which is
       prec bits below the larger term */
    outn = (prec == ARF_PREC_EXACT) ? WORD_MAX : prec - FLINT_MAX(-shift, 0);
    outn = (outn <= 0) ? 0 : (outn >= (xn + yn) * FLINT_BITS) ? xn + yn : (outn + FLINT_BITS - 1) / FLINT_BITS;

    alloc = FLINT_MAX(ARF_MUL_APPROX_ALLOC(outn, xn, yn), xn + yn);
    ARF_MUL_TMP_ALLOC(tmp, alloc)

    /* Cancellation makes the deficit of a truncated product dominate the
       rounding error; when the terms have opposite signs and comparable
       magnitude with nearly equal top limbs, use the exact product. */
    an = 0;
    if (xn + yn > ARF_MUL_SMALL_LIMBS && outn < xn + yn)
    {
        int cancel = 0;

        if ((sgnbit != ARF_SGNBIT(zm)) && shift >= -2 && shift <= 2)
        {
            ulong hi, lo, zt;
            slong s = -shift;
            umul_ppmm(hi, lo, xptr[xn - 1], yptr[yn - 1]);
            zt = zptr[zn - 1];
            if (s > 0) hi >>= s; else if (s < 0) zt >>= -s;
            hi = (hi > zt) ? hi - zt : zt - hi;
            cancel = (hi < (UWORD(1) << (FLINT_BITS - FLINT_MIN(30, FLINT_BITS - 2))));
        }

        if (!cancel)
            an = _arf_mul_approx(tmp, &e, &expfix, xptr, xn, yptr, yn, outn, xm == ym);
    }

    if (an == 0)
    {
        an = xn + yn;
        FLINT_MPN_MUL_WITH_SPECIAL_CASES(tmp, xptr, xn, yptr, yn);
        e = -1;
        expfix = 0;
    }

    /* value of the product = (A / B^an) 2^(pexp + expfix); z = (Z / B^zn) 2^zexp */
    aexp = pexp + expfix;
    shift = ARF_EXP(zm) - aexp;

    if (prec != ARF_PREC_EXACT && shift > prec + FLINT_BITS + 1 && shift > (zn + 1) * FLINT_BITS)
    {
        /* z +/- eps: the product lies more than FLINT_BITS bits below the
           precision (its deficit, if truncated, is covered by the extra
           unit as well) */
        _arb_set_round_mpn_2exp_si_eps(res, zptr, zn, ARF_SGNBIT(zm), ARF_EXP(zm) - zn * FLINT_BITS, sgnbit != ARF_SGNBIT(zm), prec, errp);
    }
    else if (prec != ARF_PREC_EXACT && -shift > prec + FLINT_BITS + 1 && -shift > (an + 1) * FLINT_BITS)
    {
        /* product +/- eps, likewise */
        _arb_set_round_mpn_2exp_si_eps(res, tmp, an, sgnbit, aexp - an * FLINT_BITS, sgnbit != ARF_SGNBIT(zm), prec, errp);
        if (e >= 0)
        {
            mag_t d;
            MAG_MAN(d) = MAG_ONE_HALF;
            MAG_EXP(d) = aexp - an * FLINT_BITS + e + 1;
            mag_fast_add(arb_radref(res), arb_radref(res), d);
        }
    }
    else
    {
        ulong sum_stack[ARF_ADD_STACK_ALLOC];
        slong frame_exp;

        /* the sum in the addition scratch (the product is in the
           multiplication scratch) */
        if (shift >= 0)
            salloc = _ARF_ADD_MPN_SIZE(zn, an, shift);
        else
            salloc = _ARF_ADD_MPN_SIZE(an, zn, -shift);
        sptr = (salloc <= ARF_ADD_STACK_ALLOC) ? sum_stack
             : (salloc <= ARF_ADD_TLS_ALLOC) ? _arf_add_tls_tmp(salloc) : flint_malloc(sizeof(ulong) * salloc);

        if (shift >= 0)
        {
            sn = _arf_add_mpn_raw(sptr, &sn_original, &ssgnbit, zptr, zn, ARF_SGNBIT(zm), tmp, an, sgnbit, shift);
            frame_exp = ARF_EXP(zm);
        }
        else
        {
            sn = _arf_add_mpn_raw(sptr, &sn_original, &ssgnbit, tmp, an, sgnbit, zptr, zn, ARF_SGNBIT(zm), -shift);
            frame_exp = aexp;
        }

        if (sn == 0)
        {
            arf_zero(arb_midref(res));
            if (errp == NULL)
                mag_zero(arb_radref(res));
            else
                mag_set(arb_radref(res), errp);
        }
        else if (e < 0)
        {
            /* value = (S / B^sn_original) 2^frame_exp */
            _arb_set_round_mpn_2exp_si_err(res, sptr, sn, ssgnbit, frame_exp - sn_original * FLINT_BITS, prec, 0, errp);
        }
        else
        {
            /* deficit of the truncated product: below 2^e units of A's
               ulp, added to the radius after the rounding error */
            mag_t d;
            _arb_set_round_mpn_2exp_si_err(res, sptr, sn, ssgnbit, frame_exp - sn_original * FLINT_BITS, prec, 0, errp);
            MAG_MAN(d) = MAG_ONE_HALF;
            MAG_EXP(d) = aexp - an * FLINT_BITS + e + 1;
            mag_fast_add(arb_radref(res), arb_radref(res), d);
        }

        ARF_ADD_TMP_FREE(sptr, salloc)
    }

    ARF_MUL_TMP_FREE(tmp, alloc)
}

void
arb_addmul(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    _arb_addmul(z, z, x, y, 0, prec);
}

void
arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    /* view y as an exact ball (shallow copy, never cleared) */
    arb_struct yb;
    *arb_midref(&yb) = *y;
    MAG_MAN(arb_radref(&yb)) = 0;
    MAG_EXP(arb_radref(&yb)) = 0;
    _arb_addmul(z, z, x, &yb, 0, prec);
}

void
arb_addmul_si(arb_t z, const arb_t x, slong y, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_addmul_arf(z, x, t, prec);
}

void
arb_addmul_ui(arb_t z, const arb_t x, ulong y, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_addmul_arf(z, x, t, prec);
}

void
arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_addmul_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_addmul_arf(z, x, t, prec);
        arf_clear(t);
    }
}
