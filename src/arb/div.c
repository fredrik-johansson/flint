/*
    Copyright (C) 2014, 2018, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb/impl.h"

/* Quotients of up to this many limbs plus guard limb are formed on the
   stack. */
#define ARB_DIV_STACK_LIMBS 64

/*
    The quotient of the finite nonzero LAGOM midpoints x and y with a
    guard limb: n + 1 limbs (n >= 1, top bit set) written to q, with the
    value q 2^e for the returned e, within one unit of the exact quotient
    (*extra = 1) or exact (*extra = 0): flint_mpn_divhigh, or the exact
    quotient when the guard limb cannot rule out an exact result (which
    is then found exactly); a divisor that is a power of two gives the
    top limbs of x. If u != NULL, it is set to an upper bound for |x / y|.
*/
static slong
_arb_div_quotient(nn_ptr q, int * extra, mag_ptr u, arf_srcptr xm, arf_srcptr ym, slong n)
{
    nn_srcptr xptr, yptr;
    slong xn, yn, exp;
    ulong s;

    ARF_GET_MPN_READONLY(xptr, xn, xm);
    ARF_GET_MPN_READONLY(yptr, yn, ym);

    if (yn == 1 && yptr[0] == LIMB_TOP)
    {
        /* x / 2^k: the top n + 1 limbs of x */
        if (xn <= n + 1)
        {
            flint_mpn_zero(q, n + 1 - xn);
            flint_mpn_copyi(q + n + 1 - xn, xptr, xn);
            *extra = 0;
        }
        else
        {
            flint_mpn_copyi(q, xptr + xn - n - 1, n + 1);
            *extra = !flint_mpn_zero_p(xptr, xn - n - 1);
        }
        exp = ARF_EXP(xm) - ARF_EXP(ym) + 1;
    }
    else
    {
        s = flint_mpn_divhigh(q, xptr, xn, yptr, yn, n);
        *extra = 1;

        if (q[0] <= FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND || q[0] >= UWORD_MAX - FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND)
        {
            /* the exact truncation: a lower approximation with deficit
               below one unit unless exact */
            *extra = flint_mpn_divhigh_exact(q, &s, xptr, xn, yptr, yn, n);
        }

        exp = ARF_EXP(xm) - ARF_EXP(ym) + (slong) s;
    }

    if (u != NULL)
    {
        MAG_MAN(u) = (q[n] >> (FLINT_BITS - MAG_BITS)) + 1;
        MAG_EXP(u) = exp;
        MAG_FAST_ADJUST_ONE_TOO_LARGE(u);
    }

    return exp - (n + 1) * FLINT_BITS;
}

/* res = x / y for LAGOM x, y (y nonzero), radius err plus the measured
   rounding error; prec >= 2 (see arb/impl.h) */
void
_arb_div_mid(arb_t res, arf_srcptr xm, arf_srcptr ym, slong prec, mag_srcptr err)
{
    ulong q_stack[ARB_DIV_STACK_LIMBS + 1];
    nn_ptr q;
    slong n, exp;
    int extra;

    n = (prec + FLINT_BITS - 1) / FLINT_BITS + ARB_GUARD_EXTRA;
    q = (n <= ARB_DIV_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));
    exp = _arb_div_quotient(q, &extra, NULL, xm, ym, n);
    _arb_set_round_guarded(res, q, n, ARF_SGNBIT(xm) ^ ARF_SGNBIT(ym), exp, prec, extra, err, NULL);
    if (q != q_stack)
        flint_free(q);
}

void
arb_div_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    mag_t zr, ym;
    int inexact;

    if (arf_is_special(y) || !arb_is_finite(x))
    {
        if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)) &&
            arf_is_finite(y) && !arf_is_zero(y))
        {
            /* inf / finite nonzero = inf */
            arf_div(arb_midref(z), arb_midref(x), y, prec, ARF_RND_DOWN);
            mag_zero(arb_radref(z));
        }
        else if (arb_is_finite(x) && arf_is_inf(y))
        {
            /* finite / inf = 0 */
            arb_zero(z);
        }
        else if (!arf_is_nan(arb_midref(x)) && mag_is_inf(arb_radref(x)) &&
            arf_is_finite(y) && !arf_is_zero(y))
        {
            /* [+/- inf] / finite nonzero = [+/- inf] */
            arb_zero_pm_inf(z);
        }
        else
        {
            arb_indeterminate(z);
        }
    }
    else if (mag_is_zero(arb_radref(x)))
    {
        if (!arf_is_zero(arb_midref(x)) && ARF_IS_LAGOM(arb_midref(x)) && ARF_IS_LAGOM(y) && prec >= 2 && prec != ARF_PREC_EXACT)
        {
            _arb_div_mid(z, arb_midref(x), y, prec, NULL);
            return;
        }

        inexact = arf_div(arb_midref(z), arb_midref(x), y, prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        else
            mag_zero(arb_radref(z));
    }
    else if (arf_is_zero(arb_midref(x)))
    {
        /* [0 +/- eps] / y  ->  [0 +/- eps / y] */
        mag_init(ym);
        arf_get_mag_lower(ym, y);
        mag_div(arb_radref(z), arb_radref(x), ym);
        arf_zero(arb_midref(z));
        mag_clear(ym);
    }
    else
    {
        mag_init(ym);
        mag_init(zr);

        arf_get_mag_lower(ym, y);
        mag_div(zr, arb_radref(x), ym);

        if (ARB_IS_LAGOM(x) && ARF_IS_LAGOM(y) && MAG_IS_LAGOM(zr) && prec >= 2 && prec != ARF_PREC_EXACT)
        {
            /* the precision clamped to the accuracy of x */
            slong acc = ARF_EXP(arb_midref(x)) - MAG_EXP(arb_radref(x));
            prec = FLINT_MAX(2, FLINT_MIN(prec, acc + MAG_BITS));
            _arb_div_mid(z, arb_midref(x), y, prec, zr);
            mag_clear(ym);
            mag_clear(zr);
            return;
        }

        inexact = arf_div(arb_midref(z), arb_midref(x), y, prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
        else
            mag_swap(arb_radref(z), zr);

        mag_clear(ym);
        mag_clear(zr);
    }
}

static void
arb_div_wide(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    mag_t a, b, t, u;

    mag_init(t);
    arb_get_mag_lower(t, y);

    if (mag_is_zero(t))
    {
        arb_indeterminate(z);
    }
    else if (arf_is_zero(arb_midref(x)))
    {
        mag_div(arb_radref(z), arb_radref(x), t);
        arf_zero(arb_midref(z));
    }
    else
    {
        if (arf_cmpabs_mag(arb_midref(x), arb_radref(x)) >= 0)
        {
            /* [a,b] /  [t,u] =  [a/u, b/t]
               [a,b] / -[t,u] = -[a/u, b/t]
              -[a,b] /  [t,u] = -[a/u, b/t]
              -[a,b] / -[t,u] =  [a/u, b/t] */
            mag_init(a);
            mag_init(b);
            mag_init(u);

            arb_get_mag_lower(a, x);
            arb_get_mag(b, x);
            arb_get_mag(u, y);

            mag_div_lower(a, a, u);
            mag_div(b, b, t);

            if ((arf_sgn(arb_midref(x)) < 0) ^ (arf_sgn(arb_midref(y)) < 0))
            {
                arb_set_interval_mag(z, a, b, prec);
                arb_neg(z, z);
            }
            else
            {
                arb_set_interval_mag(z, a, b, prec);
            }

            mag_clear(a);
            mag_clear(b);
            mag_clear(u);
        }
        else
        {
            /* [-a,b] /  [t,u] = [-a/t, b/t]
               [-a,b] / -[t,u] = [-b/t, a/t] */
            mag_init(a);
            mag_init(b);

            arb_get_mag(b, x);
            arf_get_mag_lower(a, arb_midref(x));
            mag_sub(a, arb_radref(x), a);

            if ((arf_sgn(arb_midref(x)) < 0) ^ (arf_sgn(arb_midref(y)) < 0))
                mag_swap(a, b);

            mag_div(a, a, t);
            mag_div(b, b, t);
            arb_set_interval_neg_pos_mag(z, a, b, prec);

            mag_clear(a);
            mag_clear(b);
        }
    }

    mag_clear(t);
}

void
arb_div(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    int inexact;
    slong acc, xacc, yacc;

    if (mag_is_zero(arb_radref(y)))
    {
        arb_div_arf(z, x, arb_midref(y), prec);
    }
    else if (arf_is_zero(arb_midref(y))) /* anything / 0 = nan */
    {
        arb_indeterminate(z);
    }
    else if (arf_is_zero(arb_midref(x)) && arb_is_finite(y))
    {
        arb_div_wide(z, x, y, prec);
    }
    else if (ARB_IS_LAGOM(x) && ARB_IS_LAGOM(y)) /* fast case */
    {
        /* both finite, both midpoint and radius of y not special */
        yacc = MAG_EXP(arb_midref(y)) - ARF_EXP(arb_radref(y));
        xacc = arb_rel_accuracy_bits(x);
        acc = FLINT_MIN(xacc, yacc);
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc <= 20)
        {
            arb_div_wide(z, x, y, prec);
        }
        else
        {
            mag_t t, u, v, err;
            ulong q_stack[ARB_DIV_STACK_LIMBS + 1];
            nn_ptr q;
            slong n, exp;
            int extra;

            mag_init(t);  /* no need to free */
            mag_init(u);
            mag_init(v);
            mag_init(err);

            /*               (x*yrad + y*xrad)/(y*(y-yrad))
                 <=  (1+eps) (x*yrad + y*xrad)/y^2
                 <=  (1+eps) (x*yrad/y^2 + y*xrad/y^2)
                 <=  (1+eps) ((x/y)*yrad/y + xrad/y)
                 <=  (1+eps) ((x/y)*yrad + xrad)/y
            */

            /* t = y */
            arf_get_mag_lower(t, arb_midref(y));

            /* the quotient with a guard limb, and u = x/y from it */
            n = (prec + FLINT_BITS - 1) / FLINT_BITS + ARB_GUARD_EXTRA;
            q = (n <= ARB_DIV_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));
            exp = _arb_div_quotient(q, &extra, u, arb_midref(x), arb_midref(y), n);

            /* v = (x/y)*yrad + xrad */
            MAG_MAN(v) = MAG_MAN(arb_radref(x));
            MAG_EXP(v) = MAG_EXP(arb_radref(x));
            mag_fast_addmul(v, arb_radref(y), u);
            /* v = v / y */
            mag_div(err, v, t);

            /* correct for errors */
            MAG_MAN(t) = MAG_ONE_HALF + (MAG_ONE_HALF >> 16);
            MAG_EXP(t) = 1;
            mag_fast_mul(err, err, t);

            _arb_set_round_guarded(z, q, n, ARF_SGNBIT(arb_midref(x)) ^ ARF_SGNBIT(arb_midref(y)), exp, prec, extra, err, NULL);

            if (q != q_stack)
                flint_free(q);
        }
    }
    else if (!arb_is_finite(y))
    {
        /* finite / inf = 0 */
        if (arf_is_inf(arb_midref(y)) && mag_is_finite(arb_radref(y)) &&
            arb_is_finite(x))
            arb_zero(z);
        else
            arb_indeterminate(z);
    }
    else if (!arb_is_finite(x))
    {
        if (arb_contains_zero(y) || arf_is_nan(arb_midref(x)))
        {
            arb_indeterminate(z);
        }
        else
        {
            /* +/- inf  /  finite nonzero  = +/- inf */
            if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            {
                arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_DOWN);
                mag_zero(arb_radref(z));
            }
            else if (!arf_is_nan(arb_midref(x)) && mag_is_inf(arb_radref(x)))
            {
                arb_zero_pm_inf(z);
            }
            else
            {
                arb_indeterminate(z);
            }
        }
    }
    else
    {
        yacc = arb_rel_accuracy_bits(y);
        xacc = arb_rel_accuracy_bits(x);
        acc = FLINT_MIN(xacc, yacc);
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc <= 20)
        {
            arb_div_wide(z, x, y, prec);
        }
        else
        {
            mag_t zr, xm, ym, yl, yw;

            mag_init_set_arf(xm, arb_midref(x));
            mag_init_set_arf(ym, arb_midref(y));
            mag_init(zr);
            mag_init(yl);
            mag_init(yw);

            /* (|x|*yrad + |y|*xrad)/(|y|*(|y|-yrad)) */
            mag_mul(zr, xm, arb_radref(y));
            mag_addmul(zr, ym, arb_radref(x));
            arb_get_mag_lower(yw, y);

            arf_get_mag_lower(yl, arb_midref(y));
            mag_mul_lower(yl, yl, yw);

            mag_div(zr, zr, yl);

            inexact = arf_div(arb_midref(z), arb_midref(x), arb_midref(y), prec, ARF_RND_FAST);

            if (inexact)
                arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
            else
                mag_swap(arb_radref(z), zr);

            mag_clear(xm);
            mag_clear(ym);
            mag_clear(zr);
            mag_clear(yl);
            mag_clear(yw);
        }
    }
}

void
arb_div_si(arb_t z, const arb_t x, slong y, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_ui(arb_t z, const arb_t x, ulong y, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_div_arf(z, x, t, prec);
}

void
arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_div_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_div_arf(z, x, t, prec);
        arf_clear(t);
    }
}

/* an exact arf view of a nonzero fmpz: shallow for a small one, else a
   copy of the limbs (to be cleared) */
static int
_arb_fmpz_view(arf_t t, const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        arf_init_set_si(t, *x);   /* no need to free */
        return 0;
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, x);
        return 1;
    }
}

/*
    z = x / y rounded to nearest with the measured error, through the
    guard-limb quotient of arb_div (exact quotients, e.g. 3 / 4, get
    radius 0). This is also arb_set_fmpq.
*/
void
arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, slong prec)
{
    int inexact;

    if (!fmpz_is_zero(x) && !fmpz_is_zero(y) && prec >= 2 && prec != ARF_PREC_EXACT &&
        fmpz_size(x) < ARF_MAX_LAGOM_EXP / (4 * FLINT_BITS) &&
        fmpz_size(y) < ARF_MAX_LAGOM_EXP / (4 * FLINT_BITS))
    {
        arf_t t, u;
        int ct, cu;

        ct = _arb_fmpz_view(t, x);
        cu = _arb_fmpz_view(u, y);
        _arb_div_mid(z, t, u, prec, NULL);
        if (ct)
            arf_clear(t);
        if (cu)
            arf_clear(u);
        return;
    }

    inexact = arf_fmpz_div_fmpz(arb_midref(z), x, y, prec, ARF_RND_FAST);

    if (inexact)
        arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
    else
        mag_zero(arb_radref(z));
}

void
arb_ui_div(arb_t z, ulong x, const arb_t y, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_ui(t, x);
    arb_div(z, t, y, prec);
    arb_clear(t);
}

