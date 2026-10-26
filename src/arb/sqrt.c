/*
    Copyright (C) 2014, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb/impl.h"
#include "mpn_extras/fp_small.h"

/* Roots of up to this many limbs plus guard limb are formed on the stack. */
#define ARB_SQRT_STACK_LIMBS 64

/*
    The square root or reciprocal square root of the positive LAGOM
    midpoint x with a guard limb: n + 1 limbs (n >= 1, top bit set)
    written to q, with the value q 2^e for the returned e, within one
    unit of the exact value (*extra = 1) or exact (*extra = 0):
    flint_mpn_sqrthigh / rsqrthigh, or the exact truncation when the
    guard limb cannot rule out an exact result (which is then found
    exactly). A power of two with an odd exponent is exact as a power of
    two. If u != NULL, it is set to an upper bound for the value.
*/
static slong
_arb_sqrt_root(nn_ptr q, int * extra, mag_ptr u, arf_srcptr xm, slong n, int rsqrt)
{
    nn_srcptr xptr;
    slong xn, e, exp;
    int t;

    ARF_GET_MPN_READONLY(xptr, xn, xm);
    e = ARF_EXP(xm);
    t = e & 1;

    if (xn == 1 && xptr[0] == LIMB_TOP && t)
    {
        /* x = 2^(e-1), e - 1 even: sqrt(x) = 2^((e-1)/2), 1/sqrt(x) = 2^(-(e-1)/2) */
        flint_mpn_zero(q, n);
        q[n] = LIMB_TOP;
        *extra = 0;
        exp = rsqrt ? 1 - (e - 1) / 2 : 1 + (e - 1) / 2;
    }
#if FLINT_BITS == 64
    else if (!rsqrt && xn == 1 && n == 1)
    {
        /* the exact floor s of sqrt(a B / 2^t) with its remainder r <= 2 s:
           the fraction f = sqrt - s lies in (r / (2 s + 2), r / (2 s)], so
           r / (2 y) for the double precision root y is within 2^13 units
           of the guard limb, far inside the 2^(FLINT_BITS - MAG_BITS)
           units that *extra = 1 allows */
        ulong r1, r0;
        double f, inv;

        if (_flint_mpn_isqrt_2_rem_inv(q + 1, &r1, &r0, &inv, xptr[0] >> t, (xptr[0] << (FLINT_BITS - 1)) & (-(ulong) t)))
        {
            f = ((double) r1 * 0x1p64 + (double) r0) * inv;
            q[0] = (f >= 1.0) ? UWORD_MAX : (ulong) (f * 0x1p64);
            *extra = 1;
        }
        else
        {
            q[0] = 0;
            *extra = 0;
        }
        exp = (e + t) / 2;
    }
#endif
    else if (!rsqrt)
    {
        flint_mpn_sqrthigh(q, xptr, xn, t, n);
        *extra = 1;
        if (q[0] <= FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND || q[0] >= UWORD_MAX - FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND)
            *extra = flint_mpn_sqrthigh_exact(q, xptr, xn, t, n);
        exp = (e + t) / 2;
    }
    else
    {
        flint_mpn_rsqrthigh(q, xptr, xn, t, n);
        *extra = 1;
        if (q[0] <= FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND || q[0] >= UWORD_MAX - FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND)
            *extra = flint_mpn_rsqrthigh_exact(q, xptr, xn, t, n);
        exp = (1 - t) - (e - t) / 2;
    }

    if (u != NULL)
    {
        MAG_MAN(u) = (q[n] >> (FLINT_BITS - MAG_BITS)) + 1;
        MAG_EXP(u) = exp;
        MAG_FAST_ADJUST_ONE_TOO_LARGE(u);
    }

    return exp - (n + 1) * FLINT_BITS;
}

/* res = sqrt(x) or 1/sqrt(x) for a positive LAGOM x, radius err (NULL for
   none) plus the measured rounding error; prec >= 2 */
static void
_arb_sqrt_mid(arb_t res, arf_srcptr xm, slong prec, mag_srcptr err, int rsqrt)
{
    ulong q_stack[ARB_SQRT_STACK_LIMBS + 1];
    nn_ptr q;
    slong n, exp;
    int extra;

    n = (prec + FLINT_BITS - 1) / FLINT_BITS + ARB_GUARD_EXTRA;
    q = (n <= ARB_SQRT_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));
    exp = _arb_sqrt_root(q, &extra, NULL, xm, n, rsqrt);
    _arb_set_round_guarded(res, q, n, 0, exp, prec, extra, err, NULL);
    if (q != q_stack)
        flint_free(q);
}

static void
arb_rsqrt_arf(arb_t res, const arf_t x, slong prec)
{
    if (arf_is_special(x) || arf_sgn(x) < 0)
    {
        if (arf_is_pos_inf(x))
            arb_zero(res);
        else
            arb_indeterminate(res);
        return;
    }

    if (ARF_IS_LAGOM(x) && prec >= 2 && prec != ARF_PREC_EXACT)
    {
        _arb_sqrt_mid(res, x, prec, NULL, 1);
    }
    else
    {
        int inexact;

        inexact = arf_rsqrt(arb_midref(res), x, prec, ARF_RND_FAST);

        if (inexact)
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
        else
            mag_zero(arb_radref(res));
    }
}

void
arb_rsqrt(arb_t z, const arb_t x, slong prec)
{
    if (!arb_is_finite(x) || arf_sgn(arb_midref(x)) <= 0)
    {
        if (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_zero(z);
        else
            arb_indeterminate(z);
    }
    else if (mag_is_zero(arb_radref(x)))
    {
        arb_rsqrt_arf(z, arb_midref(x), prec);
    }
    else
    {
        mag_t t, u;
        slong acc;

        mag_init(t);

        arb_get_mag_lower(t, x);

        acc = arb_rel_accuracy_bits(x);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc <= 20)
        {
            if (mag_is_zero(t))
            {
                arb_indeterminate(z);
            }
            else
            {
                mag_init(u);
                arb_get_mag(u, x);

                mag_rsqrt(t, t);
                mag_rsqrt_lower(u, u);

                arb_set_interval_mag(z, u, t, prec);

                mag_clear(u);
            }
        }
        else
        {
            /* error bound: (1/2) (x-r)^(-3/2) * r */
            mag_init(u);

            mag_rsqrt(u, t);
            mag_div(t, u, t);
            mag_mul(t, t, arb_radref(x));
            mag_mul_2exp_si(t, t, -1);

            if (ARF_IS_LAGOM(arb_midref(x)) && MAG_IS_LAGOM(t))
            {
                _arb_sqrt_mid(z, arb_midref(x), prec, t, 1);
            }
            else
            {
                arb_rsqrt_arf(z, arb_midref(x), prec);
                mag_add(arb_radref(z), arb_radref(z), t);
            }

            mag_clear(u);
        }

        mag_clear(t);
    }
}



void
arb_sqrt_arf(arb_t res, const arf_t x, slong prec)
{
    if (arf_is_special(x) || arf_sgn(x) < 0)
    {
        if (arf_is_zero(x) || arf_is_pos_inf(x))
            arb_set_arf(res, x);
        else
            arb_indeterminate(res);
        return;
    }

    if (ARF_IS_LAGOM(x) && prec >= 2 && prec != ARF_PREC_EXACT)
    {
        _arb_sqrt_mid(res, x, prec, NULL, 0);
    }
    else
    {
        if (arf_sqrt(arb_midref(res), x, prec, ARF_RND_FAST))
            arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
        else
            mag_zero(arb_radref(res));
    }
}


void
arb_sqrt(arb_t z, const arb_t x, slong prec)
{
    mag_t rx, zr;
    int inexact;

    if (mag_is_zero(arb_radref(x)))
    {
        arb_sqrt_arf(z, arb_midref(x), prec);
    }
    else if (arf_is_special(arb_midref(x)) ||
              arf_sgn(arb_midref(x)) < 0 || mag_is_inf(arb_radref(x)))
    {
        if (arf_is_pos_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
            arb_sqrt_arf(z, arb_midref(x), prec);
        else
            arb_indeterminate(z);
    }
    else  /* now both mid and rad are non-special values, mid > 0 */
    {
        slong acc;

        acc = _fmpz_sub_small(ARF_EXPREF(arb_midref(x)), MAG_EXPREF(arb_radref(x)));
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc < 0)
        {
            arb_indeterminate(z);
        }
        else if (acc <= 20)
        {
            mag_t t, u;

            mag_init(t);
            mag_init(u);

            arb_get_mag_lower(t, x);

            if (mag_is_zero(t) && arb_contains_negative(x))
            {
                arb_indeterminate(z);
            }
            else
            {
                arb_get_mag(u, x);
                mag_sqrt_lower(t, t);
                mag_sqrt(u, u);
                arb_set_interval_mag(z, t, u, prec);
            }

            mag_clear(t);
            mag_clear(u);
        }
        else if (ARB_IS_LAGOM(x)) /* small exponents, acc *and* prec >= 20 */
        {
            mag_t t, err;
            ulong q_stack[ARB_SQRT_STACK_LIMBS + 1];
            nn_ptr q;
            slong n, exp;
            int extra;

            mag_init(t); /* no need to free */
            mag_init(err);

            /* propagated error sqrt(x) - sqrt(x - r) <= r / (2 sqrt(x - r)),
               from the input alone so that it does not wait for the root:
               with acc >= 20, r < 2^-19 x, so x (1 - 2^-19) is a lower
               bound for x - r (t is a lower bound of x, so the subtraction
               keeps it one) */
            {
                nn_srcptr xptr;
                slong xn;
                ARF_GET_MPN_READONLY(xptr, xn, arb_midref(x));
                MAG_MAN(t) = (xptr[xn - 1] >> (FLINT_BITS - MAG_BITS)) - (xptr[xn - 1] >> (FLINT_BITS - MAG_BITS + 19)) - 1;
                MAG_EXP(t) = ARF_EXP(arb_midref(x));
            }
            mag_rsqrt(t, t);
            mag_fast_mul(err, arb_radref(x), t);
            MAG_EXP(err)--;

            /* the root with a guard limb */
            n = (prec + FLINT_BITS - 1) / FLINT_BITS + ARB_GUARD_EXTRA;
            q = (n <= ARB_SQRT_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));
            exp = _arb_sqrt_root(q, &extra, NULL, arb_midref(x), n, 0);

            _arb_set_round_guarded(z, q, n, 0, exp, prec, extra, err, NULL);

            if (q != q_stack)
                flint_free(q);
        }
        else
        {
            mag_init(zr);
            mag_init(rx);

            /* rx = upper bound for r / x */
            arf_get_mag_lower(rx, arb_midref(x));
            mag_div(rx, arb_radref(x), rx);

            inexact = arf_sqrt(arb_midref(z), arb_midref(x), prec, ARF_RND_FAST);

            /* zr = upper bound for sqrt(x) */
            arf_get_mag(zr, arb_midref(z));
            if (inexact)
                arf_mag_add_ulp(zr, zr, arb_midref(z), prec);

            /* propagated error:   sqrt(x) - sqrt(x-r)
                                 = sqrt(x) * [1 - sqrt(1 - r/x)]
                                <= sqrt(x) * 0.5 * (rx + rx^2)  */
            mag_addmul(rx, rx, rx);
            mag_mul(zr, zr, rx);
            mag_mul_2exp_si(zr, zr, -1);

            /* add the rounding error */
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), zr, arb_midref(z), prec);
            else
                mag_swap(arb_radref(z), zr);

            mag_clear(zr);
            mag_clear(rx);
        }
    }
}

void
arb_rsqrt_ui(arb_t z, ulong x, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    arb_rsqrt_arf(z, t, prec);
}

void
arb_sqrt_ui(arb_t z, ulong x, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, x); /* no need to free */
    arb_sqrt_arf(z, t, prec);
}

void
arb_sqrt_fmpz(arb_t z, const fmpz_t x, slong prec)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    arb_sqrt_arf(z, t, prec);
    arf_clear(t);
}
