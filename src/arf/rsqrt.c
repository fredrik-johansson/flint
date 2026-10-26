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

#define ARF_RSQRT_STACK_LIMBS 64

/*
    1/sqrt(x) = M 2^((1 - t) - (e - t) / 2) for x = A 2^e, t = e mod 2 and
    M = (2^(2-t) A)^(-1/2) in (1/2, 1), from flint_mpn_rsqrthigh (n limbs
    and a guard limb, within FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND units)
    when the guard bits certify the rounding (always for a relaxed mode),
    else from the exact truncated value with its remainder
    (flint_mpn_rsqrthigh_exact). The result is dyadic only for x a power
    of four, handled separately, so it is always inexact otherwise.
*/
int
arf_rsqrt(arf_ptr z, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    slong xn, n, fix;
    nn_srcptr xptr;
    ulong q_stack[ARF_RSQRT_STACK_LIMBS + 1];
    nn_ptr q;
    int t, inexact;
    fmpz_t e;

    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_pos_inf(z);
        else if (arf_is_pos_inf(x))
            arf_zero(z);
        else
            arf_nan(z);
        return 0;
    }

    if (ARF_SGNBIT(x))
    {
        arf_nan(z);
        return 0;
    }

    if (ARF_IS_POW2(x) && fmpz_is_odd(ARF_EXPREF(x)))
    {
        arf_set(z, x);
        fmpz_neg(ARF_EXPREF(z), ARF_EXPREF(z));
        fmpz_cdiv_q_2exp(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        fmpz_add_ui(ARF_EXPREF(z), ARF_EXPREF(z), 1);
        return 0;
    }

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_rsqrt: exact precision not supported\n");

    ARF_GET_MPN_READONLY(xptr, xn, x);
    t = fmpz_is_odd(ARF_EXPREF(x));

    n = (prec + FLINT_BITS - 1) / FLINT_BITS;
    q = (n <= ARF_RSQRT_STACK_LIMBS) ? q_stack : flint_malloc((n + 1) * sizeof(ulong));

    /* exponent (1 - t) - (e - t) / 2 (e - t is even) */
    fmpz_init(e);
    if (!COEFF_IS_MPZ(ARF_EXP(x)))
        *e = (1 - t) - ((ARF_EXP(x) - t) >> 1);
    else
    {
        fmpz_sub_ui(e, ARF_EXPREF(x), t);
        fmpz_tdiv_q_2exp(e, e, 1);
        fmpz_neg(e, e);
        fmpz_add_ui(e, e, 1 - t);
    }

    flint_mpn_rsqrthigh(q, xptr, xn, t, n);

    if (_arf_set_round_approx(z, q, n, FLINT_BIT_COUNT(FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND), e, 0, prec, rnd))
    {
        inexact = 1;
    }
    else
    {
        inexact = flint_mpn_rsqrthigh_exact(q, xptr, xn, t, n);
        FLINT_ASSERT(inexact);
        q[0] |= 1;
        _arf_set_round_mpn(z, &fix, q, n + 1, 0, prec, rnd);
        _fmpz_add_fast(ARF_EXPREF(z), e, fix);
        inexact = 1;
    }

    fmpz_clear(e);
    if (q != q_stack)
        flint_free(q);
    return inexact;
}
