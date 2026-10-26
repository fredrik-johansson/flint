/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "arb.h"
#include "arb/impl.h"

void
arb_set_round_mpn_2exp_si(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, slong prec)
{
    _arb_set_round_mpn_2exp_si_core(res, x, xn, sgnbit, exp, prec, 0, 0, NULL, NULL);
}

void
_arb_set_round_mpn_2exp_si_err(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, slong prec, int extra, mag_srcptr err)
{
    _arb_set_round_mpn_2exp_si_core(res, x, xn, sgnbit, exp, prec, 0, extra, err, NULL);
}

void
_arb_set_round_mpn_2exp_si_inplace_err(arb_t res, slong xn, int sgnbit, slong exp, slong prec, int extra, mag_srcptr err)
{
    _arb_set_round_mpn_2exp_si_core(res, ARF_PTR_D(arb_midref(res)), xn, sgnbit, exp, prec, 1, extra, err, NULL);
}

/* As above, with x the (xn-limb) buffer of the midpoint of res itself,
   obtained from ARF_GET_MPN_WRITE(x, xn, arb_midref(res)); the rounded
   mantissa replaces it. */
void
_arb_set_round_mpn_2exp_si_inplace(arb_t res, slong xn, int sgnbit, slong exp, slong prec)
{
    _arb_set_round_mpn_2exp_si_core(res, ARF_PTR_D(arb_midref(res)), xn, sgnbit, exp, prec, 1, 0, NULL, NULL);
}

/*
    res = {x, xn} 2^exp +/- eps rounded to nearest at prec bits with the
    rounding error as radius (plus err), for an unknown eps with
    0 < eps < 2^(exp + xn FLINT_BITS - prec - FLINT_BITS - 1), i.e. more
    than FLINT_BITS bits below the granularity of the measured rounding
    error: a unit added to (sub = 0) or subtracted from (sub = 1) x padded
    with zero limbs to beyond the precision stands in for eps, and one
    extra unit of 2^(e - prec - MAG_BITS) in the radius covers the
    difference. The value must be nonzero, with at most one leading zero
    bit in x.
*/
void
_arb_set_round_mpn_2exp_si_eps(arb_t res, nn_srcptr x, slong xn, int sgnbit, slong exp, int sub, slong prec, mag_srcptr err)
{
    ulong tmp_stack[ARF_ADD_STACK_ALLOC];
    nn_ptr tmp;
    slong zn;

    FLINT_ASSERT(prec != ARF_PREC_EXACT);

    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    zn = FLINT_MAX(zn, xn) + 2;
    tmp = (zn <= ARF_ADD_STACK_ALLOC) ? tmp_stack
        : (zn <= ARF_ADD_TLS_ALLOC) ? _arf_add_tls_tmp(zn) : flint_malloc(sizeof(ulong) * zn);

    flint_mpn_zero(tmp, zn - xn);
    flint_mpn_copyi(tmp + zn - xn, x, xn);
    if (sub)
        mpn_sub_1(tmp, tmp, zn, 1);
    else
        tmp[0] = 1;

    _arb_set_round_mpn_2exp_si_core(res, tmp, zn, sgnbit, exp - (zn - xn) * FLINT_BITS, prec, 0, 1, err, NULL);

    if (zn > ARF_ADD_TLS_ALLOC)
        flint_free(tmp);
}
