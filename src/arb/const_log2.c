/*
    Copyright (C) 2013, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mp_real.h"

static void
arb_const_log2_hypgeom_eval(arb_t s, slong prec)
{
    /* Zuniga's series, 11.9 bits per term (fixed/const_log2.c) */
    _mp_real_const_arb(s, mp_real_const_log2, prec);
}

_ARB_DEF_CACHED_CONSTANT(static, arb_const_log2_hypgeom, arb_const_log2_hypgeom_eval)

void
arb_const_log2(arb_t res, slong prec)
{
    if (prec < ARB_LOG_TAB2_LIMBS * FLINT_BITS - 16)
    {
        slong exp;

        /* the table is a truncation (error below one unit of its last
           limb): round it to nearest with the error measured, then add the
           table's own error */
        arb_set_round_mpn_2exp_si(res, arb_log_log2_tab, ARB_LOG_TAB2_LIMBS, 0, 0 - ARB_LOG_TAB2_LIMBS * FLINT_BITS, prec);
        mag_fast_add_2exp_si(arb_radref(res), arb_radref(res), 0 - ARB_LOG_TAB2_LIMBS * FLINT_BITS);
        (void) exp;
    }
    else
    {
        arb_const_log2_hypgeom(res, prec);
    }
}
