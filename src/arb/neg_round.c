/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb/impl.h"

void
arb_neg_round(arb_t y, const arb_t x, slong prec)
{
    if (_arb_set_round_arf_err(y, arb_midref(x), 1, prec, arb_radref(x)))
        return;

    arb_set_round(y, x, prec);
    arb_neg(y, y);
}
