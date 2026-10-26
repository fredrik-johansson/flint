/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "mpn_extras.h"

int
arf_sosq(arf_t res, const arf_t a, const arf_t b, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(a) || arf_is_special(b))
    {
        if (arf_is_zero(a))
            return arf_mul(res, b, b, prec, rnd);

        if (arf_is_zero(b))
            return arf_mul(res, a, a, prec, rnd);

        if (arf_is_nan(a) || arf_is_nan(b))
            arf_nan(res);
        else
            arf_pos_inf(res);

        return 0;
    }
    else if (a == b)
    {
        /* 2 a^2: the squaring, then one more in the exponent */
        int inexact = arf_mul(res, a, a, prec, rnd);
        fmpz_add_ui(ARF_EXPREF(res), ARF_EXPREF(res), 1);
        return inexact;
    }
    else
    {
        return _arf_fmma_mpn(res, a, a, 0, b, b, 0, prec, rnd);
    }
}
