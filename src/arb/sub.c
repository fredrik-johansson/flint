/*
    Copyright (C) 2012-2014, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb/impl.h"

void
arb_sub(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    _arb_add_sub(z, x, y, 1, prec);
}

void
arb_sub_arf(arb_t z, const arb_t x, const arf_t y, slong prec)
{
    /* view y as an exact ball (shallow copy, never cleared) */
    arb_struct yb;
    *arb_midref(&yb) = *y;
    MAG_MAN(arb_radref(&yb)) = 0;
    MAG_EXP(arb_radref(&yb)) = 0;
    _arb_add_sub(z, x, &yb, 1, prec);
}

void
arb_sub_ui(arb_t z, const arb_t x, ulong y, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_sub_arf(z, x, t, prec);
}

void
arb_sub_si(arb_t z, const arb_t x, slong y, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_sub_arf(z, x, t, prec);
}

void
arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_sub_arf(z, x, t, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_sub_arf(z, x, t, prec);
        arf_clear(t);
    }
}
