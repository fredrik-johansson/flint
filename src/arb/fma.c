/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb/impl.h"

/*
    res = x y + z: the one-pass fused multiply-add of arb_addmul (the
    product added exactly to z and rounded once to nearest, with the
    measured rounding error), which already takes a separate output.
*/
void
arb_fma_arf(arb_t res, const arb_t x, const arf_t y, const arb_t z, slong prec)
{
    /* view y as an exact ball (shallow copy, never cleared) */
    arb_struct yb;
    *arb_midref(&yb) = *y;
    MAG_MAN(arb_radref(&yb)) = 0;
    MAG_EXP(arb_radref(&yb)) = 0;
    _arb_addmul(res, z, x, &yb, 0, prec);
}

void
arb_fma(arb_t res, const arb_t x, const arb_t y, const arb_t z, slong prec)
{
    _arb_addmul(res, z, x, y, 0, prec);
}

void
arb_fma_ui(arb_t res, const arb_t x, ulong y, const arb_t z, slong prec)
{
    arf_t t;
    arf_init_set_ui(t, y); /* no need to free */
    arb_fma_arf(res, x, t, z, prec);
}

void
arb_fma_si(arb_t res, const arb_t x, slong y, const arb_t z, slong prec)
{
    arf_t t;
    arf_init_set_si(t, y); /* no need to free */
    arb_fma_arf(res, x, t, z, prec);
}

void
arb_fma_fmpz(arb_t res, const arb_t x, const fmpz_t y, const arb_t z, slong prec)
{
    arf_t t;

    if (!COEFF_IS_MPZ(*y))
    {
        arf_init_set_si(t, *y); /* no need to free */
        arb_fma_arf(res, x, t, z, prec);
    }
    else
    {
        arf_init(t);
        arf_set_fmpz(t, y);
        arb_fma_arf(res, x, t, z, prec);
        arf_clear(t);
    }
}
