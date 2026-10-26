/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_IMPL_H
#define ACB_IMPL_H

#include "acb_types.h"

void acb_gamma_stirling_eval(acb_t s, const acb_t z, slong nterms, int digamma, slong prec);

/* res = z + (-1)^neg x y (z = NULL: res = x y) in one pass per part for LAGOM balls
   with nonzero finite midpoints and finite radii; returns 0 (without
   touching res) when the operands are not of that shape */
int _acb_mul_onepass(acb_t res, const acb_t z, const acb_t x, const acb_t y, int neg, slong prec);

#endif
