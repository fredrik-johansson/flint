/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

slong
fmpz_poly_mat_find_pivot_any(const fmpz_poly_mat_t mat,
                                    slong start_row, slong end_row, slong c)
{
    slong r;

    for (r = start_row; r < end_row; r++)
    {
        if (!fmpz_poly_is_zero(fmpz_poly_mat_entry(mat, r, c)))
            return r;
    }

    return -1;
}
