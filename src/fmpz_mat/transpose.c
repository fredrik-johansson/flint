/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"

void
fmpz_mat_transpose(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong i, j;

    if (B->r != A->c || B->c != A->r)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_transpose). Incompatible dimensions.\n");
    }

    if (A == B)  /* In-place, guaranteed to be square */
    {
        for (i = 0; i < A->r - 1; i++)
            for (j = i + 1; j < A->c; j++)
                FLINT_SWAP(fmpz, *fmpz_mat_entry(B, i, j), *fmpz_mat_entry(B, j, i));
    }
    else  /* Not aliased; general case */
    {
        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                fmpz_set(fmpz_mat_entry(B, i, j), fmpz_mat_entry(A, j, i));
    }
}
