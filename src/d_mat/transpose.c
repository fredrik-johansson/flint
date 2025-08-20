/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_mat.h"

void
d_mat_transpose(d_mat_t B, const d_mat_t A)
{
    slong ii, jj, i, j, blocksize;
    blocksize = 64 / sizeof(double);

    if (B->r != A->c || B->c != A->r)
    {
        flint_throw(FLINT_ERROR, "Exception (d_mat_transpose). Incompatible dimensions.\n");
    }

    if (B == A)
    {
        /* todo: cache efficient version */
        for (i = 0; i < A->r - 1; i++)
            for (j = i + 1; j < A->c; j++)
                FLINT_SWAP(double, d_mat_entry(A, j, i), d_mat_entry(A, i, j));
    }
    else
    {
        for (ii = 0; ii < B->r; ii += blocksize)
            for (jj = 0; jj < B->c; jj += blocksize)
                for (i = ii; i < FLINT_MIN(ii + blocksize, B->r); i++)
                    for (j = jj; j < FLINT_MIN(jj + blocksize, B->c); j++)
                        d_mat_entry(B, i, j) = d_mat_entry(A, j, i);
    }
}
