/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2015 Anubhav Srivastava

******************************************************************************/


#include "fmpz_mat.h"

#define E fmpz_mat_entry

void
fmpz_mat_sqr(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong dim, n;

    n = A->r;
    
    dim = n;

    if (dim <=12)
    {
        if (dim < 3)
        {   
            fmpz_mat_sqr_bodrato(B, A);
        }
        else
        {
            fmpz_mat_mul(B, A, A);    
        }
    }
    else
    {
        slong ab;

        ab = fmpz_mat_max_bits(A);
        ab = FLINT_ABS(ab);

        if (5*(ab + ab) > dim * dim )
        {
            fmpz_mat_sqr_bodrato(B, A);
        }
        else
        {
            fmpz_mat_mul(B, A, A);
        }

    }
}
