/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "longlong.h"

void _nmod_poly_mulhigh(nn_ptr res, nn_srcptr poly1, slong len1,
                             nn_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    slong bits, bits2;

    if (len1 + len2 <= 6)
    {
        _nmod_poly_mulhigh_classical(res, poly1, len1, poly2, len2, n, mod);
        return;
    }

    bits = FLINT_BITS - (slong) mod.norm;
    bits2 = FLINT_BIT_COUNT(len1);

    if (2 * bits + bits2 <= FLINT_BITS && len1 + len2 < 16)
        _nmod_poly_mulhigh_classical(res, poly1, len1, poly2, len2, n, mod);
    else
        _nmod_poly_mul_KS(res, poly1, len1, poly2, len2, 0, mod);
}

void nmod_poly_mulhigh(nmod_poly_t res,
                   const nmod_poly_t poly1, const nmod_poly_t poly2, slong start)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;
    len_out = len1 + len2 - 1;

    if (len1 == 0 || len2 == 0 || start >= len_out)
    {
        nmod_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;

        nmod_poly_init2(temp, poly1->mod.n, len_out);

        if (len1 >= len2)
            _nmod_poly_mulhigh(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, start, poly1->mod);
        else
            _nmod_poly_mulhigh(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, start, poly1->mod);

        nmod_poly_swap(temp, res);
        nmod_poly_clear(temp);
    } else
    {
        nmod_poly_fit_length(res, len_out);

        if (len1 >= len2)
            _nmod_poly_mulhigh(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, start, poly1->mod);
        else
            _nmod_poly_mulhigh(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, start, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
