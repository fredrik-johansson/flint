/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "mpoly.h"
#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_fit_length_fit_bits(
    fq_nmod_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (d*len > A->coeffs_alloc)
    {
        A->coeffs_alloc = FLINT_MAX(d*len, 2*A->coeffs_alloc);
        A->coeffs = flint_realloc(A->coeffs, A->coeffs_alloc*sizeof(ulong));
    }

    if (bits > A->bits)
    {
        slong newN = mpoly_words_per_exp(bits, ctx->minfo);
        slong new_exps_alloc = newN*len;
        ulong * t;

        if (len < 1)
        {
            A->bits = bits;
            return;
        }

        t = (ulong *) flint_malloc(new_exps_alloc*sizeof(ulong));

        if (A->length > 0)
            mpoly_repack_monomials(t, bits, A->exps, A->bits, A->length, ctx->minfo);

        if (A->exps_alloc > 0)
            flint_free(A->exps);

        A->exps = t;
        A->exps_alloc = new_exps_alloc;
        A->bits = bits;
    }
    else
    {
        if (N*len > A->exps_alloc)
        {
            A->exps_alloc = FLINT_MAX(N*len, 2*A->exps_alloc);
            A->exps = flint_realloc(A->exps, A->exps_alloc*sizeof(ulong));
        }
    }
}
