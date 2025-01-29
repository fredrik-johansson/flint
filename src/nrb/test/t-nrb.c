/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arf.h"
#include "acf.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "nrb.h"

TEST_FUNCTION_START(nrb, state)
{
    gr_ctx_t ctx;
    slong prec;

    for (prec = NFLOAT_MIN_LIMBS * FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        nrb_ctx_init(ctx, prec);
        gr_test_ring(ctx, 100 * flint_test_multiplier(), 0 * GR_TEST_VERBOSE);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
