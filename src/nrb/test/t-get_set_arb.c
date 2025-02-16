/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "gr.h"
#include "nrb.h"

TEST_FUNCTION_START(nrb_get_set_arb, state)
{
    gr_ctx_t ctx;
    slong iter, prec, prec2;
    arb_t x, z;
    gr_ptr y;
    int big;
    int status;

    for (prec = FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        GR_MUST_SUCCEED(nrb_ctx_init(ctx, prec, 0));

        for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
        {
            arb_init(x);
            arb_init(z);
            y = gr_heap_init(ctx);

            big = n_randint(state, 2);

            prec2 = 2 + n_randint(state, 2 * prec);

            if (big)
                arb_randtest_special(x, state, prec2, 100);
            else
                arb_randtest(x, state, prec2, 20);

            status = nrb_set_arb(y, x, ctx);
            status |= nrb_get_arb(z, y, ctx);

            if (!_nrb_is_valid(y, ctx))
            {
                flint_printf("FAIL: invalid\n");
                gr_ctx_println(ctx);
                flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
                flint_printf("y = "); gr_println(y, ctx);
                flint_abort();
            }

            if (!big && status != GR_SUCCESS)
            {
                flint_printf("FAIL: unsuccessful\n");
                gr_ctx_println(ctx);
                flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
                flint_printf("y = "); gr_println(y, ctx);
                flint_abort();
            }

            if (status == GR_SUCCESS && !arb_contains(z, x))
            {
                flint_printf("FAIL: containment\n");
                gr_ctx_println(ctx);
                flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("z = "); arb_printd(z, 50); flint_printf("\n\n");
                flint_abort();
            }

            if (status == GR_SUCCESS && prec2 <= prec && arb_rel_accuracy_bits(x) > 0 &&
                    arb_rel_accuracy_bits(z) < FLINT_MIN(prec + 20, arb_rel_accuracy_bits(x) - 1))
            {
                flint_printf("FAIL: accuracy\n");
                gr_ctx_println(ctx);
                flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("z = "); arb_printd(z, 50); flint_printf("\n\n");
                flint_abort();
            }

            arb_clear(x);
            arb_clear(z);
            gr_heap_clear(y, ctx);
        }
        
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
