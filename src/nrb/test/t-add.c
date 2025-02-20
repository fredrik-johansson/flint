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

TEST_FUNCTION_START(nrb_add, state)
{
    gr_ctx_t ctx;
    slong iter, prec, ebits;
    gr_ptr x, y, z;
    arf_t xa, xb, ya, yb, za, zb, xya, xyb;
    int big;
    int status;
    int alias;

    arf_init(xa); arf_init(xb); arf_init(ya); arf_init(yb);
    arf_init(za); arf_init(zb); arf_init(xya); arf_init(xyb);

    for (prec = FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        GR_MUST_SUCCEED(nrb_ctx_init(ctx, prec, 0));

        for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
        {
            x = gr_heap_init(ctx);
            y = gr_heap_init(ctx);
            z = gr_heap_init(ctx);

            ebits = 1 + n_randint(state, 13);

            alias = n_randlimb(state);

            GR_MUST_SUCCEED(nrb_randtest_ebits(x, state, ebits, ctx));
            GR_MUST_SUCCEED(nrb_randtest_ebits(y, state, ebits, ctx));
            GR_MUST_SUCCEED(nrb_randtest_ebits(z, state, ebits, ctx));

            if (alias & 1)
            {
                GR_IGNORE(gr_set(z, x, ctx));
                GR_MUST_SUCCEED(nrb_add(z, z, y, ctx));
            }
            else if (alias & 2)
            {
                GR_IGNORE(gr_set(z, y, ctx));
                GR_MUST_SUCCEED(nrb_add(z, x, z, ctx));
            }
            else
            {
                GR_MUST_SUCCEED(nrb_add(z, x, y, ctx));
            }

            nrb_get_interval_arf(xa, xb, x, ctx, ARF_PREC_EXACT);
            nrb_get_interval_arf(ya, yb, y, ctx, ARF_PREC_EXACT);
            nrb_get_interval_arf(za, zb, z, ctx, ARF_PREC_EXACT);

            arf_add(xya, xa, ya, ARF_PREC_EXACT, ARF_RND_FLOOR);
            arf_add(xyb, xb, yb, ARF_PREC_EXACT, ARF_RND_CEIL);

            if (!(arf_cmp(za, xya) <= 0 && arf_cmp(zb, xyb) >= 0))
            {
                flint_printf("FAIL: enclosure\n");
                gr_ctx_println(ctx);
                flint_printf("alias: %d %d\n", alias & 1, alias & 2);
                flint_printf("x = "); gr_println(x, ctx); flint_printf("\n");
                flint_printf("xa = "); arf_printd(xa, 100); flint_printf("\n\n");
                flint_printf("xb = "); arf_printd(xb, 100); flint_printf("\n\n");
                flint_printf("y = "); gr_println(y, ctx); flint_printf("\n");
                flint_printf("ya = "); arf_printd(ya, 100); flint_printf("\n\n");
                flint_printf("yb = "); arf_printd(yb, 100); flint_printf("\n\n");
                flint_printf("z = "); gr_println(z, ctx); flint_printf("\n");
                flint_printf("za  = "); arf_printd(za, 100); flint_printf("\n\n");
                flint_printf("xya = "); arf_printd(xya, 100); flint_printf("\n\n");
                flint_printf("zb  = "); arf_printd(zb, 100); flint_printf("\n\n");
                flint_printf("xyb = "); arf_printd(xyb, 100); flint_printf("\n\n");

                flint_abort();
            }

            gr_heap_clear(x, ctx);
            gr_heap_clear(y, ctx);
            gr_heap_clear(z, ctx);
        }
        
        gr_ctx_clear(ctx);
    }

    arf_clear(xa); arf_clear(xb); arf_clear(ya); arf_clear(yb);
    arf_clear(za); arf_clear(zb); arf_clear(xya); arf_clear(xyb);

    TEST_FUNCTION_END(state);
}
