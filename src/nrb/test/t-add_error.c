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

TEST_FUNCTION_START(nrb_add_error, state)
{
    gr_ctx_t ctx;
    slong iter, prec, prec2;
    arf_t xa, xb, ya, yb, za, zb, t;
    gr_ptr x, y;
    int status;
    double err;
    slong err_exp;

    arf_init(xa);
    arf_init(xb);
    arf_init(ya);
    arf_init(yb);
    arf_init(za);
    arf_init(zb);
    arf_init(t);

    for (prec = FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        GR_MUST_SUCCEED(nrb_ctx_init(ctx, prec, 0));

        for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
        {
            x = gr_heap_init(ctx);
            y = gr_heap_init(ctx);

            nrb_randtest_ebits(x, state, 15, ctx);

            err = d_randtest_special(state, D_MIN_NORMAL_EXPONENT, D_MAX_NORMAL_EXPONENT);
            if (!(err >= 0.0))
                err = 0.0;

//            err_exp = n_randint(state, 4000) - 2000;
            err_exp = n_randint(state, 40) - 20;

            GR_IGNORE(gr_set(y, x, ctx));
            nrb_inplace_add_error_d_2exp_si(y, err, err_exp, ctx);

            /* y = x +/- err */
            nrb_get_interval_arf(xa, xb, x, ctx, ARF_PREC_EXACT);
            nrb_get_interval_arf(ya, yb, y, ctx, ARF_PREC_EXACT);

            arf_set_d(t, err);
            arf_mul_2exp_si(t, t, err_exp);

            arf_sub(za, xa, t, ARF_PREC_EXACT, ARF_RND_FLOOR);
            arf_add(zb, xb, t, ARF_PREC_EXACT, ARF_RND_CEIL);

            if (!_nrb_is_valid(y, ctx))
            {
                flint_printf("FAIL: invalid\n");
                gr_ctx_println(ctx);
                flint_printf("x = %{gr}\n\n", x, ctx);
                flint_printf("y = %{gr}\n\n", y, ctx);
                flint_abort();
            }

            if (!(arf_cmp(ya, za) <= 0 && arf_cmp(yb, zb) >= 0))
            {
                flint_printf("FAIL: containment\n");
                gr_ctx_println(ctx);
                flint_printf("x = %{gr}\n\n", x, ctx);
                flint_printf("err = %g * 2^%wd\n\n", err, err_exp);
                flint_printf("t = %{arf}\n\n", t);
                flint_printf("y = %{gr}\n\n", y, ctx);
                flint_printf("xa = %{arf}\n", xa);
                flint_printf("xb = %{arf}\n", xb);
                flint_printf("ya = %{arf}\n", ya);
                flint_printf("yb = %{arf}\n", yb);
                flint_printf("za = %{arf}\n", za);
                flint_printf("zb = %{arf}\n", zb);
                flint_abort();
            }

            nrb_get_rad_arf(xa, x, ctx);
            arf_add(xa, xa, t, 64, ARF_RND_UP);
            arf_set_d(xb, 1.0001);
            arf_mul(xb, xa, xb, 64, ARF_RND_UP);
            nrb_get_mid_arf(yb, x, ctx);
            arf_mul_2exp_si(yb, yb, -FLINT_BITS * NRB_N(x) - 10);
            arf_abs(yb, yb);
            arf_max(xb, xb, yb);

            nrb_get_rad_arf(ya, y, ctx);

            if (arf_cmp(ya, xb) > 0)
            {
                flint_printf("FAIL: accuracy\n");
                gr_ctx_println(ctx);
                flint_printf("x = %{gr}\n\n", x, ctx);
                flint_printf("err = %g * 2^%wd\n\n", err, err_exp);
                flint_printf("t = %{arf}\n\n", t);
                flint_printf("y = %{gr}\n\n", y, ctx);
                flint_printf("xa = %{arf}\n", xa);
                flint_printf("xb = %{arf}\n", xb);
                flint_printf("ya = %{arf}\n", ya);
                flint_abort();
            }

            gr_heap_clear(x, ctx);
            gr_heap_clear(y, ctx);
        }

        gr_ctx_clear(ctx);
    }

    arf_clear(xa);
    arf_clear(xb);
    arf_clear(ya);
    arf_clear(yb);
    arf_clear(za);
    arf_clear(zb);
    arf_clear(t);

    TEST_FUNCTION_END(state);
}
