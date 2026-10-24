/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "padic_radix.h"

TEST_FUNCTION_START(padic_radix_charpoly, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        slong n, i, j, e;
        int s1, s2;
        padic_radix_struct * c1, * c2;
        padic_radix_t d;

        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, 100));
        e = radix_limb_exponent(PADIC_RADIX_CTX_RADIX(ctx));
        n = n_randint(state, 12);

        gr_mat_init(A, n, n, ctx);
        c1 = gr_heap_init_vec(n + 1, ctx);
        c2 = gr_heap_init_vec(n + 1, ctx);
        padic_radix_init(d, ctx);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * x = gr_mat_entry_ptr(A, i, j, ctx);
                padic_radix_randtest(x, state, ctx);
                if (n_randint(state, 3) == 0)
                    GR_IGNORE(padic_radix_set_si(x, (slong) n_randint(state, 21) - 10, ctx));
                if (x->u.size != 0 && n_randint(state, 3) == 0)
                {
                    slong rel = (x->N == PADIC_RADIX_EXACT) ? 0 : x->N - x->v;
                    x->v = (slong) n_randint(state, 4 * e + 1) - 2 * e;
                    if (x->N != PADIC_RADIX_EXACT)
                        x->N = x->v + rel;
                    if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
                        GR_IGNORE(padic_radix_zero(x, ctx));
                }
            }

        s1 = _padic_radix_mat_charpoly_berkowitz(c1, A, ctx);
        s2 = _gr_mat_charpoly_berkowitz(c2, A, ctx);

        if (s1 == GR_SUCCESS && s2 == GR_SUCCESS)
            for (i = 0; i <= n; i++)
            {
                slong N = FLINT_MIN(c1[i].N, c2[i].N);
                GR_IGNORE(padic_radix_sub(d, c1 + i, c2 + i, ctx));
                if (d->u.size != 0 && d->v < N)
                {
                    flint_printf("FAIL: charpoly coefficient %wd\n", i);
                    gr_ctx_println(ctx);
                    flint_printf("c1 = "); gr_println(c1 + i, ctx);
                    flint_printf("c2 = "); gr_println(c2 + i, ctx);
                    flint_abort();
                }
            }

        gr_heap_clear_vec(c1, n + 1, ctx);
        gr_heap_clear_vec(c2, n + 1, ctx);
        padic_radix_clear(d, ctx);
        gr_mat_clear(A, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
