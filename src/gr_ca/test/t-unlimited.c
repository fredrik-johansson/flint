/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "gr_special.h"
#include "gr_tower.h"
#include "gr_ca.h"

/* The unlimited mode (GR_CA_UNLIMITED): with the precision
   limits set far too low, operations which fail in the limited mode
   succeed in the unlimited mode (the effort limits grow), the options
   stay in force (they are the starting values), and results agree. */

static void
_unl_init(gr_ctx_t K, gr_ctx_t QQ, int unlimited, slong prec)
{
    gr_ctx_init_ca(K, QQ, unlimited ? GR_CA_UNLIMITED : 0);
    if (prec != 0)
    {
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(K, GR_TOWER_OPT_PREC_LIMIT, prec));
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(K, GR_TOWER_OPT_CERTIFY_PREC_LIMIT, prec));
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(K, GR_TOWER_OPT_NUMERIC_PREC_LIMIT, prec));
    }
}

TEST_FUNCTION_START(gr_ca_unlimited, state)
{
    gr_ctx_t QQ, K;
    gr_ptr x, y, z;
    int unl;

    gr_ctx_init_fmpq(QQ);

    /* the flag sets the option; the option can be changed */
    _unl_init(K, QQ, 1, 0);
    if (gr_ca_ctx_get_option(K, GR_TOWER_OPT_UNLIMITED) != 1 ||
        gr_ca_ctx_get_option(K, GR_TOWER_OPT_PREC_LIMIT) != gr_tower_option_default(GR_TOWER_OPT_PREC_LIMIT))
    {
        flint_printf("FAIL: options of an unlimited context\n");
        flint_abort();
    }
    gr_ctx_clear(K);

    for (unl = 0; unl < 2; unl++)
    {
        truth_t t;
        int status;

        /* 1 - erf(30) = 2.6e-393 is not zero: beyond the limits of 64 bits */
        _unl_init(K, QQ, unl, 64);
        x = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_set_si(x, 30, K));
        GR_MUST_SUCCEED(gr_erf(x, x, K));
        GR_MUST_SUCCEED(gr_sub_si(x, x, 1, K));
        t = gr_is_zero(x, K);
        if ((unl && t != T_FALSE) || (!unl && t == T_TRUE))
        {
            flint_printf("FAIL: erf(30) - 1 (unlimited %d): %d\n", unl, t);
            flint_abort();
        }

        /* its sign */
        {
            int c = 0;
            GR_MUST_SUCCEED(gr_set_si(x, 30, K));
            GR_MUST_SUCCEED(gr_erf(x, x, K));
            {
                gr_ptr one = gr_heap_init(K);
                GR_MUST_SUCCEED(gr_one(one, K));
                status = gr_cmp(&c, x, one, K);
                gr_heap_clear(one, K);
            }
            if ((unl && (status != GR_SUCCESS || c != -1)) || (status == GR_SUCCESS && c != -1))
            {
                flint_printf("FAIL: sign of erf(30) - 1 (unlimited %d): %d %d\n", unl, status, c);
                flint_abort();
            }
        }
        gr_heap_clear(x, K);
        gr_ctx_clear(K);

        /* exp(-300) + exp(-301) - exp(-300) is not zero, and positive */
        _unl_init(K, QQ, unl, 64);
        x = gr_heap_init(K);
        y = gr_heap_init(K);
        z = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_set_si(x, -300, K));
        GR_MUST_SUCCEED(gr_exp(x, x, K));
        GR_MUST_SUCCEED(gr_set_si(y, -301, K));
        GR_MUST_SUCCEED(gr_exp(y, y, K));
        GR_MUST_SUCCEED(gr_add(z, x, y, K));
        GR_MUST_SUCCEED(gr_sub(z, z, x, K));
        t = gr_is_zero(z, K);
        if (t != T_FALSE)
        {
            flint_printf("FAIL: exp(-301) (unlimited %d): %d\n", unl, t);
            flint_abort();
        }
        gr_heap_clear(x, K);
        gr_heap_clear(y, K);
        gr_heap_clear(z, K);
        gr_ctx_clear(K);
    }

    /* the same results with default options, limited and unlimited, on
       random small expressions */
    {
        slong iter;
        gr_ctx_t K0, K1;
        _unl_init(K0, QQ, 0, 0);
        _unl_init(K1, QQ, 1, 0);

        for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
        {
            gr_ptr a0, a1;
            char * s0, * s1;
            int st0, st1;
            truth_t z0, z1;

            a0 = gr_heap_init(K0);
            a1 = gr_heap_init(K1);
            st0 = gr_randtest(a0, state, K0);
            /* (the same element: through its string) */
            GR_MUST_SUCCEED(gr_get_str(&s0, a0, K0));
            st1 = gr_set_str(a1, s0, K1);
            if (st0 == GR_SUCCESS && st1 == GR_SUCCESS)
            {
                st0 = gr_sqr(a0, a0, K0);
                st1 = gr_sqr(a1, a1, K1);
                st0 |= gr_add_si(a0, a0, 1, K0);
                st1 |= gr_add_si(a1, a1, 1, K1);
                z0 = gr_is_zero(a0, K0);
                z1 = gr_is_zero(a1, K1);
                if (st0 != st1 || (z0 != T_UNKNOWN && z0 != z1) || z1 == T_UNKNOWN)
                {
                    flint_printf("FAIL: limited vs unlimited (iter %wd): %s\n", iter, s0);
                    flint_printf("status %d %d, zero %d %d\n", st0, st1, z0, z1);
                    flint_abort();
                }
                GR_MUST_SUCCEED(gr_get_str(&s1, a1, K1));
                flint_free(s1);
            }
            flint_free(s0);
            gr_heap_clear(a0, K0);
            gr_heap_clear(a1, K1);
        }

        gr_ctx_clear(K0);
        gr_ctx_clear(K1);
    }

    gr_ctx_clear(QQ);
    TEST_FUNCTION_END(state);
}
