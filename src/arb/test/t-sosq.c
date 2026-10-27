/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sosq, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, z, v;
        slong prec, bits;

        arb_init(a);
        arb_init(b);
        arb_init(z);
        arb_init(v);

        bits = (n_randint(state, 20) == 0) ? 20000 : 1000;
        prec = 2 + n_randint(state, bits);
        arb_randtest_special(a, state, bits, 100);
        arb_randtest_special(b, state, bits, 100);

        arb_mul(v, a, a, prec + 100);
        arb_addmul(v, b, b, prec + 100);

        if (n_randint(state, 2))
        {
            arb_sosq(z, a, b, prec);
        }
        else
        {
            arb_set(z, a);
            arb_sosq(z, z, b, prec);
        }

        if (!arb_overlaps(z, v))
        {
            flint_printf("FAIL!\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n\n");
            flint_printf("v = "); arb_print(v); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(z);
        arb_clear(v);
    }

    TEST_FUNCTION_END(state);
}
