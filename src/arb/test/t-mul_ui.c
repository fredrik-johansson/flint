/*
    Copyright (C) 2012-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "ulong_extras.h"


/* arb_mul and its arf/si/ui/fmpz variants may round differently: the
   results must overlap and be of comparable accuracy */
#ifndef _balls_agree
#define _balls_agree _balls_agree
static int
_balls_agree(const arb_t a, const arb_t b, slong prec)
{
    slong acc1, acc2;

    if (!arb_overlaps(a, b))
        return 0;
    /* one result may be much tighter (the rounding error is measured
       rather than bounded by one ulp): fine if it is contained in the other */
    if (arb_contains(a, b) || arb_contains(b, a))
        return 1;
    if (arb_is_exact(a) || arb_is_exact(b) || !arb_is_finite(a) || !arb_is_finite(b))
        return arb_is_exact(a) == arb_is_exact(b);
    acc1 = arb_rel_accuracy_bits(a);
    acc2 = arb_rel_accuracy_bits(b);
    /* either comparable, or both essentially full precision */
    return FLINT_ABS(acc1 - acc2) <= 3 || FLINT_MIN(acc1, acc2) >= prec - 3;
}
#endif

TEST_FUNCTION_START(arb_mul_ui, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        ulong x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        x = n_randtest(state);

        prec = 2 + n_randint(state, 2000);

        arb_set_ui(b, x);
        arb_mul_ui(c, a, x, prec);
        arb_mul(d, a, b, prec);

        if (!_balls_agree(c, d, prec))
        {
            flint_printf("FAIL\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    /* aliasing */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        ulong x;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_special(a, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(b, state, 1 + n_randint(state, 2000), 100);
        arb_randtest_special(c, state, 1 + n_randint(state, 2000), 100);
        x = n_randtest(state);

        prec = 2 + n_randint(state, 2000);

        arb_set_ui(b, x);
        arb_mul_ui(c, a, x, prec);
        arb_mul_ui(a, a, x, prec);

        if (!_balls_agree(a, c, prec))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
