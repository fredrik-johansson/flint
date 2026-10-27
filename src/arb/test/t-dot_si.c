/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"


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

TEST_FUNCTION_START(arb_dot_si, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_ptr x, y;
        slong * w;
        arb_t s1, s2, z;
        slong i, len, prec;
        int initial, subtract, revx, revy;

        len = n_randint(state, 5);
        prec = 2 + n_randint(state, 200);

        initial = n_randint(state, 2);
        subtract = n_randint(state, 2);
        revx = n_randint(state, 2);
        revy = n_randint(state, 2);

        x = _arb_vec_init(len);
        y = _arb_vec_init(len);
        w = flint_malloc(sizeof(ulong) * len);
        arb_init(s1);
        arb_init(s2);
        arb_init(z);

        for (i = 0; i < len; i++)
        {
            arb_randtest(x + i, state, 2 + n_randint(state, 200), 10);
            w[i] = n_randtest(state);
            arb_set_si(y + i, w[i]);
        }

        arb_randtest(s1, state, 200, 10);
        arb_randtest(s2, state, 200, 10);
        arb_randtest(z, state, 200, 10);

        arb_dot(s1, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (y + len - 1) : y, revy ? -1 : 1,
            len, prec);

        arb_dot_si(s2, initial ? z : NULL, subtract,
            revx ? (x + len - 1) : x, revx ? -1 : 1,
            revy ? (w + len - 1) : w, revy ? -1 : 1,
            len, prec);

        if (!_balls_agree(s1, s2, prec))
        {
            flint_printf("FAIL\n\n");
            flint_printf("iter = %wd, len = %wd, prec = %wd\n\n", iter, len, prec);

            if (initial)
            {
                flint_printf("z = ", i); arb_printn(z, 100, ARB_STR_MORE); flint_printf(" (%wd)\n\n", arb_bits(z));
            }

            for (i = 0; i < len; i++)
            {
                flint_printf("x[%wd] = ", i); arb_printn(x + i, 100, ARB_STR_MORE); flint_printf(" (%wd)\n", arb_bits(x + i));
                flint_printf("y[%wd] = ", i); arb_printn(y + i, 100, ARB_STR_MORE); flint_printf(" (%wd)\n", arb_bits(y + i));
            }
            flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 100, ARB_STR_MORE); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 100, ARB_STR_MORE); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(s1);
        arb_clear(s2);
        arb_clear(z);
        _arb_vec_clear(x, len);
        _arb_vec_clear(y, len);
        flint_free(w);
    }

    TEST_FUNCTION_END(state);
}
