/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arb.h"

/* arb_fmpz_div_fmpz and arb_set_fmpq against the exact quotient: the
   midpoint is the quotient rounded to nearest (or within the quotient's
   guard limb of it), the radius at most half an ulp plus 2^-28 ulp and
   zero exactly when the quotient is exact at prec bits */
TEST_FUNCTION_START(arb_fmpz_div_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t z;
        arf_t t;
        fmpz_t x, y;
        fmpq_t q;
        mag_t w;
        slong prec;
        int ok;

        arb_init(z);
        arf_init(t);
        fmpz_init(x);
        fmpz_init(y);
        fmpq_init(q);
        mag_init(w);

        prec = 2 + n_randint(state, 400);
        fmpz_randtest(x, state, 1 + n_randint(state, n_randint(state, 4) ? 100 : 1000));
        fmpz_randtest_not_zero(y, state, 1 + n_randint(state, n_randint(state, 4) ? 100 : 1000));
        if (n_randint(state, 4) == 0)
        {
            /* exact quotients */
            fmpz_mul(x, x, y);
            if (n_randint(state, 2))
                fmpz_mul_2exp(y, y, n_randint(state, 100));
        }

        if (n_randint(state, 2))
            arb_fmpz_div_fmpz(z, x, y, prec);
        else
        {
            fmpq_set_fmpz_frac(q, x, y);
            arb_set_fmpq(z, q, prec);
        }

        fmpq_set_fmpz_frac(q, x, y);
        ok = arb_contains_fmpq(z, q);

        if (ok && !fmpz_is_zero(x))
        {
            /* exact at prec bits iff the radius is zero */
            arf_t u;
            int exact;
            arf_init(u);
            exact = (arf_set_fmpq(u, q, prec, ARF_RND_NEAR) == 0);
            ok = (exact == mag_is_zero(arb_radref(z)));
            if (ok && exact)
                ok = arf_equal(u, arb_midref(z));
            arf_clear(u);

            arf_mag_set_ulp(w, arb_midref(z), prec);
            mag_mul_ui(w, w, (UWORD(1) << 27) + 2);
            mag_mul_2exp_si(w, w, -28);
            ok = ok && mag_cmp(arb_radref(z), w) <= 0;
        }

        if (!ok)
        {
            flint_printf("FAIL\n\nprec = %wd\n", prec);
            flint_printf("x = "); fmpz_print(x); flint_printf("\n");
            flint_printf("y = "); fmpz_print(y); flint_printf("\n");
            flint_printf("z = "); arb_print(z); flint_printf("\n");
            flint_abort();
        }

        arb_clear(z);
        arf_clear(t);
        fmpz_clear(x);
        fmpz_clear(y);
        fmpq_clear(q);
        mag_clear(w);
    }

    TEST_FUNCTION_END(state);
}
