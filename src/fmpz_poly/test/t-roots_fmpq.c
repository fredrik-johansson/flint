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
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

/* Reference: rational roots from the complete factorization, sorted */
static slong
_roots_via_factor(fmpq * res, slong * exp, const fmpz_poly_t f)
{
    fmpz_poly_factor_t fac;
    slong i, j, num = 0;

    fmpz_poly_factor_init(fac);
    fmpz_poly_factor(fac, f);

    for (i = 0; i < fac->num; i++)
    {
        if (fac->p[i].length == 2)
        {
            fmpz_neg(fmpq_numref(res + num), fac->p[i].coeffs);
            fmpz_set(fmpq_denref(res + num), fac->p[i].coeffs + 1);
            fmpq_canonicalise(res + num);
            exp[num] = fac->exp[i];
            num++;
        }
    }

    /* sort */
    for (i = 1; i < num; i++)
    {
        for (j = i; j > 0 && fmpq_cmp(res + j - 1, res + j) > 0; j--)
        {
            fmpq_swap(res + j - 1, res + j);
            FLINT_SWAP(slong, exp[j - 1], exp[j]);
        }
    }

    fmpz_poly_factor_clear(fac);
    return num;
}

static void
_randtest(fmpz_poly_t f, flint_rand_t state)
{
    fmpz_poly_t g;
    fmpq_t r;
    slong i, k, m;

    fmpz_poly_init(g);
    fmpq_init(r);

    switch (n_randint(state, 4))
    {
        case 0:
            /* product of linear factors with multiplicities times random */
            fmpz_poly_randtest_not_zero(f, state, 1 + n_randint(state, 8), 1 + n_randint(state, 50));
            k = n_randint(state, 8);
            for (i = 0; i < k; i++)
            {
                if (n_randint(state, 2))
                    fmpq_randtest(r, state, 1 + n_randint(state, 60));
                else
                    fmpq_set_si(r, (slong) n_randint(state, 21) - 10, 1 + n_randint(state, n_randint(state, 2) ? 1 : 6));
                fmpz_poly_zero(g);
                fmpz_poly_set_coeff_fmpz(g, 1, fmpq_denref(r));
                fmpz_neg(fmpq_numref(r), fmpq_numref(r));
                fmpz_poly_set_coeff_fmpz(g, 0, fmpq_numref(r));
                m = 1 + (n_randint(state, 4) == 0) + (n_randint(state, 8) == 0);
                fmpz_poly_pow(g, g, m);
                fmpz_poly_mul(f, f, g);
            }
            break;

        case 1:
            /* many small integer / rational roots */
            fmpz_poly_one(f);
            k = n_randint(state, 30);
            for (i = 0; i < k; i++)
            {
                fmpz_poly_zero(g);
                fmpz_poly_set_coeff_si(g, 1, 1 + n_randint(state, 3));
                fmpz_poly_set_coeff_si(g, 0, (slong) n_randint(state, 41) - 20);
                fmpz_poly_mul(f, f, g);
            }
            if (n_randint(state, 2))
                fmpz_poly_scalar_mul_si(f, f, (slong) n_randint(state, 100) - 50);
            if (fmpz_poly_is_zero(f))
                fmpz_poly_one(f);
            break;

        case 2:
            /* Swinnerton-Dyer polynomials (no rational roots, many
               roots modulo every prime) times a linear factor */
            fmpz_poly_swinnerton_dyer(f, 1 + n_randint(state, 4));
            fmpz_poly_zero(g);
            fmpz_poly_set_coeff_si(g, 1, 1 + n_randint(state, 5));
            fmpz_poly_set_coeff_si(g, 0, (slong) n_randint(state, 21) - 10);
            fmpz_poly_mul(f, f, g);
            break;

        default:
            /* generic random polynomial */
            fmpz_poly_randtest_not_zero(f, state, 1 + n_randint(state, 30), 1 + n_randint(state, 100));
            break;
    }

    fmpz_poly_clear(g);
    fmpq_clear(r);
}

TEST_FUNCTION_START(fmpz_poly_roots_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t f;
        fmpq * r1, * r2;
        fmpz * z;
        slong * e1, * e2, * e3;
        slong n1, n2, n3, len, i, j;

        fmpz_poly_init(f);
        _randtest(f, state);

        len = FLINT_MAX(f->length, 2);
        r1 = _fmpq_vec_init(len);
        r2 = _fmpq_vec_init(len);
        z = _fmpz_vec_init(len);
        e1 = flint_malloc(sizeof(slong) * len);
        e2 = flint_malloc(sizeof(slong) * len);
        e3 = flint_malloc(sizeof(slong) * len);

        n1 = _roots_via_factor(r1, e1, f);
        n2 = fmpz_poly_roots_fmpq(r2, (n_randint(state, 4) == 0) ? NULL : e2, f);

        if (n2 == n1 && n_randint(state, 4) == 0)
            n2 = fmpz_poly_roots_fmpq(r2, e2, f);

        if (n1 != n2)
        {
            flint_printf("FAIL (fmpq, number of roots)\n");
            flint_printf("f = %{fmpz_poly}\n", f);
            flint_printf("n1 = %wd, n2 = %wd\n", n1, n2);
            flint_abort();
        }

        n2 = fmpz_poly_roots_fmpq(r2, e2, f);
        for (i = 0; i < n1; i++)
        {
            if (!fmpq_equal(r1 + i, r2 + i) || e1[i] != e2[i])
            {
                flint_printf("FAIL (fmpq, roots)\n");
                flint_printf("f = %{fmpz_poly}\n", f);
                flint_printf("i = %wd: %{fmpq} (%wd), %{fmpq} (%wd)\n", i, r1 + i, e1[i], r2 + i, e2[i]);
                flint_abort();
            }
        }

        /* integer roots */
        n3 = fmpz_poly_roots_fmpz(z, e3, f);
        for (i = j = 0; i < n1; i++)
        {
            if (fmpz_is_one(fmpq_denref(r1 + i)))
            {
                if (j >= n3 || !fmpz_equal(z + j, fmpq_numref(r1 + i)) || e3[j] != e1[i])
                {
                    flint_printf("FAIL (fmpz)\n");
                    flint_printf("f = %{fmpz_poly}\n", f);
                    flint_abort();
                }
                j++;
            }
        }

        if (j != n3)
        {
            flint_printf("FAIL (fmpz, number of roots)\n");
            flint_printf("f = %{fmpz_poly}\n", f);
            flint_abort();
        }

        fmpz_poly_clear(f);
        _fmpq_vec_clear(r1, len);
        _fmpq_vec_clear(r2, len);
        _fmpz_vec_clear(z, len);
        flint_free(e1);
        flint_free(e2);
        flint_free(e3);
    }

    TEST_FUNCTION_END(state);
}
