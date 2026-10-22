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
#include "gr_mat.h"
#include "gr_tower.h"
#include "gr_ca.h"

/* Factored denominators (GR_TOWER_OPT_FACTORED_DENOMINATORS): random
   programs over rational functions of pi, e and logarithms, run with
   the option off and on (switched at random in the middle of a run as
   well), agree; elements of the factored representation meet algebraic
   numbers, towers which grow, and the generic operations. */

#define NREG 8

/* a random generator or small rational value */
static int
_fact_atom(gr_ptr res, flint_rand_t state, slong ngen, gr_ctx_t K)
{
    slong k = n_randint(state, ngen + 2);
    int status = GR_SUCCESS;

    if (k == 0)
        status |= gr_set_si(res, (slong) n_randint(state, 7) - 3, K);
    else if (k == 1)
    {
        fmpq_t q;
        fmpq_init(q);
        fmpq_set_si(q, (slong) n_randint(state, 9) - 4, 1 + n_randint(state, 5));
        status |= gr_set_fmpq(res, q, K);
        fmpq_clear(q);
    }
    else if (k == 2)
        status |= gr_pi(res, K);
    else if (k == 3)
    {
        status |= gr_one(res, K);
        status |= gr_exp(res, res, K);
    }
    else
    {
        /* log of a prime: many generators make the towers grow */
        static const slong primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
        status |= gr_set_si(res, primes[(k - 4) % 12], K);
        status |= gr_log(res, res, K);
    }
    return status;
}

/* runs a random program (the same for the same seed), the result in res */
static int
_fact_program(gr_ptr res, ulong seed, slong len, slong ngen, int switching, int algebraic, gr_ctx_t K)
{
    flint_rand_t state;
    gr_ptr r;
    slong i, a, b, c, op;
    int status = GR_SUCCESS;

    flint_rand_init(state);
    flint_rand_set_seed(state, seed, seed + 1);
    GR_TMP_INIT_VEC(r, NREG, K);

    for (i = 0; i < NREG; i++)
        status |= _fact_atom(GR_ENTRY(r, i, K->sizeof_elem), state, ngen, K);

    for (i = 0; i < len && status == GR_SUCCESS; i++)
    {
        gr_ptr x, y, z;
        a = n_randint(state, NREG);
        b = n_randint(state, NREG);
        c = n_randint(state, NREG);
        x = GR_ENTRY(r, a, K->sizeof_elem);
        y = GR_ENTRY(r, b, K->sizeof_elem);
        z = GR_ENTRY(r, c, K->sizeof_elem);
        op = n_randint(state, 14);

        /* (the same random stream in either mode) */
        if (n_randint(state, 6) == 0 && switching)
            status |= gr_ca_ctx_set_option(K, GR_TOWER_OPT_FACTORED_DENOMINATORS,
                1 - gr_ca_ctx_get_option(K, GR_TOWER_OPT_FACTORED_DENOMINATORS));

        switch (op)
        {
            case 0: case 1: status |= gr_add(z, x, y, K); break;
            case 2: status |= gr_sub(z, x, y, K); break;
            case 3: case 4: status |= gr_mul(z, x, y, K); break;
            case 5:
                if (gr_is_zero(y, K) == T_FALSE)
                    status |= gr_div(z, x, y, K);
                break;
            case 6:
                if (gr_is_zero(x, K) == T_FALSE)
                    status |= gr_inv(z, x, K);
                break;
            case 7: status |= gr_neg(z, x, K); break;
            case 8: status |= gr_set(z, x, K); break;
            case 9: status |= _fact_atom(z, state, ngen, K); break;
            case 10:
                /* (in place) */
                status |= gr_add(x, x, y, K);
                if (gr_is_zero(y, K) == T_FALSE)
                    status |= gr_div(x, x, y, K);
                break;
            case 11:
                gr_swap(x, y, K);
                break;
            case 12:
                status |= gr_mul_si(z, x, (slong) n_randint(state, 5) - 2, K);
                status |= gr_add_si(z, z, (slong) n_randint(state, 5) - 2, K);
                break;
            default:
                /* an algebraic number joins: the towers merge */
                if (algebraic)
                {
                    status |= gr_set_si(z, 2 + n_randint(state, 2), K);
                    status |= gr_sqrt(z, z, K);
                    status |= gr_add(z, z, x, K);
                }
                else
                    status |= gr_sqr(z, x, K);
                break;
        }

        /* (keep the sizes bounded) */
        if (gr_ca_repr(z, K) != GR_CA_REPR_RATIONAL && n_randint(state, 4) == 0)
            status |= _fact_atom(z, state, ngen, K);
    }

    status |= gr_zero(res, K);
    for (i = 0; i < NREG; i++)
    {
        gr_ptr t = GR_ENTRY(r, i, K->sizeof_elem);
        status |= gr_mul_si(res, res, 3, K);
        status |= gr_add(res, res, t, K);
    }

    GR_TMP_CLEAR_VEC(r, NREG, K);
    flint_rand_clear(state);
    return status;
}

TEST_FUNCTION_START(gr_ca_factored, state)
{
    gr_ctx_t QQ, K;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    /* the representation, and the normalization of constants */
    {
        gr_ptr x, y, z;
        gr_ctx_init_ca(K, QQ, 0);
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(K, GR_TOWER_OPT_FACTORED_DENOMINATORS, 1));
        x = gr_heap_init(K);
        y = gr_heap_init(K);
        z = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_pi(x, K));
        GR_MUST_SUCCEED(gr_add_si(x, x, 1, K));
        GR_MUST_SUCCEED(gr_inv(y, x, K));
        if (gr_ca_repr(y, K) != GR_CA_REPR_FACTORED)
        {
            flint_printf("FAIL: 1/(pi+1) not factored\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_mul(z, x, y, K));
        if (gr_ca_repr(z, K) != GR_CA_REPR_RATIONAL || gr_is_one(z, K) != T_TRUE)
        {
            flint_printf("FAIL: (pi+1) / (pi+1)\n");
            flint_abort();
        }
        GR_MUST_SUCCEED(gr_mul_si(z, y, 2, K));
        GR_MUST_SUCCEED(gr_sub(z, z, y, K));
        GR_MUST_SUCCEED(gr_sub(z, z, y, K));
        if (gr_ca_repr(z, K) != GR_CA_REPR_RATIONAL || gr_is_zero(z, K) != T_TRUE)
        {
            flint_printf("FAIL: 2/(pi+1) - 2/(pi+1)\n");
            flint_abort();
        }
        /* the flat view: equality, printing, numerical value */
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(K, GR_TOWER_OPT_FACTORED_DENOMINATORS, 0));
        GR_MUST_SUCCEED(gr_pi(z, K));
        GR_MUST_SUCCEED(gr_add_si(z, z, 1, K));
        GR_MUST_SUCCEED(gr_inv(z, z, K));
        if (gr_ca_repr(z, K) != GR_CA_REPR_FLAT || gr_equal(y, z, K) != T_TRUE)
        {
            flint_printf("FAIL: flat and factored 1/(pi+1)\n");
            flint_abort();
        }
        gr_heap_clear(x, K);
        gr_heap_clear(y, K);
        gr_heap_clear(z, K);
        gr_ctx_clear(K);
    }

    /* random programs */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ptr a, b;
        ulong seed = n_randlimb(state);
        slong len = 5 + n_randint(state, 30);
        slong ngen = 1 + n_randint(state, (n_randint(state, 4) == 0) ? 14 : 4);
        int switching = n_randint(state, 2);
        int algebraic = n_randint(state, 3) == 0;
        int s1, s2;

        gr_ctx_init_ca(K, QQ, 0);
        a = gr_heap_init(K);
        b = gr_heap_init(K);

        GR_MUST_SUCCEED(gr_ca_ctx_set_option(K, GR_TOWER_OPT_FACTORED_DENOMINATORS, 0));
        s1 = _fact_program(a, seed, len, ngen, 0, algebraic, K);
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(K, GR_TOWER_OPT_FACTORED_DENOMINATORS, 1));
        s2 = _fact_program(b, seed, len, ngen, switching, algebraic, K);

        if (s1 != s2 || (s1 == GR_SUCCESS && gr_equal(a, b, K) != T_TRUE))
        {
            flint_printf("FAIL: random program (iter %wd, seed %wu, len %wd, ngen %wd, switching %d, algebraic %d)\n",
                iter, seed, len, ngen, switching, algebraic);
            flint_printf("status %d %d\n", s1, s2);
            gr_println(a, K);
            gr_println(b, K);
            flint_abort();
        }

        /* the results meet an algebraic number, and the generic matrix
           code */
        if (s1 == GR_SUCCESS && n_randint(state, 4) == 0)
        {
            gr_mat_t M;
            gr_ptr d;
            gr_mat_init(M, 2, 2, K);
            d = gr_heap_init(K);
            GR_MUST_SUCCEED(gr_set(gr_mat_entry_ptr(M, 0, 0, K), a, K));
            GR_MUST_SUCCEED(gr_set_si(gr_mat_entry_ptr(M, 0, 1, K), 2, K));
            GR_MUST_SUCCEED(gr_sqrt(gr_mat_entry_ptr(M, 0, 1, K), gr_mat_entry_ptr(M, 0, 1, K), K));
            GR_MUST_SUCCEED(gr_set(gr_mat_entry_ptr(M, 1, 0, K), b, K));
            GR_MUST_SUCCEED(gr_set(gr_mat_entry_ptr(M, 1, 1, K), gr_mat_entry_ptr(M, 0, 1, K), K));
            /* det = sqrt(2) (a - b) = 0 */
            GR_MUST_SUCCEED(gr_mat_det(d, M, K));
            if (gr_is_zero(d, K) != T_TRUE)
            {
                flint_printf("FAIL: determinant (iter %wd)\n", iter);
                gr_println(d, K);
                flint_abort();
            }
            gr_heap_clear(d, K);
            gr_mat_clear(M, K);
        }

        gr_heap_clear(a, K);
        gr_heap_clear(b, K);
        gr_ctx_clear(K);
    }

    gr_ctx_clear(QQ);
    TEST_FUNCTION_END(state);
}
