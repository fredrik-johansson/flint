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

#include "profiler.h"

#include "arb.h"

int
_gr_quaternion_mul_classical(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
int
_gr_quaternion_mul_fast(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);

TEST_FUNCTION_START(gr_quaternion, state)
{
    gr_ptr a, b;
    gr_ctx_t R, C;
    int reps = 100 * flint_test_multiplier();
    int i, test_flags;

    for (i = 0; i < 3; i++)
    {
        flint_rand_t state;
        flint_rand_init(state);

        switch (i)
        {
            case 0:
                gr_ctx_init_fmpz(R);
                a = gr_heap_init(R);
                b = gr_heap_init(R);
                fmpz_set_si(a, -1);
                fmpz_set_si(b, -1);
                break;

            case 1:
                gr_ctx_init_fmpz(R);
                a = gr_heap_init(R);
                b = gr_heap_init(R);
                fmpz_randprime(a, state, 128, 0);
                fmpz_set_si(b, -1);
                break;

            case 2:
                gr_ctx_init_fmpz(R);
                a = gr_heap_init(R);
                b = gr_heap_init(R);
                fmpz_randprime(a, state, 256, 0);
                fmpz_randprime(b, state, 256, 0);
                break;
        }

        gr_ctx_init_gr_quaternion(C, R, a, b, 0);
        gr_test_ring(C, reps, test_flags);

        /* Quick profiling */
        {
            slong bits;
            slong j;

            gr_ptr x, y, z;
            x = gr_heap_init(C);
            y = gr_heap_init(C);
            z = gr_heap_init(C);

            gr_ctx_println(C);

            for (bits = 32; bits <= 100000; bits *= 2)
            {
                for (j = 0; j < 4; j++)
                {
                    fmpz_randbits(((fmpz *) x) + j, state, bits);
                    fmpz_randbits(((fmpz *) y) + j, state, bits);
                }

                flint_printf("bits = %wd\n", bits);
                flint_printf("classical:   ");
                TIMEIT_START;
                GR_MUST_SUCCEED(_gr_quaternion_mul_classical(z, x, y, C));
                TIMEIT_STOP;
                flint_printf("fast         ");
                TIMEIT_START;
                GR_MUST_SUCCEED(_gr_quaternion_mul_fast(z, x, y, C));
                TIMEIT_STOP;
            }

            gr_heap_clear(x, C);
            gr_heap_clear(y, C);
            gr_heap_clear(z, C);
            flint_rand_clear(state);
        }

        gr_heap_clear(a, R);
        gr_heap_clear(b, R);
        gr_ctx_clear(C);
        gr_ctx_clear(R);

        flint_rand_clear(state);
    }

    {
        flint_rand_t state;
        flint_rand_init(state);

        slong bits;

        for (bits = 32; bits <= 100000; bits *= 2)
        {
            gr_ctx_init_real_arb(R, bits);
            a = gr_heap_init(R);
            b = gr_heap_init(R);
            arb_sqrt_ui(a, 3, bits);
            arb_sqrt_ui(b, 5, bits);

            gr_ctx_init_gr_quaternion(C, R, a, b, 0);

            flint_printf("\n");

            slong j;

            gr_ptr x, y, z;
            x = gr_heap_init(C);
            y = gr_heap_init(C);
            z = gr_heap_init(C);

            if (bits <= 256)
                gr_ctx_println(C);

            for (j = 0; j < 4; j++)
            {
                arb_sqrt_ui(((arb_struct *) x) + j, 11 * (j + 1), bits);
                arb_sqrt_ui(((arb_struct *) y) + j, 13 * (j + 1), bits);
            }

            flint_printf("bits = %wd\n", bits);
            flint_printf("classical:   ");
            TIMEIT_START;
            GR_MUST_SUCCEED(_gr_quaternion_mul_classical(z, x, y, C));
            TIMEIT_STOP;
            if (bits <= 256)
            {
                flint_printf("z = ");
                gr_println(z, C);
            }
            flint_printf("fast         ");
            TIMEIT_START;
            GR_MUST_SUCCEED(_gr_quaternion_mul_fast(z, x, y, C));
            TIMEIT_STOP;
            if (bits <= 256)
            {
                flint_printf("z = ");
                gr_println(z, C);
            }


            gr_heap_clear(x, C);
            gr_heap_clear(y, C);
            gr_heap_clear(z, C);
            flint_rand_clear(state);

            gr_heap_clear(a, R);
            gr_heap_clear(b, R);
            gr_ctx_clear(C);
            gr_ctx_clear(R);
        }

        flint_rand_clear(state);
    }


    for (i = 0; i < 4; i++)
    {
        // test_flags = GR_TEST_VERBOSE;
        test_flags = 0;

        switch (i)
        {
            case 0:  gr_ctx_init_fmpz(R); break;
            case 1:  gr_ctx_init_fmpq(R); break;
            case 2:  gr_ctx_init_fmpz_poly(R); break;
            case 3:  gr_ctx_init_real_arb(R, 53); break;
        }
        
        a = gr_heap_init(R);
        b = gr_heap_init(R);
        gr_set_si(a, -1-i, R);
        gr_set_si(b, -1-2*i, R);
        gr_ctx_init_gr_quaternion(C, R, a, b, 0);
        gr_test_ring(C, reps, test_flags);

        gr_heap_clear(a, R);
        gr_heap_clear(b, R);
        gr_ctx_clear(C);
        gr_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}
