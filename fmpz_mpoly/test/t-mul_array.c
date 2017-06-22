/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result, ok1, ok2;
    FLINT_TEST_INIT(state);

    flint_printf("mul_array....");
    fflush(stdout);

    /* Check f*g = g*f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       exp_bits = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k, state, len, exp_bound, coeff_bits, ctx);

          ok1 = fmpz_mpoly_mul_array(h, f, g, ctx);
             
          ok2 = fmpz_mpoly_mul_array(k, g, f, ctx);

          result = (ok1 == 0 && ok2 == 0) || fmpz_mpoly_equal(h, k, ctx);

          if (!result)
          {
             const char * vars[20];
             vars[0] = "x1", vars[1] = "x2", vars[2] = "x3", vars[3] = "x4",
             vars[4] = "x5", vars[5] = "x6", vars[6] = "x7", vars[7] = "x8",
             vars[8] = "x9", vars[9] = "x10", vars[10] = "x11", vars[11] = "x12",
             vars[12] = "x13", vars[13] = "x14", vars[14] = "x15", vars[15] = "x16",
             vars[16] = "x17", vars[17] = "x18", vars[18] = "x19", vars[19] = "x20";

             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k, vars, ctx); printf("\n\n");
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
       fmpz_mpoly_clear(k, ctx);  
    }

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k1, k2, t1, t2;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(t1, ctx);
       fmpz_mpoly_init(t2, ctx);
       fmpz_mpoly_init(k1, ctx);
       fmpz_mpoly_init(k2, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       exp_bits = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(k1, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(k2, state, len, exp_bound, coeff_bits, ctx);
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
          fmpz_mpoly_randtest(h, state, len2, exp_bound2, coeff_bits, ctx);

          fmpz_mpoly_add(t1, g, h, ctx);
          ok1 = fmpz_mpoly_mul_array(k1, f, t1, ctx);

          ok2 = fmpz_mpoly_mul_array(t1, f, g, ctx);
          if (ok2)
             ok2 = fmpz_mpoly_mul_array(t2, f, h, ctx);
          if (ok2)
             fmpz_mpoly_add(k2, t1, t2, ctx);

          result = (ok1 == 0 || ok2 == 0) || fmpz_mpoly_equal(k1, k2, ctx);

          if (!result)
          {
             const char * vars[20];
             vars[0] = "x1", vars[1] = "x2", vars[2] = "x3", vars[3] = "x4",
             vars[4] = "x5", vars[5] = "x6", vars[6] = "x7", vars[7] = "x8",
             vars[8] = "x9", vars[9] = "x10", vars[10] = "x11", vars[11] = "x12",
             vars[12] = "x13", vars[13] = "x14", vars[14] = "x15", vars[15] = "x16",
             vars[16] = "x17", vars[17] = "x18", vars[18] = "x19", vars[19] = "x20";

             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k1, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(k2, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(t1, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(t2, vars, ctx); printf("\n\n");
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
       fmpz_mpoly_clear(k1, ctx);  
       fmpz_mpoly_clear(k2, ctx);  
       fmpz_mpoly_clear(t1, ctx);  
       fmpz_mpoly_clear(t2, ctx);  
    }

    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       exp_bits = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);

          ok1 = fmpz_mpoly_mul_array(h, f, g, ctx);
             
          ok2 = fmpz_mpoly_mul_array(f, f, g, ctx);

          result = (ok1 == 0 && ok2 == 0) || fmpz_mpoly_equal(h, f, ctx);

          if (!result)
          {
             const char * vars[20];
             vars[0] = "x1", vars[1] = "x2", vars[2] = "x3", vars[3] = "x4",
             vars[4] = "x5", vars[5] = "x6", vars[6] = "x7", vars[7] = "x8",
             vars[8] = "x9", vars[9] = "x10", vars[10] = "x11", vars[11] = "x12",
             vars[12] = "x13", vars[13] = "x14", vars[14] = "x15", vars[15] = "x16",
             vars[16] = "x17", vars[17] = "x18", vars[18] = "x19", vars[19] = "x20";

             printf("FAIL\n");
             printf("Aliasing test1\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, vars, ctx); printf("\n\n");
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 10) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100);

       exp_bits = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits1 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;
       exp_bits2 = n_randint(state, 20/(nvars + mpoly_ordering_isdeg(ord)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
          fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
          fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);

          ok1 = fmpz_mpoly_mul_array(h, f, g, ctx);
             
          ok2 = fmpz_mpoly_mul_array(g, f, g, ctx);

          result = (ok1 == 0 && ok2 == 0) || fmpz_mpoly_equal(h, g, ctx);

          if (!result)
          {
             const char * vars[20];
             vars[0] = "x1", vars[1] = "x2", vars[2] = "x3", vars[3] = "x4",
             vars[4] = "x5", vars[5] = "x6", vars[6] = "x7", vars[7] = "x8",
             vars[8] = "x9", vars[9] = "x10", vars[10] = "x11", vars[11] = "x12",
             vars[12] = "x13", vars[13] = "x14", vars[14] = "x15", vars[15] = "x16",
             vars[16] = "x17", vars[17] = "x18", vars[18] = "x19", vars[19] = "x20";

             printf("FAIL\n");
             printf("Aliasing test2\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                    "len1 = %ld, exp_bits1 = %ld, exp_bound1 = %lx, "
                    "len2 = %ld, exp_bits2 = %ld, exp_bound2 = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                       len, exp_bits, exp_bound, len1, exp_bits1, exp_bound1,
                               len2, exp_bits2, exp_bound2, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, vars, ctx); printf("\n\n");
          
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

