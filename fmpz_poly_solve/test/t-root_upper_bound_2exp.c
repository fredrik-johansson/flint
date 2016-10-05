/*
  Copyright (C) 2016 Elias Tsigaridas

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_poly_solve.h"

int
main(void)
{
    int iter;
    FLINT_TEST_INIT(state);

    flint_printf("root_upper_bound_2exp ....");
    fflush(stdout);

    /* Check aliasing */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++) 
    {
        fmpz_poly_t f;
        fmpz_t b;
        slong k;
        fmpz_init(b);
        
        fmpz_poly_init(f);
             
        fmpz_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpz_poly_solve_remove_content_2exp(f);

        k = fmpz_poly_solve_root_upper_bound_2exp(f);
        if (k < 0) continue;
        
        fmpz_init_set_ui(b, 1);
        fmpz_mul_2exp(b, b, k);

        /* printf("\n Bounds: "); fmpz_print( b); printf(" \n"); */
        if ( fmpz_cmp_ui(b, 0) == 0 )
        {
            flint_printf("FAIL:\n");       
            fmpz_poly_print(f); printf("\n\n");
            printf("ERROR \n"); 
            abort();
        }


        fmpz_clear(b);
       
        fmpz_poly_clear(f);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}