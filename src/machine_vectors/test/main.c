/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-ops.c"
#include "t-mod.c"
#include "t-gemm.c"
/* must come last: it clobbers the native VEC4D_TRANSPOSE and
   vecKn_bit_shift_right_32 macros, see the comments there */
#include "t-force_generic.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(machine_vectors_ops),
    TEST_FUNCTION(machine_vectors_mod),
    TEST_FUNCTION(machine_vectors_gemm),
    TEST_FUNCTION(machine_vectors_force_generic)
};

/* main function *************************************************************/

TEST_MAIN(tests)
