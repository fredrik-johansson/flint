/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-bitwise_rs_stress.c"
#include "t-exp_bitwise_rs.c"
#include "t-exp_rs.c"
#include "t-log1p_bitwise_rs.c"
#include "t-trig_rs.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fixed_bitwise_rs_stress),
    TEST_FUNCTION(fixed_exp_bitwise_rs),
    TEST_FUNCTION(fixed_exp_rs),
    TEST_FUNCTION(fixed_log1p_bitwise_rs),
    TEST_FUNCTION(fixed_trig_rs)
};

/* main function *************************************************************/

TEST_MAIN(tests)
