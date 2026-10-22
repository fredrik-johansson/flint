/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-algebraic.c"
#include "t-api.c"
#include "t-map.c"
#include "t-flat.c"
#include "t-transcendental.c"
#include "t-modular.c"
#include "t-no_roots.c"
#include "t-zero_content.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_tower_algebraic),
    TEST_FUNCTION(gr_tower_api),
    TEST_FUNCTION(gr_tower_map),
    TEST_FUNCTION(gr_tower_flat),
    TEST_FUNCTION(gr_tower_transcendental),
    TEST_FUNCTION(gr_tower_modular),
    TEST_FUNCTION(gr_tower_no_roots),
    TEST_FUNCTION(gr_tower_zero_content)
};

/* main function *************************************************************/

TEST_MAIN(tests)
