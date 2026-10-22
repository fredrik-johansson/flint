/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

/* (before gr_special.h, which declares the Dirichlet functions only after it) */
#include "dirichlet.h"

#include "t-arith.c"
#include "t-api.c"
#include "t-richardson.c"
#include "t-elementary.c"
#include "t-roots.c"
#include "t-qqbar.c"
#include "t-numeric.c"
#include "t-gc.c"
#include "t-dense.c"
#include "t-dense_form.c"
#include "t-threads.c"
#include "t-subfields.c"
#include "t-poly.c"
#include "t-factor.c"
#include "t-trig.c"
#include "t-special.c"
#include "t-catalog.c"
#include "t-options.c"
#include "t-hypgeom.c"
#include "t-modforms.c"
#include "t-exp_torsion.c"
#include "t-factored.c"
#include "t-unlimited.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_ca_arith),
    TEST_FUNCTION(gr_ca_api),
    TEST_FUNCTION(gr_ca_richardson),
    TEST_FUNCTION(gr_ca_elementary),
    TEST_FUNCTION(gr_ca_roots),
    TEST_FUNCTION(gr_ca_qqbar),
    TEST_FUNCTION(gr_ca_numeric),
    TEST_FUNCTION(gr_ca_gc),
    TEST_FUNCTION(gr_ca_dense),
    TEST_FUNCTION(gr_ca_dense_form),
    TEST_FUNCTION(gr_ca_threads),
    TEST_FUNCTION(gr_ca_subfields),
    TEST_FUNCTION(gr_ca_poly),
    TEST_FUNCTION(gr_ca_factor),
    TEST_FUNCTION(gr_ca_trig),
    TEST_FUNCTION(gr_ca_special),
    TEST_FUNCTION(gr_ca_catalog),
    TEST_FUNCTION(gr_ca_options),
    TEST_FUNCTION(gr_ca_hypgeom),
    TEST_FUNCTION(gr_ca_modforms),
    TEST_FUNCTION(gr_ca_exp_torsion),
    TEST_FUNCTION(gr_ca_factored),
    TEST_FUNCTION(gr_ca_unlimited)
};

/* main function *************************************************************/

TEST_MAIN(tests)
