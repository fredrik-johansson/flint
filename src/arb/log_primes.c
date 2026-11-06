/*
    Copyright (C) 2013-2014, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mp_real.h"

/* The logarithms of primes and the arguments of Gaussian primes are
   computed by mp_real (machin_bsplit.c: the Machin-type sets by binary
   splitting, Zuniga's series), which also keeps the shared caches of
   their floors for its argument reductions. */

void
arb_log_primes_vec_bsplit(arb_ptr res, slong n, slong prec)
{
    if (n > MP_REAL_LOG_PRIMES_MAX)
        flint_throw(FLINT_ERROR, "(%s): at most %d primes\n", __func__, MP_REAL_LOG_PRIMES_MAX);
    _mp_real_log_primes_vec_arb(res, n, prec);
}

void
arb_atan_gauss_primes_vec_bsplit(arb_ptr res, slong n, slong prec)
{
    if (n > MP_REAL_ATAN_GAUSS_MAX)
        flint_throw(FLINT_ERROR, "(%s): at most %d Gaussian primes\n", __func__, MP_REAL_ATAN_GAUSS_MAX);
    /* (mp_real gives 2 arg(pi_j)) */
    _mp_real_atan_gauss_vec_arb(res, n, prec + 1);
    _arb_vec_scalar_mul_2exp_si(res, res, n, -1);
}
