/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("ql_a0....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with ql_a0_naive */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = (g > 1 ? 200 : 2000) + n_randint(state, 1000);
        slong bits = n_randint(state, 5);
        slong hprec = prec + 50;
        int has_t = iter % 2;
        int has_z = (iter % 4) / 2;
        slong nbt = (has_t ? 3 : 1);
        slong guard = ACB_THETA_LOW_PREC;
        slong lp = ACB_THETA_LOW_PREC;
        acb_mat_t tau;
        acb_ptr z, t, r, test;
        arb_ptr dist;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
        r = _acb_vec_init(nbt * n);
        test = _acb_vec_init(nbt * n);
        dist = _arb_vec_init(n);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (k = 0; k < g; k++)
        {
            if (has_z)
            {
                acb_urandom(&z[k], state, hprec);
            }
            if (has_t)
            {
                arb_urandom(acb_realref(&t[k]), state, hprec);
            }
        }
        acb_theta_dist_a0(dist, z, tau, lp);

        acb_theta_ql_a0(r, t, z, dist, tau, guard, prec);
        acb_theta_ql_a0_naive(test, t, z, dist, tau, guard, hprec);

            flint_printf("g = %wd, prec = %wd, has_z = %wd, has_t = %wd, tau:\n",
                g, prec, has_z, has_t);
            acb_mat_printd(tau, 5);
            flint_printf("output:\n");
            _acb_vec_printd(r, nbt * n, 5);
            flint_printf("\n");
            _acb_vec_printd(test, nbt * n, 5);
            flint_printf("\n");

        if (!_acb_vec_overlaps(r, test, nbt * n))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(r, nbt * n);
        _acb_vec_clear(test, nbt * n);
        _arb_vec_clear(dist, n);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
