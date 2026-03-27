/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

static void
_nmod_poly_mod_matrix_rows_evaluate_horner(nn_ptr res, const nmod_mat_t A, nn_srcptr h, slong n, nn_srcptr poly3, slong len3,
    nn_srcptr poly3inv, slong len3inv, nmod_t mod)
{
    nn_ptr t;
    slong len = A->r;
    slong i;

    t = _nmod_vec_init(n);

    _nmod_vec_set(res, nmod_mat_entry_ptr(A, len - 1, 0), n);

    nmod_poly_mulmod_precond_t hpre;
    _nmod_poly_mulmod_precond_init_num(hpre, h, n, poly3, len3, poly3inv, len3inv, len - 1, mod);

    for (i = len - 2; i >= 0; i--)
    {
        _nmod_poly_mulmod_precond(t, hpre, res, n, mod);
        _nmod_poly_add(res, t, n, nmod_mat_entry_ptr(A, i, 0), n, mod);
    }

    nmod_poly_mulmod_precond_clear(hpre);

    _nmod_vec_clear(t);
}

static void
_nmod_poly_mod_matrix_rows_evaluate_rectangular_precomp(nn_ptr res, const nmod_mat_t A, nn_srcptr h, slong n, nn_srcptr poly3, slong len3,
    nn_srcptr poly3inv, slong len3inv, nn_srcptr xs, slong m, const nmod_poly_mulmod_precond_t xmpre, nmod_t mod)
{
    nn_ptr s, t, q, u;
    slong len = A->r;
    slong i, j, r;

    r = (len + m - 1) / m;

    s = _nmod_vec_init(2 * n);
    t = _nmod_vec_init(2 * n);
    q = _nmod_vec_init(2 * n);
    u = _nmod_vec_init(2 * n);

#define XP(ii) (xs + (ii) * n)
#define COEFF(ii) (nmod_mat_entry_ptr(A, (ii), 0))

    _nmod_vec_set(s, COEFF((r - 1) * m), n);
    _nmod_vec_zero(s + n, n - 1);
    for (j = 0; j < len - (r - 1) * m - 1; j++)
    {
        _nmod_poly_mul(t, XP(1 + j), n, COEFF((r - 1) * m + 1 + j), n, mod);
        _nmod_vec_add(s, s, t, 2 * n - 1, mod);
    }
    _nmod_poly_divrem_newton_n_preinv(q, res, s, 2 * n - 1, poly3, len3, poly3inv, len3inv, mod);

    for (i = r - 2; i >= 0; i--)
    {
        _nmod_vec_set(s, COEFF(i * m), n);
        _nmod_vec_zero(s + n, n - 1);
        for (j = 1; j < m; j++)
        {
            _nmod_poly_mul(t, XP(j), n, COEFF(i * m + j), n, mod);
            _nmod_vec_add(s, s, t, 2 * n - 1, mod);
        }
        _nmod_poly_divrem_newton_n_preinv(q, u, s, 2 * n - 1, poly3, len3, poly3inv, len3inv, mod);
        _nmod_poly_mulmod_precond(t, xmpre, res, n, mod);
        _nmod_vec_add(res, t, u, n, mod);
    }

#undef XP
#undef COEFF

    _nmod_vec_clear(s);
    _nmod_vec_clear(t);
    _nmod_vec_clear(q);
    _nmod_vec_clear(u);
}

static void
_nmod_poly_mod_matrix_rows_evaluate_rectangular(nn_ptr res, const nmod_mat_t A, nn_srcptr h, slong n, nn_srcptr poly3, slong len3,
    nn_srcptr poly3inv, slong len3inv, nmod_t mod)
{
    nmod_poly_mulmod_precond_t xmpre;
    slong i, m;
    slong len = A->r;
    nn_ptr xs;

    m = n_sqrt(len) + 1;
    m = FLINT_MIN(m, 30);
    xs = _nmod_vec_init((m + 1) * n);
    nn_ptr xsptr[31];
    for (i = 0; i <= m; i++)
        xsptr[i] = xs + i * n;
    _nmod_poly_powers_mod_preinv_naive(xsptr, h, n, m + 1, poly3, len3, poly3inv, len3inv, mod);
    _nmod_poly_mulmod_precond_init_num(xmpre, xsptr[m], n, poly3, len3, poly3inv, len3inv, (len + m - 1) / m - 1, mod);

    _nmod_poly_mod_matrix_rows_evaluate_rectangular_precomp(res, A, h, n, poly3, len3, poly3inv, len3inv, xs, m, xmpre, mod);

    nmod_poly_mulmod_precond_clear(xmpre);
    _nmod_vec_clear(xs);
}

void
_nmod_poly_mod_matrix_rows_evaluate(nn_ptr res, const nmod_mat_t A, nn_srcptr h, slong n, nn_srcptr poly3, slong len3,
    nn_srcptr poly3inv, slong len3inv, nmod_t mod)
{
    FLINT_ASSERT(A->c == n);
    FLINT_ASSERT(n == len3 - 1);

    if (A->r <= 10)
        _nmod_poly_mod_matrix_rows_evaluate_horner(res, A, h, n, poly3, len3, poly3inv, len3inv, mod);
    else
        _nmod_poly_mod_matrix_rows_evaluate_rectangular(res, A, h, n, poly3, len3, poly3inv, len3inv, mod);
}

