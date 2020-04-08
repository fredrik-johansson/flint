/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

// Berlekamp - Massey algorithm
static slong find_min_poly(mp_limb_t *s, slong N, nmod_t mod)
{
    slong L = 0, m, n, i;
    mp_limb_t c, d_C, d_B = 1;
	
    slong deg_C = 0, deg_B = 0, deg_T = -1;
    mp_limb_t *B, *C, *T;
    B = flint_calloc(N, sizeof(*B));
    C = flint_calloc(N, sizeof(*C));
    T = flint_calloc(N, sizeof(*T));
    B[0] = C[0] = UWORD(1);

	for (n = 0, m = 1; n < N; n++, m++)
	{
		/* d_C = sum_{i = 0}^L C_i * s_{n-i} */
		d_C = s[n];
		for (i = 1; i <= L; i++)
			d_C = nmod_addmul(d_C, C[i], s[n-i], mod);
        if (d_C == 0) continue; /* C and L currently valid */

        /* C(x) = C(x) - (d_C/d_B) x^m B(x); */
        if (L <= 2*n) deg_T = deg_C, memcpy(T, C, (deg_C+1)*sizeof(*T)); /* T(X) = C(X) */
        c = nmod_neg(nmod_div(d_C, d_B, mod), mod);
        for (i = 0; i <= deg_B; ++i)
            C[m+i] = nmod_addmul(C[m+i], B[i], c, mod);
        deg_C = FLINT_MAX(deg_C, deg_B + m);
        while(C[deg_C] == UWORD(0)) --deg_C;  /* Probably unnecessary */

        if (2*L <= n) /* Increase number of errors */
        {
            L = n + 1 - L, m = 0;
            d_B = d_C, deg_B = deg_T;
            memcpy(B, T, (deg_T+1)*sizeof(*B)); /* B(x) = C(x) */
        }
	}
    /* Reverse C into s */
    for (i = 0; i <= L; ++i) s[i] = C[L-i];
    flint_free(B);
    flint_free(C);
    flint_free(T);
	return L;
}

/* Compute s_ij=(A^j y)_i for i = 0,...,ns-1, j = 0,...,num-1*/
static void make_sequences(mp_limb_t **s, slong ns, slong len, const nmod_sparse_mat_t A, mp_srcptr b) 
{
    slong i, j;
    mp_ptr y, Ay;
    y = _nmod_vec_init(A->r);
    Ay = _nmod_vec_init(A->r);
    memcpy(y, b, A->r*sizeof(*y));
    for (j = 0; j < len; ++j) 
    {
        if(j > 0) nmod_sparse_mat_mul_vec(Ay, A, y), memcpy(y, Ay, A->r*sizeof(*y));
        for (i = 0; i < ns; ++i) s[i][j] = y[i];
    }
    _nmod_vec_clear(y);
    _nmod_vec_clear(Ay);
}

/* Compute x = \Sigma_{i = 0}^{L-1} s_i * A^i * b = 0 */
static void make_sum(mp_ptr x, mp_limb_t *s, slong L, const nmod_sparse_mat_t A, mp_srcptr b)
{
    slong i;
    mp_ptr y, Ay;
    y = _nmod_vec_init(A->r);
    Ay = _nmod_vec_init(A->r);
    memcpy(y, b, A->r*sizeof(*y));
    _nmod_vec_scalar_mul_nmod(x, b, A->r, s[0], A->mod);
    for (i = 1; i < L; ++i) 
    {
        nmod_sparse_mat_mul_vec(Ay, A, y), memcpy(y, Ay, A->r*sizeof(*y));
        _nmod_vec_scalar_addmul_nmod(x, y, A->r, s[i], A->mod);
    }
    _nmod_vec_clear(y);
    _nmod_vec_clear(Ay);
}

int nmod_sparse_mat_solve_wiedemann(mp_ptr x, const nmod_sparse_mat_t A, const mp_ptr b)
{
    slong i, L, ret = 0, ns = FLINT_MIN(A->r, 2), len = 2*A->r + 1;
    mp_limb_t **s; 
    mp_ptr Ax;
    if(A->r != A->c) return 0; /* TBD: reduce to square */

    Ax = _nmod_vec_init(A->r);
    s = flint_malloc(ns * sizeof(*s));
    for (i = 0; i < ns; ++i) s[i] = flint_malloc(len*sizeof(*s[i]));
    
    make_sequences(s, ns, len, A, b);

    /* Don't have block Berlekamp yet, just try each one */
    for (i = 0; i < ns && ret == 0; ++i)
    {
        /* Get minimal polynomial */
        L = find_min_poly(s[i], len, A->mod);
        if(s[i][0]==0) continue;

        /* If \sum_{j = 0}^L s_ijA^jb = 0 => x = -1/s[0]\sum_{j = 0}^{L-1} s_i(j-1) A^jb solves Ax = b */
        make_sum(x, s[i]+1, L, A, b);
        _nmod_vec_scalar_mul_nmod(x, x, A->r, nmod_neg(nmod_inv(s[i][0], A->mod), A->mod), A->mod);

        /* Check if successful */
        nmod_sparse_mat_mul_vec(Ax, A, x);
        ret = _nmod_vec_equal(Ax, b, A->r);
    }

    _nmod_vec_clear(Ax);
    for (i = 0; i < ns; ++i) flint_free(s[i]);
    flint_free(s);
    return ret;
}

int nmod_sparse_mat_nullvector_wiedemann(mp_ptr x, const nmod_sparse_mat_t A, flint_rand_t state) 
{
    slong i, L, ret = 0, ns = FLINT_MIN(A->r, 2), len = 2*A->r + 1;
    mp_limb_t **s; 
    mp_ptr Ax, b;
    Ax = _nmod_vec_init(A->r);
    b = _nmod_vec_init(A->r);

    s = flint_malloc(ns * sizeof(*s));
    for (i = 0; i < ns; ++i) s[i] = flint_malloc(len*sizeof(*s[i]));

    _nmod_vec_randtest(x, state, A->r, A->mod);
    nmod_sparse_mat_mul_vec(b, A, x);

    if(A->r != A->c) return 0; /* TBD: reduce to square */
    make_sequences(s, ns, len, A, b);

    for (i = 0; i < ns && ret == 0; ++i)
    {
        /* Get minimal polynomial */
        L = find_min_poly(s[i], len, A->mod);

        /* \sum_{j = 0}^L s_ijA^jb = 0 => x = \sum_{j = 0}^L s_ijA^jx solves Ax = 0 */
        make_sum(x, s[i], L+1, A, x);
        nmod_sparse_mat_mul_vec(Ax, A, x);
        ret = _nmod_vec_is_zero(Ax, A->r);
    }

    _nmod_vec_clear(Ax);
    _nmod_vec_clear(b);
    for (i = 0; i < ns; ++i) flint_free(s[i]);
    flint_free(s);
    return ret;
}