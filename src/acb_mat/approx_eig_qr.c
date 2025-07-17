/*
    Copyright 2013 Timo Hartmann
    Copyright 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

/*

N=100; A = Mat(ComplexFloat_acf(384),N,N)([[randint(-10,10) + randint(-10,10)*1j for i in range(N)] for j in range(N)]) / 17; d=A.diagonalization()

*/

/* todo: gr_sosq, gr_fmma, gr_fmms */
/* todo: rearrange temp variables to avoid shallow sets */

// TODO: IMPORTANT: make sure we use nfloat_sqr when squaring, e.g. avoid gr_addmul
// TODO: squaring in dot products also

/* aliasing only allowed for a */
static int
gr_sosq(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_ctx_t ctx)
{
    int status;
    status = gr_sqr(res, a, ctx);
    status |= gr_addmul(res, b, b, ctx);
    return status;
}

int
gr_mat_approx_qr_step(gr_mat_t A, gr_mat_t Q, slong n0, slong n1, gr_srcptr shift, gr_ctx_t ctx, gr_ctx_t real_ctx)
{
    slong j, k, n;
    /* complex variables */
    gr_ptr c, s, negs, cc, cs, negcs, t;
    /* length-2 shallow vectors */
    gr_ptr v1, v1neg, v2, v2neg, v3;
    /* real variables */
    gr_ptr v, u;

    slong sz = ctx->sizeof_elem;
    slong rsz = real_ctx->sizeof_elem;

#define RE(xx) (xx)
#define IM(xx) GR_ENTRY(xx, 1, rsz)

    int status = GR_SUCCESS;

    FLINT_ASSERT(sz == 2 * rsz);

//    flint_printf("sizes %wd, %wd\n", sz, rsz);

    n = A->r;

    GR_TMP_INIT5(c, s, negs, cc, cs, ctx);
    GR_TMP_INIT2(negcs, t, ctx);
    GR_TMP_INIT2(v, u, real_ctx);

    v1 = GR_TMP_ALLOC(10 * sz);
    v1neg = GR_ENTRY(v1, 2, sz);
    v2 = GR_ENTRY(v1, 4, sz);
    v2neg = GR_ENTRY(v1, 6, sz);
    v3 = GR_ENTRY(v1, 8, sz);

    /* TODO: AVOID SOME SETS */
    /* TODO: ALLOW DOING SOME REPEATED OPS AS VECTOR OPS */
    /* TODO: Givens using length-4 dot? (c, s consecutively in vector) */

    /* Calculate Givens rotation */
    status |= gr_sub(c, GR_MAT_ENTRY(A, n0, n0, sz), shift, ctx);

//    flint_printf("gr   c = %{gr}\n\n", c, ctx);
//    flint_printf("gr   RE(c) = %{gr}\n\n", RE(c), real_ctx);
//    flint_printf("gr   IM(c) = %{gr}\n\n", IM(c), real_ctx);

    status |= gr_set(s, GR_MAT_ENTRY(A, n0 + 1, n0, sz), ctx);

//    flint_printf("gr   RE(s) = %{gr}\n\n", RE(s), real_ctx);
//    flint_printf("gr   IM(s) = %{gr}\n\n", IM(s), real_ctx);

    status |= gr_sosq(v, RE(c), IM(c), real_ctx);
    status |= gr_sosq(u, RE(s), IM(s), real_ctx);
    status |= gr_add(v, v, u, real_ctx);
    status |= gr_sqrt(v, v, real_ctx);

//    flint_printf("gr   v = %{gr}\n\n", v, real_ctx);

    /* todo: need an unknown case? */
    if (gr_is_zero(v, real_ctx) == T_TRUE)
    {
        status |= gr_one(v, real_ctx);
        status |= gr_one(c, ctx);
        status |= gr_zero(s, ctx);
    }
    else
    {
        status |= gr_inv(v, v, real_ctx);
        status |= gr_mul(RE(c), RE(c), v, real_ctx);
        status |= gr_mul(IM(c), IM(c), v, real_ctx);
        status |= gr_mul(RE(s), RE(s), v, real_ctx);
        status |= gr_mul(IM(s), IM(s), v, real_ctx);
    }

//    flint_printf("gr   c = %{gr}  s = %{gr}\n\n", c, ctx, s, ctx);

    status |= gr_conj(cc, c, ctx);
    status |= gr_conj(cs, s, ctx);

//    flint_printf("gr   cc = %{gr}  cs = %{gr}\n\n", cc, ctx, cs, ctx);

    status |= gr_neg(negs, s, ctx);
    status |= gr_neg(negcs, cs, ctx);

    gr_set_shallow(GR_ENTRY(v1, 0, sz), c, ctx);
    gr_set_shallow(GR_ENTRY(v1, 1, sz), s, ctx);
    gr_set_shallow(GR_ENTRY(v1neg, 0, sz), c, ctx);
    gr_set_shallow(GR_ENTRY(v1neg, 1, sz), negs, ctx);
    gr_set_shallow(GR_ENTRY(v2, 0, sz), cc, ctx);
    gr_set_shallow(GR_ENTRY(v2, 1, sz), cs, ctx);
    gr_set_shallow(GR_ENTRY(v2neg, 0, sz), cc, ctx);
    gr_set_shallow(GR_ENTRY(v2neg, 1, sz), negcs, ctx);

/* XXX: POSSIBLE ALIASING BUG IN ORIGINAL CODE (DEPENDING ON HOW DOT PRODUCT IS IMPLEMENTED) */
#if 1

#define BUTTERFLY(M, Mr1, Mr2, Mc1, Mc2, vv1, vv2) \
    do { \
        gr_set_shallow(GR_ENTRY(v3, 0, sz), GR_MAT_ENTRY(M, Mr1, Mc1, sz), ctx); \
        gr_set_shallow(GR_ENTRY(v3, 1, sz), GR_MAT_ENTRY(M, Mr2, Mc2, sz), ctx); \
        status |= _gr_vec_dot(t, NULL, 0, vv1, v3, 2, ctx); \
        status |= _gr_vec_dot_rev(GR_MAT_ENTRY(M, Mr2, Mc2, sz), NULL, 0, vv2, v3, 2, ctx); \
        gr_swap(GR_MAT_ENTRY(M, Mr1, Mc1, sz), t, ctx); \
    } while (0)

#else

    gr_ptr t2, t3;
    GR_TMP_INIT2(t2, t3, ctx);

#define BUTTERFLY(M, Mr1, Mr2, Mc1, Mc2, vv1, vv2) \
    do { \
        status |= gr_mul(t, GR_ENTRY(vv1, 0, sz), GR_MAT_ENTRY(M, Mr1, Mc1, sz), ctx); \
        status |= gr_mul(t2, GR_ENTRY(vv1, 1, sz), GR_MAT_ENTRY(M, Mr2, Mc2, sz), ctx); \
        status |= gr_mul(t3, GR_ENTRY(vv2, 1, sz), GR_MAT_ENTRY(M, Mr1, Mc1, sz), ctx); \
        status |= gr_add(GR_MAT_ENTRY(M, Mr1, Mc1, sz), t, t2, ctx); \
        status |= gr_mul(t, GR_ENTRY(vv2, 0, sz), GR_MAT_ENTRY(M, Mr2, Mc2, sz), ctx); \
        status |= gr_add(GR_MAT_ENTRY(M, Mr2, Mc2, sz), t3, t, ctx); \
    } while (0)

#endif

    /* Apply Givens rotation from the left */
    for (k = n0; k < n; k++)
    {
        /* x = A[n0  ,k] */
        /* y = A[n0+1,k] */
        /* A[n0,     k] = cc * x + cs * y */
        /* A[n0 + 1, k] = c  * y -  s * x */
        BUTTERFLY(A, n0, n0 + 1, k, k, v2, v1neg);
    }

    /* Apply Givens rotation from the right */
    for (k = 0; k < FLINT_MIN(n1, n0 + 3); k++)
    {
        /* x = A[k,n0  ] */
        /* y = A[k,n0+1] */
        /* A[k,n0  ] = c * x + s * y */
        /* A[k,n0+1] = cc * y - cs * x */
        BUTTERFLY(A, k, k, n0, n0 + 1, v1, v2neg);
    }

    if (Q != NULL)
    {
        for (k = 0; k < n; k++)
        {
            /* x = Q[k,n0  ] */
            /* y = Q[k,n0+1] */
            /* Q[k,n0  ] = c * x + s * y */
            /* Q[k,n0+1] = cc * y - cs * x */
            BUTTERFLY(Q, k, k, n0, n0 + 1, v1, v2neg);
        }
    }

    for (j = n0; j < n1 - 2; j++)
    {
        /* Calculate Givens rotation */
        status |= gr_set(c, GR_MAT_ENTRY(A, j + 1, j, sz), ctx);
        status |= gr_set(s, GR_MAT_ENTRY(A, j + 2, j, sz), ctx);

        status |= gr_sosq(v, RE(c), IM(c), real_ctx);
        status |= gr_sosq(u, RE(s), IM(s), real_ctx);
        status |= gr_add(v, v, u, real_ctx);
        status |= gr_sqrt(v, v, real_ctx);

        if (gr_is_zero(v, real_ctx) == T_TRUE)
        {
            status |= gr_zero(GR_MAT_ENTRY(A, j + 1, j, sz), ctx);
            status |= gr_one(v, real_ctx);
            status |= gr_one(c, ctx);
            status |= gr_zero(s, ctx);
        }
        else
        {
            status |= gr_set(RE(GR_MAT_ENTRY(A, j + 1, j, sz)), v, real_ctx);
            status |= gr_zero(IM(GR_MAT_ENTRY(A, j + 1, j, sz)), real_ctx);
            status |= gr_inv(v, v, real_ctx);
            status |= gr_mul(RE(c), RE(c), v, real_ctx);
            status |= gr_mul(IM(c), IM(c), v, real_ctx);
            status |= gr_mul(RE(s), RE(s), v, real_ctx);
            status |= gr_mul(IM(s), IM(s), v, real_ctx);
        }

//        flint_printf("gr(%wd)   c = %{gr}  s = %{gr}\n\n", j, c, ctx, s, ctx);

        status |= gr_zero(GR_MAT_ENTRY(A, j + 2, j, sz), ctx);

        status |= gr_conj(cc, c, ctx);
        status |= gr_conj(cs, s, ctx);
        status |= gr_neg(negs, s, ctx);
        status |= gr_neg(negcs, cs, ctx);

        gr_set_shallow(GR_ENTRY(v1, 0, sz), c, ctx);
        gr_set_shallow(GR_ENTRY(v1, 1, sz), s, ctx);
        gr_set_shallow(GR_ENTRY(v1neg, 0, sz), c, ctx);
        gr_set_shallow(GR_ENTRY(v1neg, 1, sz), negs, ctx);
        gr_set_shallow(GR_ENTRY(v2, 0, sz), cc, ctx);
        gr_set_shallow(GR_ENTRY(v2, 1, sz), cs, ctx);
        gr_set_shallow(GR_ENTRY(v2neg, 0, sz), cc, ctx);
        gr_set_shallow(GR_ENTRY(v2neg, 1, sz), negcs, ctx);

        /* Apply Givens rotation from the left */
        for (k = j + 1; k < n; k++)
        {
            /* x = A[j+1, k] */
            /* y = A[j+2, k] */
            /* A[j+1, k] = cc * x + cs * y */
            /* A[j+2, k] = c  * y -  s * x */
            BUTTERFLY(A, j + 1, j + 2, k, k, v2, v1neg);
        }

        /* Apply Givens rotation from the right */
        for (k = 0; k < FLINT_MIN(n1, j + 4); k++)
        {
            /* x = A[k,j+1] */
            /* y = A[k,j+2] */
            /* A[k,j+1] = c * x + s * y */
            /* A[k,j+2] = cc * y - cs * x */
            BUTTERFLY(A, k, k, j + 1, j + 2, v1, v2neg);
        }

        if (Q != NULL)
        {
            for (k = 0; k < n; k++)
            {
                /* x = Q[k,j+1] */
                /* y = Q[k,j+2] */
                /* Q[k,j+1] = c * x + s * y */
                /* Q[k,j+2] = cc * y - cs * x */
                BUTTERFLY(Q, k, k, j + 1, j + 2, v1, v2neg);
            }
        }
    }

#undef RE
#undef IM

    GR_TMP_CLEAR5(c, s, negs, cc, cs, ctx);
    GR_TMP_CLEAR2(negcs, t, ctx);
    GR_TMP_CLEAR2(v, u, real_ctx);

//    GR_TMP_CLEAR2(t2, t3, ctx);

    GR_TMP_FREE(v1, 10 * sz);

    return status;
}

int
gr_mat_approx_hessenberg_qr(gr_mat_t A, gr_mat_t Q, slong maxiter, gr_ctx_t ctx, gr_ctx_t real_ctx, gr_srcptr tol, gr_ctx_t mag_ctx)
{
    slong n, i, j, k, n0, n1, iter, total_iter;
    /* Magnitude variables */
    gr_ptr norm, tm, um, eps, ts;
    /* Complex variables */
    gr_ptr shift, s, t, a, b;
    /* Real variables */
    gr_ptr zero;

    slong prec;

    slong sz = ctx->sizeof_elem;
    slong rsz = real_ctx->sizeof_elem;

#define RE(xx) (xx)
#define IM(xx) GR_ENTRY(xx, 1, rsz)

    int status = GR_SUCCESS;

    FLINT_ASSERT(sz == 2 * rsz);

    if (gr_ctx_get_real_prec(&prec, ctx) != GR_SUCCESS)
        return GR_UNABLE;

    n = A->r;

    if (n <= 1)
        return GR_SUCCESS;

    GR_TMP_INIT5(norm, tm, um, eps, ts, mag_ctx);
    GR_TMP_INIT5(shift, s, t, a, b, ctx);
    GR_TMP_INIT(zero, real_ctx);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < FLINT_MIN(i + 2, n); j++)
        {
            status |= gr_set_other(tm, RE(GR_MAT_ENTRY(A, j, i, sz)), real_ctx, mag_ctx);
            status |= gr_set_other(um, IM(GR_MAT_ENTRY(A, j, i, sz)), real_ctx, mag_ctx);
            status |= gr_addmul(norm, tm, tm, mag_ctx);
            status |= gr_addmul(norm, um, um, mag_ctx);
        }
    }

    status |= gr_sqrt(norm, norm, mag_ctx);
    status |= gr_div_ui(norm, norm, n, mag_ctx);

    if (gr_is_zero(norm, mag_ctx) == T_TRUE)
        return GR_SUCCESS;

/*
    if (gr_is_inf(norm) == T_TRUE)
        return GR_UNABLE;
*/

    /* XXX: IN ACB, SHOULD PASS TOL TO GET BIT PRECISION TOLERANCE INSTEAD OF NFLOAT WORD PRECISION TOLERANCE */
    if (tol == NULL)
    {
        status |= gr_one(eps, mag_ctx);
        status |= gr_mul_2exp_si(eps, eps, -prec, mag_ctx);
        status |= gr_div_ui(eps, eps, 100 * n, mag_ctx);
    }
    else
    {
        status |= gr_set(eps, tol, mag_ctx);
    }

    if (maxiter <= 0)
    {
        maxiter = 14 * n;
        maxiter += prec / 10;
    }

    /* The active submatrix is A[n0:n1,n0:n1]. */
    n0 = 0;
    n1 = n;

    iter = total_iter = 0;
    status = GR_SUCCESS;

    while (1)
    {
        k = n0;

        /* flint_printf("total_iter %wd   active %wd\n", total_iter, n1 - n0); */

        while (k + 1 < n1)
        {
            status |= gr_zero(ts, mag_ctx);

            status |= gr_set_other(tm, RE(GR_MAT_ENTRY(A, k, k, sz)), real_ctx, mag_ctx);
            status |= gr_abs(tm, tm, mag_ctx);
            status |= gr_add(ts, ts, tm, mag_ctx);

            status |= gr_set_other(tm, IM(GR_MAT_ENTRY(A, k, k, sz)), real_ctx, mag_ctx);
            status |= gr_abs(tm, tm, mag_ctx);
            status |= gr_add(ts, ts, tm, mag_ctx);

            status |= gr_set_other(tm, RE(GR_MAT_ENTRY(A, k + 1, k + 1, sz)), real_ctx, mag_ctx);
            status |= gr_abs(tm, tm, mag_ctx);
            status |= gr_add(ts, ts, tm, mag_ctx);

            status |= gr_set_other(tm, IM(GR_MAT_ENTRY(A, k + 1, k + 1, sz)), real_ctx, mag_ctx);
            status |= gr_abs(tm, tm, mag_ctx);
            status |= gr_add(ts, ts, tm, mag_ctx);

            /* if s < eps * norm, s = norm */
            status |= gr_mul(tm, eps, norm, mag_ctx);
            if (gr_lt(ts, tm, mag_ctx) == T_TRUE)
                status |= gr_set(ts, norm, mag_ctx);

            /* if abs(A[k+1,k]) < eps * s, break */
            status |= gr_set_other(tm, RE(GR_MAT_ENTRY(A, k + 1, k, sz)), real_ctx, mag_ctx);
            status |= gr_abs(tm, tm, mag_ctx);
            status |= gr_set_other(um, IM(GR_MAT_ENTRY(A, k + 1, k, sz)), real_ctx, mag_ctx);
            status |= gr_abs(um, um, mag_ctx);

            /* todo: gr_hypot ... or square eps to avoid square root */
            status |= gr_sosq(tm, tm, um, mag_ctx);
            status |= gr_sqrt(tm, tm, mag_ctx);

            status |= gr_mul(um, eps, ts, mag_ctx);

            if (gr_lt(tm, um, mag_ctx) == T_TRUE)
                break;

            k++;
        }

        /* Deflation found at position (k+1, k). */
        if (k + 1 < n1)
        {
            status |= gr_zero(GR_MAT_ENTRY(A, k + 1, k, sz), ctx);
            n0 = k + 1;
            iter = 0;

            if (n0 + 1 >= n1)
            {
                /* Block of size at most two has converged. */
                n0 = 0;
                n1 = k + 1;
                if (n1 < 2)
                {
                    /* QR algorithm has converged. */
                    status = GR_SUCCESS;
                    goto cleanup;
                }
            }
        }
        else
        {
            if (iter % 30 == 10)
            {
                /* Exceptional shift */
                status |= gr_set(shift, GR_MAT_ENTRY(A, n1 - 1, n1 - 2, sz), ctx);
            }
            else if (iter % 30 == 20)
            {
                /* Exceptional shift */
                status |= gr_abs(shift, GR_MAT_ENTRY(A, n1 - 1, n1 - 2, sz), ctx);
            }
            else if (iter % 30 == 29)
            {
                /* Exceptional shift */
                status |= gr_set_other(shift, norm, mag_ctx, ctx);
            }
            else
            {
                status |= gr_add(t, GR_MAT_ENTRY(A, n1 - 2, n1 - 2, sz), GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz), ctx);
                status |= gr_sub(a, GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz), GR_MAT_ENTRY(A, n1 - 2, n1 - 2, sz), ctx);
                status |= gr_mul(a, a, a, ctx);
                status |= gr_mul(b, GR_MAT_ENTRY(A, n1 - 1, n1 - 2, sz), GR_MAT_ENTRY(A, n1 - 2, n1 - 1, sz), ctx);
                status |= gr_mul_2exp_si(b, b, 2, ctx);
                status |= gr_add(s, a, b, ctx);

                if (gr_gt(RE(s), zero, real_ctx) == T_TRUE)
                {
                    status |= gr_sqrt(s, s, ctx);
                }
                else
                {
                    status |= gr_neg(s, s, ctx);
                    status |= gr_sqrt(s, s, ctx);
                    /* multiply by i */
                    gr_swap(RE(s), IM(s), real_ctx);
                    status |= gr_neg(RE(s), RE(s), real_ctx);
                }

                status |= gr_add(a, t, s, ctx);
                status |= gr_sub(b, t, s, ctx);
                status |= gr_mul_2exp_si(a, a, -1, ctx);
                status |= gr_mul_2exp_si(b, b, -1, ctx);

                status |= gr_sub(s, GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz), a, ctx);
                status |= gr_sub(t, GR_MAT_ENTRY(A, n1 - 1, n1 - 1, sz), b, ctx);

                status |= gr_set_other(tm, RE(s), real_ctx, mag_ctx);
                status |= gr_set_other(ts, IM(s), real_ctx, mag_ctx);
                /* todo: hypot */
                status |= gr_sosq(tm, tm, ts, mag_ctx);
                status |= gr_sqrt(tm, tm, mag_ctx);

                status |= gr_set_other(um, RE(t), real_ctx, mag_ctx);
                status |= gr_set_other(ts, IM(t), real_ctx, mag_ctx);
                /* todo: hypot */
                status |= gr_sosq(um, um, ts, mag_ctx);
                status |= gr_sqrt(um, um, mag_ctx);

                if (gr_gt(tm, um, mag_ctx) == T_TRUE)
                    status |= gr_set(shift, b, ctx);
                else
                    status |= gr_set(shift, a, ctx);
            }

            iter++;
            total_iter++;

            status |= gr_mat_approx_qr_step(A, Q, n0, n1, shift, ctx, real_ctx);

            if (iter > maxiter)
            {
                status = GR_UNABLE;
                goto cleanup;
            }
        }
    }

cleanup:
    GR_TMP_CLEAR5(norm, tm, um, ts, eps, mag_ctx);
    GR_TMP_CLEAR5(shift, s, t, a, b, ctx);
    GR_TMP_CLEAR(zero, real_ctx);

    return status;
}

/*
Adapted from eigen.py in mpmath (written by Timo Hartmann)

Todo items present in the original code:

- Implement balancing
- Aggressive early deflation
*/

#include "acb.h"
#include "acb_mat.h"

static void
acb_approx_mag(mag_t res, const acb_t x)
{
    mag_t t;
    mag_init(t);
    arf_get_mag(res, arb_midref(acb_realref(x)));
    arf_get_mag(t, arb_midref(acb_imagref(x)));
    mag_hypot(res, res, t);
    mag_clear(t);
}

static void
acb_approx_mul(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_complex_mul(arb_midref(acb_realref(res)), arb_midref(acb_imagref(res)),
        arb_midref(acb_realref(x)), arb_midref(acb_imagref(x)),
        arb_midref(acb_realref(y)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static void
acb_approx_add(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_add(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_add(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static void
acb_approx_sub(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_sub(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_sub(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static void
acb_approx_set(acb_t res, const acb_t x)
{
    arf_set(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)));
}

static void
acb_approx_div_arb(acb_t res, const acb_t x, const arb_t y, slong prec)
{
    arf_div(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(y), prec, ARF_RND_DOWN);
    arf_div(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(y), prec, ARF_RND_DOWN);
}

static void
acb_approx_inv(acb_t z, const acb_t x, slong prec)
{
    arf_set(arb_midref(acb_realref(z)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(z)), arb_midref(acb_imagref(x)));

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));

    acb_inv(z, z, prec);

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));
}

static void
acb_approx_div(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_approx_inv(t, y, prec);
    acb_approx_mul(z, x, t, prec);
    acb_clear(t);
}

#include "nfloat.h"

void
acb_mat_approx_qr_step1(acb_mat_t A, acb_mat_t Q, slong n0, slong n1, const acb_t shift, slong prec)
{
    gr_mat_t T, U;
    gr_ctx_t ctx, real_ctx;
    gr_ctx_t acb_ctx;
    gr_ptr tshift;

    GR_MUST_SUCCEED(nfloat_ctx_init(real_ctx, prec, 0));
    GR_MUST_SUCCEED(nfloat_complex_ctx_init(ctx, prec, 0));

    gr_ctx_init_complex_acb(acb_ctx, prec);
    GR_TMP_INIT(tshift, ctx);

    GR_MUST_SUCCEED(gr_set_other(tshift, shift, acb_ctx, ctx));

    gr_mat_init(T, A->r, A->c, ctx);
    GR_MUST_SUCCEED(gr_mat_set_gr_mat_other(T, (gr_mat_struct *) A, acb_ctx, ctx));

    if (Q != NULL)
    {
        gr_mat_init(U, Q->r, Q->c, ctx);
        GR_MUST_SUCCEED(gr_mat_set_gr_mat_other(U, (gr_mat_struct *) Q, acb_ctx, ctx));
    }

    GR_MUST_SUCCEED(gr_mat_approx_qr_step(T, Q != NULL ? U : NULL, n0, n1, tshift, ctx, real_ctx));

    GR_MUST_SUCCEED(gr_mat_set_gr_mat_other((gr_mat_struct *) A, T, ctx, acb_ctx));
    if (Q != NULL)
        GR_MUST_SUCCEED(gr_mat_set_gr_mat_other((gr_mat_struct *) Q, U, ctx, acb_ctx));

    gr_mat_clear(T, ctx);
    if (Q != NULL)
        gr_mat_clear(U, ctx);

    GR_TMP_CLEAR(tshift, ctx);
}

int
acb_mat_approx_hessenberg_qr2(acb_mat_t A, acb_mat_t Q, const mag_t tol, slong maxiter, slong prec)
{
    gr_mat_t T, U;
    gr_ctx_t ctx, real_ctx;
    gr_ctx_t acb_ctx;
    gr_ctx_t mag_ctx;
    gr_ptr ttol;
    int result;

    GR_MUST_SUCCEED(nfloat_ctx_init(mag_ctx, MAG_BITS, 0));
    GR_MUST_SUCCEED(nfloat_ctx_init(real_ctx, prec, 0));
    GR_MUST_SUCCEED(nfloat_complex_ctx_init(ctx, prec, 0));

    gr_ctx_init_complex_acb(acb_ctx, prec);
    GR_TMP_INIT(ttol, mag_ctx);

    if (tol != NULL)
    {
        /* XXX: correct this */
        GR_MUST_SUCCEED(gr_set_ui(ttol, MAG_MAN(tol), mag_ctx));
        GR_MUST_SUCCEED(gr_mul_2exp_si(ttol, ttol, MAG_EXP(tol) - MAG_BITS, mag_ctx));
    }

    gr_mat_init(T, A->r, A->c, ctx);
    GR_MUST_SUCCEED(gr_mat_set_gr_mat_other(T, (gr_mat_struct *) A, acb_ctx, ctx));

    if (Q != NULL)
    {
        gr_mat_init(U, Q->r, Q->c, ctx);
        GR_MUST_SUCCEED(gr_mat_set_gr_mat_other(U, (gr_mat_struct *) Q, acb_ctx, ctx));
    }

    result = (GR_SUCCESS == gr_mat_approx_hessenberg_qr(T, Q != NULL ? U : NULL, maxiter, ctx, real_ctx, tol != NULL ? ttol : NULL, mag_ctx));

    GR_MUST_SUCCEED(gr_mat_set_gr_mat_other((gr_mat_struct *) A, T, ctx, acb_ctx));
    if (Q != NULL)
        GR_MUST_SUCCEED(gr_mat_set_gr_mat_other((gr_mat_struct *) Q, U, ctx, acb_ctx));

    gr_mat_clear(T, ctx);
    if (Q != NULL)
        gr_mat_clear(U, ctx);

    GR_TMP_CLEAR(ttol, mag_ctx);

    return result;
}

void
acb_mat_approx_qr_step2(acb_mat_t A, acb_mat_t Q, slong n0, slong n1, const acb_t shift, slong prec)
{
    slong j, k, n;
    acb_t c, s, negs, cc, cs, negcs, t;
    acb_struct v1[2];
    acb_struct v1neg[2];
    acb_struct v2[2];
    acb_struct v2neg[2];
    acb_struct v3[2];
    arb_t v, u;

    n = acb_mat_nrows(A);

    acb_init(c);
    acb_init(s);
    acb_init(negs);
    acb_init(cc);
    acb_init(cs);
    acb_init(negcs);
    acb_init(t);
    arb_init(v);
    arb_init(u);

    /* Calculate Givens rotation */
    acb_approx_sub(c, acb_mat_entry(A, n0, n0), shift, prec);
    acb_approx_set(s, acb_mat_entry(A, n0 + 1, n0));

    arf_sosq(arb_midref(v), arb_midref(acb_realref(c)), arb_midref(acb_imagref(c)), prec, ARF_RND_DOWN);
    arf_sosq(arb_midref(u), arb_midref(acb_realref(s)), arb_midref(acb_imagref(s)), prec, ARF_RND_DOWN);
    arf_add(arb_midref(v), arb_midref(v), arb_midref(u), prec, ARF_RND_DOWN);
    arf_sqrt(arb_midref(v), arb_midref(v), prec, ARF_RND_DOWN);

    if (arb_is_zero(v))
    {
        arb_one(v);
        acb_one(c);
        acb_zero(s);
    }
    else
    {
        acb_approx_div_arb(c, c, v, prec);
        acb_approx_div_arb(s, s, v, prec);
    }

//    flint_printf("acb  c = %{acb}  s = %{acb}\n\n", c, s);

    acb_conj(cc, c);
    acb_conj(cs, s);
    acb_neg(negs, s);
    acb_neg(negcs, cs);

    v1[0] = *c;
    v1[1] = *s;
    v1neg[0] = *c;
    v1neg[1] = *negs;
    v2[0] = *cc;
    v2[1] = *cs;
    v2neg[0] = *cc;
    v2neg[1] = *negcs;

    /* Apply Givens rotation from the left */
    for (k = n0; k < n; k++)
    {


        v3[0] = *acb_mat_entry(A, n0, k);
        v3[1] = *acb_mat_entry(A, n0 + 1, k);

        /* x = A[n0  ,k] */
        /* y = A[n0+1,k] */
        /* A[n0,     k] = cc * x + cs * y */
        /* A[n0 + 1, k] = c  * y -  s * x */

//        flint_printf("acb(%wd)  %{acb*} %{acb*}\n", k, v2, 2, v3, 2);

        acb_approx_dot(t,                           NULL, 0, v2,    1, v3,      1, 2, prec);
        acb_approx_dot(acb_mat_entry(A, n0 + 1, k), NULL, 0, v1neg, 1, v3 + 1, -1, 2, prec);
        acb_swap(t, acb_mat_entry(A, n0, k));

//        flint_printf("acb(%wd)  %{acb}  %{acb}\n", k, acb_mat_entry(A, n0, k), acb_mat_entry(A, n0 + 1, k));
    }

    /* Apply Givens rotation from the right */
    for (k = 0; k < FLINT_MIN(n1, n0 + 3); k++)
    {
        /* x = A[k,n0  ] */
        /* y = A[k,n0+1] */
        /* A[k,n0  ] = c * x + s * y */
        /* A[k,n0+1] = cc * y - cs * x */

        v3[0] = *acb_mat_entry(A, k, n0);
        v3[1] = *acb_mat_entry(A, k, n0 + 1);

        acb_approx_dot(t,                           NULL, 0, v1,    1, v3,      1, 2, prec);
        acb_approx_dot(acb_mat_entry(A, k, n0 + 1), NULL, 0, v2neg, 1, v3 + 1, -1, 2, prec);
        acb_swap(t, acb_mat_entry(A, k, n0));
    }

    if (Q != NULL)
    {
        for (k = 0; k < n; k++)
        {
            /* x = Q[k,n0  ] */
            /* y = Q[k,n0+1] */
            /* Q[k,n0  ] = c * x + s * y */
            /* Q[k,n0+1] = cc * y - cs * x */

            v3[0] = *acb_mat_entry(Q, k, n0);
            v3[1] = *acb_mat_entry(Q, k, n0 + 1);

            acb_approx_dot(t,                           NULL, 0, v1,    1, v3,      1, 2, prec);
            acb_approx_dot(acb_mat_entry(Q, k, n0 + 1), NULL, 0, v2neg, 1, v3 + 1, -1, 2, prec);
            acb_swap(t, acb_mat_entry(Q, k, n0));
        }
    }

    for (j = n0; j < n1 - 2; j++)
    {
        /* Calculate Givens rotation */
        acb_set(c, acb_mat_entry(A, j + 1, j));
        acb_set(s, acb_mat_entry(A, j + 2, j));

        arf_sosq(arb_midref(v), arb_midref(acb_realref(c)), arb_midref(acb_imagref(c)), prec, ARF_RND_DOWN);
        arf_sosq(arb_midref(u), arb_midref(acb_realref(s)), arb_midref(acb_imagref(s)), prec, ARF_RND_DOWN);
        arf_add(arb_midref(v), arb_midref(v), arb_midref(u), prec, ARF_RND_DOWN);
        arf_sqrt(arb_midref(v), arb_midref(v), prec, ARF_RND_DOWN);

        if (arb_is_zero(v))
        {
            acb_zero(acb_mat_entry(A, j + 1, j));
            arb_one(v);
            acb_one(c);
            acb_zero(s);
        }
        else
        {
            acb_set_arb(acb_mat_entry(A, j + 1, j), v);
            acb_approx_div_arb(c, c, v, prec);
            acb_approx_div_arb(s, s, v, prec);
        }

//        flint_printf("acb(%wd)  c = %{acb}  s = %{acb}\n\n", j, c, s);

        acb_zero(acb_mat_entry(A, j + 2, j));

        acb_conj(cc, c);
        acb_conj(cs, s);
        acb_neg(negs, s);
        acb_neg(negcs, cs);

        v1[0] = *c;
        v1[1] = *s;
        v1neg[0] = *c;
        v1neg[1] = *negs;
        v2[0] = *cc;
        v2[1] = *cs;
        v2neg[0] = *cc;
        v2neg[1] = *negcs;

        /* Apply Givens rotation from the left */
        for (k = j + 1; k < n; k++)
        {
            v3[0] = *acb_mat_entry(A, j + 1, k);
            v3[1] = *acb_mat_entry(A, j + 2, k);

            /* x = A[j+1, k] */
            /* y = A[j+2, k] */
            /* A[j + 1, k] = cc * x + cs * y */
            /* A[j + 2, k] = c  * y -  s * x */

            acb_approx_dot(t,                          NULL, 0, v2,    1, v3,      1, 2, prec);
            acb_approx_dot(acb_mat_entry(A, j + 2, k), NULL, 0, v1neg, 1, v3 + 1, -1, 2, prec);
            acb_swap(t, acb_mat_entry(A, j + 1, k));
        }

        /* Apply Givens rotation from the right */
        for (k = 0; k < FLINT_MIN(n1, j + 4); k++)
        {
            /* x = A[k,j+1] */
            /* y = A[k,j+2] */
            /* A[k,j+1] = c * x + s * y */
            /* A[k,j+2] = cc * y - cs * x */

            v3[0] = *acb_mat_entry(A, k, j + 1);
            v3[1] = *acb_mat_entry(A, k, j + 2);

            acb_approx_dot(t,                          NULL, 0, v1,    1, v3,      1, 2, prec);
            acb_approx_dot(acb_mat_entry(A, k, j + 2), NULL, 0, v2neg, 1, v3 + 1, -1, 2, prec);
            acb_swap(t, acb_mat_entry(A, k, j + 1));
        }

        if (Q != NULL)
        {
            for (k = 0; k < n; k++)
            {
                /* x = Q[k,j+1] */
                /* y = Q[k,j+2] */
                /* Q[k,j+1] = c * x + s * y */
                /* Q[k,j+2] = cc * y - cs * x */

                v3[0] = *acb_mat_entry(Q, k, j + 1);
                v3[1] = *acb_mat_entry(Q, k, j + 2);

                acb_approx_dot(t,                          NULL, 0, v1,    1, v3,      1, 2, prec);
                acb_approx_dot(acb_mat_entry(Q, k, j + 2), NULL, 0, v2neg, 1, v3 + 1, -1, 2, prec);
                acb_swap(t, acb_mat_entry(Q, k, j + 1));
            }
        }
    }

    acb_clear(c);
    acb_clear(s);
    acb_clear(negs);
    acb_clear(cc);
    acb_clear(cs);
    acb_clear(negcs);
    acb_clear(t);
    arb_clear(v);
    arb_clear(u);
}

#include "profiler.h"

void
acb_mat_approx_qr_step(acb_mat_t A, acb_mat_t Q, slong n0, slong n1, const acb_t shift, slong prec)
{
//    acb_mat_approx_qr_step1(A, Q, n0, n1, shift, prec);
//    return;

    if (0 && Q != NULL && prec >= 512)
    {
        flint_printf("prec %wd   %wd\n", prec, A->r);
//        acb_mat_printd(A, 15); flint_printf("\n\n");
//        acb_mat_printd(Q, 15); flint_printf("\n\n");

        acb_mat_t A2, Q2;
        acb_mat_init(A2, A->r, A->c);
        acb_mat_init(Q2, Q->r, Q->c);

        TIMEIT_START
        acb_mat_set(A2, A);
        acb_mat_set(Q2, Q);
        acb_mat_approx_qr_step2(A2, Q2, n0, n1, shift, prec);
        TIMEIT_STOP
        TIMEIT_START
        acb_mat_set(A2, A);
        acb_mat_set(Q2, Q);
        acb_mat_approx_qr_step1(A2, Q2, n0, n1, shift, prec);
        TIMEIT_STOP

        acb_mat_approx_qr_step1(A, Q, n0, n1, shift, prec);

//        flint_printf("A2: ");
//        acb_mat_printd(A2, 15); flint_printf("\n\n");
//        flint_printf("A: ");
//        acb_mat_printd(A, 15); flint_printf("\n\n");
//        flint_printf("Q2: ");
//        acb_mat_printd(Q2, 15); flint_printf("\n\n");
//        flint_printf("Q: ");
//        acb_mat_printd(Q, 15); flint_printf("\n\n");
//        flint_abort();

        acb_mat_clear(A2);
        acb_mat_clear(Q2);
    }
    else
    {
        acb_mat_approx_qr_step2(A, Q, n0, n1, shift, prec);
    }
}

void
acb_mat_approx_hessenberg_reduce_0(acb_mat_t A, acb_ptr T, slong prec)
{
    slong i, j, k, n;
    arf_t scale, scale_inv, tt, H, G, f;
    acb_ptr V1, V2;
    acb_t ff, GG, TT;
    acb_t F;

    n = acb_mat_nrows(A);
    if (n <= 2)
        return;

    arf_init(scale);
    arf_init(scale_inv);
    arf_init(tt);
    arf_init(H);
    arf_init(G);
    arf_init(f);
    acb_init(F);
    V1 = _acb_vec_init(n + 1);
    V2 = _acb_vec_init(n + 1);
    acb_init(ff);
    acb_init(GG);
    acb_init(TT);

    for (i = n - 1; i >= 2; i--)
    {
        /* Scale the vector (todo: is this needed?) */
        arf_zero(scale);
        for (k = 0; k < i; k++)
        {
            arf_abs(tt, arb_midref(acb_realref(acb_mat_entry(A, i, k))));
            arf_add(scale, scale, tt, prec, ARF_RND_DOWN);
            arf_abs(tt, arb_midref(acb_imagref(acb_mat_entry(A, i, k))));
            arf_add(scale, scale, tt, prec, ARF_RND_DOWN);
        }

        arf_ui_div(scale_inv, 1, scale, prec, ARF_RND_DOWN);

        if (arf_is_zero(scale))
        {
            acb_zero(T + i);
            acb_zero(acb_mat_entry(A, i, i - 1));
            continue;
        }

        /* Calculate parameters for Householder transformation. */
        arf_zero(H);
        for (k = 0; k < i; k++)
        {
            arf_ptr Aikr, Aiki;

            Aikr = arb_midref(acb_realref(acb_mat_entry(A, i, k)));
            Aiki = arb_midref(acb_imagref(acb_mat_entry(A, i, k)));

            arf_mul(Aikr, Aikr, scale_inv, prec, ARF_RND_DOWN);
            arf_mul(Aiki, Aiki, scale_inv, prec, ARF_RND_DOWN);

            arf_addmul(H, Aikr, Aikr, prec, ARF_RND_DOWN);
            arf_addmul(H, Aiki, Aiki, prec, ARF_RND_DOWN);
        }

        acb_set(F, acb_mat_entry(A, i, i - 1));

        /* f = abs(F) */
        arf_mul(f, arb_midref(acb_realref(F)), arb_midref(acb_realref(F)), prec, ARF_RND_DOWN);
        arf_addmul(f, arb_midref(acb_imagref(F)), arb_midref(acb_imagref(F)), prec, ARF_RND_DOWN);
        arf_sqrt(f, f, prec, ARF_RND_DOWN);

        arf_sqrt(G, H, prec, ARF_RND_DOWN);

        /* A[i,i-1] = -G scale */
        arf_mul(arb_midref(acb_realref(acb_mat_entry(A, i, i - 1))), G, scale, prec, ARF_RND_DOWN);
        arf_neg(arb_midref(acb_realref(acb_mat_entry(A, i, i - 1))), arb_midref(acb_realref(acb_mat_entry(A, i, i - 1))));
        arf_zero(arb_midref(acb_imagref(acb_mat_entry(A, i, i - 1))));

        if (arf_is_zero(f))
        {
            arb_set_arf(acb_realref(T + i), G);
            arb_zero(acb_imagref(T + i));
        }
        else
        {
            /* ff = F / f */
            arf_div(arb_midref(acb_realref(ff)), arb_midref(acb_realref(F)), f, prec, ARF_RND_DOWN);
            arf_div(arb_midref(acb_imagref(ff)), arb_midref(acb_imagref(F)), f, prec, ARF_RND_DOWN);

            /* T[i] = F + G ff */
            acb_set(T + i, F);
            arf_addmul(arb_midref(acb_realref(T + i)), arb_midref(acb_realref(ff)), G, prec, ARF_RND_DOWN);
            arf_addmul(arb_midref(acb_imagref(T + i)), arb_midref(acb_imagref(ff)), G, prec, ARF_RND_DOWN);

            /* A[i,i-1] *= ff */
            acb_approx_mul(acb_mat_entry(A, i, i - 1), acb_mat_entry(A, i, i - 1), ff, prec);
        }

        arf_addmul(H, G, f, prec, ARF_RND_DOWN);
        arf_rsqrt(H, H, prec, ARF_RND_DOWN);

        arf_mul(arb_midref(acb_realref(T + i)), arb_midref(acb_realref(T + i)), H, prec, ARF_RND_DOWN);
        arf_mul(arb_midref(acb_imagref(T + i)), arb_midref(acb_imagref(T + i)), H, prec, ARF_RND_DOWN);

        for (k = 0; k < i - 1; k++)
        {
            arf_mul(arb_midref(acb_realref(acb_mat_entry(A, i, k))),
                    arb_midref(acb_realref(acb_mat_entry(A, i, k))),
                        H, prec, ARF_RND_DOWN);
            arf_mul(arb_midref(acb_imagref(acb_mat_entry(A, i, k))),
                    arb_midref(acb_imagref(acb_mat_entry(A, i, k))),
                        H, prec, ARF_RND_DOWN);
        }

        /* todo: optimize copies below */
        /* todo: conj mid etc... */

        /* Apply Householder transformation (from the right). */
        for (j = 0; j < i; j++)
        {
            acb_conj(V1, T + i);
            acb_set(V2, acb_mat_entry(A, j, i - 1));
            for (k = 0; k < i - 1; k++)
            {
                acb_conj(V1 + k + 1, acb_mat_entry(A, i, k));
                acb_set(V2 + k + 1, acb_mat_entry(A, j, k));
            }
            acb_approx_dot(GG, NULL, 0, V1, 1, V2, 1, i, prec);

            acb_approx_mul(TT, GG, T + i, prec);
            acb_approx_sub(acb_mat_entry(A, j, i - 1),
                            acb_mat_entry(A, j, i - 1), TT, prec);
            for (k = 0; k < i - 1; k++)
            {
                acb_approx_mul(TT, GG, acb_mat_entry(A, i, k), prec);
                acb_approx_sub(acb_mat_entry(A, j, k),
                                acb_mat_entry(A, j, k), TT, prec);
            }
        }

        for (j = 0; j < n; j++)
        {
            acb_set(V1, T + i);
            acb_set(V2, acb_mat_entry(A, i - 1, j));
            for (k = 0; k < i - 1; k++)
            {
                acb_set(V1 + k + 1, acb_mat_entry(A, i, k));
                acb_set(V2 + k + 1, acb_mat_entry(A, k, j));
            }
            acb_approx_dot(GG, NULL, 0, V1, 1, V2, 1, i, prec);

            acb_conj(TT, T + i);
            acb_approx_mul(TT, GG, TT, prec);
            acb_approx_sub(acb_mat_entry(A, i - 1, j),
                            acb_mat_entry(A, i - 1, j), TT, prec);
            for (k = 0; k < i - 1; k++)
            {
                acb_conj(TT, acb_mat_entry(A, i, k));
                acb_approx_mul(TT, GG, TT, prec);
                acb_approx_sub(acb_mat_entry(A, k, j),
                                acb_mat_entry(A, k, j), TT, prec);
            }
        }
    }

    arf_clear(scale);
    arf_clear(scale_inv);
    arf_clear(tt);
    arf_clear(H);
    arf_clear(G);
    arf_clear(f);
    acb_clear(F);
    _acb_vec_clear(V1, n + 1);
    _acb_vec_clear(V2, n + 1);
    acb_clear(ff);
    acb_clear(GG);
    acb_clear(TT);
}

void
acb_mat_approx_hessenberg_reduce_1(acb_mat_t A, acb_srcptr T, slong prec)
{
    slong i, j, k, n;
    acb_t G, t;

    n = acb_mat_nrows(A);

    if (n <= 1)
    {
        if (n == 1)
            acb_one(acb_mat_entry(A, 0, 0));
        return;
    }

    acb_one(acb_mat_entry(A, 0, 0));
    acb_one(acb_mat_entry(A, 1, 1));
    acb_zero(acb_mat_entry(A, 0, 1));
    acb_zero(acb_mat_entry(A, 1, 0));

    acb_init(G);
    acb_init(t);

    for (i = 2; i < n; i++)
    {
        if (!acb_is_zero(T + i))
        {
            /* todo: rewrite using approx_dot */
            for (j = 0; j < i; j++)
            {
                acb_approx_mul(G, T + i, acb_mat_entry(A, i - 1, j), prec);
                for (k = 0; k < i - 1; k++)
                {
                    acb_approx_mul(t, acb_mat_entry(A, i, k), acb_mat_entry(A, k, j), prec);
                    acb_approx_add(G, G, t, prec);
                }

                acb_conj(t, T + i);
                acb_approx_mul(t, G, t, prec);
                acb_approx_sub(acb_mat_entry(A, i - 1, j), acb_mat_entry(A, i - 1, j), t, prec);
                for (k = 0; k < i - 1; k++)
                {
                    acb_conj(t, acb_mat_entry(A, i, k));
                    acb_approx_mul(t, G, t, prec);
                    acb_approx_sub(acb_mat_entry(A, k, j), acb_mat_entry(A, k, j), t, prec);
                }
            }
        }

        acb_one(acb_mat_entry(A, i, i));
        for (j = 0; j < i; j++)
        {
            acb_zero(acb_mat_entry(A, j, i));
            acb_zero(acb_mat_entry(A, i, j));
        }
    }

    acb_clear(G);
    acb_clear(t);
}

/* Right eigenvectors of a triu matrix. No aliasing. */
void
acb_mat_approx_eig_triu_r(acb_mat_t ER, const acb_mat_t A, slong prec)
{
    slong i, j, k, n;
    mag_t tm, smin, unfl, simin, smlnum, rmax;
    acb_t r, s, t;

    n = acb_mat_nrows(A);

    acb_mat_one(ER);

    acb_init(r);
    acb_init(s);
    acb_init(t);
    mag_init(tm);
    mag_init(smin);
    mag_init(smlnum);
    mag_init(unfl);
    mag_init(simin);
    mag_init(rmax);

    mag_set_ui_2exp_si(unfl, 1, -30 * prec);
    mag_mul_ui(smlnum, unfl, n);
    mag_mul_2exp_si(smlnum, smlnum, prec);
    mag_set_ui_2exp_si(simin, 1, prec / 2);
    mag_one(rmax);

    for (i = 1; i < n; i++)
    {
        acb_set(s, acb_mat_entry(A, i, i));

        /* smin = max(eps * abs(s), smlnum) */
        acb_approx_mag(smin, s);
        mag_mul_2exp_si(smin, smin, -prec);
        mag_max(smin, smin, smlnum);

        for (j = i - 1; j >= 0; j--)
        {
            acb_approx_dot(r, NULL, 0, acb_mat_entry(A, j, j + 1), 1, acb_mat_entry(ER, i, j + 1), 1, i - j, prec);
            acb_approx_sub(t, acb_mat_entry(A, j, j), s, prec);

            /* if abs(t) < smin: t = smin */
            acb_approx_mag(tm, t);
            if (mag_cmp(tm, smin) < 0)
            {
                acb_zero(t);
                arf_set_mag(arb_midref(acb_realref(t)), smin);
            }

            acb_approx_div(acb_mat_entry(ER, i, j), r, t, prec);
            acb_neg(acb_mat_entry(ER, i, j), acb_mat_entry(ER, i, j));

            acb_approx_mag(tm, r);
            mag_max(rmax, rmax, tm);
            if (mag_cmp(rmax, simin) > 0)
            {
                arb_t b;
                arb_init(b);
                arf_set_mag(arb_midref(b), rmax);

                for (k = j; k < i + 1; k++)
                {
                    acb_approx_div_arb(acb_mat_entry(ER, i, k),
                        acb_mat_entry(ER, i, k), b, prec);
                }

                mag_one(rmax);
                arb_clear(b);
            }
        }

        if (mag_cmp_2exp_si(rmax, 0) != 0)
        {
            arb_t b;
            arb_init(b);
            arf_set_mag(arb_midref(b), rmax);

            for (k = 0; k < i + 1; k++)
            {
                acb_approx_div_arb(acb_mat_entry(ER, i, k),
                    acb_mat_entry(ER, i, k), b, prec);
            }

            arb_clear(b);
        }
    }

    acb_mat_transpose(ER, ER);

    acb_clear(r);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tm);
    mag_clear(smin);
    mag_clear(smlnum);
    mag_clear(unfl);
    mag_clear(simin);
    mag_clear(rmax);
}

/* Left eigenvectors of a triu matrix. No aliasing. */
void
acb_mat_approx_eig_triu_l(acb_mat_t EL, const acb_mat_t A, slong prec)
{
    slong i, j, k, n;
    mag_t tm, smin, unfl, simin, smlnum, rmax;
    acb_t r, s, t;
    acb_mat_t AT;

    n = acb_mat_nrows(A);
    acb_mat_init(AT, n, n);

    acb_mat_one(EL);
    acb_mat_transpose(AT, A);

    acb_init(r);
    acb_init(s);
    acb_init(t);
    mag_init(tm);
    mag_init(smin);
    mag_init(smlnum);
    mag_init(unfl);
    mag_init(simin);
    mag_init(rmax);

    mag_set_ui_2exp_si(unfl, 1, -30 * prec);
    mag_mul_ui(smlnum, unfl, n);
    mag_mul_2exp_si(smlnum, smlnum, prec);
    mag_set_ui_2exp_si(simin, 1, prec / 2);
    mag_one(rmax);

    for (i = 0; i < n - 1; i++)
    {
        acb_set(s, acb_mat_entry(AT, i, i));

        /* smin = max(eps * abs(s), smlnum) */
        acb_approx_mag(smin, s);
        mag_mul_2exp_si(smin, smin, -prec);
        mag_max(smin, smin, smlnum);

        for (j = i + 1; j < n; j++)
        {
            acb_approx_dot(r, NULL, 0, acb_mat_entry(EL, i, i), 1, acb_mat_entry(AT, j, i), 1, j - i, prec);
            acb_approx_sub(t, acb_mat_entry(AT, j, j), s, prec);

            /* if abs(t) < smin: t = smin */
            acb_approx_mag(tm, t);
            if (mag_cmp(tm, smin) < 0)
            {
                acb_zero(t);
                arf_set_mag(arb_midref(acb_realref(t)), smin);
            }

            acb_approx_div(acb_mat_entry(EL, i, j), r, t, prec);
            acb_neg(acb_mat_entry(EL, i, j), acb_mat_entry(EL, i, j));

            acb_approx_mag(tm, r);
            mag_max(rmax, rmax, tm);
            if (mag_cmp(rmax, simin) > 0)
            {
                arb_t b;
                arb_init(b);
                arf_set_mag(arb_midref(b), rmax);

                for (k = i; k < j + 1; k++)
                {
                    acb_approx_div_arb(acb_mat_entry(EL, i, k),
                        acb_mat_entry(EL, i, k), b, prec);
                }

                mag_one(rmax);
                arb_clear(b);
            }
        }

        if (mag_cmp_2exp_si(rmax, 0) != 0)
        {
            arb_t b;
            arb_init(b);
            arf_set_mag(arb_midref(b), rmax);

            for (k = i; k < n; k++)
            {
                acb_approx_div_arb(acb_mat_entry(EL, i, k),
                    acb_mat_entry(EL, i, k), b, prec);
            }

            arb_clear(b);
        }
    }

    acb_mat_clear(AT);

    acb_clear(r);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tm);
    mag_clear(smin);
    mag_clear(smlnum);
    mag_clear(unfl);
    mag_clear(simin);
    mag_clear(rmax);
}

int
acb_mat_approx_hessenberg_qr(acb_mat_t A, acb_mat_t Q, const mag_t tol, slong maxiter, slong prec)
{
    slong n, i, j, k, n0, n1, iter, total_iter;
    mag_t norm, tm, um, eps, ts;
    acb_t shift, s, t, a, b;
    int result;

    n = acb_mat_nrows(A);

    if (n <= 1)
        return 1;

    mag_init(norm);
    mag_init(tm);
    mag_init(um);
    mag_init(eps);
    mag_init(ts);
    acb_init(shift);
    acb_init(s);
    acb_init(t);
    acb_init(a);
    acb_init(b);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < FLINT_MIN(i + 2, n); j++)
        {
            arf_get_mag(tm, arb_midref(acb_realref(acb_mat_entry(A, j, i))));
            arf_get_mag(um, arb_midref(acb_imagref(acb_mat_entry(A, j, i))));
            mag_addmul(norm, tm, tm);
            mag_addmul(norm, um, um);
        }
    }

    mag_sqrt(norm, norm);
    mag_div_ui(norm, norm, n);

    if (mag_is_zero(norm))
        return 1;

    if (mag_is_inf(norm))
        return 0;

    if (tol == NULL)
    {
        mag_one(eps);
        mag_mul_2exp_si(eps, eps, -prec);
        mag_div_ui(eps, eps, 100 * n);
    }
    else
    {
        mag_set(eps, tol);
    }

    if (maxiter <= 0)
    {
        maxiter = 14 * n;
        maxiter += prec / 10;
    }

    /* The active submatrix is A[n0:n1,n0:n1]. */
    n0 = 0;
    n1 = n;

    iter = total_iter = 0;
    result = 0;

    while (1)
    {
        k = n0;

        /* flint_printf("total_iter %wd   active %wd\n", total_iter, n1 - n0); */

        while (k + 1 < n1)
        {
            mag_zero(ts);
            arf_get_mag(tm, arb_midref(acb_realref(acb_mat_entry(A, k, k))));
            mag_add(ts, ts, tm);
            arf_get_mag(tm, arb_midref(acb_imagref(acb_mat_entry(A, k, k))));
            mag_add(ts, ts, tm);
            arf_get_mag(tm, arb_midref(acb_realref(acb_mat_entry(A, k + 1, k + 1))));
            mag_add(ts, ts, tm);
            arf_get_mag(tm, arb_midref(acb_imagref(acb_mat_entry(A, k + 1, k + 1))));
            mag_add(ts, ts, tm);

            /* if s < eps * norm, s = norm */
            mag_mul(tm, eps, norm);
            if (mag_cmp(ts, tm) < 0)
                mag_set(ts, norm);

            /* if abs(A[k+1,k]) < eps * s, break */
            arf_get_mag(tm, arb_midref(acb_realref(acb_mat_entry(A, k + 1, k))));
            arf_get_mag(um, arb_midref(acb_imagref(acb_mat_entry(A, k + 1, k))));
            mag_hypot(tm, tm, um);
            mag_mul(um, eps, ts);
            if (mag_cmp(tm, um) < 0)
                break;

            k++;
        }

        /* Deflation found at position (k+1, k). */
        if (k + 1 < n1)
        {
            acb_zero(acb_mat_entry(A, k + 1, k));
            n0 = k + 1;
            iter = 0;

            if (n0 + 1 >= n1)
            {
                /* Block of size at most two has converged. */
                n0 = 0;
                n1 = k + 1;
                if (n1 < 2)
                {
                    /* QR algorithm has converged. */
                    result = 1;
                    goto cleanup;
                }
            }
        }
        else
        {
            if (iter % 30 == 10)
            {
                /* Exceptional shift */
                acb_set(shift, acb_mat_entry(A, n1 - 1, n1 - 2));
            }
            else if (iter % 30 == 20)
            {
                /* Exceptional shift */
                acb_abs(acb_realref(shift), acb_mat_entry(A, n1 - 1, n1 - 2), prec);
                arb_zero(acb_imagref(shift));
            }
            else if (iter % 30 == 29)
            {
                /* Exceptional shift */
                acb_zero(shift);
                arf_set_mag(arb_midref(acb_realref(shift)), norm);
            }
            else
            {
                acb_approx_add(t, acb_mat_entry(A, n1 - 2, n1 - 2), acb_mat_entry(A, n1 - 1, n1 - 1), prec);
                acb_approx_sub(a, acb_mat_entry(A, n1 - 1, n1 - 1), acb_mat_entry(A, n1 - 2, n1 - 2), prec);
                acb_approx_mul(a, a, a, prec);
                acb_approx_mul(b, acb_mat_entry(A, n1 - 1, n1 - 2), acb_mat_entry(A, n1 - 2, n1 - 1), prec);
                acb_mul_2exp_si(b, b, 2);
                acb_approx_add(s, a, b, prec);

                if (arb_is_positive(acb_realref(s)))
                {
                    acb_sqrt(s, s, prec);
                    acb_get_mid(s, s);
                }
                else
                {
                    acb_neg(s, s);
                    acb_sqrt(s, s, prec);
                    acb_get_mid(s, s);
                    acb_mul_onei(s, s);
                }

                acb_approx_add(a, t, s, prec);
                acb_approx_sub(b, t, s, prec);
                acb_mul_2exp_si(a, a, -1);
                acb_mul_2exp_si(b, b, -1);

                acb_approx_sub(s, acb_mat_entry(A, n1 - 1, n1 - 1), a, prec);
                acb_approx_sub(t, acb_mat_entry(A, n1 - 1, n1 - 1), b, prec);
                acb_get_mag(tm, s);
                acb_get_mag(um, t);

                if (mag_cmp(tm, um) > 0)
                    acb_set(shift, b);
                else
                    acb_set(shift, a);
            }

            mag_zero(arb_radref(acb_realref(shift)));
            mag_zero(arb_radref(acb_imagref(shift)));

            iter++;
            total_iter++;

            acb_mat_approx_qr_step(A, Q, n0, n1, shift, prec);

            if (iter > maxiter)
            {
                result = 0;
                goto cleanup;
            }
        }
    }

cleanup:
    mag_clear(norm);
    mag_clear(tm);
    mag_clear(um);
    mag_clear(ts);
    mag_clear(eps);
    acb_clear(shift);
    acb_clear(s);
    acb_clear(t);
    acb_clear(a);
    acb_clear(b);

    return result;
}

int
acb_mat_approx_eig_qr(acb_ptr E, acb_mat_t L, acb_mat_t R, const acb_mat_t A, const mag_t tol, slong maxiter, slong prec)
{
    slong n, i, j;
    acb_ptr T;
    acb_mat_t Acopy, Q;
    int result;

    n = acb_mat_nrows(A);

    T = _acb_vec_init(n);
    acb_mat_init(Acopy, n, n);
    acb_mat_get_mid(Acopy, A);

    acb_mat_approx_hessenberg_reduce_0(Acopy, T, prec);

    if (L != NULL || R != NULL)
    {
        acb_mat_init(Q, n, n);
        acb_mat_set(Q, Acopy);
        acb_mat_approx_hessenberg_reduce_1(Q, T, prec);
    }

    for (i = 0; i < n; i++)
        for (j = i + 2; j < n; j++)
            acb_zero(acb_mat_entry(Acopy, j, i));

if (0)
{
    result = acb_mat_approx_hessenberg_qr(Acopy,
        (L != NULL || R != NULL) ? Q : NULL, tol, maxiter, prec);
}
else
{
    acb_mat_t Acopy2, Q2;
    acb_mat_init(Acopy2, A->r, A->c);

    flint_printf("n %wd   prec %wd\n", n, prec);

    TIMEIT_START
    acb_mat_set(Acopy2, Acopy);
    if (L != NULL || R != NULL)
    {
        acb_mat_init(Q2, Q->r, Q->c);
        acb_mat_set(Q2, Q);
    }
    result = acb_mat_approx_hessenberg_qr(Acopy2,
        (L != NULL || R != NULL) ? Q2 : NULL, tol, maxiter, prec);
    TIMEIT_STOP

    TIMEIT_START
    acb_mat_set(Acopy2, Acopy);
    if (L != NULL || R != NULL)
    {
        acb_mat_init(Q2, Q->r, Q->c);
        acb_mat_set(Q2, Q);
    }
    result = acb_mat_approx_hessenberg_qr2(Acopy2,
        (L != NULL || R != NULL) ? Q2 : NULL, tol, maxiter, prec);
    TIMEIT_STOP

    result = acb_mat_approx_hessenberg_qr(Acopy,
        (L != NULL || R != NULL) ? Q : NULL, tol, maxiter, prec);
}

    for (i = 0; i < n; i++)
        acb_get_mid(E + i, acb_mat_entry(Acopy, i, i));

    if (R != NULL)
    {
        acb_mat_t ER;
        acb_mat_init(ER, n, n);
        acb_mat_approx_eig_triu_r(ER, Acopy, prec);
        acb_mat_approx_mul(R, Q, ER, prec);
        acb_mat_get_mid(R, R);  /* zero radii */
        acb_mat_clear(ER);
    }

    if (L != NULL)
    {
        acb_mat_t EL;
        acb_mat_init(EL, n, n);
        acb_mat_approx_eig_triu_l(EL, Acopy, prec);
        acb_mat_conjugate_transpose(Q, Q);
        acb_mat_approx_mul(L, EL, Q, prec);
        acb_mat_get_mid(L, L);  /* zero radii */
        acb_mat_clear(EL);
    }

    if (L != NULL || R != NULL)
        acb_mat_clear(Q);

    _acb_vec_clear(T, n);
    acb_mat_clear(Acopy);

    return result;
}
