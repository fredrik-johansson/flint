/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "gr_vec.h"
#include "gr_poly.h"

/* evaluate given precomputed powers up to x^m inclusive */
int
_gr_poly_evaluate_rectangular_precomp(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr xpow, slong m, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong i, r;
    gr_ptr s;
    int status = GR_SUCCESS;

    if (len <= 1)
    {
        if (len == 0)
            return gr_zero(res, ctx);
        else
            return gr_set(res, poly, ctx);
    }

    if (len <= m)
        return _gr_vec_dot(res, poly, 0, GR_ENTRY(poly, 1, sz), GR_ENTRY(xpow, 1, sz), m - 1, ctx);

    r = (len + m - 1) / m;

    GR_TMP_INIT(s, ctx);

    status |= _gr_vec_dot(res, GR_ENTRY(poly, (r - 1) * m, sz), 0, GR_ENTRY(xpow, 1, sz), GR_ENTRY(poly, (r - 1) * m + 1, sz), len - (r - 1) * m - 1, ctx);

    for (i = r - 2; i >= 0; i--)
    {
        status |= gr_mul(res, res, GR_ENTRY(xpow, m, sz), ctx);

        /* todo: can be better to do dot + add as a single operation, allowing
           a multiplication by 1 in the dot product */
        status |= _gr_vec_dot(s, GR_ENTRY(poly, i * m, sz), 0, GR_ENTRY(xpow, 1, sz), GR_ENTRY(poly, i * m + 1, sz), m - 1, ctx);
        status |= gr_add(res, res, s, ctx);
    }

    GR_TMP_CLEAR(s, ctx);
    return status;
}

typedef struct
{
    gr_srcptr f;
    gr_srcptr f_prime;
    gr_ptr w;
    gr_srcptr z;
    slong deg;
    gr_ctx_struct * ctx;
    slong m;
    slong a;
    slong b;
    int status;
}
work_chunk_t;

static void
aberth_worker(void * wwork)
{
    work_chunk_t * work = wwork;

    gr_ctx_struct * ctx = work->ctx;
    slong a = work->a;
    slong b = work->b;
    slong m = work->m;
    slong deg = work->deg;
    gr_srcptr z = work->z;
    gr_ptr w = work->w;
    gr_srcptr f = work->f;
    gr_srcptr f_prime = work->f_prime;

    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr u, T, U, P, Q;
    gr_srcptr zk, zj;
    gr_ptr wk;
    gr_ptr zkpow;
    slong j, k;

    GR_TMP_INIT5(u, T, U, P, Q, ctx);
    GR_TMP_INIT_VEC(zkpow, m + 1, ctx);

    for (k = a; k < b; k++)
    {
        zk = GR_ENTRY(z, k, sz);
        wk = GR_ENTRY(w, k, sz);

        status |= gr_zero(P, ctx);
        status |= gr_one(Q, ctx);

        for (j = 0; j < deg; j++)
        {
            if (j != k)
            {
                zj = GR_ENTRY(z, j, sz);
                status |= gr_sub(u, zk, zj, ctx);
                status |= gr_mul(P, P, u, ctx);
                status |= gr_add(P, P, Q, ctx);
                status |= gr_mul(Q, Q, u, ctx);
            }
        }

        status |= _gr_vec_set_powers(zkpow, zk, m + 1, ctx);
        status |= _gr_poly_evaluate_rectangular_precomp(T, f, deg + 1, zkpow, m, ctx);
        status |= _gr_poly_evaluate_rectangular_precomp(U, f_prime, deg, zkpow, m, ctx);

        status |= gr_mul(wk, T, Q, ctx);
        status |= gr_mul(u, Q, U, ctx);
        status |= gr_submul(u, T, P, ctx);
        status |= gr_div(wk, wk, u, ctx);
        if (status != GR_SUCCESS)
            status = gr_zero(wk, ctx);
    }

    GR_TMP_CLEAR5(u, T, U, P, Q, ctx);
    GR_TMP_CLEAR_VEC(zkpow, m + 1, ctx);

    work->status = status;
}

int
_gr_poly_refine_roots_aberth_threaded(gr_ptr w, gr_srcptr f, gr_srcptr f_prime, slong deg, gr_srcptr z, slong thread_limit, gr_ctx_t ctx)
{
    slong i, num_threads, num_workers;
    thread_pool_handle * handles;
    slong chunk_size;
    work_chunk_t * work;
    slong m;
    int status = GR_SUCCESS;
    TMP_INIT;
    TMP_START;

    m = n_sqrt(deg);
    m = FLINT_MAX(m, 1);

    num_workers = flint_request_threads(&handles, thread_limit);
    num_threads = num_workers + 1;

    work = TMP_ALLOC(num_threads * sizeof(work_chunk_t));

    chunk_size = (deg + num_threads - 1) / num_threads;

    for (i = 0; i < num_threads; i++)
    {
        work[i].f = f;
        work[i].f_prime = f_prime;
        work[i].w = w;
        work[i].z = z;
        work[i].deg = deg;
        work[i].ctx = ctx;
        work[i].m = m;
        work[i].a = i * chunk_size;
        work[i].b = FLINT_MIN((i + 1) * chunk_size, deg);
    }

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, aberth_worker, &work[i]);

    aberth_worker(&work[num_workers]);

    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    for (i = 0; i < num_workers; i++)
        status |= work[i].status;

    flint_give_back_threads(handles, num_workers);

    TMP_END;
    return status;
}

int
_gr_poly_refine_roots_aberth(gr_ptr w, gr_srcptr f, gr_srcptr f_prime, slong deg, gr_srcptr z, int progressive, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr u, T, U, P, Q;
    gr_srcptr zk, zj;
    gr_ptr wk;
    gr_ptr zkpow;
    gr_ptr zprog = NULL;
    slong j, k, m;

    if (progressive == 2)
    {
        return _gr_poly_refine_roots_aberth_threaded(w, f, f_prime, deg, z, flint_get_num_available_threads(), ctx);
    }

    m = n_sqrt(deg);
    m = FLINT_MAX(m, 1);

    GR_TMP_INIT5(u, T, U, P, Q, ctx);
    GR_TMP_INIT_VEC(zkpow, m + 1, ctx);

    if (progressive)
    {
        GR_TMP_INIT_VEC(zprog, deg, ctx);
        status |= _gr_vec_set(zprog, z, deg, ctx);
    }

    for (k = 0; k < deg; k++)
    {
        zk = progressive ? GR_ENTRY(zprog, k, sz) : GR_ENTRY(z, k, sz);
        wk = GR_ENTRY(w, k, sz);

        status |= gr_zero(P, ctx);
        status |= gr_one(Q, ctx);

        /* todo: vectorize subtractions (want _gr_vec_add_scalar) */
        /* todo: have some kind of _gr_vec_sum_reciprocals method? */
        for (j = 0; j < deg; j++)
        {
            if (j != k)
            {
                zj = progressive ? GR_ENTRY(zprog, j, sz) : GR_ENTRY(z, j, sz);

                status |= gr_sub(u, zk, zj, ctx);
                status |= gr_mul(P, P, u, ctx);
                status |= gr_add(P, P, Q, ctx);
                status |= gr_mul(Q, Q, u, ctx);
            }
        }

#if 0
        status |= _gr_poly_evaluate_horner(T, f, deg + 1, zk, ctx);
        status |= _gr_poly_evaluate_horner(U, f_prime, deg, zk, ctx);
#else
        status |= _gr_vec_set_powers(zkpow, zk, m + 1, ctx);
        status |= _gr_poly_evaluate_rectangular_precomp(T, f, deg + 1, zkpow, m, ctx);
        status |= _gr_poly_evaluate_rectangular_precomp(U, f_prime, deg, zkpow, m, ctx);
#endif

        status |= gr_mul(wk, T, Q, ctx);
        status |= gr_mul(u, Q, U, ctx);
        status |= gr_submul(u, T, P, ctx);
        status |= gr_div(wk, wk, u, ctx);
        if (status != GR_SUCCESS)
            status = gr_zero(wk, ctx);

        if (progressive)
            status |= gr_sub(GR_ENTRY(zprog, k, sz), GR_ENTRY(zprog, k, sz), wk, ctx);
    }

    GR_TMP_CLEAR5(u, T, U, P, Q, ctx);
    GR_TMP_CLEAR_VEC(zkpow, m + 1, ctx);

    if (progressive)
        GR_TMP_CLEAR_VEC(zprog, deg, ctx);

    return status;
}

int
_gr_poly_refine_roots_wdk(gr_ptr w, gr_srcptr f, slong deg, gr_srcptr z, int progressive, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_ptr u, T, Q;
    gr_srcptr zk, zj;
    gr_ptr wk;
    gr_ptr zkpow, zprog = NULL;
    slong j, k, m;
    int monic;

    m = n_sqrt(deg);
    m = FLINT_MAX(m, 1);

    GR_TMP_INIT3(u, T, Q, ctx);
    GR_TMP_INIT_VEC(zkpow, m + 1, ctx);

    if (progressive)
    {
        GR_TMP_INIT_VEC(zprog, deg, ctx);
        status |= _gr_vec_set(zprog, z, deg, ctx);
    }

    monic = (gr_is_one(GR_ENTRY(f, deg, sz), ctx) == T_TRUE);

    for (k = 0; k < deg; k++)
    {
        zk = progressive ? GR_ENTRY(zprog, k, sz) : GR_ENTRY(z, k, sz);
        wk = GR_ENTRY(w, k, sz);

        status |= gr_one(Q, ctx);

        /* todo: vectorize subtractions (want _gr_vec_add_scalar) */
        for (j = 0; j < deg; j++)
        {
            if (j != k)
            {
                zj = progressive ? GR_ENTRY(zprog, j, sz) : GR_ENTRY(z, j, sz);

                status |= gr_sub(u, zk, zj, ctx);
                status |= gr_mul(Q, Q, u, ctx);
            }
        }

#if 0
        status |= _gr_poly_evaluate_horner(T, f, deg + 1, zk, ctx);
#else
        status |= _gr_vec_set_powers(zkpow, zk, m + 1, ctx);
        status |= _gr_poly_evaluate_rectangular_precomp(T, f, deg + 1, zkpow, m, ctx);
#endif

        if (!monic)
            status |= gr_mul(Q, Q, GR_ENTRY(f, deg, sz), ctx);

        status |= gr_div(wk, T, Q, ctx);
        if (status != GR_SUCCESS)
            status = gr_zero(wk, ctx);

        if (progressive)
            status |= gr_sub(GR_ENTRY(zprog, k, sz), GR_ENTRY(zprog, k, sz), wk, ctx);
    }

    GR_TMP_CLEAR3(u, T, Q, ctx);
    GR_TMP_CLEAR_VEC(zkpow, m + 1, ctx);

    if (progressive)
        GR_TMP_CLEAR_VEC(zprog, deg, ctx);

    return status;
}
