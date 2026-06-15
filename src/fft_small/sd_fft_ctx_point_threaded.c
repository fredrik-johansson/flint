/*
    Copyright (C) 2024 The FLINT authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"
#include "machine_vectors.h"
#include "thread_pool.h"
#include "fft_small.h"

/*
    The pointwise multiply/square iterate over independent BLK_SZ blocks
    (block I lives at sd_fft_ctx_blk_offset(I) = I << LG_BLK_SZ).  Splitting on
    the block index is exact (each block's arithmetic is independent) and needs
    only ONE fork-join with no internal barrier, since nothing consumes the
    result until the following inverse transform.  BLK_SZ = 256 doubles = 2 KB,
    so thread chunks never share a cache line.

    The _range bodies below are the serial sd_fft_ctx_point_mul / _point_sqr
    kernels (mpn_mul.c) with the outer block loop bounded to [Istart, Istop).
*/

void
sd_fft_ctx_point_mul_range(
    const sd_fft_ctx_t Q, double * a, const double * b, ulong m_, ulong depth,
    ulong Istart, ulong Istop)
{
    vec8d m    = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) m_, Q->p));
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);

    FLINT_ASSERT(depth >= LG_BLK_SZ);

    for (ulong I = Istart; I < Istop; I++)
    {
        double * ax = a + sd_fft_ctx_blk_offset(I);
        const double * bx = b + sd_fft_ctx_blk_offset(I);
        ulong jj = 0; do {
            vec8d x0, x1, b0, b1;
            x0 = vec8d_load(ax + jj + 0);
            x1 = vec8d_load(ax + jj + 8);
            b0 = vec8d_load(bx + jj + 0);
            b1 = vec8d_load(bx + jj + 8);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            x0 = vec8d_mulmod(x0, b0, n, ninv);
            x1 = vec8d_mulmod(x1, b1, n, ninv);
            vec8d_store(ax + jj + 0, x0);
            vec8d_store(ax + jj + 8, x1);
        } while (jj += 16, jj < BLK_SZ);
    }
}

void
sd_fft_ctx_point_sqr_range(
    const sd_fft_ctx_t Q, double * a, ulong m_, ulong depth,
    ulong Istart, ulong Istop)
{
    vec8d m    = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) m_, Q->p));
    vec8d n    = vec8d_set_d(Q->p);
    vec8d ninv = vec8d_set_d(Q->pinv);

    FLINT_ASSERT(depth >= LG_BLK_SZ);

    for (ulong I = Istart; I < Istop; I++)
    {
        double * ax = a + sd_fft_ctx_blk_offset(I);
        ulong jj = 0; do {
            vec8d x0, x1;
            x0 = vec8d_load(ax + jj + 0);
            x1 = vec8d_load(ax + jj + 8);
            x0 = vec8d_mulmod(x0, x0, n, ninv);
            x1 = vec8d_mulmod(x1, x1, n, ninv);
            x0 = vec8d_mulmod(x0, m, n, ninv);
            x1 = vec8d_mulmod(x1, m, n, ninv);
            vec8d_store(ax + jj + 0, x0);
            vec8d_store(ax + jj + 8, x1);
        } while (jj += 16, jj < BLK_SZ);
    }
}

typedef struct
{
    const sd_fft_ctx_struct * Q;
    double * a;
    const double * b;       /* NULL => square */
    ulong m_, depth, Istart, Istop;
}
sd_pw_par_t;

static void
sd_pw_par_worker(void * v)
{
    sd_pw_par_t * W = (sd_pw_par_t *) v;
    if (W->b == NULL)
        sd_fft_ctx_point_sqr_range(W->Q, W->a, W->m_, W->depth, W->Istart, W->Istop);
    else
        sd_fft_ctx_point_mul_range(W->Q, W->a, W->b, W->m_, W->depth, W->Istart, W->Istop);
}

static void
sd_pw_threaded(
    const sd_fft_ctx_t Q, double * a, const double * b, ulong m_, ulong depth,
    thread_pool_handle * handles, slong nhandles)
{
    ulong nblk = n_pow2(depth - LG_BLK_SZ);
    ulong nthreads = (ulong)(nhandles + 1);
    ulong t;
    sd_pw_par_t * W;

    if (nhandles < 1 || nblk < PW_PAR_MIN_BLK)
    {
        if (b == NULL) sd_fft_ctx_point_sqr(Q, a, m_, depth);
        else           sd_fft_ctx_point_mul(Q, a, b, m_, depth);
        return;
    }

    if (nthreads > nblk)
        nthreads = nblk;

    W = FLINT_ARRAY_ALLOC(nthreads, sd_pw_par_t);
    for (t = 0; t < nthreads; t++)
    {
        W[t].Q = Q; W[t].a = a; W[t].b = b; W[t].m_ = m_; W[t].depth = depth;
        W[t].Istart = (t * nblk) / nthreads;
        W[t].Istop  = ((t + 1) * nblk) / nthreads;
    }

    for (t = nthreads - 1; t >= 1; t--)
        thread_pool_wake(global_thread_pool, handles[t - 1], 0,
                         sd_pw_par_worker, W + t);
    sd_pw_par_worker(W + 0);
    for (t = nthreads - 1; t >= 1; t--)
        thread_pool_wait(global_thread_pool, handles[t - 1]);

    flint_free(W);
}

void
sd_fft_ctx_point_mul_threaded(
    const sd_fft_ctx_t Q, double * a, const double * b, ulong m_, ulong depth,
    thread_pool_handle * handles, slong nhandles)
{
    sd_pw_threaded(Q, a, b, m_, depth, handles, nhandles);
}

void
sd_fft_ctx_point_sqr_threaded(
    const sd_fft_ctx_t Q, double * a, ulong m_, ulong depth,
    thread_pool_handle * handles, slong nhandles)
{
    sd_pw_threaded(Q, a, NULL, m_, depth, handles, nhandles);
}
