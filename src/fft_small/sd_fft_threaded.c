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
    Threaded sd_fft_trunc.  The serial sd_fft_trunc_internal (k > 2) splits the
    2^k blocks into an l1 x l2 matrix (k1 = k/2, k2 = k - k1, l1 = 2^k1,
    l2 = 2^k2, row-major) and runs two phases:

        phase 0  columns : l2 strided block transforms      (a-indexed)
        ---- barrier ----
        phase 1  rows    : l1(+partial) contiguous transforms (b-indexed)

    Distinct columns touch disjoint strided blocks; distinct rows touch disjoint
    contiguous block ranges.  We therefore fan each phase out over its tasks and
    use the wait at the end of each phase as the barrier.  Sub-transforms call
    the unchanged serial leaf routines, so the result is bit-identical.
*/

#define FFT_PHASE_COLS 0
#define FFT_PHASE_ROWS 1

typedef struct
{
    const sd_fft_ctx_struct * Q;
    double * x;
    ulong j, k1, k2, l2;
    ulong z1, z2, n1p;      /* columns: itrunc = z1 + (a < z2), otrunc = n1p */
    ulong z2p, n1, n2;      /* rows:    itrunc = z2p, otrunc = (b < n1) ? l2 : n2 */
    ulong start, stop;      /* this worker's task slice */
    int phase;
    int full;               /* use the no_trunc leaf routines */
}
sd_fft_par_t;

static void
sd_fft_par_worker(void * v)
{
    sd_fft_par_t * W = (sd_fft_par_t *) v;
    const sd_fft_ctx_struct * Q = W->Q;
    double * x = W->x;

    if (W->phase == FFT_PHASE_COLS)
    {
        if (W->full)
            for (ulong a = W->start; a < W->stop; a++)
                sd_fft_no_trunc_block(Q, x + BLK_SZ*a, W->l2, W->k1, W->j);
        else
            for (ulong a = W->start; a < W->stop; a++)
                sd_fft_trunc_block(Q, x + BLK_SZ*a, W->l2, W->k1, W->j,
                                   W->z1 + (a < W->z2), W->n1p);
    }
    else
    {
        if (W->full)
            for (ulong b = W->start; b < W->stop; b++)
                sd_fft_no_trunc_internal(Q, x + BLK_SZ*(b << W->k2), W->k2,
                                         (W->j << W->k1) + b);
        else
            for (ulong b = W->start; b < W->stop; b++)
                sd_fft_trunc_internal(Q, x + BLK_SZ*(b << W->k2), W->k2,
                                      (W->j << W->k1) + b,
                                      W->z2p, (b < W->n1) ? W->l2 : W->n2);
    }
}

/* One fork-join over [0, ntasks): wake n-1 workers, run worker 0 inline, wait. */
static void
sd_fft_par_run(const sd_fft_par_t * base, sd_fft_par_t * W,
               ulong ntasks, int phase,
               thread_pool_handle * handles, slong nhandles)
{
    ulong nthreads, t;

    if (ntasks < 1)
        return;

    nthreads = (ulong)(nhandles + 1);
    if (nthreads > ntasks)
        nthreads = ntasks;          /* never wake idle workers */

    for (t = 0; t < nthreads; t++)
    {
        W[t]       = *base;
        W[t].phase = phase;
        W[t].start = (t * ntasks) / nthreads;
        W[t].stop  = ((t + 1) * ntasks) / nthreads;
    }

    for (t = nthreads - 1; t >= 1; t--)
        thread_pool_wake(global_thread_pool, handles[t - 1], 0,
                         sd_fft_par_worker, W + t);
    sd_fft_par_worker(W + 0);
    for (t = nthreads - 1; t >= 1; t--)
        thread_pool_wait(global_thread_pool, handles[t - 1]);
}

static void
sd_fft_trunc_internal_threaded(
    const sd_fft_ctx_t Q, double * x, ulong k, ulong j,
    ulong itrunc, ulong otrunc,
    thread_pool_handle * handles, slong nhandles)
{
    ulong k1, k2, l2;
    sd_fft_par_t base, * W;

    /* Serial fallback for: no budget, too small to amortize two fork-joins,
       leaf sizes, or the degenerate itrunc/otrunc cases that the serial
       routine handles cheaply (zero-fill / no-op). */
    if (nhandles < 1 || k <= 2 || k < PAR_MIN_K || itrunc < 1 || otrunc < 1)
    {
        sd_fft_trunc_internal(Q, x, k, j, itrunc, otrunc);
        return;
    }

    k1 = k/2;
    k2 = k - k1;
    l2 = n_pow2(k2);

    W = FLINT_ARRAY_ALLOC(nhandles + 1, sd_fft_par_t);

    base.Q = Q; base.x = x; base.j = j;
    base.k1 = k1; base.k2 = k2; base.l2 = l2;
    base.z1 = 0; base.z2 = 0; base.n1p = 0;
    base.z2p = 0; base.n1 = 0; base.n2 = 0;

    if (itrunc == otrunc && otrunc == n_pow2(k))
    {
        /* full == sd_fft_no_trunc_internal: columns, barrier, rows */
        base.full = 1;
        sd_fft_par_run(&base, W, l2,         FFT_PHASE_COLS, handles, nhandles);
        sd_fft_par_run(&base, W, n_pow2(k1), FFT_PHASE_ROWS, handles, nhandles);
    }
    else
    {
        /* truncated == sd_fft_trunc_internal (k > 2): columns, barrier, rows */
        ulong n1  = otrunc >> k2;
        ulong n2  = otrunc & (l2 - 1);
        ulong z1  = itrunc >> k2;
        ulong z2  = itrunc & (l2 - 1);
        ulong n1p = n1 + (n2 != 0);
        ulong z2p = n_min(l2, itrunc);

        base.full = 0;
        base.z1 = z1; base.z2 = z2; base.n1p = n1p;
        base.z2p = z2p; base.n1 = n1; base.n2 = n2;

        sd_fft_par_run(&base, W, z2p, FFT_PHASE_COLS, handles, nhandles);
        sd_fft_par_run(&base, W, n1p, FFT_PHASE_ROWS, handles, nhandles);
    }

    flint_free(W);
}

void
sd_fft_trunc_threaded(
    sd_fft_ctx_t Q, double * d, ulong L, ulong itrunc, ulong otrunc,
    thread_pool_handle * handles, slong nhandles)
{
    FLINT_ASSERT(itrunc <= n_pow2(L));
    FLINT_ASSERT(otrunc <= n_pow2(L));

    if (nhandles < 1)
    {
        sd_fft_trunc(Q, d, L, itrunc, otrunc);
        return;
    }

    sd_fft_ctx_fit_depth(Q, L);     /* the only Q mutation; before any worker */

    if (L > LG_BLK_SZ)
    {
        ulong new_itrunc = n_cdiv(itrunc, BLK_SZ);
        ulong new_otrunc = n_cdiv(otrunc, BLK_SZ);

        for (int i = 0; i < ((-(int) itrunc) & (BLK_SZ - 1)); i++)
            d[itrunc + i] = 0.0;

        sd_fft_trunc_internal_threaded(Q, d, L - LG_BLK_SZ, 0,
                                       new_itrunc, new_otrunc, handles, nhandles);
        return;
    }

    /* L <= LG_BLK_SZ: too small to thread; serial basecases. */
    sd_fft_trunc(Q, d, L, itrunc, otrunc);
}
