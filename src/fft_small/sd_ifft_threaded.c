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
    Threaded sd_ifft_trunc.  The inverse is the mirror of the forward transform:
    rows first, then columns.

    Full case (== sd_ifft_no_trunc_internal):
        rows    : b in [0, l1)   contiguous, disjoint
        ---- barrier ----
        columns : a in [0, l2)   strided,    disjoint

    Truncated case (== sd_ifft_trunc_internal, k > 2): four ordered groups.
    The serial routine relies on the ORDER between groups (columns consume row
    outputs; the partial row is sandwiched between the column groups), so we
    parallelize only WITHIN each group -- whose iterations touch provably
    disjoint blocks -- and keep a barrier between groups, never reordering:

        (1) complete rows   b in [0, n1)
        (2) rightmost cols  a in [n2, z2p)
        (3) last partial row (single task)
        (4) leftmost cols   a in [0, n2)

    Group (3) runs between barriers, so the handles are idle at that instant and
    it is safe to recurse into the threaded driver for it.  This is NOT true for
    the multi-task groups: their workers run concurrently and must call the
    serial leaf routines, otherwise a worker would wake handles its siblings are
    already running on.
*/

typedef enum { IFFT_ROWS_FULL, IFFT_COLS_FULL, IFFT_COLS_TRUNC } sd_ifft_phase;

typedef struct
{
    const sd_fft_ctx_struct * Q;
    double * x;
    ulong j, k1, k2, l2;
    ulong z1, col_thr, col_n;   /* COLS_TRUNC: z = z1 + (a < col_thr), n = col_n */
    int   col_f;                /* COLS_TRUNC: f flag */
    ulong start, stop;
    sd_ifft_phase phase;
}
sd_ifft_par_t;

static void
sd_ifft_par_worker(void * v)
{
    sd_ifft_par_t * W = (sd_ifft_par_t *) v;
    const sd_fft_ctx_struct * Q = W->Q;
    double * x = W->x;

    switch (W->phase)
    {
    case IFFT_ROWS_FULL:
        for (ulong b = W->start; b < W->stop; b++)
            sd_ifft_no_trunc_internal(Q, x + BLK_SZ*(b << W->k2), W->k2,
                                      (W->j << W->k1) + b);
        break;
    case IFFT_COLS_FULL:
        for (ulong a = W->start; a < W->stop; a++)
            sd_ifft_no_trunc_block(Q, x + BLK_SZ*a, W->l2, W->k1, W->j);
        break;
    case IFFT_COLS_TRUNC:
        for (ulong a = W->start; a < W->stop; a++)
            sd_ifft_trunc_block(Q, x + BLK_SZ*a, W->l2, W->k1, W->j,
                                W->z1 + (a < W->col_thr), W->col_n, W->col_f);
        break;
    }
}

/* One fork-join over the true-index interval [lo, hi). */
static void
sd_ifft_par_run(const sd_ifft_par_t * base, sd_ifft_par_t * W,
                ulong lo, ulong hi, sd_ifft_phase phase,
                thread_pool_handle * handles, slong nhandles)
{
    ulong ntasks, nthreads, t;

    if (hi <= lo)
        return;

    ntasks = hi - lo;
    nthreads = (ulong)(nhandles + 1);
    if (nthreads > ntasks)
        nthreads = ntasks;

    for (t = 0; t < nthreads; t++)
    {
        W[t]       = *base;
        W[t].phase = phase;
        W[t].start = lo + (t * ntasks) / nthreads;
        W[t].stop  = lo + ((t + 1) * ntasks) / nthreads;
    }

    for (t = nthreads - 1; t >= 1; t--)
        thread_pool_wake(global_thread_pool, handles[t - 1], 0,
                         sd_ifft_par_worker, W + t);
    sd_ifft_par_worker(W + 0);
    for (t = nthreads - 1; t >= 1; t--)
        thread_pool_wait(global_thread_pool, handles[t - 1]);
}

static void
sd_ifft_trunc_internal_threaded(
    const sd_fft_ctx_t Q, double * x, ulong k, ulong j,
    ulong z, ulong n, int f,
    thread_pool_handle * handles, slong nhandles)
{
    ulong k1, k2, l2;
    sd_ifft_par_t base, * W;

    if (nhandles < 1 || k <= 2 || k < PAR_MIN_K)
    {
        sd_ifft_trunc_internal(Q, x, k, j, z, n, f);
        return;
    }

    k1 = k/2;
    k2 = k - k1;
    l2 = n_pow2(k2);

    W = FLINT_ARRAY_ALLOC(nhandles + 1, sd_ifft_par_t);
    base.Q = Q; base.x = x; base.j = j;
    base.k1 = k1; base.k2 = k2; base.l2 = l2;
    base.z1 = 0; base.col_thr = 0; base.col_n = 0; base.col_f = 0;

    if (!f && z == n && n == n_pow2(k))
    {
        /* full == sd_ifft_no_trunc_internal: rows, barrier, columns */
        sd_ifft_par_run(&base, W, 0, n_pow2(k1), IFFT_ROWS_FULL, handles, nhandles);
        sd_ifft_par_run(&base, W, 0, l2,         IFFT_COLS_FULL, handles, nhandles);
    }
    else
    {
        /* truncated == sd_ifft_trunc_internal (k > 2): four ordered groups */
        ulong n1  = n >> k2;
        ulong n2  = n & (l2 - 1);
        ulong z1  = z >> k2;
        ulong z2  = z & (l2 - 1);
        int   fp  = (n2 + f) > 0;
        ulong z2p = n_min(l2, z);
        ulong mm  = n_min(n2, z2);
        ulong mp  = n_max(n2, z2);

        base.z1 = z1;

        /* (1) complete rows */
        sd_ifft_par_run(&base, W, 0, n1, IFFT_ROWS_FULL, handles, nhandles);

        /* (2) rightmost columns: z = z1 + (a < mp), (n, f) = (n1, fp) */
        base.col_thr = mp; base.col_n = n1; base.col_f = fp;
        sd_ifft_par_run(&base, W, n2, z2p, IFFT_COLS_TRUNC, handles, nhandles);

        /* (3) last partial row: handles are idle here, so thread it recursively */
        if (fp)
            sd_ifft_trunc_internal_threaded(Q, x + BLK_SZ*(n1 << k2), k2,
                                            (j << k1) + n1, z2p, n2, f,
                                            handles, nhandles);

        /* (4) leftmost columns: z = z1 + (a < mm), (n, f) = (n1 + 1, 0) */
        base.col_thr = mm; base.col_n = n1 + 1; base.col_f = 0;
        sd_ifft_par_run(&base, W, 0, n2, IFFT_COLS_TRUNC, handles, nhandles);
    }

    flint_free(W);
}

void
sd_ifft_trunc_threaded(
    sd_fft_ctx_t Q, double * d, ulong L, ulong trunc,
    thread_pool_handle * handles, slong nhandles)
{
    FLINT_ASSERT(trunc <= n_pow2(L));

    if (nhandles < 1)
    {
        sd_ifft_trunc(Q, d, L, trunc);
        return;
    }

    sd_fft_ctx_fit_depth(Q, L);

    if (L > LG_BLK_SZ)
    {
        ulong new_trunc = n_cdiv(trunc, BLK_SZ);

        sd_ifft_trunc_internal_threaded(Q, d, L - LG_BLK_SZ, 0,
                                        new_trunc, new_trunc, 0, handles, nhandles);
        return;
    }

    sd_ifft_trunc(Q, d, L, trunc);
}
