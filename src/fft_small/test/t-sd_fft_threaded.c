/*
    Copyright (C) 2024 The FLINT authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "ulong_extras.h"
#include "fft_small.h"
#include "thread_pool.h"
#include "fft_small.h"

/*
    The threaded transforms must be BIT-IDENTICAL to the serial ones for any
    handle count, because they run the same leaf operations on the same blocks.
    So the strongest, simplest check is memcmp of serial vs threaded output.

    Thread counts are chosen to exercise: the nh == 0 serial fallback (1), the
    (t*ntasks)/nthreads boundary arithmetic for counts that do not divide the
    task counts (3, 5), and the nthreads > ntasks clamp at small L (8).
*/

static slong
_grab_handles(slong want, thread_pool_handle ** handles_out)
{
    thread_pool_handle * h = NULL;
    slong nh = 0;

    flint_set_num_threads(want);
    if (want > 1)
    {
        h = FLINT_ARRAY_ALLOC(want - 1, thread_pool_handle);
        nh = thread_pool_request(global_thread_pool, h, want - 1);
    }
    *handles_out = h;
    return nh;
}

static void
_give_back_handles(thread_pool_handle * h, slong nh)
{
    slong i;
    for (i = 0; i < nh; i++)
        thread_pool_give_back(global_thread_pool, h[i]);
    flint_free(h);
}

TEST_FUNCTION_START(sd_fft_threaded, state)
{
    static const slong tcounts[] = { 1, 2, 3, 5, 8 };
    const int nT = sizeof(tcounts) / sizeof(tcounts[0]);

    sd_fft_ctx_t Q;
    sd_fft_ctx_init_prime(Q, UWORD(0x0003f00000000001));

    for (ulong L = LG_BLK_SZ + 1; L <= 18; L++)
    {
        ulong N = n_pow2(L);
        double * in   = FLINT_ARRAY_ALLOC(N, double);
        double * fref = flint_aligned_alloc(64, FLINT_MAX(64, N * sizeof(double)));
        double * rref = flint_aligned_alloc(64, FLINT_MAX(64, N * sizeof(double)));
        double * got  = flint_aligned_alloc(64, FLINT_MAX(64, N * sizeof(double)));

        for (ulong rep = 0; rep < 10; rep++)
        {
            ulong itrunc, otrunc, trunc;
            int ti;

            if (rep == 0)
            {
                itrunc = otrunc = trunc = N;        /* full transform */
            }
            else
            {
                itrunc = 1 + n_randint(state, N);
                otrunc = 1 + n_randint(state, N);
                trunc  = 1 + n_randint(state, N);
            }

            for (ulong i = 0; i < N; i++)
                in[i] = (double) n_randint(state, Q->mod.n);

            /* serial references (independent of thread count) */
            memcpy(fref, in, N * sizeof(double));
            sd_fft_trunc(Q, fref, L, itrunc, otrunc);

            memcpy(rref, in, N * sizeof(double));
            sd_fft_trunc(Q, rref, L, trunc, trunc);
            sd_ifft_trunc(Q, rref, L, trunc);

            for (ti = 0; ti < nT; ti++)
            {
                thread_pool_handle * h;
                slong nh = _grab_handles(tcounts[ti], &h);

                /* forward: threaded == serial bit-for-bit on [0, otrunc) */
                memcpy(got, in, N * sizeof(double));
                sd_fft_trunc_threaded(Q, got, L, itrunc, otrunc, h, nh);
                if (memcmp(fref, got, otrunc * sizeof(double)) != 0)
                {
                    flint_printf("FAIL forward: L=%wu itrunc=%wu otrunc=%wu "
                                 "nthreads=%wd nh=%wd\n",
                                 L, itrunc, otrunc, tcounts[ti], nh);
                    fflush(stdout);
                    flint_abort();
                }

                /* fft o ifft: threaded == serial bit-for-bit on [0, trunc) */
                memcpy(got, in, N * sizeof(double));
                sd_fft_trunc_threaded (Q, got, L, trunc, trunc, h, nh);
                sd_ifft_trunc_threaded(Q, got, L, trunc, h, nh);
                if (memcmp(rref, got, trunc * sizeof(double)) != 0)
                {
                    flint_printf("FAIL round-trip: L=%wu trunc=%wu "
                                 "nthreads=%wd nh=%wd\n",
                                 L, trunc, tcounts[ti], nh);
                    fflush(stdout);
                    flint_abort();
                }

                _give_back_handles(h, nh);
            }
        }

        flint_free(in);
        flint_aligned_free(fref);
        flint_aligned_free(rref);
        flint_aligned_free(got);
    }

    sd_fft_ctx_clear(Q);

    TEST_FUNCTION_END(state);
}
