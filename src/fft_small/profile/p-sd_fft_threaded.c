/*
    Copyright (C) 2024 The FLINT authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#if !defined(_POSIX_C_SOURCE) || _POSIX_C_SOURCE < 199309L
# define _POSIX_C_SOURCE 199309L   /* clock_gettime / CLOCK_MONOTONIC */
#endif

#include <math.h>
#include <string.h>
#include <time.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fft_small.h"
#include "thread_pool.h"
#include "fft_small.h"

/*
    Sweep transform depth L x thread count, timing each pipeline stage
    (forward x2, pointwise, inverse) separately and reporting speedup against
    the serial (nh == 0) path.  Use this to set PAR_MIN_K, PW_PAR_MIN_BLK, and
    the 2*np regime cutoff in mpn_mul.c:

      * PAR_MIN_K      : smallest L at which the fwd2/ifft columns clear ~1.0x
                         speedup at T = 2  (then k = L - LG_BLK_SZ).
      * PW_PAR_MIN_BLK : same crossover read off the pw column (in blocks,
                         2^(L-8)); pays off sooner than the FFT.
      * 2*np cutoff    : the largest np for which serialize+thread beats
                         concurrent-across-primes is where the total speedup
                         column exceeds np.

    On NUMA hardware, also run under `numactl --interleave=all`; the strided
    column phase makes the forward transform the NUMA-sensitive stage.
*/

static double
now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
}

int main(void)
{
    const ulong Lmin = 18, Lmax = 30;
    const slong threadset[] = { 1, 2, 4, 8, 16, 32 };
    const int nT = sizeof(threadset) / sizeof(threadset[0]);

    sd_fft_ctx_t Q;
    flint_rand_t state;

    sd_fft_ctx_init_prime(Q, UWORD(0x0003f00000000001));
    flint_rand_init(state);

    flint_printf("%4s %4s %10s %10s %10s %11s %8s\n",
                 "L", "T", "fwd2(ms)", "pw(ms)", "ifft(ms)", "total(ms)", "speedup");

    for (ulong L = Lmin; L <= Lmax; L++)
    {
        ulong N  = n_pow2(L);
        ulong m_ = n_randint(state, Q->mod.n);
        double serial_total = 0.0;
        double *a0, *b0, *a, *b;
        int ti;

        sd_fft_ctx_fit_depth(Q, L);     /* grow tables once, outside timing */

        a0 = flint_aligned_alloc(64, N * sizeof(double));
        b0 = flint_aligned_alloc(64, N * sizeof(double));
        a  = flint_aligned_alloc(64, N * sizeof(double));
        b  = flint_aligned_alloc(64, N * sizeof(double));

        for (ulong i = 0; i < N; i++)
        {
            a0[i] = (double) n_randint(state, Q->mod.n);
            b0[i] = (double) n_randint(state, Q->mod.n);
        }

        for (ti = 0; ti < nT; ti++)
        {
            slong T = threadset[ti];
            thread_pool_handle * h = NULL;
            slong nh = 0;
            int reps = (L <= 22) ? 20 : (L <= 26 ? 5 : 2);
            double tf = 1e30, tp = 1e30, tv = 1e30, total;
            int r;

            flint_set_num_threads(T);
            if (T > 1)
            {
                h  = FLINT_ARRAY_ALLOC(T - 1, thread_pool_handle);
                nh = thread_pool_request(global_thread_pool, h, T - 1);
            }

            for (r = 0; r < reps; r++)
            {
                double s;

                memcpy(a, a0, N * sizeof(double));
                memcpy(b, b0, N * sizeof(double));

                s = now_sec();                              /* forward x2 */
                sd_fft_trunc_threaded(Q, b, L, N, N, h, nh);
                sd_fft_trunc_threaded(Q, a, L, N, N, h, nh);
                tf = fmin(tf, now_sec() - s);

                s = now_sec();                              /* pointwise */
                sd_fft_ctx_point_mul_threaded(Q, a, b, m_, L, h, nh);
                tp = fmin(tp, now_sec() - s);

                s = now_sec();                              /* inverse */
                sd_ifft_trunc_threaded(Q, a, L, N, h, nh);
                tv = fmin(tv, now_sec() - s);
            }

            total = tf + tp + tv;
            if (T == 1)
                serial_total = total;

            flint_printf("%4wu %4wd %10.3f %10.3f %10.3f %11.3f %8.2f\n",
                         L, T, 1e3 * tf, 1e3 * tp, 1e3 * tv, 1e3 * total,
                         serial_total / total);

            if (nh > 0)
            {
                slong i;
                for (i = 0; i < nh; i++)
                    thread_pool_give_back(global_thread_pool, h[i]);
            }
            flint_free(h);
        }

        flint_aligned_free(a0);
        flint_aligned_free(b0);
        flint_aligned_free(a);
        flint_aligned_free(b);
        flint_printf("\n");
    }

    flint_rand_clear(state);
    sd_fft_ctx_clear(Q);
    return 0;
}
