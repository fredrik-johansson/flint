/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "impl.h"

/*
    Powers |x|^y = exp(t), t = y log |x|, of floating-point mantissas
    (mp_real.h): with |t| < 2^tb, the logarithm is computed to a relative
    accuracy of about 2^-(pb + tb + 12) and multiplied exactly by the
    mantissa of y, so that t is known to an absolute accuracy below
    2^-(pb + 10), and exp(t) is evaluated by the exponential
    kernels (the huge one beyond 2^(FLINT_BITS - 3)).  The callers handle
    the exact cases (|x| = 1, y = 0, x a power of two with an integer y,
    y = +-1/2, ...) and the special values.
*/

int
_mp_real_pow_fp(nn_ptr y, ulong * err, nn_ptr q, slong * qn, int * tsgn,
    nn_srcptr xd, slong xn, slong xe, nn_srcptr yd, slong yn, slong ye, int ysgn, slong pb, slong tbmax)
{
    nn_ptr L, P, Tm;
    slong w, tb, len, ue, pl, te, texp, sh, dd, extra, nk = _mp_real_elem_limbs(pb);
    ulong errL, et, eprop;
    int lneg, sgn;
    TMP_INIT;

    /* |t| < 2^tb with tb = e_y + bits(|E| + 1) (|log |x|| < |E| + 1);
       log |x| to tb + 12 more bits of relative accuracy */
    {
        ulong aE = (xe < 0) ? -(ulong) xe : (ulong) xe;

        tb = FLINT_MAX(ye, -FLINT_BITS * nk - 8) + FLINT_BIT_COUNT(aE + 1);
        if (tb > tbmax + 2)
            tb = tbmax + 2;
        tb = FLINT_MAX(tb, 0);
    }

    /* t = y L with L the logarithm at lb = pb + tb + 12 bits (the
       kernel adds the bits lost to cancellation near 1), truncated to
       w >= nk limbs, tb + 12 more bits than the result: the propagated
       error below takes a right shift */
    w = nk + (tb + 12 + FLINT_BITS - 1) / FLINT_BITS;

    TMP_START;

    for (extra = 0; ; extra += FLINT_BITS)
    {
        L = TMP_ALLOC(sizeof(ulong) * (_mp_real_elem_limbs(pb + tb + 12 + extra) + xn + 5));

        /* log |x| = (L, len) 2^ue within errL units */
        lneg = _mp_real_log_fp(L, &len, &ue, &errL, xd, xn, xe, 0, pb + tb + 12 + extra);
        if (len == 0)
        {
            TMP_END;
            return MP_REAL_POW_ONE;
        }

        P = TMP_ALLOC(sizeof(ulong) * (len + yn + w + 1));
        Tm = P + len + yn;

        /* t = y L exactly from the mantissas: P 2^(ue + e_y - FLINT_BITS yn),
           within errL |Y| 2^(ue + e_y - FLINT_BITS yn) < errL 2^(ue + e_y) */
        if (len >= yn)
            flint_mpn_mul(P, L, len, yd, yn);
        else
            flint_mpn_mul(P, yd, yn, L, len);
        pl = len + yn;
        while (P[pl - 1] == 0)
            pl--;

        /* the top w limbs, normalised: t = Tm 2^te within the units for
           the dropped bits plus the error of L, errL 2^(ue + e_y - te) */
        et = _mp_real_fp_normalize(Tm, &te, P, pl, ue + ye - FLINT_BITS * yn, 0, w);

        texp = te + FLINT_BITS * w;
        *tsgn = lneg ^ ysgn;

        /* tiny t: 1 + delta, |delta| < 2|t| */
        if (texp < -FLINT_BITS * nk - 4)
        {
            TMP_END;
            return MP_REAL_POW_TINY;
        }

        if (texp > tbmax)
        {
            TMP_END;
            return MP_REAL_POW_HUGE;
        }

        /* exp(t +- delta) = exp(t)(1 +- 2 delta) and exp(t) < 2^(FLINT_BITS nk + 1)
           units: delta 2^(FLINT_BITS nk + 2) units, for the truncation
           delta = et 2^te and the error of L, errL 2^(ue + e_y) */
        sh = te + FLINT_BITS * nk + 2;
        dd = ue + ye + FLINT_BITS * nk + 2;
        eprop = 0;

        FLINT_ASSERT(sh <= 0);
        if (sh <= -FLINT_BITS)
            eprop += 1;
        else
            eprop += (et >> (-sh)) + 1;

        if (dd <= -FLINT_BITS)
            eprop += 1;
        else if (dd <= 0)
            eprop += (errL >> (-dd)) + 1;
        else if (dd + (slong) FLINT_BIT_COUNT(errL) <= FLINT_BITS - 4)
            eprop += errL << dd;     /* (below 2^(FLINT_BITS nk - pb + 12): the
                                       output units are finer than pb bits) */
        else
            continue;   /* (rare: the logarithm's units too coarse for the count) */

        break;
    }

    if (texp <= MP_REAL_EXP_FP_MAX_EXP)
    {
        slong k = _mp_real_exp_fp(y, err, Tm, w, texp, *tsgn, pb);
        sgn = (k < 0);
        q[0] = (k < 0) ? -(ulong) k : (ulong) k;
        *qn = (k != 0);
    }
    else
        sgn = _mp_real_exp_fp_huge(y, err, q, qn, Tm, w, texp, *tsgn, pb);

    *err += eprop;
    TMP_END;
    return sgn;
}
