/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NFLOAT_IMPL_H
#define NFLOAT_IMPL_H

#include "nfloat.h"

/* Whether a result with sign bit sgnbit must be rounded away from zero
   (as opposed to truncated) to respect the directed rounding mode. */
FLINT_FORCE_INLINE int
nfloat_should_round_up(int sgnbit, gr_ctx_t ctx)
{
    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        return 0;
    else
        return sgnbit ^ ((NFLOAT_CTX_FLAGS(ctx) & NFLOAT_RND_CEIL) != 0);
}

/* Checks the exponent of a nonzero finite result. */
#define NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx) \
    do { \
        if (FLINT_UNLIKELY((exp) < NFLOAT_MIN_EXP)) \
            return _nfloat_underflow((res), (sgnbit), (ctx)); \
        if (FLINT_UNLIKELY((exp) > NFLOAT_MAX_EXP)) \
            return _nfloat_overflow((res), (sgnbit), (ctx)); \
    } while (0)

/* Increments the normalized n-limb mantissa d of magnitude-truncated
   result, adjusting the exponent on carry out. */
#define NFLOAT_MANT_INCREMENT(d, n, exp) \
    do { \
        if (mpn_add_1((d), (d), (n), 1)) \
        { \
            (d)[(n) - 1] = UWORD(1) << (FLINT_BITS - 1); \
            (exp)++; \
        } \
    } while (0)

/* Copying and zeroing for short operands. Plain loops are turned into
   calls to memcpy / memset by the compiler, which costs more than the
   operation itself for a few limbs. */
FLINT_FORCE_INLINE void
_nfloat_copy_limbs(nn_ptr d, nn_srcptr s, slong n)
{
    switch (n)
    {
        case 8: d[7] = s[7]; /* fallthrough */
        case 7: d[6] = s[6]; /* fallthrough */
        case 6: d[5] = s[5]; /* fallthrough */
        case 5: d[4] = s[4]; /* fallthrough */
        case 4: d[3] = s[3]; /* fallthrough */
        case 3: d[2] = s[2]; /* fallthrough */
        case 2: d[1] = s[1]; /* fallthrough */
        case 1: d[0] = s[0]; /* fallthrough */
        case 0: break;
        default: flint_mpn_copyi(d, s, n);
    }
}

FLINT_FORCE_INLINE void
_nfloat_zero_limbs(nn_ptr d, slong n)
{
    switch (n)
    {
        case 8: d[7] = 0; /* fallthrough */
        case 7: d[6] = 0; /* fallthrough */
        case 6: d[5] = 0; /* fallthrough */
        case 5: d[4] = 0; /* fallthrough */
        case 4: d[3] = 0; /* fallthrough */
        case 3: d[2] = 0; /* fallthrough */
        case 2: d[1] = 0; /* fallthrough */
        case 1: d[0] = 0; /* fallthrough */
        case 0: break;
        default: flint_mpn_zero(d, n);
    }
}

#endif
