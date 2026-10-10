/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"
#include "nfloat.h"

/*
    Complex elementary functions, as compositions of the real functions
    using formulas without catastrophic cancellation except where the
    function itself is ill-conditioned (no internal guard bits; each
    component typically has an error of a few ulp relative to the
    magnitude of the result, but small components of large results can
    have large relative errors, and directed rounding is not respected).
    Real arguments in the real domain of the function use the real version.
*/

#define RE(z) NFLOAT_COMPLEX_RE(z, ctx)
#define IM(z) NFLOAT_COMPLEX_IM(z, ctx)
#define CTMP(T) ulong T[2 * NFLOAT_MAX_ALLOC]
#define RTMP(T) ulong T[NFLOAT_MAX_ALLOC]

/* (the operations leave garbage on over- or underflow) */
#define CHK(op) do { int __st = (op); if (__st != GR_SUCCESS) return __st; } while (0)

typedef enum
{
    C_EXP, C_EXPM1, C_EXP_PI_I, C_EXP2, C_EXP10,
    C_LOG, C_LOG1P, C_LOG_PI_I, C_LOG2, C_LOG10,
    C_SIN, C_COS, C_TAN, C_COT, C_SEC, C_CSC,
    C_SIN_PI, C_COS_PI, C_TAN_PI, C_COT_PI, C_SEC_PI, C_CSC_PI,
    C_SINC, C_SINC_PI,
    C_SINH, C_COSH, C_TANH, C_COTH, C_SECH, C_CSCH,
    C_ASIN, C_ACOS, C_ATAN, C_ACOT, C_ASEC, C_ACSC,
    C_ASINH, C_ACOSH, C_ATANH, C_ACOTH, C_ASECH, C_ACSCH,
    C_ASIN_PI, C_ACOS_PI, C_ATAN_PI, C_ACOT_PI, C_ASEC_PI, C_ACSC_PI
}
cfunc_t;

/* real domains */
enum { D_NONE, D_ALL, D_POS, D_GT_M1, D_ABS_LE1, D_ABS_LT1, D_ABS_GE1, D_ABS_GT1, D_GE1, D_ASECH };

typedef int (*_real_func_t)(nfloat_ptr, nfloat_srcptr, gr_ctx_t);

static const struct
{
    _real_func_t real;
    unsigned char domain;
}
_cfunc_real[] =
{
    { nfloat_exp, D_ALL }, { nfloat_expm1, D_ALL }, { NULL, D_NONE }, { nfloat_exp2, D_ALL }, { nfloat_exp10, D_ALL },
    { nfloat_log, D_POS }, { nfloat_log1p, D_GT_M1 }, { NULL, D_NONE }, { nfloat_log2, D_POS }, { nfloat_log10, D_POS },
    { nfloat_sin, D_ALL }, { nfloat_cos, D_ALL }, { nfloat_tan, D_ALL }, { nfloat_cot, D_ALL }, { nfloat_sec, D_ALL }, { nfloat_csc, D_ALL },
    { nfloat_sin_pi, D_ALL }, { nfloat_cos_pi, D_ALL }, { nfloat_tan_pi, D_ALL }, { nfloat_cot_pi, D_ALL }, { nfloat_sec_pi, D_ALL }, { nfloat_csc_pi, D_ALL },
    { nfloat_sinc, D_ALL }, { nfloat_sinc_pi, D_ALL },
    { nfloat_sinh, D_ALL }, { nfloat_cosh, D_ALL }, { nfloat_tanh, D_ALL }, { nfloat_coth, D_ALL }, { nfloat_sech, D_ALL }, { nfloat_csch, D_ALL },
    { nfloat_asin, D_ABS_LE1 }, { nfloat_acos, D_ABS_LE1 }, { nfloat_atan, D_ALL }, { nfloat_acot, D_ALL }, { nfloat_asec, D_ABS_GE1 }, { nfloat_acsc, D_ABS_GE1 },
    { nfloat_asinh, D_ALL }, { nfloat_acosh, D_GE1 }, { nfloat_atanh, D_ABS_LT1 }, { nfloat_acoth, D_ABS_GT1 }, { nfloat_asech, D_ASECH }, { nfloat_acsch, D_ALL },
    { nfloat_asin_pi, D_ABS_LE1 }, { nfloat_acos_pi, D_ABS_LE1 }, { nfloat_atan_pi, D_ALL }, { nfloat_acot_pi, D_ALL }, { nfloat_asec_pi, D_ABS_GE1 }, { nfloat_acsc_pi, D_ABS_GE1 },
};

/* whether the real number a lies in the real domain d */
static int
_in_domain(nfloat_srcptr a, int d, gr_ctx_t ctx)
{
    int c, s;
    RTMP(one);

    if (d == D_ALL)
        return 1;
    if (d == D_NONE)
        return 0;

    s = NFLOAT_IS_ZERO(a) ? 0 : (NFLOAT_SGNBIT(a) ? -1 : 1);

    if (d == D_POS)
        return s > 0;

    nfloat_one(one, ctx);
    if (d == D_GT_M1)
    {
        if (s >= 0)
            return 1;
        nfloat_cmpabs(&c, a, one, ctx);
        return c < 0;
    }

    nfloat_cmpabs(&c, a, one, ctx);

    switch (d)
    {
        case D_ABS_LE1: return c <= 0;
        case D_ABS_LT1: return c < 0;
        case D_ABS_GE1: return c >= 0;
        case D_ABS_GT1: return c > 0;
        case D_GE1: return s > 0 && c >= 0;
        default: return s > 0 && c <= 0;
    }
}

/* res = i^k z */
static void
_c_mul_i_pow(nfloat_complex_ptr res, nfloat_complex_srcptr z, int k, gr_ctx_t ctx)
{
    CTMP(t);

    k &= 3;
    if (k == 0)
    {
        nfloat_complex_set(res, (nfloat_complex_ptr) z, ctx);
        return;
    }
    if (k == 2)
    {
        nfloat_complex_neg(res, z, ctx);
        return;
    }

    /* i (a + bi) = -b + ai, -i (a + bi) = b - ai */
    nfloat_set(RE(t), IM(z), ctx);
    nfloat_set(IM(t), RE(z), ctx);
    if (k == 1)
        nfloat_neg(RE(t), RE(t), ctx);
    else
        nfloat_neg(IM(t), IM(t), ctx);
    nfloat_complex_set(res, t, ctx);
}

/* res = x c (c real) */
static int
_c_mul_real(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_srcptr c, gr_ctx_t ctx)
{
    CHK(nfloat_mul(RE(res), RE(x), c, ctx));
    return nfloat_mul(IM(res), IM(x), c, ctx);
}

static int
_c_div_real(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_srcptr c, gr_ctx_t ctx)
{
    CHK(nfloat_div(RE(res), RE(x), c, ctx));
    return nfloat_div(IM(res), IM(x), c, ctx);
}

/* sin and cos of a + bi (of pi (a + bi) if pi); s or c may be NULL */
static int
_c_sin_cos(nfloat_complex_ptr s, nfloat_complex_ptr c, nfloat_srcptr a, nfloat_srcptr b, int pi, gr_ctx_t ctx)
{
    RTMP(sa); RTMP(ca); RTMP(sb); RTMP(cb); RTMP(t);

    /* sin(a + bi) = sin a cosh b + i cos a sinh b,
       cos(a + bi) = cos a cosh b - i sin a sinh b */
    if (pi)
    {
        CHK(nfloat_sin_cos_pi(sa, ca, a, ctx));
        CHK(nfloat_pi(t, ctx));
        CHK(nfloat_mul(t, t, b, ctx));
        CHK(nfloat_sinh_cosh(sb, cb, t, ctx));
    }
    else
    {
        CHK(nfloat_sin_cos(sa, ca, a, ctx));
        CHK(nfloat_sinh_cosh(sb, cb, b, ctx));
    }

    if (s != NULL)
    {
        CHK(nfloat_mul(RE(s), sa, cb, ctx));
        CHK(nfloat_mul(IM(s), ca, sb, ctx));
    }

    if (c != NULL)
    {
        CHK(nfloat_mul(RE(c), ca, cb, ctx));
        CHK(nfloat_mul(IM(c), sa, sb, ctx));
        nfloat_neg(IM(c), IM(c), ctx);
    }

    return GR_SUCCESS;
}

/* tan(a + bi) (of pi (a + bi) if pi), b != 0 */
static int
_c_tan(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, int pi, gr_ctx_t ctx)
{
    RTMP(sa); RTMP(ca); RTMP(T); RTMP(S); RTMP(D); RTMP(u);

    /* tan(a + bi) = (sin a cos a + i sinh b cosh b) / (cos^2 a + sinh^2 b);
       dividing by cosh^2 b, with T = tanh b, S = sech b:
       (sin a cos a S^2 + i T) / (cos^2 a S^2 + T^2), a sum of
       nonnegative terms in the denominator */
    if (pi)
    {
        CHK(nfloat_pi(u, ctx));
        CHK(nfloat_mul(u, u, b, ctx));
        b = u;
    }

    CHK(nfloat_tanh(T, b, ctx));

    /* tan(bi) = i tanh(b) */
    if (NFLOAT_IS_ZERO(a))
    {
        nfloat_zero(RE(res), ctx);
        return nfloat_set(IM(res), T, ctx);
    }

    if (pi)
        CHK(nfloat_sin_cos_pi(sa, ca, a, ctx));
    else
        CHK(nfloat_sin_cos(sa, ca, a, ctx));

    /* sech b underflows only for huge |b|, where its square does not
       affect the denominator (the real part underflows) */
    if (nfloat_sech(S, b, ctx) != GR_SUCCESS)
    {
        if (NFLOAT_EXP(b) < 8)
            return GR_UNABLE;
        nfloat_zero(RE(res), ctx);
        return nfloat_inv(IM(res), T, ctx);
    }

    CHK(nfloat_sqr(S, S, ctx));
    CHK(nfloat_mul(sa, sa, S, ctx));
    CHK(nfloat_mul(sa, sa, ca, ctx));
    CHK(nfloat_sqr(ca, ca, ctx));
    CHK(nfloat_mul(ca, ca, S, ctx));
    CHK(nfloat_sqr(u, T, ctx));
    CHK(nfloat_add(D, ca, u, ctx));
    CHK(nfloat_div(RE(res), sa, D, ctx));
    return nfloat_div(IM(res), T, D, ctx);
}

/* exp(a + bi) */
static int
_c_exp(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    RTMP(E); RTMP(s); RTMP(c);

    CHK(nfloat_sin_cos(s, c, b, ctx));
    if (NFLOAT_IS_ZERO(a))
    {
        nfloat_set(RE(res), c, ctx);
        nfloat_set(IM(res), s, ctx);
        return GR_SUCCESS;
    }
    CHK(nfloat_exp(E, a, ctx));
    CHK(nfloat_mul(RE(res), E, c, ctx));
    return nfloat_mul(IM(res), E, s, ctx);
}

/* log(a + bi), not both zero */
static int
_c_log(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    RTMP(t); RTMP(u); RTMP(one);
    nfloat_srcptr p, q;
    slong e;

    CHK(nfloat_atan2(IM(res), b, a, ctx));

    /* p = the larger of |a|, |b| */
    p = a; q = b;
    if (NFLOAT_IS_ZERO(a) || (!NFLOAT_IS_ZERO(b) && NFLOAT_EXP(b) > NFLOAT_EXP(a)))
    {
        p = b; q = a;
    }
    e = NFLOAT_EXP(p);

    if (e == 0 || e == 1)
    {
        /* |z| close to 1: log(|z|) = log1p((p - 1)(p + 1) + q^2) / 2, where
           p - 1 is exact for 1/2 <= p <= 2 */
        nfloat_one(one, ctx);
        nfloat_abs(t, p, ctx);
        CHK(nfloat_add(u, t, one, ctx));
        CHK(nfloat_sub(t, t, one, ctx));
        CHK(nfloat_mul(t, t, u, ctx));
        if (!NFLOAT_IS_ZERO(q))
        {
            CHK(nfloat_sqr(u, q, ctx));
            CHK(nfloat_add(t, t, u, ctx));
        }
        if (NFLOAT_IS_ZERO(t))
            return nfloat_zero(RE(res), ctx);
        CHK(nfloat_log1p(RE(res), t, ctx));
        NFLOAT_EXP(RE(res)) -= 1;
        return GR_SUCCESS;
    }

    CHK(nfloat_hypot(t, a, b, ctx));
    return nfloat_log(RE(res), t, ctx);
}

/* atan(a + bi), b != 0, a + bi != +-i */
static int
_c_atan(nfloat_complex_ptr res, nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    RTMP(t); RTMP(u); RTMP(v); RTMP(one);

    /* atan is odd: b > 0 avoids cancellation in 1 + q below */
    if (NFLOAT_SGNBIT(b))
    {
        RTMP(na); RTMP(nb);
        nfloat_neg(na, a, ctx);
        nfloat_neg(nb, b, ctx);
        CHK(_c_atan(res, na, nb, ctx));
        return nfloat_complex_neg(res, res, ctx);
    }

    /* Re = atan2(2a, (1 - b)(1 + b) - a^2) / 2,
       Im = log1p(4b / (a^2 + (1 - b)^2)) / 4, the argument of log1p
       positive */
    nfloat_one(one, ctx);
    CHK(nfloat_sub(t, one, b, ctx));
    CHK(nfloat_add(u, one, b, ctx));
    CHK(nfloat_mul(u, t, u, ctx));
    CHK(nfloat_sqr(t, t, ctx));
    if (!NFLOAT_IS_ZERO(a))
    {
        CHK(nfloat_sqr(v, a, ctx));
        CHK(nfloat_sub(u, u, v, ctx));
        CHK(nfloat_add(t, t, v, ctx));
    }

    CHK(nfloat_div(t, b, t, ctx));
    NFLOAT_EXP(t) += 2;
    CHK(nfloat_log1p(IM(res), t, ctx));
    NFLOAT_EXP(IM(res)) -= 2;

    if (NFLOAT_IS_ZERO(a))
    {
        /* atan(bi) for |b| > 1 has real part sgn(b) pi/2 */
        if (NFLOAT_SGNBIT(u))
        {
            CHK(nfloat_pi(RE(res), ctx));
            NFLOAT_EXP(RE(res)) -= 1;
            NFLOAT_SGNBIT(RE(res)) = NFLOAT_SGNBIT(b);
        }
        else
            nfloat_zero(RE(res), ctx);
        return GR_SUCCESS;
    }

    nfloat_mul_2exp_si(v, a, 1, ctx);
    CHK(nfloat_atan2(RE(res), v, u, ctx));
    NFLOAT_EXP(RE(res)) -= 1;
    return GR_SUCCESS;
}

/* square roots of 1 - z and 1 + z (or z - 1 and z + 1 if minus) */
static int
_c_sqrt_1pm(nfloat_complex_ptr s1, nfloat_complex_ptr s2, nfloat_complex_srcptr z, int minus, gr_ctx_t ctx)
{
    CTMP(t);

    nfloat_complex_one(t, ctx);
    if (minus)
        CHK(nfloat_complex_sub(s1, z, t, ctx));
    else
        CHK(nfloat_complex_sub(s1, t, z, ctx));
    CHK(nfloat_complex_add(s2, t, z, ctx));
    CHK(nfloat_complex_sqrt(s1, s1, ctx));
    return nfloat_complex_sqrt(s2, s2, ctx);
}

/* Kahan's formulas for asin, acos and acosh */
static int
_c_asin_acos(nfloat_complex_ptr res, nfloat_complex_srcptr z, cfunc_t f, gr_ctx_t ctx)
{
    CTMP(s1); CTMP(s2);
    RTMP(t); RTMP(u);

    if (f == C_ACOSH)
    {
        /* s1 = sqrt(z - 1), s2 = sqrt(z + 1):
           Re = asinh(Re(conj(s1) s2)), Im = 2 atan2(Im s1, Re s2) */
        CHK(_c_sqrt_1pm(s1, s2, z, 1, ctx));
        CHK(nfloat_mul(t, RE(s1), RE(s2), ctx));
        CHK(nfloat_mul(u, IM(s1), IM(s2), ctx));
        CHK(nfloat_add(t, t, u, ctx));
        CHK(nfloat_atan2(IM(res), IM(s1), RE(s2), ctx));
        NFLOAT_EXP(IM(res)) += !NFLOAT_IS_ZERO(IM(res));
        return nfloat_asinh(RE(res), t, ctx);
    }

    /* s1 = sqrt(1 - z), s2 = sqrt(1 + z) */
    CHK(_c_sqrt_1pm(s1, s2, z, 0, ctx));

    if (f == C_ASIN)
    {
        /* Re = atan2(a, Re(s1 s2)), Im = asinh(Im(conj(s1) s2)) */
        CHK(nfloat_mul(t, RE(s1), RE(s2), ctx));
        CHK(nfloat_mul(u, IM(s1), IM(s2), ctx));
        CHK(nfloat_sub(t, t, u, ctx));
        CHK(nfloat_mul(u, RE(s1), IM(s2), ctx));
        CHK(nfloat_submul(u, IM(s1), RE(s2), ctx));
        CHK(nfloat_atan2(RE(res), RE(z), t, ctx));
        return nfloat_asinh(IM(res), u, ctx);
    }
    else
    {
        /* Re = 2 atan2(Re s1, Re s2), Im = asinh(Im(conj(s2) s1)) */
        CHK(nfloat_mul(u, RE(s2), IM(s1), ctx));
        CHK(nfloat_submul(u, IM(s2), RE(s1), ctx));
        CHK(nfloat_atan2(RE(res), RE(s1), RE(s2), ctx));
        NFLOAT_EXP(RE(res)) += !NFLOAT_IS_ZERO(RE(res));
        return nfloat_asinh(IM(res), u, ctx);
    }
}

static int
_nfloat_complex_func(nfloat_complex_ptr res, nfloat_complex_srcptr x, cfunc_t f, gr_ctx_t ctx)
{
    CTMP(z); CTMP(w);
    RTMP(t);
    nfloat_srcptr a = RE(x), b = IM(x);
    int recip_arg = 0, recip_res = 0, rot = 0, div_pi = 0;

    if (NFLOAT_COMPLEX_IS_SPECIAL(x, ctx) && !(NFLOAT_IS_ZERO(a) || NFLOAT_IS_ZERO(b)))
        return GR_UNABLE;
    if (NFLOAT_IS_SPECIAL(a) && !NFLOAT_IS_ZERO(a))
        return GR_UNABLE;
    if (NFLOAT_IS_SPECIAL(b) && !NFLOAT_IS_ZERO(b))
        return GR_UNABLE;

    /* real arguments in the real domain */
    if (NFLOAT_IS_ZERO(b) && _in_domain(a, _cfunc_real[f].domain, ctx))
    {
        CHK(_cfunc_real[f].real(RE(res), a, ctx));
        return nfloat_zero(IM(res), ctx);
    }

    /* reduce to fewer cases: functions of 1/z, reciprocals of functions,
       hyperbolic functions as trigonometric functions of iz */
    if (f >= C_ASIN_PI)
    {
        div_pi = 1;
        f = (cfunc_t) (f - C_ASIN_PI + C_ASIN);
    }

    switch (f)
    {
        case C_ACOT: f = C_ATAN; recip_arg = 1; break;
        case C_ASEC: f = C_ACOS; recip_arg = 1; break;
        case C_ACSC: f = C_ASIN; recip_arg = 1; break;
        case C_ACOTH: f = C_ATANH; recip_arg = 1; break;
        case C_ASECH: f = C_ACOSH; recip_arg = 1; break;
        case C_ACSCH: f = C_ASINH; recip_arg = 1; break;
        case C_COT: f = C_TAN; recip_res = 1; break;
        case C_SEC: f = C_COS; recip_res = 1; break;
        case C_CSC: f = C_SIN; recip_res = 1; break;
        case C_COT_PI: f = C_TAN_PI; recip_res = 1; break;
        case C_SEC_PI: f = C_COS_PI; recip_res = 1; break;
        case C_CSC_PI: f = C_SIN_PI; recip_res = 1; break;
        case C_COTH: f = C_TANH; recip_res = 1; break;
        case C_SECH: f = C_COSH; recip_res = 1; break;
        case C_CSCH: f = C_SINH; recip_res = 1; break;
        default: break;
    }

    if (recip_arg)
        CHK(nfloat_complex_inv(z, x, ctx));
    else
        nfloat_complex_set(z, (nfloat_complex_ptr) x, ctx);

    /* sinh z = -i sin(iz), cosh z = cos(iz), tanh z = -i tan(iz),
       asinh z = -i asin(iz), atanh z = -i atan(iz) */
    switch (f)
    {
        case C_SINH: f = C_SIN; rot = 3; break;
        case C_COSH: f = C_COS; rot = 4; break;
        case C_TANH: f = C_TAN; rot = 3; break;
        case C_ASINH: f = C_ASIN; rot = 3; break;
        case C_ATANH: f = C_ATAN; rot = 3; break;
        default: break;
    }

    if (rot)
        _c_mul_i_pow(z, z, 1, ctx);

    a = RE(z);
    b = IM(z);

    if (NFLOAT_IS_ZERO(b) && f != C_EXP_PI_I && f != C_LOG_PI_I &&
        _in_domain(a, _cfunc_real[f].domain, ctx))
    {
        CHK(_cfunc_real[f].real(RE(w), a, ctx));
        nfloat_zero(IM(w), ctx);
    }
    else switch (f)
    {
        case C_EXP:
            CHK(_c_exp(w, a, b, ctx));
            break;

        case C_EXPM1:
            if (NFLOAT_EXP(a) <= -1 && NFLOAT_EXP(b) <= 0)
            {
                /* Re = expm1(a) cos b - 2 sin^2(b/2), Im = exp(a) sin b */
                RTMP(s); RTMP(c);
                CHK(nfloat_mul_2exp_si(t, b, -1, ctx));
                CHK(nfloat_sin(s, t, ctx));
                CHK(nfloat_sqr(s, s, ctx));
                NFLOAT_EXP(s) += !NFLOAT_IS_ZERO(s);
                CHK(nfloat_sin_cos(IM(w), c, b, ctx));
                CHK(nfloat_expm1(t, a, ctx));
                CHK(nfloat_mul(RE(w), t, c, ctx));
                CHK(nfloat_sub(RE(w), RE(w), s, ctx));
                CHK(nfloat_exp(t, a, ctx));
                CHK(nfloat_mul(IM(w), IM(w), t, ctx));
            }
            else
            {
                CHK(_c_exp(w, a, b, ctx));
                nfloat_one(t, ctx);
                CHK(nfloat_sub(RE(w), RE(w), t, ctx));
            }
            break;

        case C_EXP_PI_I:
            /* exp(pi i z) = exp(-pi b) (cos(pi a) + i sin(pi a)) */
            CHK(nfloat_sin_cos_pi(IM(w), RE(w), a, ctx));
            if (!NFLOAT_IS_ZERO(b))
            {
                CHK(nfloat_pi(t, ctx));
                CHK(nfloat_mul(t, t, b, ctx));
                nfloat_neg(t, t, ctx);
                CHK(nfloat_exp(t, t, ctx));
                CHK(_c_mul_real(w, w, t, ctx));
            }
            break;

        case C_EXP2:
        case C_EXP10:
            CHK(nfloat_set_ui(t, (f == C_EXP2) ? 2 : 10, ctx));
            CHK(nfloat_log(t, t, ctx));
            CHK(_c_mul_real(w, z, t, ctx));
            CHK(_c_exp(w, RE(w), IM(w), ctx));
            break;

        case C_LOG:
        case C_LOG2:
        case C_LOG10:
        case C_LOG_PI_I:
            if (NFLOAT_IS_ZERO(a) && NFLOAT_IS_ZERO(b))
                return GR_DOMAIN;
            CHK(_c_log(w, a, b, ctx));
            if (f == C_LOG2 || f == C_LOG10)
            {
                CHK(nfloat_set_ui(t, (f == C_LOG2) ? 2 : 10, ctx));
                CHK(nfloat_log(t, t, ctx));
                CHK(_c_div_real(w, w, t, ctx));
            }
            else if (f == C_LOG_PI_I)
            {
                /* log(z) / (pi i) */
                CHK(nfloat_pi(t, ctx));
                CHK(_c_div_real(w, w, t, ctx));
                _c_mul_i_pow(w, w, 3, ctx);
            }
            break;

        case C_LOG1P:
            if (NFLOAT_EXP(a) <= -1 && NFLOAT_EXP(b) <= -1)
            {
                /* Re = log1p(a (2 + a) + b^2) / 2, Im = atan2(b, 1 + a) */
                RTMP(u);
                nfloat_one(t, ctx);
                CHK(nfloat_add(u, t, a, ctx));
                CHK(nfloat_atan2(IM(w), b, u, ctx));
                CHK(nfloat_add(u, u, t, ctx));
                CHK(nfloat_mul(u, u, a, ctx));
                CHK(nfloat_sqr(t, b, ctx));
                CHK(nfloat_add(u, u, t, ctx));
                CHK(nfloat_log1p(RE(w), u, ctx));
                NFLOAT_EXP(RE(w)) -= !NFLOAT_IS_ZERO(RE(w));
            }
            else
            {
                nfloat_one(t, ctx);
                CHK(nfloat_add(t, t, a, ctx));
                if (NFLOAT_IS_ZERO(t) && NFLOAT_IS_ZERO(b))
                    return GR_DOMAIN;
                CHK(_c_log(w, t, b, ctx));
            }
            break;

        case C_SIN:
        case C_SIN_PI:
            CHK(_c_sin_cos(w, NULL, a, b, f == C_SIN_PI, ctx));
            break;

        case C_COS:
        case C_COS_PI:
            CHK(_c_sin_cos(NULL, w, a, b, f == C_COS_PI, ctx));
            break;

        case C_TAN:
        case C_TAN_PI:
            CHK(_c_tan(w, a, b, f == C_TAN_PI, ctx));
            break;

        case C_SINC:
        case C_SINC_PI:
            CHK(_c_sin_cos(w, NULL, a, b, f == C_SINC_PI, ctx));
            CHK(nfloat_complex_div(w, w, z, ctx));
            if (f == C_SINC_PI)
            {
                CHK(nfloat_pi(t, ctx));
                CHK(_c_div_real(w, w, t, ctx));
            }
            break;

        case C_ASIN:
        case C_ACOS:
        case C_ACOSH:
            CHK(_c_asin_acos(w, z, f, ctx));
            break;

        case C_ATAN:
            CHK(_c_atan(w, a, b, ctx));
            break;

        default:
            return GR_UNABLE;
    }

    if (rot)
        _c_mul_i_pow(w, w, rot, ctx);

    if (recip_res)
        CHK(nfloat_complex_inv(w, w, ctx));

    if (div_pi)
    {
        CHK(nfloat_pi(t, ctx));
        CHK(_c_div_real(w, w, t, ctx));
    }

    nfloat_complex_set(res, w, ctx);
    return GR_SUCCESS;
}

#define DEF(name, F) \
    int nfloat_complex_ ## name(nfloat_complex_ptr res, nfloat_complex_srcptr x, gr_ctx_t ctx) \
    { return _nfloat_complex_func(res, x, F, ctx); }

DEF(exp, C_EXP)
DEF(expm1, C_EXPM1)
DEF(exp_pi_i, C_EXP_PI_I)
DEF(exp2, C_EXP2)
DEF(exp10, C_EXP10)
DEF(log, C_LOG)
DEF(log1p, C_LOG1P)
DEF(log_pi_i, C_LOG_PI_I)
DEF(log2, C_LOG2)
DEF(log10, C_LOG10)
DEF(sin, C_SIN)
DEF(cos, C_COS)
DEF(tan, C_TAN)
DEF(cot, C_COT)
DEF(sec, C_SEC)
DEF(csc, C_CSC)
DEF(sin_pi, C_SIN_PI)
DEF(cos_pi, C_COS_PI)
DEF(tan_pi, C_TAN_PI)
DEF(cot_pi, C_COT_PI)
DEF(sec_pi, C_SEC_PI)
DEF(csc_pi, C_CSC_PI)
DEF(sinc, C_SINC)
DEF(sinc_pi, C_SINC_PI)
DEF(sinh, C_SINH)
DEF(cosh, C_COSH)
DEF(tanh, C_TANH)
DEF(coth, C_COTH)
DEF(sech, C_SECH)
DEF(csch, C_CSCH)
DEF(asin, C_ASIN)
DEF(acos, C_ACOS)
DEF(atan, C_ATAN)
DEF(acot, C_ACOT)
DEF(asec, C_ASEC)
DEF(acsc, C_ACSC)
DEF(asinh, C_ASINH)
DEF(acosh, C_ACOSH)
DEF(atanh, C_ATANH)
DEF(acoth, C_ACOTH)
DEF(asech, C_ASECH)
DEF(acsch, C_ACSCH)
DEF(asin_pi, C_ASIN_PI)
DEF(acos_pi, C_ACOS_PI)
DEF(atan_pi, C_ATAN_PI)
DEF(acot_pi, C_ACOT_PI)
DEF(asec_pi, C_ASEC_PI)
DEF(acsc_pi, C_ACSC_PI)

int
nfloat_complex_sin_cos(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    CTMP(s); CTMP(c);

    if (NFLOAT_COMPLEX_IS_SPECIAL(x, ctx) && !(NFLOAT_IS_ZERO(RE(x)) && NFLOAT_IS_ZERO(IM(x))) &&
        ((NFLOAT_IS_SPECIAL(RE(x)) && !NFLOAT_IS_ZERO(RE(x))) || (NFLOAT_IS_SPECIAL(IM(x)) && !NFLOAT_IS_ZERO(IM(x)))))
        return GR_UNABLE;

    CHK(_c_sin_cos(s, c, RE(x), IM(x), 0, ctx));
    nfloat_complex_set(res1, s, ctx);
    nfloat_complex_set(res2, c, ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_sin_cos_pi(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    CTMP(s); CTMP(c);

    if (NFLOAT_COMPLEX_IS_SPECIAL(x, ctx) &&
        ((NFLOAT_IS_SPECIAL(RE(x)) && !NFLOAT_IS_ZERO(RE(x))) || (NFLOAT_IS_SPECIAL(IM(x)) && !NFLOAT_IS_ZERO(IM(x)))))
        return GR_UNABLE;

    CHK(_c_sin_cos(s, c, RE(x), IM(x), 1, ctx));
    nfloat_complex_set(res1, s, ctx);
    nfloat_complex_set(res2, c, ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_sinh_cosh(nfloat_complex_ptr res1, nfloat_complex_ptr res2, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    CTMP(z);

    /* sinh z = -i sin(iz), cosh z = cos(iz) */
    _c_mul_i_pow(z, x, 1, ctx);
    CHK(nfloat_complex_sin_cos(res1, res2, z, ctx));
    _c_mul_i_pow(res1, res1, 3, ctx);
    return GR_SUCCESS;
}

int
nfloat_complex_pow(nfloat_complex_ptr res, nfloat_complex_srcptr x, nfloat_complex_srcptr y, gr_ctx_t ctx)
{
    CTMP(t);

    if (NFLOAT_COMPLEX_IS_SPECIAL(x, ctx) || NFLOAT_COMPLEX_IS_SPECIAL(y, ctx))
    {
        if (NFLOAT_COMPLEX_IS_ZERO(y, ctx))
            return nfloat_complex_one(res, ctx);

        if (NFLOAT_COMPLEX_IS_ZERO(x, ctx))
        {
            /* 0^y = 0 for Re(y) > 0 */
            if (!NFLOAT_IS_SPECIAL(RE(y)) && !NFLOAT_SGNBIT(RE(y)))
                return nfloat_complex_zero(res, ctx);
            return GR_DOMAIN;
        }
    }

    /* real powers of positive reals */
    if (NFLOAT_IS_ZERO(IM(x)) && NFLOAT_IS_ZERO(IM(y)) && !NFLOAT_SGNBIT(RE(x)))
    {
        CHK(nfloat_pow(RE(res), RE(x), RE(y), ctx));
        return nfloat_zero(IM(res), ctx);
    }

    /* integer powers (|y| < 2^16) by repeated squaring */
    if (NFLOAT_IS_ZERO(IM(y)) && !NFLOAT_IS_SPECIAL(RE(y)) && NFLOAT_EXP(RE(y)) <= 16)
    {
        RTMP(u);

        if (nfloat_nint(u, RE(y), ctx) == GR_SUCCESS && nfloat_equal(u, RE(y), ctx) == T_TRUE)
        {
            fmpz_t e;
            int status;

            fmpz_init(e);
            status = nfloat_get_fmpz(e, u, ctx);
            if (status == GR_SUCCESS)
                status = gr_generic_pow_fmpz(res, x, e, ctx);
            fmpz_clear(e);
            return status;
        }
    }

    /* exp(y log(x)) */
    CHK(nfloat_complex_log(t, x, ctx));
    CHK(nfloat_complex_mul(t, t, y, ctx));
    return nfloat_complex_exp(res, t, ctx);
}
