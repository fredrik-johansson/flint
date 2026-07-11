/* Hand-optimized replacements for the n = 1 and n = 2 specialized exp
   series (see WINNERS_OPT in dev/gen_fixed_exp_rs_hard.py, which skips
   these two sizes and defers to this file).  At these sizes the RS
   code generator's machinery -- the shift/split framework, the calls
   to flint_mpn_sqrhigh / flint_mpn_mulhigh_n / mpn_addmul_1, the
   fac20_inv rescaling -- dwarfs the handful of word multiplies the
   series actually needs, so the evaluators are written directly.

   Both compute exp(x) for a reduced argument x with x < 2^-r (r = 12
   for n = 1, r = 16 for n = 2), returning the unit limb (always 1,
   since exp(x) < 1 + 2x < 2 here) and n fraction limbs.  Horner is
   used on the tail coefficients c_k = 1/k!, each a hardcoded n-limb
   fixed-point value (the binary point above the top limb).  The
   truncated tail x^(N+1)/(N+1)! is far below the output ulp, and the
   products are "high" products keeping only the top n limbs; at n = 2
   the a0*b0 partial is dropped entirely (MULH2), which the error
   budget (thousands of ulp) easily absorbs. */

/* --- 1-limb high product: top limb of a*b for fractions a, b < 1 --- */
#define FIXED_EXP_MULH1(a, b) \
    ({ ulong _h, _l; umul_ppmm(_h, _l, (a), (b)); (void) _l; _h; })

/* --- 2-limb "sloppy" high product: top 2 limbs of (a1:a0)*(b1:b0),
   dropping the a0*b0 partial (its contribution is below the retained
   precision).  Keeps the a1*b1 high part and both cross products. --- */
#define FIXED_EXP_MULH2(rh, rl, a1, a0, b1, b0) \
    do { \
        ulong __c1h, __c1l, __c2h, __c2l, __sh, __sl; \
        umul_ppmm(rh, rl, a1, b1); \
        umul_ppmm(__c1h, __c1l, a1, b0); \
        umul_ppmm(__c2h, __c2l, a0, b1); \
        add_ssaaaa(rh, rl, rh, rl, UWORD(0), __c1h); \
        add_ssaaaa(rh, rl, rh, rl, UWORD(0), __c2h); \
        add_ssaaaa(__sh, __sl, UWORD(0), __c1l, UWORD(0), __c2l); \
        add_ssaaaa(rh, rl, rh, rl, UWORD(0), __sh); \
    } while (0)

/* c_k = floor(2^64 / k!): the 1-limb Horner coefficients */
#define FIXED_EXP1_C2 UWORD(0x8000000000000000)
#define FIXED_EXP1_C3 UWORD(0x2aaaaaaaaaaaaaaa)
#define FIXED_EXP1_C4 UWORD(0x0aaaaaaaaaaaaaaa)

/* exp(x), x < 2^-12, one fraction limb.  N = 4 terms suffice: the
   dropped x^5/5! is about 2^-67, far below the 2^-64 output ulp.
   exp = 1 + x + x^2 (C2 + x (C3 + x C4)). */
void
_fixed_exp_rs_opt_1(nn_ptr res, nn_srcptr x)
{
    ulong X = x[0], p, x2, h, l;

    p = FIXED_EXP1_C3 + FIXED_EXP_MULH1(X, FIXED_EXP1_C4);
    p = FIXED_EXP1_C2 + FIXED_EXP_MULH1(X, p);
    x2 = FIXED_EXP_MULH1(X, X);
    p = FIXED_EXP_MULH1(x2, p);

    /* fraction = x + x^2 p (< 2^-11, so the unit is exactly 1) */
    add_ssaaaa(h, l, UWORD(0), X, UWORD(0), p);
    res[0] = l;
    res[1] = 1 + h;
}

/* c_k = floor(2^128 / k!) as (hi, lo): the 2-limb Horner coefficients */
#define FIXED_EXP2_C2_HI UWORD(0x8000000000000000)
#define FIXED_EXP2_C2_LO UWORD(0x0000000000000000)
#define FIXED_EXP2_C3_HI UWORD(0x2aaaaaaaaaaaaaaa)
#define FIXED_EXP2_C3_LO UWORD(0xaaaaaaaaaaaaaaaa)
#define FIXED_EXP2_C4_HI UWORD(0x0aaaaaaaaaaaaaaa)
#define FIXED_EXP2_C4_LO UWORD(0xaaaaaaaaaaaaaaaa)
#define FIXED_EXP2_C5_HI UWORD(0x0222222222222222)
#define FIXED_EXP2_C5_LO UWORD(0x2222222222222222)
#define FIXED_EXP2_C6_HI UWORD(0x005b05b05b05b05b)
#define FIXED_EXP2_C6_LO UWORD(0x05b05b05b05b05b0)
#define FIXED_EXP2_C7_HI UWORD(0x000d00d00d00d00d)
#define FIXED_EXP2_C7_LO UWORD(0x00d00d00d00d00d0)

/* exp(x), x < 2^-16, two fraction limbs.  N = 7 terms: the dropped
   x^8/8! is about 2^-143, below the 2^-128 output ulp.
   exp = 1 + x + x^2 (C2 + x (C3 + ... + x C7)). */
void
_fixed_exp_rs_opt_2(nn_ptr res, nn_srcptr x)
{
    ulong X1 = x[1], X0 = x[0];
    ulong p1, p0, h, l, x2h, x2l;

    p1 = FIXED_EXP2_C7_HI; p0 = FIXED_EXP2_C7_LO;
    FIXED_EXP_MULH2(h, l, X1, X0, p1, p0);
    add_ssaaaa(p1, p0, FIXED_EXP2_C6_HI, FIXED_EXP2_C6_LO, h, l);
    FIXED_EXP_MULH2(h, l, X1, X0, p1, p0);
    add_ssaaaa(p1, p0, FIXED_EXP2_C5_HI, FIXED_EXP2_C5_LO, h, l);
    FIXED_EXP_MULH2(h, l, X1, X0, p1, p0);
    add_ssaaaa(p1, p0, FIXED_EXP2_C4_HI, FIXED_EXP2_C4_LO, h, l);
    FIXED_EXP_MULH2(h, l, X1, X0, p1, p0);
    add_ssaaaa(p1, p0, FIXED_EXP2_C3_HI, FIXED_EXP2_C3_LO, h, l);
    FIXED_EXP_MULH2(h, l, X1, X0, p1, p0);
    add_ssaaaa(p1, p0, FIXED_EXP2_C2_HI, FIXED_EXP2_C2_LO, h, l);

    FIXED_EXP_MULH2(x2h, x2l, X1, X0, X1, X0);
    FIXED_EXP_MULH2(h, l, x2h, x2l, p1, p0);

    /* fraction = x + x^2 p (< 2^-15, so the unit is exactly 1) */
    add_ssaaaa(h, l, h, l, X1, X0);
    res[0] = l;
    res[1] = h;
    res[2] = 1;
}
