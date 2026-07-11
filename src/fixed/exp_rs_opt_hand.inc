/* Hand-optimized replacements for the n = 1 and n = 2 specialized exp
   series (see WINNERS_OPT in dev/gen_fixed_exp_rs_hard.py, which skips
   these two sizes and defers to this file).  At these sizes the RS
   code generator's machinery -- the shift/split framework, the calls
   to flint_mpn_sqrhigh / flint_mpn_mulhigh_n / mpn_addmul_1, the
   fac20_inv rescaling -- dwarfs the handful of word multiplies the
   series actually needs, so the evaluators are written directly.

   Both compute exp(x) for a reduced argument x with x < 2^-r (r = 12
   for n = 1, r = 16 for n = 2), returning n fraction limbs and the
   unit limb, which is always 1 (exp(x) < 1 + 2x < 2 here).

   The work is done at MIXED precision.  Writing the series as

       exp(x) = 1 + x + x^2 P,   P = C2 + x Q,   Q = C3 + x R,   ...

   with C_k = 1/k!, an error in an inner Horner value is attenuated by
   the powers of x multiplying it: err(P) enters the result scaled by
   x^2 < 2^-2r, and each level below that by a further x < 2^-r.  So at
   n = 2 (output ulp 2^-128, x < 2^-16), holding the series error near
   2^-122,

       P needs 90 bits, Q 74, R 58, S 42, U 26

   -- only P and Q need two limbs, while R, S and U (three of the five
   Horner steps) fit in ONE and are evaluated on x's top limb alone.
   The same argument at n = 1 leaves P needing just 34 bits.

   The products keep only the top limbs: the 2x2 form drops the a0*b0
   partial entirely, multiplying by a one-limb value drops a column,
   and the squaring uses the symmetry of its two cross products.  The
   error budget -- FIXED_EXP_BITWISE_RS_MAX_ERR is thousands of ulp
   against an output of hundreds of bits -- absorbs all of this many
   times over; the measured error is a handful of ulp. */

/* top limb of a * b for fractions a, b < 1 (statement form: ISO C has
   no braced-group expressions) */
#define FIXED_EXP_MULH1(r, a, b) \
    do { \
        ulong __h, __l; \
        umul_ppmm(__h, __l, (a), (b)); \
        (void) __l; \
        (r) = __h; \
    } while (0)

/* top 2 limbs of (a1:a0) * (b1:b0), dropping the a0*b0 partial */
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

/* top 2 limbs of (a1:a0) * b with b a ONE-limb fraction: the b0 column
   of the product above vanishes, leaving two multiplies */
#define FIXED_EXP_MULH2_1(rh, rl, a1, a0, b) \
    do { \
        ulong __ch, __cl; \
        umul_ppmm(rh, rl, a1, b); \
        umul_ppmm(__ch, __cl, a0, b); \
        add_ssaaaa(rh, rl, rh, rl, UWORD(0), __ch); \
        (void) __cl; \
    } while (0)

/* top 2 limbs of (a1:a0)^2: the two cross products coincide */
#define FIXED_EXP_SQRH2(rh, rl, a1, a0) \
    do { \
        ulong __ch, __cl; \
        umul_ppmm(rh, rl, a1, a1); \
        umul_ppmm(__ch, __cl, a1, a0); \
        add_ssaaaa(rh, rl, rh, rl, UWORD(0), __ch); \
        add_ssaaaa(rh, rl, rh, rl, UWORD(0), __ch); \
        (void) __cl; \
    } while (0)

/* c_k = floor(2^64 / k!), one limb */
#define FIXED_EXP_E2 UWORD(0x8000000000000000)
#define FIXED_EXP_E3 UWORD(0x2aaaaaaaaaaaaaaa)
#define FIXED_EXP_E4 UWORD(0x0aaaaaaaaaaaaaaa)
#define FIXED_EXP_E5 UWORD(0x0222222222222222)
#define FIXED_EXP_E6 UWORD(0x005b05b05b05b05b)
#define FIXED_EXP_E7 UWORD(0x000d00d00d00d00d)

/* c_k = floor(2^128 / k!), two limbs */
#define FIXED_EXP_D2_HI UWORD(0x8000000000000000)
#define FIXED_EXP_D2_LO UWORD(0x0000000000000000)
#define FIXED_EXP_D3_HI UWORD(0x2aaaaaaaaaaaaaaa)
#define FIXED_EXP_D3_LO UWORD(0xaaaaaaaaaaaaaaaa)

/* exp(x), x < 2^-12, one fraction limb.  N = 4 terms: the dropped
   x^5/5! is about 2^-67, below the 2^-64 output ulp.  Written as
   1 + x + x^2 (C2 + x C3 + x^2 C4), reusing x^2 so that the two small
   terms are formed independently rather than in a serial chain. */
void
_fixed_exp_rs_opt_1(nn_ptr res, nn_srcptr x)
{
    ulong X = x[0], x2, p, a, b, h, l;

    FIXED_EXP_MULH1(x2, X, X);
    FIXED_EXP_MULH1(a, X, FIXED_EXP_E3);
    FIXED_EXP_MULH1(b, x2, FIXED_EXP_E4);
    p = FIXED_EXP_E2 + a + b;
    FIXED_EXP_MULH1(p, x2, p);

    add_ssaaaa(h, l, UWORD(0), X, UWORD(0), p);
    res[0] = l;
    res[1] = 1 + h;
}

/* exp(x), x < 2^-16, two fraction limbs.  N = 7 terms: the dropped
   x^8/8! is about 2^-143, below the 2^-128 output ulp.

   R = C4 + x (C5 + x (C6 + x C7)) needs only 58 bits, so it and the
   two Horner steps below it run at a single limb, on x's top limb
   alone; Q = C3 + x R and P = C2 + x Q then need two. */
void
_fixed_exp_rs_opt_2(nn_ptr res, nn_srcptr x)
{
    ulong X1 = x[1], X0 = x[0];
    ulong u, s, r, q1, q0, p1, p0, h, l, x2h, x2l;

    /* inner tail at one limb */
    FIXED_EXP_MULH1(u, X1, FIXED_EXP_E7);
    u += FIXED_EXP_E6;
    FIXED_EXP_MULH1(s, X1, u);
    s += FIXED_EXP_E5;
    FIXED_EXP_MULH1(r, X1, s);
    r += FIXED_EXP_E4;

    /* Q = C3 + x R (R is a single limb) */
    FIXED_EXP_MULH2_1(h, l, X1, X0, r);
    add_ssaaaa(q1, q0, FIXED_EXP_D3_HI, FIXED_EXP_D3_LO, h, l);

    /* P = C2 + x Q */
    FIXED_EXP_MULH2(h, l, X1, X0, q1, q0);
    add_ssaaaa(p1, p0, FIXED_EXP_D2_HI, FIXED_EXP_D2_LO, h, l);

    /* x^2 P, then the fraction x + x^2 P */
    FIXED_EXP_SQRH2(x2h, x2l, X1, X0);
    FIXED_EXP_MULH2(h, l, x2h, x2l, p1, p0);
    add_ssaaaa(h, l, h, l, X1, X0);

    res[0] = l;
    res[1] = h;
    res[2] = 1;
}
