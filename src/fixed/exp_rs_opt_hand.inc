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


/* n_mulhi is newer than FLINT 3.3; define it when building against an
   older release so that this file stands alone.  (This block can go
   once the module only targets releases that ship it.) */
#if __FLINT_RELEASE < __FLINT_RELEASE_NUM(3, 4, 0)
static inline ulong
n_mulhi(ulong a, ulong b)
{
#if FLINT_BITS == 64 && defined(__GNUC__)
    return ((__uint128_t) a * (__uint128_t) b) >> 64;
#else
    ulong hi, lo;
    umul_ppmm(hi, lo, a, b);
    return hi;
#endif
}
#endif

/* Top 2 limbs of (a1:a0) * (b1:b0), discarding the low words of both
   cross products and the a0*b0 partial entirely: three multiplies and
   two adds.  This mirrors the helper in ulong_extras/ll_is_prime.c,
   which is file-local there and would be worth promoting to a header. */
FLINT_FORCE_INLINE void
_fixed_mulhi_2x2_sloppy(ulong * r1, ulong * r0,
    ulong a1, ulong a0, ulong b1, ulong b0)
{
    ulong s2, u2, v3, v2, r3, r2;

    s2 = n_mulhi(a0, b1);
    u2 = n_mulhi(a1, b0);
    add_ssaaaa(r3, r2, UWORD(0), s2, UWORD(0), u2);
    umul_ppmm(v3, v2, a1, b1);
    add_ssaaaa(r3, r2, r3, r2, v3, v2);
    *r0 = r2;
    *r1 = r3;
}

/* Top 2 limbs of (a1:a0) * b with b a ONE-limb fraction: the b0 column
   of the product above vanishes, leaving two multiplies. */
FLINT_FORCE_INLINE void
_fixed_mulhi_2x1(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b)
{
    ulong rh, rl, ch, cl;

    umul_ppmm(rh, rl, a1, b);
    umul_ppmm(ch, cl, a0, b);
    add_ssaaaa(rh, rl, rh, rl, UWORD(0), ch);
    (void) cl;
    *r0 = rl;
    *r1 = rh;
}

/* Top 2 limbs of (a1:a0)^2: the two cross products coincide, so this
   takes two multiplies rather than three -- measurably better than
   routing the squaring through the generic 2x2 form above. */
FLINT_FORCE_INLINE void
_fixed_sqrhi_2x2(ulong * r1, ulong * r0, ulong a1, ulong a0)
{
    ulong rh, rl, ch, cl;

    umul_ppmm(rh, rl, a1, a1);
    umul_ppmm(ch, cl, a1, a0);
    add_ssaaaa(rh, rl, rh, rl, UWORD(0), ch);
    add_ssaaaa(rh, rl, rh, rl, UWORD(0), ch);
    (void) cl;
    *r0 = rl;
    *r1 = rh;
}

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
    ulong X = x[0], x2, p, h, l;

    x2 = n_mulhi(X, X);
    p = FIXED_EXP_E2 + n_mulhi(X, FIXED_EXP_E3)
                     + n_mulhi(x2, FIXED_EXP_E4);
    p = n_mulhi(x2, p);

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
    u = FIXED_EXP_E6 + n_mulhi(X1, FIXED_EXP_E7);
    s = FIXED_EXP_E5 + n_mulhi(X1, u);
    r = FIXED_EXP_E4 + n_mulhi(X1, s);

    /* Q = C3 + x R (R is a single limb) */
    _fixed_mulhi_2x1(&h, &l, X1, X0, r);
    add_ssaaaa(q1, q0, FIXED_EXP_D3_HI, FIXED_EXP_D3_LO, h, l);

    /* P = C2 + x Q */
    _fixed_mulhi_2x2_sloppy(&h, &l, X1, X0, q1, q0);
    add_ssaaaa(p1, p0, FIXED_EXP_D2_HI, FIXED_EXP_D2_LO, h, l);

    /* x^2 P, then the fraction x + x^2 P */
    _fixed_sqrhi_2x2(&x2h, &x2l, X1, X0);
    _fixed_mulhi_2x2_sloppy(&h, &l, x2h, x2l, p1, p0);
    add_ssaaaa(h, l, h, l, X1, X0);

    res[0] = l;
    res[1] = h;
    res[2] = 1;
}
