#!/usr/bin/env python3
"""Regenerate src/fixed/sin_cos_bitwise_rs_small.inc (run from the
top-level FLINT directory:
python3 dev/gen_fixed_sin_cos_bitwise_small.py
> src/fixed/sin_cos_bitwise_rs_small.inc): per-size register
implementations of the fixed_sin_cos_bitwise_rs reduction for
n = 1..7.

This combines the two earlier small paths.  The DECISIONS are the
table-driven ones of the exp reduction (compare the angle against the
tabulated A_i = atan(2^-i)), while the STATE that has to be carried
along is a rotating vector, as in the atan vectoring: the deferred
rotation W = prod (1 + i 2^-i) is applied on the fly, so no used[]
list is needed at all.

Registers: the reduced angle t0..t{n-1} (which holds x/2, the table
being the undoubled angles), the rotation W = (wu.p{n-1}..p0,
q{n-1}..q0), and nothing else.  W's real part needs an explicit units
limb wu: unlike atan's X, which only grows from 1, wx DECREASES here
(wx -= wy 2^-i), down to about 0.72, so wu falls from 1 to 0.  The
imaginary part stays in [0, 0.92) and fits n fraction limbs with no
carry out.  Both updates read the OLD components, so u = trunc(wx >>
i) and v = trunc(wy >> i) are materialized first; the accept mask then
gates wx -= v, wy += u and the subtraction of A_i from t.

Because the table holds the UNDOUBLED angles, entries obey
A_i < 2^-i, so after the boundary step the residual is strictly below
the window and the single-limb model needs no stray-bit check (atan's
geometric residual, bounded only by X 2^-i < 2^(1-i), does).  The
boundary step and the step at i = R are each emitted TWICE: the table
entries are floors, so truncation creep can leave the remainder above
the entry after one subtraction (see exp_bitwise_rs.c), and two
passes restore the invariant.

The reduction runs to R = r + 1 and the residual is doubled before
the series, since the rotation angle of a factor is 2 A_i.
"""

SUB = {2: "sub_ddmmss", 3: "sub_dddmmmsss", 4: "sub_ddddmmmmssss",
       5: "sub_dddddmmmmmsssss", 6: "sub_ddddddmmmmmmssssss",
       7: "sub_dddddddmmmmmmmsssssss",
       8: "sub_ddddddddmmmmmmmmssssssss",
       9: "sub_dddddddddmmmmmmmmmsssssssss"}
ADD = {2: "add_ssaaaa", 3: "add_sssaaaaaa", 4: "add_ssssaaaaaaaa",
       5: "add_sssssaaaaaaaaaa", 6: "add_ssssssaaaaaaaaaaaa",
       7: "add_sssssssaaaaaaaaaaaaaa",
       8: "add_ssssssssaaaaaaaaaaaaaaaa"}

NMIN, NMAX = 1, 7


def add_into(o, ind, dst, src, n):
    if n == 1:
        o("%s%s0 += %s;" % (ind, dst, src[0]))
        return
    dl = ["%s%d" % (dst, j) for j in range(n - 1, -1, -1)]
    o("%s%s(%s," % (ind, ADD[n], ", ".join(dl)))
    o("%s    %s," % (ind, ", ".join(dl)))
    o("%s    %s);" % (ind, ", ".join(src[j] for j in range(n - 1, -1, -1))))


def sub_from(o, ind, dsts, srcs, k):
    """dsts[0..k-1] -= srcs[0..k-1] (dsts listed low to high)"""
    if k == 1:
        o("%s%s -= %s;" % (ind, dsts[0], srcs[0]))
        return
    dl = [dsts[j] for j in range(k - 1, -1, -1)]
    sl = [srcs[j] for j in range(k - 1, -1, -1)]
    o("%s%s(%s," % (ind, SUB[k], ", ".join(dl)))
    o("%s    %s," % (ind, ", ".join(dl)))
    o("%s    %s);" % (ind, ", ".join(sl)))


def masked_step(o, n, c, ind="        "):
    """one masked accept of the factor at index i.  Assumes Ap points
    at the tabulated angle and that u0.. / v0.. hold trunc(wx >> i)
    and trunc(wy >> i) (n limbs each), taken from the OLD W."""
    tl = ["t%d" % j for j in range(n - 1, -1, -1)]
    al = ["Ap[%d]" % j for j in range(n - 1, -1, -1)]
    el = ["e%d" % j for j in range(n - 1, -1, -1)]
    # borrow-capturing chain: bw = all-ones exactly when A_i > t
    o("%s%s(bw, %s," % (ind, SUB[n + 1], ", ".join(el)))
    o("%s    UWORD(0), %s," % (ind, ", ".join(tl)))
    o("%s    UWORD(0), %s);" % (ind, ", ".join(al)))
    o("%sm = ~bw;                /* accept iff no borrow */" % ind)
    for j in range(n):
        o("%st%d = (t%d & bw) | (e%d & m);" % (ind, j, j, j))
    # wx -= v & m over the n + 1 limbs (the units limb falls 1 -> 0)
    sub_from(o, ind, ["p%d" % j for j in range(n)] + ["wu"],
             ["v%d & m" % j for j in range(n)] + ["UWORD(0)"], n + 1)
    # wy += u & m (stays below 0.92: no carry out)
    add_into(o, ind, "q", ["u%d & m" % j for j in range(n)], n)


def emit_uv(o, n, c, ind, b):
    """materialize u = trunc(wx >> i) and v = trunc(wy >> i), n limbs
    each, from the OLD W.  wx's limbs are (p0..p{n-1}, wu) and wy's
    are (q0..q{n-1}, 0).  b is the C expression of the runtime shift,
    or None for the b == 0 boundary form."""
    if b is None:
        for j in range(n):
            k = j + c
            src = ("p%d" % k) if k < n else ("wu" if k == n else "UWORD(0)")
            o("%su%d = %s;" % (ind, j, src))
        for j in range(n):
            k = j + c
            src = ("q%d" % k) if k < n else "UWORD(0)"
            o("%sv%d = %s;" % (ind, j, src))
        return
    # b > 0: one funnel shift per limb; the limbs above n - 1 - c are
    # zero (wu >> b vanishes for b >= 1)
    for j in range(n - 1 - c):
        o("%su%d = MPN_RIGHT_SHIFT_LOW(p%d, p%d, %s);"
          % (ind, j, j + c + 1, j + c, b))
    o("%su%d = MPN_RIGHT_SHIFT_LOW(wu, p%d, %s);"
      % (ind, n - 1 - c, n - 1, b))
    for j in range(n - c, n):
        o("%su%d = UWORD(0);" % (ind, j))
    for j in range(n - 1 - c):
        o("%sv%d = MPN_RIGHT_SHIFT_LOW(q%d, q%d, %s);"
          % (ind, j, j + c + 1, j + c, b))
    o("%sv%d = MPN_RIGHT_SHIFT_LOW(UWORD(0), q%d, %s);"
      % (ind, n - 1 - c, n - 1, b))
    for j in range(n - c, n):
        o("%sv%d = UWORD(0);" % (ind, j))


def gen_one(o, n):
    o("static void")
    o("_fixed_sin_cos_bitwise_rs_%d(nn_ptr ysin, nn_ptr ycos," % n)
    o("    nn_srcptr x, int r)")
    o("{")
    o("    ulong " + ", ".join("t%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("p%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("q%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("e%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("u%d" % j for j in range(n)) + ";")
    o("    ulong " + ", ".join("v%d" % j for j in range(n)) + ";")
    o("    ulong wu, bw, m, h, lt;")
    o("    ulong tt[%d], sr[%d], cr[%d];" % (n, n + 1, n + 1))
    o("    ulong wx[%d], wy[%d];" % (n + 1, n + 1))
    o("    slong i, nc, R;")
    o("    int pass;")
    o("")
    o("    r = FLINT_MAX(r, 16);")
    o("    r = FLINT_MIN(r, FLINT_BITS * %d - 16);" % n)
    o("    R = (slong) r + 1;")
    o("")
    o("    _fixed_atans_ensure(%d, R);" % n)
    o("    nc = _fixed_atans_n;")
    o("")
    o("    /* t = x / 2 (the table holds the undoubled angles) */")
    o("    mpn_rshift(tt, x, %d, 1);" % n)
    for j in range(n):
        o("    t%d = tt[%d];" % (j, j))
    o("")
    o("    wu = 1;")
    for j in range(n):
        o("    p%d = 0;" % j)
    for j in range(n):
        o("    q%d = 0;" % j)
    o("")
    o("#define AP(ii) (_fixed_atans + (ii) * nc + (nc - %d))" % n)
    o("")

    for c in range(n):
        i0 = 1 if c == 0 else 64 * c
        o("    /* window %d */" % c)
        guard = "" if c == 0 else "if (R >= %d)\n    " % (64 * c)
        o("    %s{" % guard)
        o("        nn_srcptr Ap;")
        o("")
        o("        /* boundary step i = %d, twice (truncation creep) */"
          % i0)
        o("        Ap = AP(%d);" % i0)
        o("        for (pass = 0; pass < 2; pass++)")
        o("        {")
        if c == 0:
            emit_uv(o, n, 0, "            ", "1")
        else:
            emit_uv(o, n, c, "            ", None)
        masked_step(o, n, c, ind="            ")
        o("        }")
        o("")
        o("        h = t%d;" % (n - 1 - c))
        o("")
        o("        for (i = %d; i <= FLINT_MIN(R, %d); i++)"
          % (i0 + 1, 64 * c + 63))
        o("        {")
        o("            int b = (int) (i - %d);" % (64 * c))
        o("")
        o("            Ap = AP(i);")
        o("            lt = Ap[%d];" % (n - 1 - c))
        o("            if (h < lt)")
        o("                continue;    /* certain reject */")
        o("")
        emit_uv(o, n, c, "            ", "b")
        masked_step(o, n, c, ind="            ")
        o("            h = t%d;" % (n - 1 - c))
        o("        }")
        o("    }")
        o("")

    o("    /* the step at i = R, twice (truncation creep) */")
    o("    {")
    o("        nn_srcptr Ap = AP(R);")
    o("        int b = (int) (R & (FLINT_BITS - 1));")
    o("")
    o("        for (pass = 0; pass < 2; pass++)")
    o("        {")
    o("            switch (R / FLINT_BITS)")
    o("            {")
    for c in range(n):
        o("            case %d:" % c)
        if c >= 1:
            o("                if (b == 0)")
            o("                {")
            emit_uv(o, n, c, "                    ", None)
            masked_step(o, n, c, ind="                    ")
            o("                    break;")
            o("                }")
        emit_uv(o, n, c, "                ", "b")
        masked_step(o, n, c, ind="                ")
        o("                break;")
    o("            }")
    o("        }")
    o("    }")
    o("")
    o("#undef AP")
    o("")
    o("    /* the residual angle is twice the reduced remainder */")
    for j in range(n):
        o("    tt[%d] = t%d;" % (j, j))
    o("    mpn_lshift(tt, tt, %d, 1);" % n)
    o("")
    o("    if (r < 32)")
    o("        _fixed_sin_cos_rs16(sr, cr, tt, %d);" % n)
    o("    else")
    o("        fixed_sin_cos_rs(sr, cr, tt, %d);" % n)
    o("")
    for j in range(n):
        o("    wx[%d] = p%d;" % (j, j))
    o("    wx[%d] = wu;" % n)
    for j in range(n):
        o("    wy[%d] = q%d;" % (j, j))
    o("    wy[%d] = 0;" % n)
    o("")
    o("    _fixed_sin_cos_recon(ysin, ycos, sr, cr, wx, wy, %d);" % n)
    o("}")
    o("")


def main():
    out = []
    o = out.append
    o("/* Generated by dev/gen_fixed_sin_cos_bitwise_small.py --")
    o("   register implementations of the fixed_sin_cos_bitwise_rs")
    o("   reduction for small n -- do not edit by hand. */")
    o("")
    for n in range(NMIN, NMAX + 1):
        gen_one(o, n)
    o("static void (* const _fixed_sin_cos_bitwise_rs_small_tab[])")
    o("    (nn_ptr, nn_ptr, nn_srcptr, int) = {")
    o("    NULL,")
    for n in range(NMIN, NMAX + 1):
        o("    _fixed_sin_cos_bitwise_rs_%d," % n)
    o("};")
    print("\n".join(out))


if __name__ == "__main__":
    main()
