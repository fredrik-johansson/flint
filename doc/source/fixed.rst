.. _fixed:

**fixed.h** -- efficient low-level fixed-point real arithmetic
===============================================================================

This module provides low-level, low-overhead functions for
fixed-point real arithmetic, intended as kernels for
arbitrary-precision floating-point arithmetic and ball arithmetic.
The scope currently includes routines for elementary function
evaluation following [HJ2024]_ and [Joh2014c]_.

A fixed-point number `(x, n)` is an unsigned `n`-limb fraction
``x[0], ..., x[n-1]`` representing
`\sum_i x_i \, 2^{\mathrm{FLINT\_BITS} (i - n)}`, i.e. `0 \le x < 1`
with unit in the last place (ulp) `2^{-\mathrm{FLINT\_BITS}\, n}`.
Outputs of size `n + 1` additionally carry an integer (units) limb at
index `n`.

On 32-bit limbs, all generated straight-line and register
implementations are disabled and evaluation goes through the generic
and fallback code paths, which compile cleanly but have not been
exercised on 32-bit hardware; :func:`fixed_exp_bitwise_rs`
additionally requires `n \ge 2` there.

Taylor series evaluation
-------------------------------------------------------------------------------

The series evaluation functions require an argument reduced below
`2^{-32}` (checked with a ``FLINT_ASSERT``) and dispatch internally on
the top limb of `x`: nonzero selects hardcoded straight-line routines
for the 32-bit reduction range, zero selects hardcoded or windowed
general routines exploiting a whole leading zero limb (or more).  All
evaluation uses rectangular splitting.  Error bounds are given in ulp
by the following macros; for :func:`fixed_exp_rs` and the hyperbolic
functions they are one-sided (the computed result never exceeds the
true value), for the alternating functions two-sided.

.. macro:: FIXED_EXP_RS_MAX_ERR(n)
           FIXED_SIN_RS_MAX_ERR(n)
           FIXED_COS_RS_MAX_ERR(n)
           FIXED_SIN_COS_RS_MAX_ERR(n)
           FIXED_SINH_RS_MAX_ERR(n)
           FIXED_COSH_RS_MAX_ERR(n)
           FIXED_SINH_COSH_RS_MAX_ERR(n)
           FIXED_ATAN_RS_MAX_ERR(n)
           FIXED_ATANH_RS_MAX_ERR(n)

    Bounds, in ulp, for the error of the corresponding function
    below.  These are small constants.

.. function:: void fixed_exp_rs(nn_ptr res, nn_srcptr x, slong n)

    Sets `(res, n + 1)` to an approximation of `\exp((x, n))`,
    requiring `x < 2^{-32}`.

.. function:: void fixed_sin_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_cos_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n)
              void fixed_sinh_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_cosh_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, nn_srcptr x, slong n)

    Set `(res, n + 1)` to an approximation of the respective function
    of `(x, n)`, requiring `x < 2^{-32}`.  The combined versions allow
    either output pointer to be *NULL*.

.. function:: void fixed_atan_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_atanh_rs(nn_ptr res, nn_srcptr x, slong n)

    Sets `(res, n)` to an approximation of `\operatorname{atan}((x, n))`
    resp. `\operatorname{atanh}((x, n))`, requiring `x < 2^{-32}`.

Fallbacks
-------------------------------------------------------------------------------

The following fallbacks run at constant full precision with
coefficients generated on the fly, so they support arbitrary `n`; they
serve out-of-table sizes, 32-bit machines, and the test code (which
compares the optimized functions against the fallbacks called with one
extra limb of precision).

.. function:: void _fixed_exp_rs_fallback(nn_ptr res, nn_srcptr x, slong n)
              void _fixed_sin_cos_rs_fallback(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n, int alternating)
              void _fixed_atan_rs_fallback(nn_ptr res, nn_srcptr x, slong n, int alternating)

    As the corresponding functions above, requiring only
    `x < 2^{-32}` (`x < 2^{-16}` for the atan/atanh fallback; the
    number of terms is chosen from the actual leading zero bits of
    `x`).

.. function:: void _fixed_atanh_rs16(nn_ptr res, nn_srcptr x, slong n)

    Internal: as :func:`fixed_atanh_rs` but for the wider range
    `x < 2^{-16}`, backed by generated straight-line code for
    `n \le 4`; used by the small reduction parameters of
    :func:`fixed_log1p_bitwise_rs`.

Exponential with bitwise argument reduction
-------------------------------------------------------------------------------

.. macro:: FIXED_EXP_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_exp_bitwise_rs`,
    currently `9 r + 100`.  The bound grows linearly with `r` because
    all work happens at the output precision; callers wanting sub-ulp
    accuracy should pad the precision by one limb themselves, which
    absorbs the bound with more than 50 bits to spare.

.. function:: void fixed_exp_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n + 1)` to an approximation of `\exp((x, n))` for any
    `0 \le x < 1` (`n \ge 2` on 32-bit limbs, where `n = 1` would
    clamp `r` below the series contract).  Passing `r = 0` selects a
    tuned default from a built-in table of crossovers (the largest
    tabulated value serving all larger `n`); the table can be
    regenerated for a given machine with
    ``src/fixed/tune/tune-bitwise-r.c``.  The argument is reduced below `2^{-r}` (where
    `r \ge 32` is a tuning parameter, internally clamped to
    `\mathrm{FLINT\_BITS} \, n - 16`) by subtracting in turn each
    logarithm `L_i = \log(1 + 2^{-i})`, `i = 0, 1, \ldots, r`, for
    which `L_i \le x`; the Taylor series is evaluated on the reduced
    argument; and the used factors `1 + 2^{-i}` are multiplied back
    in, each by a single shift-and-add.  When the reduced series is
    long, `\sinh` is evaluated instead (half the terms) and the
    exponential reconstructed as
    `\exp(t) = \sinh(t) + \sqrt{1 + \sinh(t)^2}`.

    The logarithm table is generated at runtime and cached per thread
    at the largest precision and index range requested so far.

    Small `r` minimizes table and reduction work, large `r` shortens
    the series; as a rule of thumb on 64-bit machines, `r = 32` is
    best up to about 8 limbs, `r = 64` up to about 32 limbs, and
    `r = 128` or `192` beyond.  See ``profile/p-exp_bitwise_rs.c``.

Logarithm with bitwise argument reduction
-------------------------------------------------------------------------------

.. macro:: FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_log1p_bitwise_rs`,
    currently `3 r + 64`; as for the exponential, all work happens at
    the output precision, so callers wanting sub-ulp accuracy should
    pad the precision by one limb themselves.

.. function:: void fixed_log1p_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of `\log(1 + (x, n))` for any
    `0 \le x < 1`, by the dual of the bitwise exponential reduction
    (the L-mode BKM recurrence): a product `P` of factors
    `1 + 2^{-i}`, `i = 1, \ldots, r`, is built up greedily below
    `X = 1 + x`, each accepted factor costing one shift-and-add; the
    residual is evaluated as
    `\log(X/P) = 2 \operatorname{atanh}((X - P)/(X + P))`, where the
    numerator is the exact deficit maintained through the reduction
    and the single division fuses the normalization by `P` with the
    atanh transformation, so that the quotient satisfies the
    :func:`fixed_atanh_rs` contract directly (and the odd series
    needs half the terms of `\log(1 + u)`); finally the tabulated
    logarithms of the used factors are added.  The same cached table
    serves :func:`fixed_exp_bitwise_rs` and this function.

    Passing `r = 0` selects a tuned default, as for
    :func:`fixed_exp_bitwise_rs`.  Otherwise requires `r \ge 16`:
    values below 32 shorten the reduction
    further and take effect for `n \le 4`, where the atanh series
    for the wider range `t < 2^{-16}` is available; the general path
    clamps `r` up to 32, matching the :func:`fixed_atanh_rs`
    contract.  Generated register implementations (decisions and
    updates in straight-line masked carry chains) serve
    `n \le 7`.

    The significant length of `P` grows gradually (by `i` bits per
    accepted factor), so the per-factor work starts out at a single
    limb and `P` is only ever truncated once its exact length exceeds
    `n` limbs.

Sine and cosine with bitwise argument reduction
-------------------------------------------------------------------------------

.. macro:: FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_sin_cos_bitwise_rs`,
    currently `6 r + 128`; as elsewhere in this module all work
    happens at the output precision, so callers wanting sub-ulp
    accuracy should pad the precision by one limb themselves.

.. function:: void fixed_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n, int r)

    Sets `(ysin, n+1)` and `(ycos, n+1)` to approximations of
    `\sin((x, n))` and `\cos((x, n))` for any `0 \le x < 1`; either
    output may be ``NULL``.  The reduction is the rotation analog of
    :func:`fixed_exp_bitwise_rs`, based on the identity of [Joh2022]:
    with `\alpha = 2 \operatorname{atan}(b/a)`,

    .. math ::

        e^{ix} = e^{i(x - k\alpha)}
            \left( \frac{a + bi}{a - bi} \right)^k,

    whose rotation unit is a ratio of complex conjugates and is
    therefore *exactly* unimodular, so that -- unlike a CORDIC
    rotation -- no normalization by the growing modulus is needed.
    Taking `a = 2^i`, `b = 1` gives the rotation angles
    `2 \operatorname{atan}(2^{-i})`.  The cached table, however,
    holds the *undoubled* angles
    `A_i = \operatorname{atan}(2^{-i})`, and the greedy reduction is
    applied to `x/2` for `i = 1, \ldots, r+1`.  This is not merely
    an economy of one table shared with
    :func:`fixed_atan_bitwise_rs`: the windowed decision model of
    the shared reduction requires table entries below `2^{-i}`,
    which `A_i` satisfies -- as `\log(1 + 2^{-i})` does -- but
    `2 \operatorname{atan}(2^{-i}) \approx 2^{1-i}` does not.
    Concavity of `\operatorname{atan}` gives `A_{i-1} < 2 A_i`, so
    each index is used at most once; writing
    `x/2 = \sum A_i + t'`, the residual angle is `2 t'`, whence the
    reduction to `r + 1`.  The sine and cosine of
    the residual `t < 2^{-r}` are then evaluated by
    :func:`fixed_sin_cos_rs`, or, once the series is long enough to
    amortize it, by :func:`fixed_sin_rs` together with
    `\cos t = \sqrt{1 - \sin^2 t}` at the cost of one squaring and
    one integer square root -- the same tradeoff as the `\sinh` path
    of :func:`fixed_exp_bitwise_rs`.

    The rotation is applied in deferred form: with
    `W = \prod (1 + i 2^{-i})` over the used indices (the real
    denominators `2^i` cancel between `W` and its conjugate), each
    factor costs two shifts and two add/subtracts, and the result is
    read off from

    .. math ::

        (\cos x, \sin x) = (c + is) \, W^2 / |W|^2,

    where two squarings supply `\operatorname{Re} W^2` and `|W|^2`
    from the same pair.  Since the *same* truncated `W` is used
    against its own conjugate, the reconstruction is self
    normalizing: truncations perturb only the angle of the
    correction, never its modulus.

    Passing `r = 0` selects a tuned default.  Otherwise `r \ge 16`:
    values below 32 shorten the reduction further and are served by
    the wider-range series for `t < 2^{-16}`, which pays up to about
    `n = 4`.  Requires `n \ge 2` when ``FLINT_BITS == 32``.  The angle
    table is shared with :func:`fixed_atan_bitwise_rs`.

Inverse tangent with bitwise argument reduction
-------------------------------------------------------------------------------

.. macro:: FIXED_ATAN_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_atan_bitwise_rs`,
    currently `4 r + 64`.

.. function:: void fixed_atan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of
    `\operatorname{atan}((x, n))` for any `0 \le x < 1`, by the
    vectoring dual of the rotation above (as
    :func:`fixed_log1p_bitwise_rs` is the dual of
    :func:`fixed_exp_bitwise_rs`).  The vector `(X, Y) = (1, x)` is
    rotated towards the real axis by the factors `1 - i 2^{-i}`,
    `i = 1, \ldots, r`, applying a factor -- two shifts and two
    add/subtracts on the old components -- whenever it keeps
    `Y \ge 0`, that is whenever `Y \ge X 2^{-i}`.  This is a greedy
    on the angle with steps `A_i = \operatorname{atan}(2^{-i})`, so
    afterwards `\operatorname{atan}(Y/X) < 2^{-r}` and a single
    division yields a residual meeting the :func:`fixed_atan_rs`
    contract.  Because the angle is scale invariant, the growth of
    `|Z|` needs no compensation at all.  The tabulated angles of the
    used factors -- the same `A_i` cached for
    :func:`fixed_sin_cos_bitwise_rs`, which uses them against `x/2`
    -- are simply summed at the end; the total is below
    `\sum_{i \ge 1} A_i \approx 0.898 < 1`, so no rescaling is
    needed.

    Note that the decisions of the vectoring loop never consult the
    table (they compare `Y` against `X 2^{-i}`), so only the
    rotation direction constrains how the shared angles are scaled.

    Passing `r = 0` selects a tuned default.  Otherwise `r \ge 16`,
    values below 32 being served by the wider-range series for
    `t < 2^{-16}`; the shorter reduction pays up to about `n = 6`.
    Generated register implementations (the vectoring in straight-line
    masked borrow chains, bit-for-bit identical to the generic code)
    serve `n \le 7`.
