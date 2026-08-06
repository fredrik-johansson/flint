.. _fft-small:

**fft_small.h** -- FFT modulo word-size primes
===============================================================================

This module currently requires building FLINT with support for
AVX2 or NEON instructions.

The truncated FFT routines used internally by this module follow the
approach described by van der Hoeven in [vdH2004]_.

Integer multiplication
--------------------------------------------------------------------------------

.. type:: mpn_ctx_struct
          mpn_ctx_t

    Context object for multiplications allowing non-FFT moduli.
    The structure contains FFT context objects for multiple FFT primes
    (currently 8) together with tables for Chinese remaindering.

.. function:: void mpn_ctx_init(mpn_ctx_t R, ulong p)

    Initialize multiplication context object with initial prime ``p``.
    Usually ``p`` is a constant provided by :func:`get_default_mpn_ctx`.

.. function:: void mpn_ctx_clear(mpn_ctx_t R)

    Free memory allocated by the context object.

.. function:: mpn_ctx_struct * get_default_mpn_ctx(void)

    Return a pointer to a cached thread-local context object used by default
    for multiplications. Calling :func:`flint_cleanup` or :func:`flint_cleanup_master`
    frees the cache.

.. function:: void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong * r1, const ulong * i1, ulong n1, const ulong * i2, ulong n2)
              void mpn_mul_default_mpn_ctx(nn_ptr r1, nn_srcptr i1, slong n1, nn_srcptr i2, slong n2)

    Writes to ``r1`` the product of the integers ``(i1, n1)`` and ``(i2, n2)``.
    Assumes that `n_1 \ge n_2 \ge 1`, respectively using a given context
    object ``R`` or the default thread-local object.

Polynomial arithmetic
---------------------------------------------------------------------------------

.. function:: void _nmod_poly_mul_mid_mpn_ctx(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)
              void _nmod_poly_mul_mid_default_mpn_ctx(nn_ptr res, slong zl, slong zh, nn_srcptr a, slong an, nn_srcptr b, slong bn, nmod_t mod)

    Writes to ``z`` the middle product containing coefficients in the
    range `[zl, zh)` of the product of the polynomials  ``(a, an)`` and ``(b, bn)``,
    respectively using a given context object ``R`` or the default thread-local object.
    Assumes that `an \ge bn \ge 1`.

.. function:: int _fmpz_poly_mul_mid_mpn_ctx(fmpz * z, ulong zl, ulong zh, const fmpz * a, ulong an, const fmpz * b, ulong bn, mpn_ctx_t R)
              int _fmpz_poly_mul_mid_default_mpn_ctx(fmpz * z, ulong zl, ulong zh, const fmpz * a, ulong an, const fmpz * b, ulong bn)

    Like the ``nmod`` functions. Performs the multiplication and returns 1
    if there are sufficiently many primes ``R`` to compute the result;
    otherwise returns 0 without touching the output.

.. function:: void _nmod_poly_divrem_mpn_ctx(ulong * q, ulong * r, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)

    Polynomial division with remainder.

Preconditioned polynomial arithmetic
---------------------------------------------------------------------------------

.. type:: mul_precomp_struct

.. function:: void _mul_precomp_init(mul_precomp_struct * M, const ulong * b, ulong bn, ulong btrunc, ulong depth, nmod_t mod, mpn_ctx_t R)
              void _mul_precomp_clear(mul_precomp_struct * M)

    Represents ``(b, bn)`` in transformed form for preconditioned multiplication.

.. function:: int _nmod_poly_mul_mid_precomp(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, mul_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

    Polynomial multiplication given a precomputed transform ``M``.
    Returns 1 if successful, 0 if the precomputed transform is too short.

.. type:: nmod_poly_divrem_precomp_struct

.. function:: void _nmod_poly_divrem_precomp_init(nmod_poly_divrem_precomp_struct * M, const ulong * b, ulong bn, ulong Bn, nmod_t mod, mpn_ctx_t R)
              void _nmod_poly_divrem_precomp_clear(nmod_poly_divrem_precomp_struct * M)

    Represents ``(b, bn)`` and its inverse in transformed form for preconditioned multiplication.

.. function:: int _nmod_poly_divrem_precomp(ulong * q, ulong * r, const ulong * a, ulong an, nmod_poly_divrem_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

    Polynomial multiplication given a precomputed transform ``M``.
    Returns 1 if successful, 0 if the precomputed transform is too short.


Transform plans and operands
--------------------------------------------------------------------------------

.. type:: fft_small_plan_struct
          fft_small_plan_t

    Describes one convolution shape: the primes used, transform depth,
    truncation lengths and coefficient windows. A plan is computed once
    and reused for any number of transforms of that shape.

.. type:: fft_small_op_struct
          fft_small_op_t

    One operand in transformed representation: a buffer of evaluations
    for each prime, together with the plan parameters it was built for
    and its current domain (primal data or pointwise products).

.. function:: void fft_small_op_init(fft_small_op_t X, const fft_small_plan_t P)
              void fft_small_op_init_borrowed(fft_small_op_t X, const fft_small_plan_t P, double * data)
              ulong fft_small_op_sizeof_data(const fft_small_plan_t P)
              void fft_small_op_clear(fft_small_op_t X)

    Operand storage management. The *borrowed* variant places the
    operand on caller-provided storage of ``fft_small_op_sizeof_data``
    bytes, 4096-aligned, which must outlive the operand and is not freed
    by *clear*.

.. function:: void fft_small_export_mpn_signed(ulong * z, ulong zn, int * sign, const fft_small_op_t X, ulong nslots, const fft_small_plan_t P)
              void fft_small_export_mpn_signed_trunc(ulong * z, ulong zn, int * sign, const fft_small_op_t X, ulong nslots, ulong lo_limbs, const fft_small_plan_t P)

    Chinese-remainder reconstruction of a signed integer result: the
    value is assembled in two's complement and the sign resolved from
    the top limb, with the magnitude written to *z*. The *trunc* variant
    writes the limbs of the value starting at limb *lo_limbs*, with an error against the
    exact value within `(-1.5, +0.5)` ulp of the lowest returned limb --
    equivalently, at most 1 from the floor-truncated value: every slot
    whose coefficient span can reach the first returned limb is included
    with carries propagated, so the error is the truncation itself
    (one-sided, below 1 ulp) plus the total mass of the wholly dropped
    slots, under half an ulp either way. With ``lo_limbs = 0`` the
    export is exact.

.. function:: void fft_small_export_nmod_range(ulong * z, const fft_small_op_t X, ulong zl, ulong zh, nmod_t mod, const fft_small_plan_t P)

    Chinese remaindering restricted to the coefficient window
    `[zl, zh)`, which must lie inside the plan's window; *z* receives
    `zh - zl` reduced coefficients. Used for polynomial middle products.

.. macro:: FLINT_FFT_SMALL_ALIGNMENT

    Alignment, in bytes, of transform data buffers.

.. var:: ulong flint_fft_small_max_transformed_ring_size

    Upper bound, in bytes, on the transform storage a single transformed
    ring or context may plan for: its expected number of simultaneously
    live elements times the per-element data size. Constructors decline
    above the bound, so drivers fall back to slower algorithms rather
    than exhaust memory; the matrix multiplication drivers respond by
    multiplying in blocks instead. Mutable for tuning; the default is
    4 GiB on 64-bit machines.


Negacyclic pointwise multiplication
-------------------------------------------------------------------------------

These functions multiply residues modulo `2^N + 1` by splitting the
operands into `m` digits of `b = N/m` bits, so that the product is a
length-`m` negacyclic convolution of the digit sequences, computed via
the small-prime FFT machinery of this module and Chinese remaindering.
They serve as the pointwise stage of Schoenhage-Strassen convolutions:
see :func:`fft_convolution_negmul` and the ``fmpz_poly`` multiplication
functions, which engage this engine for large coefficient sizes.

Operands and results occupy `n + 1 = N/\mathtt{FLINT\_BITS} + 1` limbs.
An *operand* may have top limb `0` or `1` with arbitrary remaining
limbs; its value is interpreted modulo `2^N + 1` (in particular the
representation need not be reduced, and both representatives of a
residue with top limb set are accepted). The *result* is returned
fully reduced: the unique representative in `[0, 2^N]`, with top limb
`1` exactly when the value is `2^N` and the remaining limbs then
zero.

.. function:: void mpn_negmul_ctx_init(mpn_negmul_ctx_struct * C, mpn_ctx_struct * R, slong N, slong m)

    Initializes a context for multiplication modulo `2^N + 1` with
    digit count *m*. Requires `N = mb` for an integer digit size *b*,
    *m* a power of two at least 16, and *N* a multiple of
    ``FLINT_BITS``. The number of FFT primes is chosen internally as
    the smallest count, between 3 and 8, whose exact prime product
    exceeds the coefficient bound `2^{2b + \log_2 m + 1}` and for
    which a digit-splitting input pass is instantiated; fewer primes
    means proportionally less work in every stage. Raises an exception
    if no supported configuration exists for the given shape. The
    primes and their transform contexts are taken from *R*, normally
    ``get_default_mpn_ctx()``.

.. function:: void mpn_negmul_ctx_init_auto(mpn_negmul_ctx_struct * C, mpn_ctx_struct * R, slong N)

    As :func:`mpn_negmul_ctx_init`, choosing the digit count *m*
    internally over the instantiated digit sizes dividing *N* by an
    estimated-cost comparison. For power-of-two *N* the selection
    moves from `(np, b) = (3, 64)` to `(6, 128)` at half the
    transform length once three-prime capacity runs out, flattening
    the large-*N* cost staircase. Raises an exception if no
    supported splitting of *N* exists.

.. function:: void mpn_negmul2(mpn_negmul_ctx_struct * C, nn_ptr z, nn_srcptr x, nn_srcptr y, mpn_negmul_scratch_struct * S)
              void mpn_negmul_scratch_init(mpn_negmul_scratch_struct * S, const mpn_negmul_ctx_struct * C)
              void mpn_negmul_scratch_clear(mpn_negmul_scratch_struct * S)

    The context is read-only during multiplication; all mutable state
    lives in a scratch object, so a single context may serve several
    threads concurrently, each with its own scratch.
    :func:`mpn_negmul` is equivalent to :func:`mpn_negmul2` with the
    scratch owned by the context. The pointwise stage of the
    ``fmpz_poly`` Schoenhage-Strassen convolutions multiplies its
    coefficients in parallel over the available threads this way.

.. function:: void mpn_negmul_ctx_clear(mpn_negmul_ctx_struct * C)

    Releases the memory held by the context.

.. function:: void mpn_negmul(mpn_negmul_ctx_struct * C, nn_ptr z, nn_srcptr x, nn_srcptr y)

    Sets `z = xy \bmod 2^N + 1` for semi-reduced operands as described
    above. Aliasing of *z* with *x* or *y* is permitted. The digit
    residues are computed by the templated input passes of this
    module, transformed, multiplied pointwise, inverse transformed,
    recombined by Chinese remaindering, and folded back into the
    semi-reduced representation.

