.. _machine-vectors:

**machine_vectors.h** -- SIMD-accelerated operations on fixed-length vectors
===============================================================================

Vector types and operations mapping onto the target's SIMD instructions,
together with the ``flint_sgemm`` and ``flint_dgemm`` matrix
multiplication kernels built on top of them.

Backends are selected automatically: AVX2 (and AVX512 for the ``vec8dz``
family) on x86, NEON on ARM, and otherwise a generic backend, using GNU
vector extensions where the compiler supports them and plain ISO C
structs elsewhere. The backends live in ``machine_vectors_avx2.h``,
``machine_vectors_neon.h`` and ``machine_vectors_generic.h``, which are
included by ``machine_vectors.h`` and must not be included directly.
Exactly one of ``FLINT_MACHINE_VECTORS_AVX2``,
``FLINT_MACHINE_VECTORS_NEON`` and ``FLINT_MACHINE_VECTORS_GENERIC``
ends up defined and identifies the backend in use, which matters
because each backend provides a slightly different superset of the
common interface; the availability of each operation is listed with
the operation below.

The generic backend implements the union of the AVX2 and NEON
interfaces (everything below except the AVX512-only ``z`` types), with
the same semantics; in particular the modular arithmetic is bit for bit
identical to the other backends, see below, so ``fft_small`` builds
wherever the word size is 64 bits.

It is not, however, enabled everywhere it builds. Without a hardware
fused multiply-add the generic backend has to reach the exact remainder
of a modular product through 64-bit integer arithmetic, which leaves the
transforms several times slower, far enough that the ``fft_small`` based
algorithms lose to the ones they would otherwise replace. ``configure``
therefore enables the module only where the target has AVX2, NEON, or a
fused multiply-add, which on x86-64 means that a build for the bare
baseline architecture does without it. Tuning the crossover thresholds
per backend, so that the module could be enabled everywhere and used
only where it wins, would be the better answer.

The AVX2 and NEON backends, and the
integer and modular operations of the generic backend, require 64-bit
words since the integer vectors have :type:`ulong` lanes; a 32-bit
build gets only the floating point part of the generic backend, which
is what ``flint_sgemm``/``flint_dgemm`` need.

The entire header assumes IEEE 754 double precision arithmetic in
round-to-nearest mode. It must not be compiled with ``-ffast-math`` or
anything else that licenses value-changing floating point
transformations; the modular arithmetic depends on exact rounding.

For the vector operations to use the target's instructions, FLINT must
be built with appropriate compiler flags. ``configure`` chooses these
from the detected CPU; ``--enable-avx2`` and ``--enable-avx512`` force
them on.

Defining ``FLINT_MACHINE_VECTORS_FORCE_GENERIC`` before including this
header selects the generic backend even on a target with AVX2 or NEON,
and ``FLINT_MACHINE_VECTORS_STRICT_C`` additionally selects the ISO C
tier over GNU vector extensions. These are intended for testing and
profiling the portable code paths.

The strict ISO C tier has no way to express a vector operation, so
whether its operations become SIMD instructions is entirely up to the
compiler's SLP vectorizer; its performance therefore varies between
compilers and compiler versions, unlike the other backends. It exists
so that the header works on compilers without GNU vector extensions,
and is not the fallback used on GCC or clang.

Both tiers compute ``mulmod`` in one of two ways, chosen by whether the
target defines ``__FP_FAST_FMA``. With a hardware fma the exact
remainder is obtained the way the AVX2 and NEON backends obtain it, from
the error of the product; without one, fma() is a software routine and
the remainder goes through an exact 64 bit integer product instead. Both
compute the same integer, so results do not depend on the choice, but
the cost does: on x86-64 with AVX2 the fma path makes ``vec4d_mulmod``
2.6x faster, and the strict tier then matches the AVX2 backend. This is
what the generic backend is for, since a target with vectors and fma but
no hand written backend here is the usual case.

The two tiers are still worth measuring against each other on a new
target, since which one a compiler handles better is not obvious.
With GCC 13 on x86-64 the GNU vector tier is the faster of the two for
``flint_dgemm`` and ``flint_sgemm``, by about 1.2x and 1.6x with AVX2
enabled and by a little at the x86-64 baseline, which is why it is the
default. On a target that takes the integer path, the strict tier is the
faster one for ``mulmod``: the remainder is then a 64 bit integer
product, one instruction in a scalar register but synthesised from 32
bit multiplies in a vector one, since x86 has no vector 64 bit multiply
below AVX512DQ.
``src/machine_vectors/profile/p-backends.c`` measures both.

The generic backends express a fused multiply-add as ``a * b + c``,
which a compiler fuses into an FMA instruction only when floating-point
contraction is enabled. GCC in a strict ISO mode, which is how FLINT is
built, does not contract by default; code using these operations in a
performance-critical loop should request contraction, as
``machine_vectors/gemm.c`` does with ``#pragma GCC optimize
("fp-contract=fast")``. The AVX2, AVX512 and NEON backends use fused
intrinsics and are unaffected.

Some functions may require that vectors are aligned in memory.

For testing, ``src/machine_vectors/test`` checks each operation of the
active backend against scalar references (``t-ops.c``), checks the
modular arithmetic contracts stated below (``t-mod.c``), and
additionally instantiates both generic tiers under renamed identifiers
so that they are exercised, and compared bit for bit against the native
backend, even in an AVX2 or NEON build (``t-force_generic.c``).

Types
-------------------------------------------------------------------------------

.. type:: vec1n
          vec2n
          vec4n
          vec8n

    Vector with 1, 2, 4, or 8 :type:`ulong` entries.

.. type:: vec1d
          vec2d
          vec4d
          vec8d

    Vector with 1, 2, 4, or 8 ``double`` entries.

.. type:: vec1f
          vec4f
          vec8f
          vec16f

    Vector with 1, 4, 8, or 16 ``float`` entries.

.. type:: vec8dz
          vec16dz
          vec16fz
          vec32fz
          vec8nz

    Vectors backed by AVX512 registers, available only when building
    with AVX512F support: 8 or 16 ``double`` entries, 16 or 32 ``float``
    entries, and 8 :type:`ulong` entries respectively. The ``z`` suffix
    distinguishes these from the equally-named types built from pairs of
    narrower registers; ``vec8d``, for instance, remains a pair of AVX2
    registers, which gives more instruction level parallelism in
    existing code.

Printing
-------------------------------------------------------------------------------

.. function:: void vec4d_print(vec4d a)
              void vec4n_print(vec4n a)

Access and conversions
-------------------------------------------------------------------------------

.. function:: vec1d vec1d_load(const double * a)
              vec4d vec4d_load(const double * a)
              vec8d vec8d_load(const double * a)

.. function:: vec1d vec1d_load_aligned(const double * a)
              vec4d vec4d_load_aligned(const double * a)
              vec8d vec8d_load_aligned(const double * a)

.. function:: vec1d vec1d_load_unaligned(const double * a)
              vec4d vec4d_load_unaligned(const double * a)
              vec8d vec8d_load_unaligned(const double * a)
              vec4n vec4n_load_unaligned(const ulong * a)
              vec8n vec8n_load_unaligned(const ulong * a)

.. function:: void vec1d_store(double * z, vec1d a)
              void vec4d_store(double * z, vec4d a)
              void vec8d_store(double * z, vec8d a)

.. function:: void vec1d_store_aligned(double * z, vec1d a)
              void vec4d_store_aligned(double * z, vec4d a)
              void vec8d_store_aligned(double * z, vec8d a)

.. function:: void vec1d_store_unaligned(double * z, vec1d a)
              void vec4d_store_unaligned(double * z, vec4d a)
              void vec4n_store_unaligned(ulong * z, vec4n a)
              void vec8d_store_unaligned(double * z, vec8d a)

.. function:: double vec4d_get_index(vec4d a, const int i)
              double vec8d_get_index(vec8d a, int i)

    Extract the entry at index `i`.

.. function:: vec1d vec1d_set_d(double a)
              vec4d vec4d_set_d(double a)
              vec4n vec4n_set_n(ulong a)
              vec8d vec8d_set_d(double a)
              vec8n vec8n_set_n(ulong a)

    Set all entries to the same value.

.. function:: vec4d vec4d_set_d4(double a0, double a1, double a2, double a3)
              vec4n vec4n_set_n4(ulong a0, ulong a1, ulong a2, ulong a3)
              vec8d vec8d_set_d8(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7)

    Create vector from distinct entries.

.. function:: vec1n vec1d_convert_limited_vec1n(vec1d a)
              vec2n vec2d_convert_limited_vec2n(vec2d a)
              vec4n vec4d_convert_limited_vec4n(vec4d a)

    Given that each entry in the input vector is an exact integer in
    ``[0, 2^{52})``, convert it to :type:`ulong`.
    Note that ``vec1d`` and ``vec2d`` functions are only available on NEON.

.. function:: vec2d vec2n_convert_limited_vec2d(vec2n a)
              vec4d vec4n_convert_limited_vec4d(vec4n a)
              vec8d vec8n_convert_limited_vec8d(vec8n a)

    The inverse conversion, from :type:`ulong` entries to exact
    ``double`` entries. Requires the same assumption.

Permutations
-------------------------------------------------------------------------------

.. function:: vec4d vec4d_unpacklo(vec4d a, vec4d b)
              vec4d vec4d_unpackhi(vec4d a, vec4d b)
              vec4d vec4d_permute_0_2_1_3(vec4d a)
              vec4d vec4d_permute_3_1_2_0(vec4d a)
              vec4d vec4d_permute_3_2_1_0(vec4d a)
              vec4d vec4d_permute2_0_2(vec4d a, vec4d b)
              vec4d vec4d_permute2_1_3(vec4d a, vec4d b)
              vec4d vec4d_unpack_lo_permute_0_2_1_3(vec4d u, vec4d v)
              vec4d vec4d_unpack_hi_permute_0_2_1_3(vec4d u, vec4d v)
              vec4d vec4d_unpackhi_permute_3_1_2_0(vec4d u, vec4d v)
              vec4d vec4d_unpacklo_permute_3_1_2_0(vec4d u, vec4d v)

.. macro:: VEC4D_TRANSPOSE(z0, z1, z2, z3, a0, a1, a2, a3)

    Sets the rows ``z`` to the transpose of the 4x4 matrix
    given by rows ``a``.

Comparisons
-------------------------------------------------------------------------------

.. function:: int vec1d_same(double a, double b)
              int vec4d_same(vec4d a, vec4d b)
              int vec8d_same(vec8d a, vec8d b)

    Check whether the vectors are equal.

.. function:: vec4d vec4d_cmp_ge(vec4d a, vec4d b)
              vec4d vec4d_cmp_gt(vec4d a, vec4d b)

    Entrywise comparisons.

Arithmetic and basic operations
-------------------------------------------------------------------------------

.. function:: vec1d vec1d_round(vec1d a)
              vec4d vec4d_round(vec4d a)
              vec8d vec8d_round(vec8d a)

.. function:: vec1d vec1d_zero()
              vec4d vec4d_zero()
              vec8d vec8d_zero()

.. function:: vec1d vec1d_one()
              vec4d vec4d_one()
              vec8d vec8d_one()

.. function:: vec1d vec1d_add(vec1d a, vec1d b)
              vec1d vec1d_sub(vec1d a, vec1d b)
              vec4d vec4d_add(vec4d a, vec4d b)
              vec4d vec4d_sub(vec4d a, vec4d b)
              vec4n vec4n_add(vec4n a, vec4n b)
              vec4n vec4n_sub(vec4n a, vec4n b)
              vec8d vec8d_add(vec8d a, vec8d b)
              vec8d vec8d_sub(vec8d a, vec8d b)

.. function:: vec1d vec1d_addsub(vec1d a, vec1d b)
              vec4d vec4d_addsub(vec4d a, vec4d b)

.. function:: vec1d vec1d_neg(vec1d a)
              vec4d vec4d_neg(vec4d a)
              vec8d vec8d_neg(vec8d a)

.. function:: vec1d vec1d_abs(vec1d a)
              vec4d vec4d_abs(vec4d a)

.. function:: vec1d vec1d_max(vec1d a, vec1d b)
              vec1d vec1d_min(vec1d a, vec1d b)
              vec4d vec4d_max(vec4d a, vec4d b)
              vec4d vec4d_min(vec4d a, vec4d b)
              vec8d vec8d_max(vec8d a, vec8d b)
              vec8d vec8d_min(vec8d a, vec8d b)

.. function:: vec1d vec1d_mul(vec1d a, vec1d b)
              vec4d vec4d_mul(vec4d a, vec4d b)
              vec8d vec8d_mul(vec8d a, vec8d b)

.. function:: vec1d vec1d_half(vec1d a)
              vec4d vec4d_half(vec4d a)

.. function:: vec1d vec1d_div(vec1d a, vec1d b)
              vec4d vec4d_div(vec4d a, vec4d b)
              vec8d vec8d_div(vec8d a, vec8d b)

.. function:: vec1d vec1d_fmadd(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fmadd(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fmadd(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_fmsub(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fmsub(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fmsub(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_fnmadd(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fnmadd(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fnmadd(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_fnmsub(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fnmsub(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fnmsub(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_blendv(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_blendv(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_blendv(vec8d a, vec8d b, vec8d c)

.. function:: vec4n vec4n_bit_shift_right(vec4n a, ulong b)
              vec8n vec8n_bit_shift_right(vec8n a, ulong b)

.. function:: vec4n vec4n_bit_and(vec4n a, vec4n b)
              vec8n vec8n_bit_and(vec8n a, vec8n b)


Modular arithmetic
-------------------------------------------------------------------------------

These functions are used internally by the small-prime FFT. The
``double`` variants represent residues as integer valued doubles and
assume an odd modulus `n < 2^{50}` together with the precomputed
``ninv``, which must be exactly the correctly rounded double ``1.0/n``.
The moduli ``mpn_ctx`` picks are just below `2^{50}`, but the modulus
is not always one of those: :func:`nmod_poly_mul` transforms over the
caller's own modulus whenever that is prime, below `2^{50}` and
2-adically deep enough, from 20 bits upwards, and nothing in these
operations depends on where the modulus came from. The following
contracts, verified by ``t-mod.c`` over the whole range of sizes, are
derived in ``src/fft_small/mulmod_satisfies_bounds.c``:

- ``mulmod(a, b, n, ninv)`` computes the exact integer `a b - q n` with
  `q = \operatorname{round}(\operatorname{fl}(\operatorname{fl}(a b)
  \cdot ninv))`, for any integer valued operands with
  `|a b| < 4 n^2` (and `|a|, |b| < 2^{62}` so that the products of the
  generic backend do not overflow, which is no restriction in
  practice). The result is congruent to `a b` modulo `n` and lies in
  `(-9/8\, n, 9/8\, n)` when `|a b| < 2 n^2` and in
  `(-7/4\, n, 7/4\, n)` when `|a b| < 4 n^2`. If
  :func:`fft_small_mulmod_satisfies_bounds` holds for `n`, which is the
  case for the moduli ``fft_small`` selects, these ranges tighten to
  `(-n, n)` and `(-3/2\, n, 3/2\, n)` respectively.

- ``nmulmod`` is exactly ``-mulmod`` (both give `+0.0` when the result
  is zero).

- The reductions ``reduce_to_pm1n``, ``reduce_to_pm1no`` and
  ``reduce_to_0n`` compute `a - \operatorname{round}(a \cdot ninv)
  \cdot n` (the last one followed by ``reduce_pm1no_to_0n``) exactly,
  for every integer valued `a` with `|a| \le 2^{53} - n`, which is
  every value a ``double`` represents exactly with the quotient times
  the modulus still representable. The result is congruent to `a` and lies in
  `(-n, n)` (respectively `[0, n)`); `-0.0` is never produced. Note
  that the quotient is not restricted to `|q| \le 8`, which is as far
  as a naive floating point `a - q n` stays exact; ``fft_small`` does
  call these with larger quotients.

Since none of the operations above involve a fused multiply-add in the
computation of `q`, and the remainder is exact on every backend, their
results are identical bit for bit across the AVX2, NEON and generic
backends; ``fft_small`` computes the same transforms everywhere. This
is checked by ``t-force_generic.c``.

The remaining reductions are pure range maps whose exact per lane
semantics (relevant only on the boundary of their domains and for
signed zeros) are: ``reduce_pm1no_to_0n`` adds `n` exactly when the
sign bit of the lane is set, so `-0.0` maps to `+n`;
``reduce_2n_to_n`` subtracts `n` exactly when `a - n \ge +0.0`;
``reduce_0n_to_pmhn`` subtracts `n` exactly when `a > n/2`; and
``reduce_pm1n_to_pmhn`` adds `\mp n` exactly when `|a| > n/2`.

.. function:: int vec1d_same_mod(vec1d a, vec1d b, vec1d n, vec1d ninv)
              int vec4d_same_mod(vec4d a, vec4d b, vec4d n, vec4d ninv)

    Return whether `a` and `b` are the same mod `n`.

.. function:: vec1d vec1d_reduce_pm1no_to_0n(vec1d a, vec1d n)
              vec1d vec4d_reduce_pm1no_to_0n(vec4d a, vec4d n)
              vec8d vec8d_reduce_pm1no_to_0n(vec8d a, vec8d n)

    Return `a \bmod n` reduced to `[0,n)` assuming `a \in (-n,n)`.

.. function:: vec1d vec1d_reduce_to_pm1n(vec1d a, vec1d n, vec1d ninv)
              vec4d vec4d_reduce_to_pm1n(vec4d a, vec4d n, vec4d ninv)
              vec8d vec8d_reduce_to_pm1n(vec8d a, vec8d n, vec8d ninv)

    Return `a \bmod n` reduced to `(-n,n)` for integer valued
    `|a| \le 2^{53} - n`; see the contract above. Also available as
    ``vec2d_reduce_to_pm1n`` on the NEON and generic backends.

.. function:: vec1d vec1d_reduce_to_pm1no(vec1d a, vec1d n, vec1d ninv)
              vec4d vec4d_reduce_to_pm1no(vec4d a, vec4d n, vec4d ninv)
              vec8d vec8d_reduce_to_pm1no(vec8d a, vec8d n, vec8d ninv)

    The same operation under its historical name (the trailing ``o``
    is for the open interval `(-n, n)`). Also available as
    ``vec2d_reduce_to_pm1no`` on the NEON and generic backends.

.. function:: vec1d vec1d_reduce_0n_to_pmhn(vec1d a, vec1d n)
              vec4d vec4d_reduce_0n_to_pmhn(vec4d a, vec4d n)

    Return `a \bmod n` reduced to `[-n/2, n/2]` given `a \in [0,n]`.

.. function:: vec1d vec1d_reduce_pm1n_to_pmhn(vec1d a, vec1d n)
              vec4d vec4d_reduce_pm1n_to_pmhn(vec4d a, vec4d n)
              vec8d vec8d_reduce_pm1n_to_pmhn(vec8d a, vec8d n)

    Return `a \bmod n` reduced to `[-n/2, n/2]` given `a \in [-n,n]`.

.. function:: vec1d vec1d_reduce_2n_to_n(vec1d a, vec1d n)
              vec4d vec4d_reduce_2n_to_n(vec4d a, vec4d n)
              vec8d vec8d_reduce_2n_to_n(vec8d a, vec8d n)

    Return `a \bmod n` reduced to `[0,n)` given `a \in [0,2n)`.

.. function:: vec1d vec1d_reduce_to_0n(vec1d a, vec1d n, vec1d ninv)
              vec4d vec4d_reduce_to_0n(vec4d a, vec4d n, vec4d ninv)
              vec8d vec8d_reduce_to_0n(vec8d a, vec8d n, vec8d ninv)

    Return `a \bmod n` reduced to `[0,n)`.

.. function:: vec1d vec1d_mulmod(vec1d a, vec1d b, vec1d n, vec1d ninv)
              vec4d vec4d_mulmod(vec4d a, vec4d b, vec4d n, vec4d ninv)
              vec8d vec8d_mulmod(vec8d a, vec8d b, vec8d n, vec8d ninv)

    Return an integer congruent to `ab` modulo `n`, in the ranges
    stated in the contract above. Also available as ``vec2d_mulmod`` on
    the NEON and generic backends.

.. function:: vec1d vec1d_nmulmod(vec1d a, vec1d b, vec1d n, vec1d ninv)
              vec4d vec4d_nmulmod(vec4d a, vec4d b, vec4d n, vec4d ninv)
              vec8d vec8d_nmulmod(vec8d a, vec8d b, vec8d n, vec8d ninv)

    Exactly the negation of ``mulmod``, so an integer congruent to
    `-ab` modulo `n` in the mirrored ranges. Also available as
    ``vec2d_nmulmod`` on the NEON and generic backends.

.. function:: vec4n vec4n_addmod(vec4n a, vec4n b, vec4n n)
              vec8n vec8n_addmod(vec8n a, vec8n b, vec8n n)

    Return `a + b \bmod n` in `[0,n)` given `a, b \in [0, n)`; any
    `n \ne 0`, including `n > 2^{63}`, is allowed. Also available as
    ``vec1n_addmod`` and ``vec2n_addmod`` on the NEON and generic
    backends.

.. function:: vec4n vec4n_addmod_limited(vec4n a, vec4n b, vec4n n)
              vec8n vec8n_addmod_limited(vec8n a, vec8n b, vec8n n)

    The same with the assumption `n < 2^{63}`, which saves work. Also
    available as ``vec2n_addmod_limited`` on the NEON and generic
    backends.

Matrix multiplication
-------------------------------------------------------------------------------

These functions compute a matrix product in single or double precision.
They are always available: when FLINT is built with BLAS, the default is
to call it, and otherwise FLINT's own kernels are used. The intended use
is as a building block for exact linear algebra over `\mathbb{Z}` and
`\mathbb{Z}/n\mathbb{Z}`, for example in :func:`nmod_mat_mul_blas` and
:func:`fmpz_mat_mul_blas`.

All of these functions compute `C = AB` for row-major matrices, where
*C* is *m* by *n*, *A* is *m* by *k* and *B* is *k* by *n*, with
*ldc*, *lda* and *ldb* the respective leading dimensions (the number of
entries between the start of consecutive rows, which must be at least
the number of columns). No transposition or accumulation is performed:
the previous contents of *C* are overwritten, and `k = 0` sets *C* to
zero. This is equivalent to ``cblas_sgemm`` or ``cblas_dgemm`` called
with ``CblasRowMajor``, ``CblasNoTrans``, ``CblasNoTrans``, ``alpha``
equal to 1 and ``beta`` equal to 0. Aliasing of *C* with *A* or *B* is
not allowed.

The FLINT kernels are multithreaded internally according to
:func:`flint_get_num_threads`, using FLINT's thread pool. They handle
arbitrary dimensions, including thin and unbalanced shapes, without
requiring any padding or alignment of the input.

.. function:: void flint_sgemm(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
              void flint_dgemm(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

    Sets `C = AB`, calling either the BLAS or the FLINT implementation
    according to :var:`flint_gemm_use_blas`.

.. function:: void flint_sgemm_blas(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
              void flint_dgemm_blas(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

    Sets `C = AB` using ``cblas_sgemm`` or ``cblas_dgemm``. These raise
    an exception if FLINT was built without BLAS support.

.. function:: void flint_sgemm_fallback(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
              void flint_dgemm_fallback(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

    Sets `C = AB` using FLINT's own kernels, which are always available.

.. var:: int flint_gemm_use_blas

    Selects the implementation used by :func:`flint_sgemm` and
    :func:`flint_dgemm`. It is initialized to 1 if FLINT was built with
    BLAS support and 0 otherwise, and may be set to either value at
    runtime, for example to compare the two implementations. Setting it
    to 1 in a build without BLAS support will result in an exception
    when a gemm is attempted.
