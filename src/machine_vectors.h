/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023 Mathieu Gouttenoire

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MACHINE_VECTORS_H
#define MACHINE_VECTORS_H

#define ALIGN_STRUCT(x) __attribute__((aligned(x)))

#include <math.h>
#include <string.h>

#include "flint.h"

/*
    The AVX2 and NEON backends provide integer vectors with ulong lanes
    and use intrinsics that exist only in 64-bit mode, so they require a
    64-bit word size; a 32-bit build uses the generic backends. This
    matters because machine_vectors.h is included by machine_vectors/
    gemm.c on every target, whereas it was previously reached only from
    fft_small, which is built only in 64-bit configurations.

    FLINT_MACHINE_VECTORS_FORCE_GENERIC selects the generic backends even
    on a target that has AVX2 or NEON; FLINT_MACHINE_VECTORS_STRICT_C
    additionally selects the strict ISO C tier over GNU vector
    extensions. Both are intended for testing and profiling.
*/
#if defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        || defined(FLINT_MACHINE_VECTORS_STRICT_C) \
        || FLINT_BITS != 64
# define FLINT_MACHINE_VECTORS_GENERIC 1
#endif

#if !defined(FLINT_MACHINE_VECTORS_GENERIC)
# if defined(__GNUC__)
#  if defined(__AVX2__)
#   include <immintrin.h>
#  elif defined(__ARM_NEON)
#   include <arm_neon.h>
#  endif
# elif defined(_MSC_VER)
#  if defined(__AVX2__)
#   include <intrin.h>
#  elif defined(_M_ARM64)
#   include <arm_neon.h>
#  endif
# endif
#endif

#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Exactly one of FLINT_MACHINE_VECTORS_AVX2, FLINT_MACHINE_VECTORS_NEON,
    FLINT_MACHINE_VECTORS_GENERIC ends up defined; it identifies the backend
    in use and can be tested by code that uses backend specific operations.
*/
#if defined(__AVX2__) && !defined(FLINT_MACHINE_VECTORS_GENERIC)
# define FLINT_MACHINE_VECTORS_AVX2 1
# include "machine_vectors_avx2.h"
#elif (defined(__ARM_NEON) || defined(_M_ARM64)) \
        && !defined(FLINT_MACHINE_VECTORS_GENERIC)
# define FLINT_MACHINE_VECTORS_NEON 1
# include "machine_vectors_neon.h"
#else
# ifndef FLINT_MACHINE_VECTORS_GENERIC
#  define FLINT_MACHINE_VECTORS_GENERIC 1
# endif
# include "machine_vectors_generic.h"
#endif



/* gemm ********************************************************************/

/*
    C = A * B for row-major matrices with no transposes and no
    accumulation: C is m x n, A is m x k, B is k x n, with leading
    dimensions ldc, lda, ldb. Equivalent to cblas_sgemm/cblas_dgemm with
    CblasRowMajor, CblasNoTrans, CblasNoTrans, alpha = 1, beta = 0.
    These are always available; they use FLINT's thread pool according to
    flint_get_num_threads().
*/
/*
    Three interchangeable implementations of each. The _blas versions
    call cblas (and abort if FLINT was built without BLAS); the
    _fallback versions are FLINT's own kernels, always available. The
    unsuffixed versions dispatch on flint_gemm_use_blas, which is
    initialized to FLINT_USES_BLAS and may be set at runtime.
*/
FLINT_DLL extern int flint_gemm_use_blas;

void flint_sgemm(slong m, slong k, slong n,
                 const float * A, slong lda,
                 const float * B, slong ldb,
                 float * C, slong ldc);

void flint_dgemm(slong m, slong k, slong n,
                 const double * A, slong lda,
                 const double * B, slong ldb,
                 double * C, slong ldc);

void flint_sgemm_blas(slong m, slong k, slong n,
                      const float * A, slong lda,
                      const float * B, slong ldb,
                      float * C, slong ldc);

void flint_dgemm_blas(slong m, slong k, slong n,
                      const double * A, slong lda,
                      const double * B, slong ldb,
                      double * C, slong ldc);

void flint_sgemm_fallback(slong m, slong k, slong n,
                          const float * A, slong lda,
                          const float * B, slong ldb,
                          float * C, slong ldc);

void flint_dgemm_fallback(slong m, slong k, slong n,
                          const double * A, slong lda,
                          const double * B, slong ldb,
                          double * C, slong ldc);

#ifdef __cplusplus
}
#endif

#endif
