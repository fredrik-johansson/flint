/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nrb.h"
#include "gr.h"
#include "gr_mat.h"

int _nrb_methods_initialized = 0;

gr_static_method_table _nrb_methods;

gr_method_tab_input _nrb_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) nrb_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL, (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_CTX_HAS_REAL_PREC, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC, (gr_funcptr) _nfloat_ctx_set_real_prec},
    {GR_METHOD_CTX_GET_REAL_PREC, (gr_funcptr) _nfloat_ctx_get_real_prec},

    {GR_METHOD_INIT,            (gr_funcptr) nrb_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) nrb_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) nrb_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) nrb_set},
    {GR_METHOD_RANDTEST,        (gr_funcptr) nrb_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) nrb_write},
    {GR_METHOD_ZERO,            (gr_funcptr) nrb_zero},
    {GR_METHOD_ONE,             (gr_funcptr) nrb_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) nrb_neg_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) nrb_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) nrb_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) nrb_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) nrb_equal},
    {GR_METHOD_SET,             (gr_funcptr) nrb_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) nrb_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) nrb_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) nrb_set_fmpz},

    {GR_METHOD_NEG,             (gr_funcptr) nrb_neg},
    {GR_METHOD_ADD,             (gr_funcptr) nrb_add},
    {GR_METHOD_SUB,             (gr_funcptr) nrb_sub},
    {GR_METHOD_MUL,             (gr_funcptr) nrb_mul},

    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) nrb_mul_2exp_si},

    {GR_METHOD_ABS,             (gr_funcptr) nrb_abs},

    {GR_METHOD_SIN,             (gr_funcptr) nrb_sin},
    {GR_METHOD_COS,             (gr_funcptr) nrb_cos},

    {0,                         (gr_funcptr) NULL},
};

int
nrb_ctx_init(gr_ctx_t ctx, slong prec)
{
    slong nlimbs;

    if (prec <= 0 || prec > NFLOAT_MAX_LIMBS * FLINT_BITS)
        return GR_UNABLE;

    nlimbs = (prec + FLINT_BITS - 1) / FLINT_BITS;

    ctx->which_ring = GR_CTX_RR_NRB;
    ctx->sizeof_elem = sizeof(ulong) * (nlimbs + NRB_HEADER_LIMBS);
    ctx->size_limit = WORD_MAX;

    NFLOAT_CTX_NLIMBS(ctx) = nlimbs;
    /* No NaNs or Infs; disallow underflow */
    NFLOAT_CTX_FLAGS(ctx) = 0;
    NFLOAT_CTX_RND(ctx) = 0;

    ctx->methods = _nrb_methods;

    if (!_nrb_methods_initialized)
    {
        gr_method_tab_init(_nrb_methods, _nrb_methods_input);
        _nrb_methods_initialized = 1;
    }

    return GR_SUCCESS;
}

int
nrb_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Real numbers (nrb, prec = ");
    gr_stream_write_si(out, NFLOAT_CTX_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

