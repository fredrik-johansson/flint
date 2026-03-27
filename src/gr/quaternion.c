/*
    Copyright (C) 2026 Fredrik Johansson
    Copyright (C) 2026 Kerem Kelleboz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "gr_generic.h"

#define QUATERNION_K_AS_GEN 1

typedef struct
{
	gr_ctx_struct * basef_ctx;
	int flags;
} _gr_quaternion_ctx_struct;

typedef gr_ctx_struct gr_quaternion_ctx_struct;
typedef gr_quaternion_ctx_struct gr_quaternion_ctx_t[1];

#define GR_QUATERNION_CTX(ring_ctx) ((_gr_quaternion_ctx_struct *)((ring_ctx)))
#define GR_QUATERNION_BASEF_CTX(ring_ctx) (GR_QUATERNION_CTX(ring_ctx)->basef_ctx)
#define GR_QUATERNION_FLAGS(ring_ctx) (GR_QUATERNION_CTX(ring_ctx)->flags)

#define BASEF_CTX GR_QUATERNION_BASEF_CTX
#define RE(x, ctx) (x)
#define IM(x, ctx) GR_ENTRY((x), 1, GR_QUATERNION_BASEF_CTX(ctx)->sizeof_elem)
#define JM(x, ctx) GR_ENTRY((x), 2, GR_QUATERNION_BASEF_CTX(ctx)->sizeof_elem)
#define KM(x, ctx) GR_ENTRY((x), 3, GR_QUATERNION_BASEF_CTX(ctx)->sizeof_elem)


static int
_gr_quaternion_ctx_write(gr_stream_t out, gr_quaternion_ctx_t ctx)
{
	gr_stream_write(out, "Quaternion algebra over ");
	gr_ctx_write(out, BASEF_CTX(ctx));
	return GR_SUCCESS;
}

static void
_gr_quaternion_ctx_clear(gr_ctx_t ctx)
{
}

static truth_t _gr_quaternion_ctx_is_ring(gr_ctx_t ctx) { return gr_ctx_is_ring(BASEF_CTX(ctx)); }
static truth_t _gr_quaternion_ctx_is_integral_domain(gr_ctx_t ctx) { return T_FALSE; }
static truth_t _gr_quaternion_ctx_is_field(gr_ctx_t ctx) { return T_FALSE; }
static truth_t _gr_quaternion_ctx_is_commutative_ring(gr_ctx_t ctx) { return gr_ctx_is_zero_ring(BASEF_CTX(ctx)); }
static truth_t _gr_quaternion_ctx_is_rational_vector_space(gr_ctx_t ctx) { return gr_ctx_is_rational_vector_space(BASEF_CTX(ctx)); }
static truth_t _gr_quaternion_ctx_is_real_vector_space(gr_ctx_t ctx) { return gr_ctx_is_real_vector_space(BASEF_CTX(ctx)); }

static truth_t _gr_quaternion_ctx_is_threadsafe(gr_ctx_t ctx) { return gr_ctx_is_threadsafe(BASEF_CTX(ctx)); }
static truth_t _gr_quaternion_ctx_is_finite(gr_ctx_t ctx) { return gr_ctx_is_finite(BASEF_CTX(ctx)); }
static truth_t _gr_quaternion_ctx_is_finite_characteristic(gr_ctx_t ctx) { return gr_ctx_is_finite_characteristic(BASEF_CTX(ctx)); }
static truth_t _gr_quaternion_ctx_is_exact(gr_ctx_t ctx) { return gr_ctx_is_exact(BASEF_CTX(ctx)); }

static int
_gr_quaternion_ctx_ngens(slong *n, gr_ctx_t ctx)
{
    *n = (GR_QUATERNION_FLAGS(ctx) & QUATERNION_K_AS_GEN) ? 3 : 2;
    return GR_SUCCESS;
}

static int
_gr_quaternion_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    slong num_gens;
    if (GR_QUATERNION_FLAGS(ctx) & QUATERNION_K_AS_GEN)
        num_gens = 3;
    else
        num_gens = 2;
    if (i < 0 || i >= num_gens)    
        return GR_DOMAIN;

    *name = flint_malloc(2);
    if (*name == NULL)
        return GR_UNABLE;
    if( i == 0 )
        (*name)[0] = 'I';
    else if( i == 1 )
        (*name)[0] = 'J';
    else 
        (*name)[0] = 'K';    
    (*name)[1] = '\0';
    return GR_SUCCESS;
}

static void
_gr_quaternion_init(gr_ptr x, gr_quaternion_ctx_t ctx)
{
	gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
	gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

	gr_init(a1, basef_ctx);
	gr_init(a2, basef_ctx);
	gr_init(a3, basef_ctx);
	gr_init(a4, basef_ctx);
}

static void
_gr_quaternion_clear(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

    gr_clear(a1, basef_ctx);
    gr_clear(a2, basef_ctx);
    gr_clear(a3, basef_ctx);
    gr_clear(a4, basef_ctx);
}

static int
_gr_quaternion_write(gr_stream_t out, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;
    slong i;

    int a1zero = (gr_is_zero(a1, basef_ctx) == T_TRUE);
    int a2zero = (gr_is_zero(a2, basef_ctx) == T_TRUE);
    int a3zero = (gr_is_zero(a3, basef_ctx) == T_TRUE);
    int a4zero = (gr_is_zero(a4, basef_ctx) == T_TRUE);

    if (a2zero && a3zero && a4zero)
    {
        status |= gr_write(out, a1, basef_ctx);
    }
    else if (a1zero && a3zero && a4zero)
    {
        if (gr_is_one(a2, basef_ctx) != T_TRUE)
        {
            status |= gr_stream_write(out, "(");
            status |= gr_write(out, a2, basef_ctx);
            status |= gr_stream_write(out, ") * ");
        }

        status |= gr_stream_write(out, "I");
    }
    else if (a1zero && a2zero && a4zero)
    {
        if (gr_is_one(a3, basef_ctx) != T_TRUE)
        {
            status |= gr_stream_write(out, "(");
            status |= gr_write(out, a3, basef_ctx);
            status |= gr_stream_write(out, ") * ");
        }

        status |= gr_stream_write(out, "J");
    }
    else if (a1zero && a2zero && a3zero)
    {
        if (gr_is_one(a4, basef_ctx) != T_TRUE)
        {
            status |= gr_stream_write(out, "(");
            status |= gr_write(out, a4, basef_ctx);
            status |= gr_stream_write(out, ") * ");
        }

        if (GR_QUATERNION_FLAGS(ctx) & QUATERNION_K_AS_GEN)
            status |= gr_stream_write(out, "K");
        else
            status |= gr_stream_write(out, "I * J");
    }
    else
    {
        status |= gr_stream_write(out, "(");
        status |= gr_write(out, a1, basef_ctx);
        status |= gr_stream_write(out, ") + (");
        status |= gr_write(out, a2, basef_ctx);
        status |= gr_stream_write(out, ") * I + (");
        status |= gr_write(out, a3, basef_ctx);
        status |= gr_stream_write(out, ") * J + (");
        status |= gr_write(out, a4, basef_ctx);

        if (GR_QUATERNION_FLAGS(ctx) & QUATERNION_K_AS_GEN)
            status |= gr_stream_write(out, ") * K");
        else
            status |= gr_stream_write(out, ") * I * J");
    }

    return status;
}

static void
_gr_quaternion_set_shallow(gr_ptr x, gr_srcptr y, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_srcptr b1 = RE(y, ctx), b2 = IM(y, ctx), b3 = JM(y, ctx), b4 = KM(y, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

    gr_set_shallow(a1, b1, basef_ctx);
    gr_set_shallow(a2, b2, basef_ctx);
    gr_set_shallow(a3, b3, basef_ctx);
    gr_set_shallow(a4, b4, basef_ctx);
}

static void
_gr_quaternion_swap(gr_ptr x, gr_ptr y, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ptr b1 = RE(y, ctx), b2 = IM(y, ctx), b3 = JM(y, ctx), b4 = KM(y, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

    gr_swap(a1, b1, basef_ctx);
    gr_swap(a2, b2, basef_ctx);
    gr_swap(a3, b3, basef_ctx);
    gr_swap(a4, b4, basef_ctx);
}

static int
_gr_quaternion_randtest(gr_ptr x, flint_rand_t state, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_randtest(a1, state, basef_ctx);
    status |= gr_randtest(a2, state, basef_ctx);
    status |= gr_randtest(a3, state, basef_ctx);
    status |= gr_randtest(a4, state, basef_ctx);

    return status;
}

static int
_gr_quaternion_zero(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_zero(a1, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_one(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_one(a1, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_neg_one(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_neg_one(a1, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_i(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;
    
    status |= gr_zero(a1, basef_ctx);
    status |= gr_one(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_j(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;
    
    status |= gr_zero(a1, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_one(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_k(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;
    
    status |= gr_zero(a1, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_one(a4, basef_ctx);
    return status;
}

static int _gr_quaternion_set_real(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(res, ctx), a2 = IM(res, ctx), a3 = JM(res, ctx), a4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(a1, x, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    
    return status;
}

static int
_gr_quaternion_gens(gr_vec_t vec, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong num_gens;

    if (GR_QUATERNION_FLAGS(ctx) & QUATERNION_K_AS_GEN)
        num_gens = 3;
    else
        num_gens = 2;

    gr_vec_set_length(vec, num_gens, ctx);

    status |= _gr_quaternion_i(gr_vec_entry_ptr(vec, 0, ctx), ctx);
    status |= _gr_quaternion_j(gr_vec_entry_ptr(vec, 1, ctx), ctx);
    if (num_gens == 3)
        status |= _gr_quaternion_k(gr_vec_entry_ptr(vec, 2, ctx), ctx);

    return status;
}

static int
_gr_quaternion_gens_recursive(gr_vec_t vec, gr_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n;
    slong num_gens;

    if (GR_QUATERNION_FLAGS(ctx) & QUATERNION_K_AS_GEN)
        num_gens = 3;
    else
        num_gens = 2;

    gr_vec_init(vec1, 0, BASEF_CTX(ctx));
    status = gr_gens_recursive(vec1, BASEF_CTX(ctx));
    n = vec1->length;

    gr_vec_set_length(vec, n + num_gens, ctx);

    for (i = 0; i < n; i++)
    {
        gr_ptr x = gr_vec_entry_ptr(vec, i, ctx);
        gr_srcptr y = gr_vec_entry_srcptr(vec1, i, BASEF_CTX(ctx));

        status |= _gr_quaternion_set_real(x, y, ctx);
    }

    status |= _gr_quaternion_i(gr_vec_entry_ptr(vec, n, ctx), ctx);
    status |= _gr_quaternion_j(gr_vec_entry_ptr(vec, n+1, ctx), ctx);
    if (num_gens == 3)
        status |= _gr_quaternion_k(gr_vec_entry_ptr(vec, n+2, ctx), ctx);

    gr_vec_clear(vec1, BASEF_CTX(ctx));

    return status;
}

static truth_t
_gr_quaternion_equal(gr_srcptr x, gr_srcptr y, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1, a2, a3, a4, b1, b2, b3, b4;
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

    a1 = RE(x, ctx);
    a2 = IM(x, ctx);
    a3 = JM(x, ctx);
    a4 = KM(x, ctx);
    b1 = RE(y, ctx);
    b2 = IM(y, ctx);
    b3 = JM(y, ctx);
    b4 = KM(y, ctx);

    truth_t eq1, eq2, eq3, eq4;

    eq1 = gr_equal(a1, b1, basef_ctx);
    if (eq1 == T_FALSE)
        return eq1;
    eq2 = gr_equal(a2, b2, basef_ctx);
    if (eq2 == T_FALSE)
        return eq2;
    eq3 = gr_equal(a3, b3, basef_ctx);
    if (eq3 == T_FALSE)
        return eq3;
    eq4 = gr_equal(a4, b4, basef_ctx);
    if (eq4 == T_FALSE)
        return eq4;
    return eq4;    
}

static truth_t
_gr_quaternion_is_zero(gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1, a2, a3, a4;
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

    a1 = RE(x, ctx);
    a2 = IM(x, ctx);
    a3 = JM(x, ctx);
    a4 = KM(x, ctx);

    truth_t eq1, eq2, eq3, eq4;
    eq1 = gr_is_zero(a1, basef_ctx);
    if (eq1 == T_FALSE)
        return eq1;

    eq2 = gr_is_zero(a2, basef_ctx);
    if (eq2 == T_FALSE)
        return eq2;

    eq3 = gr_is_zero(a3, basef_ctx);
    if (eq3 == T_FALSE)
        return eq3;

    eq4 = gr_is_zero(a4, basef_ctx);
    if (eq4 == T_FALSE)
        return eq4;

    return eq4;
}

static truth_t
_gr_quaternion_is_one(gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1, a2, a3, a4;
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

    a1 = RE(x, ctx);
    a2 = IM(x, ctx);
    a3 = JM(x, ctx);
    a4 = KM(x, ctx);

    truth_t eq1, eq2, eq3, eq4;
    eq1 = gr_is_one(a1, basef_ctx);
    if (eq1 == T_FALSE)
        return eq1;

    eq2 = gr_is_zero(a2, basef_ctx);
    if (eq2 == T_FALSE)
        return eq2;

    eq3 = gr_is_zero(a3, basef_ctx);
    if (eq3 == T_FALSE)
        return eq3;

    eq4 = gr_is_zero(a4, basef_ctx);
    if (eq4 == T_FALSE)
        return eq4;

    return eq4;
}

static truth_t
_gr_quaternion_is_neg_one(gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1, a2, a3, a4;
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);

    a1 = RE(x, ctx);
    a2 = IM(x, ctx);
    a3 = JM(x, ctx);
    a4 = KM(x, ctx);

    truth_t eq1, eq2, eq3, eq4;
    eq1 = gr_is_neg_one(a1, basef_ctx);
    if (eq1 == T_FALSE)
        return eq1;

    eq2 = gr_is_zero(a2, basef_ctx);
    if (eq2 == T_FALSE)
        return eq2;

    eq3 = gr_is_zero(a3, basef_ctx);
    if (eq3 == T_FALSE)
        return eq3;

    eq4 = gr_is_zero(a4, basef_ctx);
    if (eq4 == T_FALSE)
        return eq4;

    return eq4;
}

static int
_gr_quaternion_neg(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ptr b1 = RE(res, ctx), b2 = IM(res, ctx), b3 = JM(res, ctx), b4 = KM(res, ctx);

    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_neg(b1, a1, basef_ctx);
    status |= gr_neg(b2, a2, basef_ctx);
    status |= gr_neg(b3, a3, basef_ctx);
    status |= gr_neg(b4, a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_set(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ptr b1 = RE(res, ctx), b2 = IM(res, ctx), b3 = JM(res, ctx), b4 = KM(res, ctx);

    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(b1, a1, basef_ctx);
    status |= gr_set(b2, a2, basef_ctx);
    status |= gr_set(b3, a3, basef_ctx);
    status |= gr_set(b4, a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_set_si(gr_ptr res, slong x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(res, ctx), a2 = IM(res, ctx), a3 = JM(res, ctx), a4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_si(a1, x, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_set_ui(gr_ptr res, ulong x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(res, ctx), a2 = IM(res, ctx), a3 = JM(res, ctx), a4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_ui(a1, x, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_set_fmpz(gr_ptr res, const fmpz_t x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(res, ctx), a2 = IM(res, ctx), a3 = JM(res, ctx), a4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set_fmpz(a1, x, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    return status;
}

static int
_gr_quaternion_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_srcptr b1 = RE(y, ctx), b2 = IM(y, ctx), b3 = JM(y, ctx), b4 = KM(y, ctx);
    gr_ptr c1 = RE(res, ctx), c2 = IM(res, ctx), c3 = JM(res, ctx), c4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_add(c1, a1, b1, basef_ctx);
    status |= gr_add(c2, a2, b2, basef_ctx);
    status |= gr_add(c3, a3, b3, basef_ctx);
    status |= gr_add(c4, a4, b4, basef_ctx);

    return status;
}

static int
_gr_quaternion_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_srcptr b1 = RE(y, ctx), b2 = IM(y, ctx), b3 = JM(y, ctx), b4 = KM(y, ctx);
    gr_ptr c1 = RE(res, ctx), c2 = IM(res, ctx), c3 = JM(res, ctx), c4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_sub(c1, a1, b1, basef_ctx);
    status |= gr_sub(c2, a2, b2, basef_ctx);
    status |= gr_sub(c3, a3, b3, basef_ctx);
    status |= gr_sub(c4, a4, b4, basef_ctx);

    return status;
}

static int
_gr_quaternion_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), b1 = IM(x, ctx), c1 = JM(x, ctx), d1 = KM(x, ctx);
    gr_srcptr a2 = RE(y, ctx), b2 = IM(y, ctx), c2 = JM(y, ctx), d2 = KM(y, ctx);
    gr_ptr a3 = RE(res, ctx), b3 = IM(res, ctx), c3 = JM(res, ctx), d3 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    gr_ptr a1a2, b1b2, c1c2, d1d2;
    GR_TMP_INIT4(a1a2, b1b2, c1c2, d1d2, basef_ctx);
    gr_ptr a1b2, b1a2, c1d2, d1c2;
    GR_TMP_INIT4(a1b2, b1a2, c1d2, d1c2, basef_ctx);
    gr_ptr a1c2, b1d2, c1a2, d1b2;
    GR_TMP_INIT4(a1c2, b1d2, c1a2, d1b2, basef_ctx);
    gr_ptr a1d2, b1c2, c1b2, d1a2;
    GR_TMP_INIT4(a1d2, b1c2, c1b2, d1a2, basef_ctx);

    gr_ptr re_st1, re_st2;
    GR_TMP_INIT2(re_st1, re_st2, basef_ctx);
    gr_ptr im_st1, im_st2;
    GR_TMP_INIT2(im_st1, im_st2, basef_ctx);
    gr_ptr jm_st1, jm_st2;
    GR_TMP_INIT2(jm_st1, jm_st2, basef_ctx);
    gr_ptr km_st1, km_st2;
    GR_TMP_INIT2(km_st1, km_st2, basef_ctx);

    status |= gr_mul(a1a2, a1, a2, basef_ctx);   
    status |= gr_mul(b1b2, b1, b2, basef_ctx);   
    status |= gr_mul(c1c2, c1, c2, basef_ctx);   
    status |= gr_mul(d1d2, d1, d2, basef_ctx);

    status |= gr_mul(a1b2, a1, b2, basef_ctx);   
    status |= gr_mul(b1a2, b1, a2, basef_ctx);   
    status |= gr_mul(c1d2, c1, d2, basef_ctx);   
    status |= gr_mul(d1c2, d1, c2, basef_ctx);

    status |= gr_mul(a1c2, a1, c2, basef_ctx);   
    status |= gr_mul(b1d2, b1, d2, basef_ctx);   
    status |= gr_mul(c1a2, c1, a2, basef_ctx);   
    status |= gr_mul(d1b2, d1, b2, basef_ctx);

    status |= gr_mul(a1d2, a1, d2, basef_ctx);   
    status |= gr_mul(b1c2, b1, c2, basef_ctx);   
    status |= gr_mul(c1b2, c1, b2, basef_ctx);   
    status |= gr_mul(d1a2, d1, a2, basef_ctx);
    
    status |= gr_sub(re_st1, a1a2, b1b2, basef_ctx);
    status |= gr_sub(re_st2, re_st1, c1c2, basef_ctx);
    status |= gr_sub(a3, re_st2, d1d2, basef_ctx);

    status |= gr_add(im_st1, a1b2, b1a2, basef_ctx);
    status |= gr_add(im_st2, im_st1, c1d2, basef_ctx);
    status |= gr_sub(b3, im_st2, d1c2, basef_ctx);

    status |= gr_sub(jm_st1, a1c2, b1d2, basef_ctx);
    status |= gr_add(jm_st2, jm_st1, c1a2, basef_ctx);
    status |= gr_add(c3, jm_st2, d1b2, basef_ctx);

    status |= gr_add(km_st1, a1d2, b1c2, basef_ctx);
    status |= gr_sub(km_st2, km_st1, c1b2, basef_ctx);
    status |= gr_add(d3, km_st2, d1a2, basef_ctx);

    GR_TMP_CLEAR4(a1a2, b1b2, c1c2, d1d2, basef_ctx);
    GR_TMP_CLEAR4(a1b2, b1a2, c1d2, d1c2, basef_ctx);
    GR_TMP_CLEAR4(a1c2, b1d2, c1a2, d1b2, basef_ctx);
    GR_TMP_CLEAR4(a1d2, b1c2, c1b2, d1a2, basef_ctx);
    GR_TMP_CLEAR2(re_st1, re_st2, basef_ctx);
    GR_TMP_CLEAR2(im_st1, im_st2, basef_ctx);
    GR_TMP_CLEAR2(jm_st1, jm_st2, basef_ctx);
    GR_TMP_CLEAR2(km_st1, km_st2, basef_ctx);

    return status;
}

static int
_gr_quaternion_inv(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ptr b1 = RE(res, ctx), b2 = IM(res, ctx), b3 = JM(res, ctx), b4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    gr_ptr s1, s2, s3, s4, m_st1, m_st2, m;
    GR_TMP_INIT4(s1, s2, s3, s4, basef_ctx);
    GR_TMP_INIT3(m_st1, m_st2, m, basef_ctx);

    status |= gr_sqr(s1, a1, basef_ctx);
    status |= gr_sqr(s2, a2, basef_ctx);
    status |= gr_sqr(s3, a3, basef_ctx);
    status |= gr_sqr(s4, a4, basef_ctx);
    status |= gr_add(m_st1, s1, s2, basef_ctx);
    status |= gr_add(m_st2, m_st1, s3, basef_ctx);
    status |= gr_add(m, m_st2, s4, basef_ctx);

    status |= _gr_vec_div_scalar(b1, x, 4, m, basef_ctx);
    status |= gr_neg(b2, b2, basef_ctx);
    status |= gr_neg(b3, b3, basef_ctx);
    status |= gr_neg(b4, b4, basef_ctx);

    GR_TMP_CLEAR4(s1, s2, s3, s4, basef_ctx);
    GR_TMP_CLEAR3(m_st1, m_st2, m, basef_ctx);

    return status;
}

static int
_gr_quaternion_abs(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_srcptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ptr b1 = RE(res, ctx), b2 = IM(res, ctx), b3 = JM(res, ctx), b4 = KM(res, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    gr_ptr s1, s2, s3, s4, m_st1, m_st2, m;
    GR_TMP_INIT4(s1, s2, s3, s4, basef_ctx);
    GR_TMP_INIT3(m_st1, m_st2, m, basef_ctx);

    status |= gr_sqr(s1, a1, basef_ctx);
    status |= gr_sqr(s2, a2, basef_ctx);
    status |= gr_sqr(s3, a3, basef_ctx);
    status |= gr_sqr(s4, a4, basef_ctx);
    status |= gr_add(m_st1, s1, s2, basef_ctx);
    status |= gr_add(m_st2, m_st1, s3, basef_ctx);
    status |= gr_add(m, m_st2, s4, basef_ctx);

    status |= gr_sqrt(b1, m, basef_ctx);
    status |= gr_zero(b2, basef_ctx);
    status |= gr_zero(b3, basef_ctx);
    status |= gr_zero(b4, basef_ctx);

    GR_TMP_CLEAR4(s1, s2, s3, s4, basef_ctx);
    GR_TMP_CLEAR3(m_st1, m_st2, m, basef_ctx);

    return status;
}

static int
_gr_quaternion_re(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(RE(res, ctx), RE(x, ctx), basef_ctx);
    status |= gr_zero(IM(res, ctx), basef_ctx);
    status |= gr_zero(JM(res, ctx), basef_ctx);
    status |= gr_zero(KM(res, ctx), basef_ctx);
    return status;
}

static int
_gr_quaternion_im(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(IM(res, ctx), IM(x, ctx), basef_ctx);
    status |= gr_zero(RE(res, ctx), basef_ctx);
    status |= gr_zero(JM(res, ctx), basef_ctx);
    status |= gr_zero(KM(res, ctx), basef_ctx);
    return status;
}

static int
_gr_quaternion_jm(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(JM(res, ctx), JM(x, ctx), basef_ctx);
    status |= gr_zero(RE(res, ctx), basef_ctx);
    status |= gr_zero(IM(res, ctx), basef_ctx);
    status |= gr_zero(KM(res, ctx), basef_ctx);
    return status;
}

static int
_gr_quaternion_km(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(KM(res, ctx), KM(x, ctx), basef_ctx);
    status |= gr_zero(RE(res, ctx), basef_ctx);
    status |= gr_zero(IM(res, ctx), basef_ctx);
    status |= gr_zero(JM(res, ctx), basef_ctx);
    return status;
}

static int
_gr_quaternion_conj(gr_ptr res, gr_srcptr x, gr_quaternion_ctx_t ctx)
{
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_set(RE(res, ctx), RE(x, ctx), basef_ctx);
    status |= gr_neg(IM(res, ctx), IM(x, ctx), basef_ctx);
    status |= gr_neg(JM(res, ctx), IM(x, ctx), basef_ctx);
    status |= gr_neg(KM(res, ctx), IM(x, ctx), basef_ctx);

    return status;
}

static int
_gr_quaternion_pi(gr_ptr x, gr_quaternion_ctx_t ctx)
{
    gr_ptr a1 = RE(x, ctx), a2 = IM(x, ctx), a3 = JM(x, ctx), a4 = KM(x, ctx);
    gr_ctx_struct * basef_ctx = BASEF_CTX(ctx);
    int status = GR_SUCCESS;

    status |= gr_pi(a1, basef_ctx);
    status |= gr_zero(a2, basef_ctx);
    status |= gr_zero(a3, basef_ctx);
    status |= gr_zero(a4, basef_ctx);
    
    return status;
}

int _gr_quaternion_methods_initialized = 0;

gr_static_method_table _gr_quaternion_methods;

gr_method_tab_input _gr_quaternion_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_quaternion_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_quaternion_ctx_clear},

    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) _gr_quaternion_ctx_is_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr)  _gr_quaternion_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr)  _gr_quaternion_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr)  _gr_quaternion_ctx_is_field},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _gr_quaternion_ctx_is_rational_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _gr_quaternion_ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) _gr_quaternion_ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr)  _gr_quaternion_ctx_is_threadsafe},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) _gr_quaternion_ctx_is_finite},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) _gr_quaternion_ctx_is_finite_characteristic},

    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) _gr_quaternion_ctx_is_exact},
    {GR_METHOD_CTX_NGENS,       (gr_funcptr) _gr_quaternion_ctx_ngens},
    {GR_METHOD_CTX_GEN_NAME,    (gr_funcptr) _gr_quaternion_ctx_gen_name},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_quaternion_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_quaternion_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_quaternion_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_quaternion_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_quaternion_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_quaternion_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_quaternion_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_quaternion_one},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) _gr_quaternion_neg_one},
    {GR_METHOD_GENS,           (gr_funcptr) _gr_quaternion_gens},
    {GR_METHOD_GENS_RECURSIVE,  (gr_funcptr) _gr_quaternion_gens_recursive},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_quaternion_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_quaternion_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_quaternion_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_quaternion_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_quaternion_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_quaternion_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_quaternion_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_quaternion_set_fmpz},
    /*{GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_quaternion_set_other},*/
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_quaternion_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_quaternion_add},
/*
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_quaternion_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_quaternion_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_quaternion_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_quaternion_add_fmpq},
*/
    {GR_METHOD_SUB,             (gr_funcptr) _gr_quaternion_sub},
/*
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_quaternion_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_quaternion_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_quaternion_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_quaternion_sub_fmpq},
*/
    {GR_METHOD_MUL,             (gr_funcptr) _gr_quaternion_mul},
/*
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_quaternion_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_quaternion_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_quaternion_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_quaternion_mul_fmpq},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_quaternion_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_quaternion_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_quaternion_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_quaternion_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_quaternion_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_quaternion_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_quaternion_div_fmpq},
*/
/*
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_quaternion_is_invertible},
*/
    {GR_METHOD_INV,             (gr_funcptr) _gr_quaternion_inv},
/*
    {GR_METHOD_SQRT,             (gr_funcptr) _gr_quaternion_sqrt},
    {GR_METHOD_NUMERATOR,       (gr_funcptr) _gr_quaternion_numerator},
    {GR_METHOD_DENOMINATOR,     (gr_funcptr) _gr_quaternion_denominator},
*/

    {GR_METHOD_ABS,               (gr_funcptr) _gr_quaternion_abs},
    {GR_METHOD_RE,               (gr_funcptr) _gr_quaternion_re},
    {GR_METHOD_IM,               (gr_funcptr) _gr_quaternion_im},
    {GR_METHOD_CONJ,             (gr_funcptr) _gr_quaternion_conj},

    {GR_METHOD_I,               (gr_funcptr) _gr_quaternion_i},
    {GR_METHOD_PI,              (gr_funcptr) _gr_quaternion_pi},

    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_gr_quaternion(gr_ctx_t ctx, gr_ctx_t basef_ctx, int flags)
{
    ctx->which_ring = GR_CTX_GR_QUATERNION;
    ctx->sizeof_elem = 4 * basef_ctx->sizeof_elem;
    ctx->size_limit = WORD_MAX;

    GR_QUATERNION_BASEF_CTX(ctx) = basef_ctx;
    GR_QUATERNION_FLAGS(ctx) = flags;

    ctx->methods = _gr_quaternion_methods;

    if (!_gr_quaternion_methods_initialized)
    {
        gr_method_tab_init(_gr_quaternion_methods, _gr_quaternion_methods_input);
        _gr_quaternion_methods_initialized = 1;
    }
}


