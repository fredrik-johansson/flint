/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"
#include "gr_ca/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/*
    Factored denominators (GR_TOWER_OPT_FACTORED_DENOMINATORS).

    In a tower without algebraic generators -- a field Q(t_1, ..., t_r)
    of rational functions of transcendental generators (pi, e,
    logarithms, values of special functions) -- the flat elements need no reduction, and their
    arithmetic is that of fmpz_mpoly_q, dominated by the polynomial gcds
    which keep the fractions canonical. With the option, the arithmetic
    of such elements uses fmpz_mpoly_q_factored instead: the denominator
    is kept as a product of irreducible factors, and the gcds become
    trial divisions by the known factors. A flat operand enters the
    factored form at its first operation (its denominator is factored
    once), and the results stay factored; every other routine reads
    elements through their flat view, which multiplies the denominator
    out. Results which are rational numbers move to the trivial tower as
    usual.

    The arithmetic below applies when the operands are rational numbers
    or elements of one tower F without algebraic generators (an
    operation which would merge towers goes through the flat
    representation); _gr_ca_factored_binary returns -1 when it
    does not apply.

    The option is off by default. Measured against the flat arithmetic
    (towers over Q(pi), Q(pi, e), Q(pi, e, log 2)): the factored form
    wins in some linear algebra over these fields (Berkowitz
    determinants of 5 x 5 matrices of reciprocals of linear forms in
    three generators, 1.9 times faster; LU solving, 2.5 times), but
    loses in sums of many terms with distinct denominators (sum of
    1/(pi + k), k <= 200: 15 times slower), in iterated rational maps,
    whose divisions factor the growing numerators of the divisors
    (x -> (x^2 + pi)/(x + e): 14 times slower), and in other
    determinants (8 x 8 reciprocals in two generators: 1.3 times
    slower). Neither cost guards (the sizes of the polynomials to
    factor, the number of factors of a sum) nor restricting the option
    to the linear algebra routines gave a combination which avoids the
    losses while keeping the gains. The catalog of src/python/
    gr_tower_profile.py, where about half of the cases perform some
    arithmetic of this kind (mostly multiples of pi before an
    exponential), takes the same time either way.
*/

int
_gr_ca_factored_tower_ok(gr_tower_flat_struct * F, gr_ctx_t ctx)
{
    return CA_CTX(ctx)->options[GR_TOWER_OPT_FACTORED_DENOMINATORS] != 0 &&
           F != CA_CTX(ctx)->trivial && F->T->length == 0 && F->T->num_trans > 0;
}

/* the flat data of x (FACTORED) in the current context of F */
void
_gr_ca_factored_get_flat(fmpz_mpoly_q_t res, const gr_ca_t x, gr_tower_flat_struct * F)
{
    const fmpz_mpoly_ctx_struct * old = x->elem.factored.mctx;
    if (old == F->mctx)
        fmpz_mpoly_q_factored_get_fmpz_mpoly_q(res, &x->elem.factored.data, F->mctx);
    else
    {
        fmpz_mpoly_q_t t;
        fmpz_mpoly_q_init(t, old);
        fmpz_mpoly_q_factored_get_fmpz_mpoly_q(t, &x->elem.factored.data, old);
        gr_tower_flat_convert(res, t, old, F);
        fmpz_mpoly_q_clear(t, old);
    }
}

/* the level (prefix length) of the factored data r of F: that of its
   numerator and of its denominator factors */
static slong
_factored_level(const fmpz_mpoly_q_factored_t r, gr_tower_flat_struct * F)
{
    const fmpz_mpoly_factor_struct * d = fmpz_mpoly_q_factored_denref(r);
    fmpz_mpoly_q_struct t;
    slong i, level;

    /* (shallow fractions: a numerator, a denominator factor) */
    t.num = *fmpz_mpoly_q_factored_numref(r);
    t.den.length = 0;
    level = gr_tower_flat_level(&t, F);
    t.num.length = 0;
    for (i = 0; i < d->num; i++)
    {
        t.den = d->poly[i];
        level = FLINT_MAX(level, gr_tower_flat_level(&t, F));
    }
    return level;
}

/*
    The operand x as factored data in the current context of F (x is a
    rational number or an element of F): *ptr points to x's own data when
    it is factored and up to date, otherwise to tmp, which is then
    initialized (*owned = 1) and must be cleared by the caller.
*/
static void
_operand(const fmpz_mpoly_q_factored_struct ** ptr, fmpz_mpoly_q_factored_t tmp, int * owned,
    const gr_ca_t x, gr_tower_flat_struct * F, gr_ctx_t ctx)
{
    gr_tower_flat_ensure(F);
    *owned = 0;

    if (x->repr == CA_REPR_FACTORED && x->elem.factored.mctx == F->mctx)
    {
        *ptr = &x->elem.factored.data;
        return;
    }

    fmpz_mpoly_q_factored_init(tmp, F->mctx);
    *owned = 1;
    *ptr = tmp;

    if (x->repr == CA_REPR_RATIONAL)
        fmpz_mpoly_q_factored_set_fmpq(tmp, &x->elem.q, F->mctx);
    else
    {
        /* (flat, or factored in a superseded context: the flat view, up
           to date, whose denominator is factored here, once) */
        gr_ca_struct * v = _gr_ca_flat_view(x);
        fmpz_mpoly_q_factored_set_fmpz_mpoly_q(tmp, &v->elem.flat.data, F->mctx);
    }
}

/* res = the factored data r (consumed) of F, or a rational number of the
   trivial tower if r is constant */
static void
_install(gr_ca_t res, gr_tower_flat_struct * F, fmpz_mpoly_q_factored_t r, gr_ctx_t ctx)
{
    if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_factored_numref(r), F->mctx) &&
        fmpz_mpoly_q_factored_denref(r)->num == 0)
    {
        fmpq_t q;
        fmpq_init(q);
        if (fmpz_mpoly_is_zero(fmpz_mpoly_q_factored_numref(r), F->mctx))
            fmpq_zero(q);
        else
        {
            const fmpz_mpoly_factor_struct * d = fmpz_mpoly_q_factored_denref(r);
            fmpz_mpoly_get_fmpz(fmpq_numref(q), fmpz_mpoly_q_factored_numref(r), F->mctx);
            fmpz_mul(fmpq_numref(q), fmpq_numref(q), d->constant_den);
            fmpz_set(fmpq_denref(q), d->constant);
            fmpq_canonicalise(q);
        }
        GR_MUST_SUCCEED(_gr_ca_set_fmpq(res, q, ctx));
        fmpq_clear(q);
        fmpz_mpoly_q_factored_clear(r, F->mctx);
        return;
    }

    _gr_ca_fresh(res, F, _factored_level(r, F), ctx);
    fmpz_mpoly_q_clear(&res->elem.flat.data, res->elem.flat.mctx);
    CA_POISON_FLAT(res);
    res->repr = CA_REPR_FACTORED;
    res->elem.factored.mctx = F->mctx;
    res->elem.factored.data = *r;       /* (moved) */
    res->reduced_version = F->ideal_version;
}

/* the tower of the factored arithmetic of x and y, or NULL */
static gr_tower_flat_struct *
_common(const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
{
    gr_tower_flat_struct * T0 = CA_CTX(ctx)->trivial;
    gr_tower_flat_struct * F;

    if (x->F != T0 && y->F != T0 && x->F != y->F)
        return NULL;
    F = (x->F != T0) ? x->F : y->F;
    if (F == T0 || !_gr_ca_factored_tower_ok(F, ctx))
        return NULL;
    /* (rational operands are rational, and a dense form needs algebraic
       generators: the other representations are flat or factored) */
    if ((x->repr != CA_REPR_RATIONAL && x->repr != CA_REPR_FLAT && x->repr != CA_REPR_FACTORED) ||
        (y->repr != CA_REPR_RATIONAL && y->repr != CA_REPR_FLAT && y->repr != CA_REPR_FACTORED))
        return NULL;
    /* (two rational numbers: not here) */
    if (x->repr == CA_REPR_RATIONAL && y->repr == CA_REPR_RATIONAL)
        return NULL;
    return F;
}

/*
    res = x op y (CA_FACTORED_ADD, _SUB, _MUL, _DIV) in factored form,
    or -1 if the factored arithmetic does not apply (then nothing is
    done).
*/
int
_gr_ca_factored_binary(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, int op, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F = _common(x, y, ctx);
    const fmpz_mpoly_q_factored_struct * px, * py;
    fmpz_mpoly_q_factored_t tx, ty, r;
    int ox, oy, status = GR_SUCCESS;

    if (F == NULL)
        return -1;

    _operand(&px, tx, &ox, x, F, ctx);
    _operand(&py, ty, &oy, y, F, ctx);
    fmpz_mpoly_q_factored_init(r, F->mctx);

    if (op == CA_FACTORED_ADD)
        fmpz_mpoly_q_factored_add(r, px, py, F->mctx);
    else if (op == CA_FACTORED_SUB)
        fmpz_mpoly_q_factored_sub(r, px, py, F->mctx);
    else if (op == CA_FACTORED_MUL)
        fmpz_mpoly_q_factored_mul(r, px, py, F->mctx);
    else if (fmpz_mpoly_q_factored_is_zero(py, F->mctx))
        status = GR_DOMAIN;
    else
        fmpz_mpoly_q_factored_div(r, px, py, F->mctx);

    if (ox) fmpz_mpoly_q_factored_clear(tx, F->mctx);
    if (oy) fmpz_mpoly_q_factored_clear(ty, F->mctx);

    if (status == GR_SUCCESS)
        _install(res, F, r, ctx);
    else
        fmpz_mpoly_q_factored_clear(r, F->mctx);
    return status;
}

/* res = 1/x (x nonzero) in factored form, or -1 if it does not apply */
int
_gr_ca_factored_inv(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F = x->F;
    const fmpz_mpoly_q_factored_struct * px;
    fmpz_mpoly_q_factored_t tx, r;
    int ox, status = GR_SUCCESS;

    if (F == CA_CTX(ctx)->trivial || !_gr_ca_factored_tower_ok(F, ctx) ||
        (x->repr != CA_REPR_FLAT && x->repr != CA_REPR_FACTORED))
        return -1;

    _operand(&px, tx, &ox, x, F, ctx);
    fmpz_mpoly_q_factored_init(r, F->mctx);

    if (fmpz_mpoly_q_factored_is_zero(px, F->mctx))
        status = GR_DOMAIN;
    else
        fmpz_mpoly_q_factored_inv(r, px, F->mctx);

    if (ox) fmpz_mpoly_q_factored_clear(tx, F->mctx);

    if (status == GR_SUCCESS)
        _install(res, F, r, ctx);
    else
        fmpz_mpoly_q_factored_clear(r, F->mctx);
    return status;
}

POP_OPTIONS
