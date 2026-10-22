/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_CA_H
#define GR_CA_H

/*
    Castle arithmetic: exact real and complex numbers as elements of
    towers of field extensions of Q (gr_tower.h), built on a common
    ground field and evaluated dynamically. A castle field is the
    complex field of the numbers the towers can represent, or its real
    or algebraic subfield, or a view of one. Elements carry a pointer to
    a shared tower which is created and extended on demand, and their
    data is converted and reduced lazily; the zero test is complete
    (under Schanuel's conjecture for termination). Castle arithmetic
    succeeds the Calcium design (ca.h), from which it inherits the
    "ca" mnemonic.
*/

#include "fmpq_poly.h"
#include "fmpz_mpoly_q.h"
#include "dirichlet.h"
#include "gr_tower.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Elements of the castle fields. The fields may be read (they are
    documented with the module) but not written: an element refers to
    a tower shared with other elements, and its data is converted and
    reduced lazily under the context's lock. The accessors below bring
    an element up to date before exposing its data. The layout may
    change between versions of FLINT.
*/
#define GR_CA_REPR_FLAT 0       /* elem.flat: fmpz_mpoly_q in the polynomial context elem.flat.mctx */
#define GR_CA_REPR_RATIONAL 1   /* elem.q: a rational number (element of the trivial tower) */
#define GR_CA_REPR_DENSE 2      /* elem.dense: a polynomial in the first generator of the tower */
#define GR_CA_REPR_FACTORED 3   /* elem.factored: fmpz_mpoly_q_factored (a factored denominator) in elem.factored.mctx */

struct _gr_tower_dense_field_struct;

typedef struct
{
    gr_tower_flat_struct * F;       /* the tower (via F->T) and its flat machinery */
    slong level;                    /* prefix length, valid for structure_version of the tower */
    ulong structure_version;
    ulong reduced_version;          /* ideal version of F for which data is known to be reduced (0: unknown) */
    slong shallow;                  /* the data may be shared with another element (internal) */
    slong repr;                     /* GR_CA_REPR_* */
#ifdef CA_UNION_DEBUG
    struct
#else
    union
#endif
    {
        fmpq q;                                    /* RATIONAL */
        struct
        {
            fmpz_mpoly_q_struct data;
            fmpz_mpoly_ctx_struct * mctx;          /* context in which data is stored */
        }
        flat;                                      /* FLAT */
        struct
        {
            fmpq_poly_struct poly;                 /* canonical, of length < the degree of the field */
            const struct _gr_tower_dense_field_struct * nf;   /* the field: generator and modulus (internal) */
        }
        dense;                                     /* DENSE */
        struct
        {
            fmpz_mpoly_q_factored_struct data;
            fmpz_mpoly_ctx_struct * mctx;          /* context in which data is stored */
        }
        factored;                                  /* FACTORED */
    }
    elem;
}
gr_ca_struct;

typedef gr_ca_struct gr_ca_t[1];

/* element printing (the GR_TOWER_PRINT_* flags of gr_tower.h) */
void gr_ca_ctx_set_print(gr_ctx_t ctx, int flags, slong digits);

/* choice of generators in castle fields (gr_ca_ctx_set_gen_flags):
   by default, the square root of a negative rational number -A B^2 (A
   squarefree) is B sqrt(-A) with sqrt(-A) a single generator (a root of
   x^2 + A; sqrt(-1) = i and sqrt(-3) = 2 zeta_3 + 1 are roots of unity);
   with GR_TOWER_GENS_SPLIT_IMAGINARY it is i B sqrt(A), with sqrt(A) the
   product of the roots of the primes dividing A. Roots of unity of
   composite orders are products of roots of prime power orders, and
   square roots of positive rationals products of square roots of primes,
   unless GR_TOWER_GENS_COMPOSITE_ROOTS (zeta_N as one generator, a root
   of Phi_N, of the least common multiple of the orders requested, as in
   Calcium) or GR_TOWER_GENS_COMPOSITE_RADICALS (sqrt(A) with A the
   squarefree part as one generator) is set */
#define GR_TOWER_GENS_SPLIT_IMAGINARY 1
#define GR_TOWER_GENS_COMPOSITE_ROOTS 2
#define GR_TOWER_GENS_COMPOSITE_RADICALS 4
#define GR_TOWER_COMPOSITE_ROOTS_DEGREE 128   /* the cyclotomic degree limit set by GR_TOWER_GENS_COMPOSITE_ROOTS when none is set */
void gr_ca_ctx_set_gen_flags(gr_ctx_t ctx, int flags);
/* tuning options (GR_TOWER_OPT_*) of a castle field, shared by its towers
   and read whenever they apply (a change takes effect at once, also for
   existing elements); GR_DOMAIN for a value out of range */
int gr_ca_ctx_set_option(gr_ctx_t ctx, slong option, slong value);
slong gr_ca_ctx_get_option(gr_ctx_t ctx, slong option);
int gr_ca_ctx_gen_flags(gr_ctx_t ctx);
int gr_ca_ctx_print_flags(gr_ctx_t ctx);
int gr_ca_ctx_field_flags(gr_ctx_t ctx);

/* Castle field of algebraic and transcendental numbers over the base field
   (the rational numbers): elements
   carry a pointer to a shared tower, created and extended on demand.
   flags: GR_CA_REAL, GR_CA_ALGEBRAIC (subfields), and
   GR_CA_UNLIMITED: no operation fails because a precision,
   degree or number of tries was exhausted -- an operation which would
   fail for that reason is repeated with all such effort limits doubled
   (the tuning options are their starting values), so that it may take
   arbitrarily long, and does not terminate when no method decides it
   (an unknown identity between transcendental constants). */
void gr_ctx_init_ca(gr_ctx_t ctx, gr_ctx_t base, int flags);
/* A view of the castle context parent sharing its towers and elements,
   restricted to a subfield (GR_CA_REAL, GR_CA_ALGEBRAIC
   added to the restrictions of parent); parent may be cleared first. */
void gr_ctx_init_ca_view(gr_ctx_t ctx, gr_ctx_t parent, int field_flags);
int gr_ca_ctx_same_state(gr_ctx_t ctx1, gr_ctx_t ctx2);

/* The methods of the gr interface of the castle fields, callable directly
   (the element functions of gr.h dispatch to them). Rational elements,
   and dense elements of the same field, take lock-free paths; other
   operations hold the context's lock. */
void gr_ca_init(gr_ca_t x, gr_ctx_t ctx);
void gr_ca_clear(gr_ca_t x, gr_ctx_t ctx);
void gr_ca_swap(gr_ca_t x, gr_ca_t y, gr_ctx_t ctx);
void gr_ca_set_shallow(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_randtest(gr_ca_t res, flint_rand_t state, gr_ctx_t ctx);
int gr_ca_write(gr_stream_t out, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_set_str(gr_ca_t res, const char * s, gr_ctx_t ctx);
#ifdef FEXPR_H
int gr_ca_get_fexpr(fexpr_t res, const gr_ca_t x, gr_ctx_t ctx);
#endif
int gr_ca_zero(gr_ca_t res, gr_ctx_t ctx);
int gr_ca_one(gr_ca_t res, gr_ctx_t ctx);
int gr_ca_pi(gr_ca_t res, gr_ctx_t ctx);
int gr_ca_i(gr_ca_t res, gr_ctx_t ctx);
int gr_ca_gens(gr_vec_t vec, gr_ctx_t ctx);
int gr_ca_set(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_set_si(gr_ca_t res, slong c, gr_ctx_t ctx);
int gr_ca_set_ui(gr_ca_t res, ulong c, gr_ctx_t ctx);
int gr_ca_set_fmpz(gr_ca_t res, const fmpz_t c, gr_ctx_t ctx);
int gr_ca_set_fmpq(gr_ca_t res, const fmpq_t c, gr_ctx_t ctx);
int gr_ca_set_other(gr_ca_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
int gr_ca_get_si(slong * res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_get_ui(ulong * res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_get_fmpz(fmpz_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_get_fmpq(fmpq_t res, const gr_ca_t x, gr_ctx_t ctx);
truth_t gr_ca_is_zero(const gr_ca_t x, gr_ctx_t ctx);
truth_t gr_ca_is_one(const gr_ca_t x, gr_ctx_t ctx);
truth_t gr_ca_equal(const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_neg(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_add(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_sub(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_mul(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_div(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_inv(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_add_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx);
int gr_ca_add_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx);
int gr_ca_add_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_ca_add_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_ca_sub_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx);
int gr_ca_sub_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx);
int gr_ca_sub_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_ca_sub_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_ca_mul_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx);
int gr_ca_mul_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx);
int gr_ca_mul_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_ca_mul_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_ca_div_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx);
int gr_ca_div_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx);
int gr_ca_div_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx);
int gr_ca_div_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx);
int gr_ca_pow(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_sqrt(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_abs(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_conj(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_exp(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_log(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_poly_mullow(gr_ca_struct * res, const gr_ca_struct * p1, slong len1, const gr_ca_struct * p2, slong len2, slong n, gr_ctx_t ctx);
int gr_ca_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);
int gr_ca_poly_factor(gr_poly_t c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);
int gr_ca_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int gr_ca_mat_det(gr_ca_t res, const gr_mat_t A, gr_ctx_t ctx);
int gr_ca_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);

int gr_ca_get_qqbar(qqbar_t res, const gr_ca_t x, gr_ctx_t ctx);

/* Elementary functions of the castle field, through exp, log and sqrt;
   inverse functions are checked against the principal branch
   numerically (GR_UNABLE on the branch cuts). Rounding and comparisons
   apply to real elements (GR_DOMAIN otherwise). */
int gr_ca_sin(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_cos(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_tan(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_sinh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_cosh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_tanh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_asin(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_acos(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_atan(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_asinh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_acosh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_atanh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_re(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_im(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_arg(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_sgn(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_csgn(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
truth_t gr_ca_is_real(const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_cmp(int * res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_cmpabs(int * res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_floor(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_ceil(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_trunc(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_nint(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_get_d(double * res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_root_ui(gr_ca_t res, const gr_ca_t x, ulong n, gr_ctx_t ctx);
int gr_ca_get_acb(acb_t res, const gr_ca_t x, slong prec, gr_ctx_t ctx);

/* Special functions of the castle field (see special.c): values at
   canonical arguments become generators, the functional equations and
   special values being applied first. */
int gr_ca_gamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_rgamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_lgamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_dirichlet_l(gr_ca_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_ca_t s, gr_ctx_t ctx);
int gr_ca_lerch_phi(gr_ca_t res, const gr_ca_t z, const gr_ca_t s, const gr_ca_t a, gr_ctx_t ctx);
int gr_ca_beta(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx);
int gr_ca_digamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_polygamma(gr_ca_t res, const gr_ca_t s, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_erf(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_erfc(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_erfi(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_lambertw(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_lambertw_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t k, gr_ctx_t ctx);
int gr_ca_zeta(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_hurwitz_zeta(gr_ca_t res, const gr_ca_t s, const gr_ca_t a, gr_ctx_t ctx);
int gr_ca_polylog(gr_ca_t res, const gr_ca_t s, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_dilog(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_elliptic_k(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_modular_lambda(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_modular_j(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_modular_delta(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_dedekind_eta(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_eisenstein_e(gr_ca_t res, ulong k, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_eisenstein_g(gr_ca_t res, ulong k, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_jacobi_theta(gr_ca_t res1, gr_ca_t res2, gr_ca_t res3, gr_ca_t res4, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_jacobi_theta_1(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_jacobi_theta_2(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_jacobi_theta_3(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_jacobi_theta_4(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_weierstrass_p(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_weierstrass_p_prime(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_weierstrass_sigma(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_elliptic_roots(gr_ca_t e1, gr_ca_t e2, gr_ca_t e3, const gr_ca_t tau, gr_ctx_t ctx);
int gr_ca_elliptic_invariants(gr_ca_t g2, gr_ca_t g3, const gr_ca_t tau, gr_ctx_t ctx);

/* the table of hypergeometric identities (for the tests) */
slong _gr_tower_hypgeom_table_length(void);
const char * _gr_tower_hypgeom_table_entry(slong i);
int _gr_tower_hypgeom_table_shape(slong * p, slong * q, int * used, int * zfree, fmpq_t z0, slong i);
int _gr_tower_hypgeom_table_eval(gr_ptr params, int * cond, gr_ptr value, slong i, gr_srcptr * vals, gr_ctx_t ctx);
int gr_ca_hypgeom_pfq_vec(gr_ca_t res, const gr_ca_struct * a, slong p, const gr_ca_struct * b, slong q, const gr_ca_t z, int flags, gr_ctx_t ctx);
int gr_ca_hypgeom_pfq(gr_ca_t res, const gr_vec_t a, const gr_vec_t b, const gr_ca_t z, int flags, gr_ctx_t ctx);
int gr_ca_hypgeom_0f1(gr_ca_t res, const gr_ca_t b, const gr_ca_t z, int flags, gr_ctx_t ctx);
int gr_ca_hypgeom_1f1(gr_ca_t res, const gr_ca_t a, const gr_ca_t b, const gr_ca_t z, int flags, gr_ctx_t ctx);
int gr_ca_hypgeom_2f1(gr_ca_t res, const gr_ca_t a, const gr_ca_t b, const gr_ca_t c, const gr_ca_t z, int flags, gr_ctx_t ctx);
int gr_ca_elliptic_e(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_euler(gr_ca_t res, gr_ctx_t ctx);
int gr_ca_catalan(gr_ca_t res, gr_ctx_t ctx);
void gr_ca_ctx_stats(gr_ctx_t ctx);
slong gr_ca_ctx_num_towers(gr_ctx_t ctx);
gr_tower_struct * gr_ca_get_tower(slong * level, const gr_ca_t x, gr_ctx_t ctx);
const fmpz_mpoly_q_struct * gr_ca_get_data(const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_repr(const gr_ca_t x, gr_ctx_t ctx);
int gr_ca_get_fmpq_poly(fmpq_poly_t res, fmpz_poly_t modulus, const gr_ca_t x, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
