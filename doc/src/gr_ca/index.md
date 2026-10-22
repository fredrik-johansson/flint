---
title: gr_ca.h
header: gr_ca.h
description: Castle arithmetic: exact real and complex numbers in towers of extension fields built on demand over the rationals, with complete zero tests for algebraic numbers, relation searches for transcendental ones, and the gr interface.
category: Exact real and complex numbers
order: 0.1
---

*Castle arithmetic* computes exactly with real and complex numbers: a
*castle field* keeps a collection of towers of field extensions
([](gr_tower)) built on the common ground field $\mathbb{Q}$, and every
element is a number of one of the towers, which are created, extended,
refined and merged on demand as numbers are computed. Field arithmetic
is formal arithmetic in a tower; the zero test is complete for algebraic
numbers (dynamic evaluation) and decides numbers involving exponentials
and logarithms by searching for relations between the generators
(complete under Schanuel's conjecture), and the special functions add
their known identities. The design, its guarantees and its relation to
the earlier design of [](ca) are described in
[](guide/types/castle-arithmetic); [](gr_ca/using) shows how to use
it; the algorithms are in [](algorithms/exact).

A castle field is a [](gr) context: elements have type `gr_ca_t` and
are used through the generic functions of `gr.h` (`gr_add`,
`gr_is_zero`, `gr_sqrt`, `gr_exp`, `gr_poly_roots`, `gr_mat_det`, ...),
which dispatch to the functions below. There are four fields: the
complex field of all the numbers the towers can represent, its real
subfield, and the complex and real algebraic numbers. The contexts are
thread-safe.

```c example name=gr_ca_first
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_ca.h>

int main(void)
{
    gr_ctx_t QQ, K;
    gr_ptr x, y;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_ca(K, QQ, 0);
    x = gr_heap_init(K);
    y = gr_heap_init(K);

    /* sqrt(5 + 2 sqrt(6)) - sqrt(2) - sqrt(3) */
    status |= gr_set_str(x, "sqrt(5 + 2*sqrt(6)) - sqrt(2) - sqrt(3)", K);
    flint_printf("zero: %d\n", gr_is_zero(x, K) == T_TRUE);

    /* exp(pi*i/3) */
    status |= gr_set_str(y, "exp(pi*i/3)", K);
    status |= gr_println(y, K);

    gr_heap_clear(x, K);
    gr_heap_clear(y, K);
    gr_ctx_clear(K);
    gr_ctx_clear(QQ);
    flint_cleanup();
    return status;
}
```
```output
zero: 1
a4+1 {a4 = exp(2*pi*i/3)}
```

The first number is zero: $6 = 2 \cdot 3$ makes $\sqrt 6 = \sqrt 2
\sqrt 3$, and the square root of $5 + 2\sqrt 6$ is found in
$\mathbb{Q}(\sqrt 2, \sqrt 3)$ to be $\sqrt 2 + \sqrt 3$. The second
number is printed in terms of the generators of its tower, followed by
their definitions: $e^{\pi i/3} = 1 + \zeta_3$, since the roots of unity
of a castle field are products of powers of roots of unity of prime
power orders (the names $a_1, a_2, \ldots$ are given in order of
creation, so that `a4` follows $\sqrt 2$, $\sqrt 3$ and auxiliary
generators of the computation).

## Contexts {#contexts}

::: function
void gr_ctx_init_ca(gr_ctx_t ctx, gr_ctx_t base, int flags)

Initializes `ctx` to a castle field over `base` (currently the
rational numbers, `gr_ctx_init_fmpq`). `flags` selects the field: the
complex field of all the numbers the towers can represent by default,
its real subfield with `GR_CA_REAL`, and the algebraic numbers with
`GR_CA_ALGEBRAIC`; the four contexts print as "Complex field (castle
arithmetic)", "Real field (castle arithmetic)", "Complex algebraic
field (castle arithmetic)" and "Real algebraic field (castle
arithmetic)". With `GR_CA_UNLIMITED`, no operation fails because a
precision, degree or number of tries was exhausted
([](#unlimited)). `GR_TOWER_MERGE_EXPRESS` enables low-effort
expression searches when towers are merged.

::: details
In a real context, an operation whose result would not be real fails
with `GR_DOMAIN`: $i$, roots and logarithms of negative numbers,
inverse trigonometric functions outside their real domains,
conversions of nonreal numbers; `gr_poly_roots` returns the real roots
only (for rational polynomials, only the real roots are isolated), and
`gr_factor` gives quadratic factors for the pairs of nonreal roots.
Since all elements are real, the tests are sign tests. In an algebraic
context, $\pi$ is not available and the elementary functions only at
the arguments where their values are algebraic ($\exp(0)$,
$\log(1)$, $\cos(0)$, $1^y$, $x^{p/q}$, ..., decided by the
Lindemann–Weierstrass and Gelfond–Schneider theorems), the others
failing with `GR_DOMAIN`. The restrictions apply to the operations
called by the user, not to the intermediate values of the
implementation ($\arccos(1/2) = \pi/3$ in the real field goes through
complex logarithms). The real and algebraic fields can also be created
as views of a complex field sharing its elements
(`gr_ctx_init_ca_view`).

Every operation on a context holds a recursive mutex of the context
(operations may refine or merge the shared towers and update the
cached representation of elements), so operations on one context are
serialized, except for the lock-free arithmetic of rational elements
and of dense elements of the same number field; distinct contexts are
independent. How the towers are managed (sharing, collection, naming
of generators) is described in [](gr_ca/using).
:::
:::

::: macro
GR_CA_REAL
GR_CA_ALGEBRAIC
GR_CA_UNLIMITED

Flags of `gr_ctx_init_ca`: the real subfield, the algebraic subfield,
and the unlimited mode (`GR_TOWER_OPT_UNLIMITED`). The flags share their
space with the merge flags of [](gr_tower) (`GR_TOWER_MERGE_EXPRESS`).
:::

::: function
void gr_ctx_init_ca_view(gr_ctx_t ctx, gr_ctx_t parent, int field_flags)
int gr_ca_ctx_same_state(gr_ctx_t ctx1, gr_ctx_t ctx2)

Initializes `ctx` to a *view* of the castle field `parent`: a context
sharing its towers, registry, generator names and lock (and hence
its elements), restricted further by `field_flags`
(`GR_CA_REAL`, `GR_CA_ALGEBRAIC`, added
to those of `parent`). The shared state is reference counted:
`parent` and its views may be cleared in any order.
`gr_set_other` between contexts of the same state
(`gr_ca_ctx_same_state`) copies the element and checks
that it lies in the target field; between different states it
transfers the element structurally (the generators it involves are
recreated in the target state from their definitions, with the
roots matched by their enclosures, so that even roots closer than
their printed approximations keep their identity). The real field
is thus the real subfield of a complex field whose computations it
shares: its values are computed through the complex numbers where
that is the natural way, and only the values entering it are
checked. As everywhere in `gr`, the operations of a view assume
that their operands are elements of the view: an element of the
complex field that is used directly in its real view without
going through `gr_set_other` (which is possible since the
elements have the same representation) is not checked for realness,
and the sign tests, comparisons and rounding functions of the view
then act on its real part.

::: details
Membership of the subfield is enforced where values enter it from
outside the operations: conversions, `gr_set_str` (whose definitions
may introduce auxiliary nonreal generators, as long as the value is
real) and `gr_gens`, which lists only the definitions lying in the
subfield (in a real field, not $i$ or $\exp(i)$). The realness test is
numerical when the imaginary part is separated from zero, and otherwise
the exact test $\overline{x} = x$. The results of the elementary
functions, of root finding and of factorization in a real field are
rewritten in real generators (the tangent normal form of the
trigonometric constants, the half-angle substitution for the
trigonometric functions): see [](algorithms/exact/castle-real-fields).
:::
:::

## The unlimited mode {#unlimited}

A castle field gives up on an operation (`GR_UNABLE`, or `T_UNKNOWN`
for a predicate) when the effort it is allowed to spend is exhausted:
the precision of the numerical separations and relation searches, the
degrees up to which polynomials are factored and irreducibility is
proved, the numbers of places tried, the sizes of special values
computed exactly. These bounds are tuning options
(`gr_ca_ctx_set_option`, [](gr_tower#tuning-options)). In the
*unlimited mode* (the flag `GR_CA_UNLIMITED`, or the option
`GR_TOWER_OPT_UNLIMITED`), an operation which would fail for this
reason is repeated with all such *effort limits* doubled, and again,
until it succeeds: an unlimited castle field never fails prematurely,
and can replace [](qqbar) where a computation may take as long as it
needs. The options keep their meaning as the starting values, and the
options which only choose between algorithms are unaffected; each
round tries all the methods of an operation again with the larger
limits, so that a fast heuristic is never skipped for a slow one.

An operation which fails for other reasons (an unsupported function, a
point on a branch cut) fails at once. An operation which no method can
decide -- a zero test of a number which is zero by an identity between
transcendental constants that the relation searches do not know --
does not terminate in the unlimited mode; in the limited mode it gives
`T_UNKNOWN` after the effort limits. See
[](algorithms/exact/effort-limits).

## Choice of generators, options and printing {#options}

::: function
void gr_ca_ctx_set_gen_flags(gr_ctx_t ctx, int flags)
int gr_ca_ctx_gen_flags(gr_ctx_t ctx)

How new generators are chosen. With
`GR_TOWER_GENS_SPLIT_IMAGINARY`, the square root of a negative
rational number is $i \sqrt{A}$ (the generators $i$ and
$\sqrt{p}$ for the primes $p$ dividing $A$, shared with the real
square roots) rather than the single generator $\sqrt{-A}$, which
is the default (0). By default, roots of unity of composite orders
are products of roots of unity of prime power orders (linearly
disjoint fields: canonical, sparse representations, monomials for
the powers of roots of unity), and square roots of positive
rationals are products of square roots of primes. With
`GR_TOWER_GENS_COMPOSITE_ROOTS` (which sets the option
`GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT` to
`GR_TOWER_COMPOSITE_ROOTS_DEGREE`, 128, when it is 0;
clearing it sets the option to 0), a root of unity is a power of
one generator $\zeta_N$ (a root of $\Phi_N$) for $N$ the least common
multiple of the orders requested, as in Calcium: arithmetic within
one cyclotomic field is dense univariate (3–4 times faster for
$\mathbb{Q}(\zeta_{60})$ and $\mathbb{Q}(\zeta_{255})$), but
computations involving many orders build large fields (the DFT
benchmark is several times slower); towers merged in this mode adjoin
a missing root of unity as $\zeta_L$ for $L$ the least common multiple
of its order and those present, which become powers of it, up to
the degree limit and `GR_TOWER_CYCLOTOMIC_ORDER_LIMIT`, beyond
which roots of prime power orders are used. The flags are a view of
the options `GR_TOWER_OPT_SPLIT_IMAGINARY`,
`GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT` and
`GR_TOWER_OPT_COMPOSITE_RADICALS`. With
`GR_TOWER_GENS_COMPOSITE_RADICALS`, the square root of a
positive rational $A B^2 / D^2$ ($A$ squarefree) is $B \sqrt{A} / D$
with $\sqrt{A}$ one generator (relations such as
$\sqrt{6} = \sqrt{2}\sqrt{3}$ are then found by the merges rather
than by construction). Elements created before a change remain
valid and compatible with the new ones.
:::

::: function
int gr_ca_ctx_set_option(gr_ctx_t ctx, slong option, slong value)
slong gr_ca_ctx_get_option(gr_ctx_t ctx, slong option)

Sets or gets a tuning option ([](gr_tower#tuning-options)) of the
castle field, shared by its views and its towers; a change takes effect
at once, also for the existing elements. Returns `GR_DOMAIN` for
a value out of the range of the option. The Python interface takes
the names (`gr_tower_option_name`), also as keyword arguments
of the constructors (`ComplexField_castle(cyclotomic_degree_limit=64)`).
:::

::: function
void gr_ca_ctx_set_print(gr_ctx_t ctx, int flags, slong digits)
int gr_ca_ctx_print_flags(gr_ctx_t ctx)
int gr_ca_ctx_field_flags(gr_ctx_t ctx)

Printing of elements, by any combination of
`GR_TOWER_PRINT_NUMERIC` (a numerical value with `digits`
digits), `GR_TOWER_PRINT_SYMBOLIC` (the expression in the
generators) and `GR_TOWER_PRINT_DEFS` (the definitions of the
generators the expression involves, transitively, in definition
order); the default is the expression with the definitions. For
instance $2 + e^2 + 3 e^{e^2}$ prints as
`3*t2+t1+2 {t1 = exp(2); t2 = exp(t1)}`, as `3*t2+t1+2`
without the definitions, as `4863.92` numerically, and as
`4863.92 {3*t2+t1+2 where t1 = exp(2); t2 = exp(t1)}` with
everything. The definitions are written as `sqrt(x)`,
`root(x, n)`, `exp(x)`, `log(x)`, `exp(2*pi*i/n)` and, for an
algebraic generator without a closed form,
`root(m(a), approximation)` with the minimal polynomial `m` over the
field below, in the order of creation of the generators (a definition
only involves earlier ones). A rational number prints as such. The
field flags are those of `gr_ctx_init_ca`.

::: details
The numerical value omits a real or imaginary part which is zero
to within the working precision (as for $\sqrt{-163}$, printed as
`12.7671*i`).

The printed form with the definitions is read back by
`gr_set_str`: the definitions are evaluated in order, each in
the names defined before it -- the names are local to the string,
so a string printed in one context or session is read in another --
then the expression in all of them (a `root(m, approximation)`
definition takes the unique root of `m` within the precision of the
approximation -- the digits printed for a root are enough to
separate it from the other roots of its polynomial -- and the string
is rejected with `GR_UNABLE` when several roots agree to the given
digits). A string without definitions is read in the generators of
the context by their names, with `pi`, `i` and the elementary
functions (`sqrt(2)`, `exp(a1)`). Conversions between castle
contexts (`gr_set_other`), for instance from the complex field
into the real or algebraic subfield, do not go through strings but
transfer the generators structurally; the conversion fails if the
number does not belong to the target field.
:::
:::

::: macro
GR_TOWER_PRINT_NUMERIC
GR_TOWER_PRINT_SYMBOLIC
GR_TOWER_PRINT_DEFS

:::

The symbolic expression of an element (`gr_get_fexpr`) is
`Where(expr, Def(a_1, Sqrt(2)), Def(t_1, Exp(a_1)), ...)` with the
generators as the symbols `a_1`, `t_1`, ... (`Pi` and
`NumberI` for $\pi$ and $i$), the definitions as `Exp`,
`Log`, `Sqrt`, `Pow(x, Div(1, n))`, `Exp(Div(Mul(2, Pi, NumberI), n))` and `PolynomialRootNearest(List(...), Decimal(...))` (with the origin polynomial over $\mathbb{Q}$ when
the generator has one, otherwise the coefficients of the minimal
polynomial over the field below as expressions). The LaTeX form
follows from `fexpr_get_str_latex`.

## Conversions and access to the towers {#towers}

::: function
int gr_ca_get_qqbar(qqbar_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_get_acb(acb_t res, const gr_ca_t x, slong prec, gr_ctx_t ctx)
int gr_ca_get_d(double * res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_root_ui(gr_ca_t res, const gr_ca_t x, ulong n, gr_ctx_t ctx)
gr_tower_struct * gr_ca_get_tower(slong * level, const gr_ca_t x, gr_ctx_t ctx)
const fmpz_mpoly_q_struct * gr_ca_get_data(const gr_ca_t x, gr_ctx_t ctx)
void gr_ca_ctx_stats(gr_ctx_t ctx)
slong gr_ca_ctx_num_towers(gr_ctx_t ctx)

Conversions and access to the underlying tower representation
of an element (the tower stays valid while `x` lives in it: towers
without elements are collected). `gr_ca_ctx_num_towers`
returns the number of towers of the context after collecting the
unused ones. The `qqbar` conversion rebuilds the algebraic part
of the tower over $\mathbb{Q}$ when the tower has transcendental
generators which neither the element nor the moduli involve.
`gr_ca_get_data` returns the element as a
`fmpz_mpoly_q` over the flat context of its tower, converting
the element to that representation in place if it is currently
stored as a rational number, in dense form or with a factored denominator (see
`gr_ca_struct`); the pointer stays valid until the
next operation involving `x`.
`gr_ca_ctx_stats` prints the number of towers, their
largest degree and the number of steps; with the environment
variable `GR_TOWER_STATS_VERBOSE` set, also the sizes of the
moduli of every tower with at least two steps (flat terms and
coefficient bits, or "deferred", and nested coefficients).
:::

## Elements and their representations {#elements}

::: type
gr_ca_struct
gr_ca_t

The elements of the castle fields (`elem_size` of the context); a
`gr_ca_t` is an array of length one of
`gr_ca_struct`, permitting it to be passed by
reference. The structure is public so that its data can be read
without copies, but it must not be written: an element refers to a
tower shared with other elements, and the context converts, reduces
and moves the data of elements lazily, under its lock. The fields:

::: details
* `slong repr`: which member of the union `elem` holds the
  element, one of the constants below.
* `fmpq elem.q` (`GR_CA_REPR_RATIONAL`): a rational
  number.
* `fmpz_mpoly_q_struct elem.flat.data` and
  `fmpz_mpoly_ctx_struct * elem.flat.mctx`
  (`GR_CA_REPR_FLAT`): a fraction over the flat
  context of the tower. The context and the reduction of the data
  may lag behind the tower (they are brought up to date by the next
  operation or accessor), so that the data should be read through
  `gr_ca_get_data`.
* `fmpq_poly_struct elem.dense.poly`
  (`GR_CA_REPR_DENSE`): a canonical polynomial of
  degree less than $d$ in the first generator $a$ of the tower, of
  degree $d$ (the modulus is obtained with
  `gr_ca_get_fmpq_poly`). The other member of this
  branch, `elem.dense.nf`, is internal.
* `fmpz_mpoly_q_factored_struct elem.factored.data` and
  `fmpz_mpoly_ctx_struct * elem.factored.mctx`
  (`GR_CA_REPR_FACTORED`): a rational function with a factored
  denominator ([](fmpz_mpoly_q#factored)) in a tower without algebraic
  generators, with the option `GR_TOWER_OPT_FACTORED_DENOMINATORS`.
* `F`, `level`, `structure_version`, `reduced_version`,
  `shallow`: the tower (through `F->T`) and the bookkeeping of
  the lazy updates (internal).

The layout may change between versions of FLINT; the accessors
below are the stable interface. A rational or dense element may be
read while other threads read it; reading a flat element races with
the lock-holding thread, which may convert it to the dense form, so
that it should be read through the locked accessors.
:::
:::

::: macro
GR_CA_REPR_FLAT
GR_CA_REPR_RATIONAL
GR_CA_REPR_DENSE
GR_CA_REPR_FACTORED

The representations of `gr_ca_struct`.
:::

::: function
int gr_ca_repr(const gr_ca_t x, gr_ctx_t ctx)

The current representation of `x` (one of the constants above). It
is a property of the storage, not of the value: the same number may
be flat in one element and dense in another, and an element may
change representation after an operation that involves it.
:::

::: function
int gr_ca_get_fmpq_poly(fmpq_poly_t res, fmpz_poly_t modulus, const gr_ca_t x, gr_ctx_t ctx)

Sets `res` to `x` as a polynomial in the first generator $a$ of
its tower, of degree less than the degree of $a$, and `modulus`
(unless `NULL`) to the minimal polynomial of $a$ (monic, with
integer coefficients), whatever the representation of `x`. A rational
number is a constant polynomial with modulus $x$. Returns
`GR_DOMAIN` if `x` involves another generator or $a$ is not
algebraic, and `GR_UNABLE` if the modulus of $a$ is not monic
integral or `x` has an algebraic denominator. The polynomial is
canonical when the step of $a$ is proven: equal elements give equal
polynomials. For example, $(1 + \sqrt{2})^3$ gives $7 + 5x$ with
modulus $x^2 - 2$.
:::

## The gr methods {#methods}

::: function
void gr_ca_init(gr_ca_t x, gr_ctx_t ctx)
void gr_ca_clear(gr_ca_t x, gr_ctx_t ctx)
void gr_ca_swap(gr_ca_t x, gr_ca_t y, gr_ctx_t ctx)
void gr_ca_set_shallow(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_randtest(gr_ca_t res, flint_rand_t state, gr_ctx_t ctx)
int gr_ca_write(gr_stream_t out, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_set_str(gr_ca_t res, const char * s, gr_ctx_t ctx)
int gr_ca_get_fexpr(fexpr_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_zero(gr_ca_t res, gr_ctx_t ctx)
int gr_ca_one(gr_ca_t res, gr_ctx_t ctx)
int gr_ca_pi(gr_ca_t res, gr_ctx_t ctx)
int gr_ca_i(gr_ca_t res, gr_ctx_t ctx)
int gr_ca_gens(gr_vec_t vec, gr_ctx_t ctx)
int gr_ca_set(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_set_si(gr_ca_t res, slong c, gr_ctx_t ctx)
int gr_ca_set_ui(gr_ca_t res, ulong c, gr_ctx_t ctx)
int gr_ca_set_fmpz(gr_ca_t res, const fmpz_t c, gr_ctx_t ctx)
int gr_ca_set_fmpq(gr_ca_t res, const fmpq_t c, gr_ctx_t ctx)
int gr_ca_set_other(gr_ca_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
int gr_ca_get_si(slong * res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_get_ui(ulong * res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_get_fmpz(fmpz_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_get_fmpq(fmpq_t res, const gr_ca_t x, gr_ctx_t ctx)
truth_t gr_ca_is_zero(const gr_ca_t x, gr_ctx_t ctx)
truth_t gr_ca_is_one(const gr_ca_t x, gr_ctx_t ctx)
truth_t gr_ca_equal(const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_neg(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_add(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_sub(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_mul(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_div(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_inv(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_add_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx)
int gr_ca_add_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx)
int gr_ca_add_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx)
int gr_ca_add_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx)
int gr_ca_sub_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx)
int gr_ca_sub_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx)
int gr_ca_sub_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx)
int gr_ca_sub_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx)
int gr_ca_mul_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx)
int gr_ca_mul_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx)
int gr_ca_mul_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx)
int gr_ca_mul_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx)
int gr_ca_div_si(gr_ca_t res, const gr_ca_t x, slong c, gr_ctx_t ctx)
int gr_ca_div_ui(gr_ca_t res, const gr_ca_t x, ulong c, gr_ctx_t ctx)
int gr_ca_div_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t c, gr_ctx_t ctx)
int gr_ca_div_fmpq(gr_ca_t res, const gr_ca_t x, const fmpq_t c, gr_ctx_t ctx)
int gr_ca_pow(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_sqrt(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_abs(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_conj(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_exp(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_log(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_poly_mullow(gr_ca_struct * res, const gr_ca_struct * p1, slong len1, const gr_ca_struct * p2, slong len2, slong n, gr_ctx_t ctx)
int gr_ca_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
int gr_ca_poly_factor(gr_poly_t c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
int gr_ca_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
int gr_ca_mat_det(gr_ca_t res, const gr_mat_t A, gr_ctx_t ctx)
int gr_ca_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)

The methods of the `gr` interface of the castle fields, which the
generic element functions of `gr.h` (`gr_add` and so on)
dispatch to, and which may also be called directly; the semantics
are those of the generic functions. Operations on rational elements,
and on dense elements of the same field (with rational scalars),
take lock-free paths; the others hold the context's lock.
`gr_ca_set_other` converts an element `x` of the
context `x_ctx`. In the real and algebraic views (see
`gr_ctx_init_ca_view`), operations whose result
leaves the subfield return `GR_DOMAIN`.
:::

## Elementary functions, parts and rounding {#elementary}

::: function
int gr_ca_sin(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_cos(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_tan(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_sinh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_cosh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_tanh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_asin(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_acos(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_atan(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_asinh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_acosh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_atanh(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_re(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_im(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_arg(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_sgn(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_csgn(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
truth_t gr_ca_is_real(const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_cmp(int * res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_cmpabs(int * res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_floor(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_ceil(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_trunc(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_nint(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)

The elementary functions and real-number operations of the lazy
field, as installed in its method table. The inverse functions
return `GR_UNABLE` on their branch cuts (where the logarithmic
formulas may differ from the principal branch). As for
`qqbar_t`, the rounding functions act on the real part
(`floor(2 + i sqrt(2)) = 2`), `csgn` is the sign of the real part,
or of the imaginary part when the real part is zero, and `cmpabs`
compares $x \bar x$ and $y \bar y$; `cmp` returns `GR_DOMAIN` for
elements known to be nonreal (`GR_UNABLE` when realness cannot be
decided). All of these are exact: enclosures decide when they
separate, and the exact zero test decides otherwise.
:::

## Special functions {#special}

::: function
int gr_ca_gamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_rgamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_lgamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_beta(gr_ca_t res, const gr_ca_t x, const gr_ca_t y, gr_ctx_t ctx)
int gr_ca_digamma(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_polygamma(gr_ca_t res, const gr_ca_t s, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_erf(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_erfc(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_erfi(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_lambertw(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_lambertw_fmpz(gr_ca_t res, const gr_ca_t x, const fmpz_t k, gr_ctx_t ctx)
int gr_ca_zeta(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_hurwitz_zeta(gr_ca_t res, const gr_ca_t s, const gr_ca_t a, gr_ctx_t ctx)
int gr_ca_polylog(gr_ca_t res, const gr_ca_t s, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_lerch_phi(gr_ca_t res, const gr_ca_t z, const gr_ca_t s, const gr_ca_t a, gr_ctx_t ctx)
int gr_ca_dirichlet_l(gr_ca_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const gr_ca_t s, gr_ctx_t ctx)
int gr_ca_dilog(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_elliptic_k(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_elliptic_e(gr_ca_t res, const gr_ca_t x, gr_ctx_t ctx)
int gr_ca_euler(gr_ca_t res, gr_ctx_t ctx)
int gr_ca_catalan(gr_ca_t res, gr_ctx_t ctx)
int gr_ca_hypgeom_0f1(gr_ca_t res, const gr_ca_t b, const gr_ca_t z, int flags, gr_ctx_t ctx)
int gr_ca_hypgeom_1f1(gr_ca_t res, const gr_ca_t a, const gr_ca_t b, const gr_ca_t z, int flags, gr_ctx_t ctx)
int gr_ca_hypgeom_2f1(gr_ca_t res, const gr_ca_t a, const gr_ca_t b, const gr_ca_t c, const gr_ca_t z, int flags, gr_ctx_t ctx)
int gr_ca_hypgeom_pfq(gr_ca_t res, const gr_vec_t a, const gr_vec_t b, const gr_ca_t z, int flags, gr_ctx_t ctx)
int gr_ca_hypgeom_pfq_vec(gr_ca_t res, const gr_ca_struct * a, slong p, const gr_ca_struct * b, slong q, const gr_ca_t z, int flags, gr_ctx_t ctx)

The special functions of the castle field, as installed in its method
table, with canonical arguments (see *Special functions* above).
Poles give `GR_DOMAIN`. The order of `polygamma` must be an
integer; `polylog` with a non-integer order is evaluated at the
roots of unity (through the Lerch transcendent), and `lerch_phi` at
$z = 1$, at the roots of unity, and at $a = 1$ with an integer $s$
(`GR_UNABLE` otherwise). `gr_dirichlet_l` dispatches to
`gr_ca_dirichlet_l` for the castle fields. In the real fields, a value
which is not real gives `GR_DOMAIN`; in the algebraic fields, only
algebraic values (such as $\Gamma(n)$ and $\zeta(-n)$) are returned,
and `GR_UNABLE` otherwise, since the transcendence of special
function values is not known in general. The printed forms
(`gamma(...)`, `erf(...)`, `erfi(...)`, `lambertw(x, k)`,
`digamma(...)`, `polygamma(m, x)`, `polylog(s, x)`,
`zeta(...)`, `elliptic_k(...)`, `elliptic_e(...)`, `euler`,
`catalan`) are read back by the parser (the hypergeometric
functions, of several arguments, are not). The hypergeometric
functions take the flag 1 for the regularized functions (as
`gr_hypgeom_pfq`); their normal form is described in
*Hypergeometric functions* above.
:::

## Modular forms, theta and elliptic functions {#modular}

::: function
int gr_ca_modular_lambda(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_modular_j(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_modular_delta(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_dedekind_eta(gr_ca_t res, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_eisenstein_e(gr_ca_t res, ulong k, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_eisenstein_g(gr_ca_t res, ulong k, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_jacobi_theta(gr_ca_t res1, gr_ca_t res2, gr_ca_t res3, gr_ca_t res4, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_jacobi_theta_1(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_jacobi_theta_2(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_jacobi_theta_3(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_jacobi_theta_4(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_weierstrass_p(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_weierstrass_p_prime(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_weierstrass_sigma(gr_ca_t res, const gr_ca_t z, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_elliptic_roots(gr_ca_t e1, gr_ca_t e2, gr_ca_t e3, const gr_ca_t tau, gr_ctx_t ctx)
int gr_ca_elliptic_invariants(gr_ca_t g2, gr_ca_t g3, const gr_ca_t tau, gr_ctx_t ctx)

The modular lambda function, the $j$-invariant, the discriminant
$\Delta = \eta^{24}$, the Dedekind eta function, the Eisenstein series
$E_k$ (normalized to constant term 1) and $G_k = 2\zeta(k) E_k$ for even
$k \ge 2$ (`GR_DOMAIN` for other $k$, as for the complex numbers;
`GR_UNABLE` for $k > 1000$), the Jacobi theta functions (with the conventions of
`acb_modular_theta`, $q = e^{\pi i \tau}$), and the Weierstrass
functions $\wp$, $\wp'$, $\sigma$, the roots $e_1, e_2, e_3$ and the
invariants $g_2, g_3$ of the lattice $\mathbb{Z} + \tau \mathbb{Z}$
(as `acb_elliptic_p`; $\wp$ at a lattice point gives
`GR_DOMAIN`), in the normal form
described in *Modular forms and theta functions* above. They are
installed in the method table of the castle field (`gr_modular_j` and
so on). A point with $\operatorname{Im}(\tau) \le 0$ gives
`GR_DOMAIN`; `GR_UNABLE` is returned when $\tau$ cannot be reduced
exactly (for instance when the sign of $\operatorname{Im}(\tau)$ is
undecided). In the algebraic fields, only algebraic values ($j$ and
$\lambda$ at CM points) are returned. `modular_lambda`,
`modular_j`, `modular_delta` and `dedekind_eta` are read back by
the parser.
:::
