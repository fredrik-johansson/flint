---
title: Castle arithmetic
summary: Exact real and complex arithmetic with castle fields: numbers in towers of extension fields built on demand over the rationals, proved predicates and the meaning of unknown answers, dynamic evaluation, the unlimited mode, the real and algebraic subfields, and how castle arithmetic relates to qqbar and to the earlier Calcium design.
order: 3
---

FLINT computes exactly with real and complex numbers in *castle
arithmetic*. A number is represented symbolically, as an element of a
*tower* of field extensions of $\mathbb{Q}$, such as
$\mathbb{Q}(\sqrt 2, \sqrt 3)$, $\mathbb{Q}(\pi, e^{i \pi/5})$ or
$\mathbb{Q}(\log 2, \sqrt[3]{\log 2 + 1})$; arithmetic is formal arithmetic
in the tower, and numerical values are available to any precision. A
*castle field* is the collection of all the towers of a computation,
built on the common ground field $\mathbb{Q}$: the towers are created,
extended, refined and merged on demand as numbers are computed, and the
relations found between their generators are recorded in them, so that
every later computation benefits. The name follows *ball arithmetic*
([](guide/types/ball-arithmetic)): where a ball is an approximation that
carries its error, a castle number carries the structure that makes it
exact.

The point of exact numbers is equality. Ball arithmetic can prove that
two numbers differ, by computing them accurately enough that their
balls are disjoint, but never that they are equal; castle arithmetic
decides equality of algebraic numbers always, of numbers built with
exponentials and logarithms whenever Schanuel's conjecture holds (every
answer is proved, the conjecture being needed only for termination),
and of values of special functions through their known identities.

Castle arithmetic is implemented by the modules [](gr_ca) (the castle
fields, as [](gr) contexts) and [](gr_tower) (the towers). It succeeds
the design of the [](ca) module, formerly the library Calcium
[@Joh2021calcium], from which it inherits the "ca" mnemonic; see
[](#legacy) for how the two differ. For individual algebraic numbers in
canonical form, see [](qqbar). The algorithms are described in
[](algorithms/exact).

## Castle semantics {#semantics}

A castle field $K$ is a field of complex numbers: the complex field of
all the numbers its towers can represent, or its real or algebraic
subfield. An element of $K$ is an element $x$ of the top field $F$ of
some tower $\mathbb{Q} \subset \cdots \subset F$, where every generator
of the tower is a definite complex number, fixed by its *definition*
($\sqrt 2$, $\exp(u)$ for an element $u$ of a lower field, the root of
a polynomial nearest to a given approximation, $\Gamma(1/3)$, ...) and
by a numerical enclosure, which selects the conjugate or the branch.
The value of $x$ is thus a definite complex number, and the field
operations of $K$ are exact.

Every answer of castle arithmetic is proved. An operation either
returns `GR_SUCCESS` with the exact result, or fails: `GR_DOMAIN` when
the result does not exist (division by zero, a pole, a number outside
the real field of a real context), `GR_UNABLE` when it could not be
computed. A predicate returns `T_TRUE` or `T_FALSE` only when it is
proved, and `T_UNKNOWN` otherwise. The *zero test* behind all
predicates is

- complete for algebraic numbers: the generators of a tower need not be
  known to be irreducible, since an operation that exposes a
  factorization refines the tower (*dynamic evaluation*,
  [](algorithms/exact/dynamic-evaluation));
- complete under Schanuel's conjecture for numbers built with $\exp$,
  $\log$ and the trigonometric functions: a number which is zero must be
  so because of a relation between the generators, which is found by a
  lattice search and proved exactly (Richardson's algorithm,
  [](algorithms/exact/zero-testing#castle));
- for the values of special functions, as complete as the identities it
  knows: shifts, reflection and multiplication formulas, special values
  and transformations of the gamma, zeta, polylogarithm, hypergeometric,
  elliptic and modular functions
  ([](algorithms/exact/castle-special-functions)).

Within the default *effort limits*, a zero test which needs more work
than allowed returns `T_UNKNOWN`; in the *unlimited mode* the work grows
until the answer is found ([](#effort)).

## A first program {#first}

A castle field is created with `gr_ctx_init_ca` over the rational
numbers, and its elements are used through the functions of
[](gr):

```c example name=castle_first
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_special.h>
#include <flint/gr_ca.h>

int main(void)
{
    gr_ctx_t QQ, C;
    gr_ptr x, y;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_ca(C, QQ, 0);
    x = gr_heap_init(C);
    y = gr_heap_init(C);

    /* x = exp(pi*i/4) */
    status |= gr_set_str(x, "exp(pi*i/4)", C);
    status |= gr_println(x, C);

    /* x^2 = i, a power of the generator zeta_8 */
    status |= gr_sqr(y, x, C);
    status |= gr_println(y, C);

    /* log(x) = pi*i/4 */
    status |= gr_log(y, x, C);
    status |= gr_println(y, C);

    /* the identity sqrt(2) + sqrt(3) = sqrt(5 + 2 sqrt(6)) */
    status |= gr_set_str(x, "sqrt(2) + sqrt(3) - sqrt(5 + 2*sqrt(6))", C);
    flint_printf("%d\n", gr_is_zero(x, C) == T_TRUE);

    gr_heap_clear(x, C);
    gr_heap_clear(y, C);
    gr_ctx_clear(C);
    gr_ctx_clear(QQ);
    flint_cleanup();
    return status;
}
```
```output
a1 {a1 = exp(2*pi*i/8)}
a1^2 {a1 = exp(2*pi*i/8)}
pi*i/4
1
```

The number $e^{\pi i/4}$ is the root of unity $\zeta_8$, a generator of
its own tower; its square $i$ is a power of it, and its logarithm is
read off from the definition. The last line shows that the identity
$\sqrt 2 + \sqrt 3 = \sqrt{5 + 2\sqrt 6}$ is proved.

Numbers print as an expression in the generators of their tower,
followed by the definitions of the generators involved; the names
$a_1, a_2, \ldots$ (algebraic generators) and $t_1, t_2, \ldots$
(transcendental generators) are given in order of creation, and $\pi$
and $i$ keep their names. The printed form is read back by
`gr_set_str` in any castle field. `gr_ca_ctx_set_print` selects a
numerical display instead ([](gr_ca#options)).

## Towers {#towers}

The towers of a castle field are built as numbers are computed. The
square root of an element is first sought in the tower of the element,
and adjoined to it if it is not found there; an exponential or a
logarithm is adjoined as a transcendental generator, unless a relation
with the generators present expresses it ($\exp(\log x) = x$); roots of
rational numbers and of unity are products of canonical generators of
prime power orders ($\sqrt 6 = \sqrt 2 \sqrt 3$, $\zeta_{15} = \zeta_3^a
\zeta_5^b$); and an operation on numbers of different towers takes place
in a tower which contains both, found among the existing ones by the
definitions of the generators or built by merging them
([](algorithms/exact/field-towers)). Relations found by zero tests are
recorded in the towers: a generator which turns out to be expressible
through the others is eliminated, and all the numbers built on the tower
see the simplification.

A tower is shared by all the numbers built in it, and a castle field
collects the towers in which no number lives any more, so that a long
computation with few live numbers keeps few towers. How a castle field
manages its towers, and how to inspect them, is described in
[](gr_ca/using).

## Predicates {#predicates}

A predicate of castle arithmetic returns a `truth_t`, and `T_UNKNOWN`
is a legitimate answer in the default mode: the code that calls it
has three branches.

```c
truth_t z = gr_is_zero(x, C);

if (z == T_TRUE)
{
    ...  /* x = 0, proved */
}
else if (z == T_FALSE)
{
    ...  /* x != 0, proved */
}
else
{
    ...  /* not decided within the effort limits */
}
```

The functions which need a predicate internally (division, inversion,
the generic polynomial and matrix algorithms of [](gr_poly) and
[](gr_mat)) return `GR_UNABLE` when it is unknown. For algebraic
numbers, `T_UNKNOWN` can only come from the effort limits; for numbers
with transcendental generators, also from a zero which no known
relation explains (an identity that the relation searches do not know).

## Effort and the unlimited mode {#effort}

Deciding a predicate may require numerical evaluation at a high
precision, factoring polynomials over a tower of high degree, or
proving irreducibility at many places. The work spent is bounded by the
*tuning options* of the castle field (`gr_ca_ctx_set_option`, the
options of [](gr_tower#tuning-options)), which select algorithms and
limit their effort; the defaults are chosen so that typical
computations are fast and hopeless ones fail quickly.

A castle field created with the flag `GR_CA_UNLIMITED` never gives up
because of these limits: an operation which would fail after
exhausting one of them is repeated with all the effort limits doubled,
until it succeeds. The options remain in force as the starting values.
This makes an unlimited castle field a replacement for [](qqbar)
where a computation must not fail, at the price of possibly unbounded
time: a zero test that no method can decide does not terminate.

```c example name=castle_unlimited
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_special.h>
#include <flint/gr_ca.h>

int main(void)
{
    gr_ctx_t QQ, C;
    gr_ptr x;
    int unlimited;

    gr_ctx_init_fmpq(QQ);

    for (unlimited = 0; unlimited <= 1; unlimited++)
    {
        gr_ctx_init_ca(C, QQ, unlimited ? GR_CA_UNLIMITED : 0);
        /* a deliberately low limit: 64 bits of numerical evaluation */
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(C, GR_TOWER_OPT_NUMERIC_PREC_LIMIT, 64));
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(C, GR_TOWER_OPT_PREC_LIMIT, 64));
        GR_MUST_SUCCEED(gr_ca_ctx_set_option(C, GR_TOWER_OPT_CERTIFY_PREC_LIMIT, 64));

        /* erf(30) - 1 = -2.6e-393 */
        x = gr_heap_init(C);
        GR_MUST_SUCCEED(gr_set_si(x, 30, C));
        GR_MUST_SUCCEED(gr_erf(x, x, C));
        GR_MUST_SUCCEED(gr_sub_si(x, x, 1, C));
        truth_t z = gr_is_zero(x, C);
        flint_printf("unlimited = %d: %s\n", unlimited,
            z == T_TRUE ? "zero" : z == T_FALSE ? "nonzero" : "unknown");
        gr_heap_clear(x, C);
        gr_ctx_clear(C);
    }

    gr_ctx_clear(QQ);
    flint_cleanup();
    return 0;
}
```
```output
unlimited = 0: unknown
unlimited = 1: nonzero
```

## Real and algebraic fields {#subfields}

The flags `GR_CA_REAL` and `GR_CA_ALGEBRAIC` of `gr_ctx_init_ca` give
the real subfield and the algebraic subfields of the complex castle
field. In a real field, an operation whose result is not real fails
with `GR_DOMAIN` ($\sqrt{-1}$, $\log(-1)$, a nonreal polynomial root),
comparisons and rounding are available, and results are written in
real generators ($\sin(\pi/5)$ as a polynomial in $\tan(\pi/5)$, $\cos
1$ as a rational function of $\tan(1/2)$). In an algebraic field, $\pi$
is not available and the elementary functions only at the arguments
where their values are algebraic, by the Lindemann–Weierstrass and
Gelfond–Schneider theorems. The subfields can share the towers of a
complex field as *views* (`gr_ctx_init_ca_view`), so that conversions
between them are cheap ([](gr_ca#contexts)).

## Castle arithmetic, qqbar and Calcium {#legacy}

FLINT has had two other designs for exact numbers, which remain
available.

A [](qqbar) number is an algebraic number in canonical form: its
minimal polynomial over $\mathbb{Q}$ and an isolating ball. Equality is
structural and every predicate is decided, but an arithmetic operation
computes the minimal polynomial of the result, whose degree is that of
the compositum of the operands: sums of many square roots, or the
roots of a polynomial combined with each other, quickly become
expensive ([](algorithms/exact/qqbar-arithmetic)). Castle arithmetic
keeps such numbers in towers with many small steps, where arithmetic is
cheap, and computes a minimal polynomial only when asked
(`gr_ca_get_qqbar`).

The [](ca) module (Calcium) represents a number as an element of a
*formal field* $\mathbb{Q}(a_1, \ldots, a_n)$, the fraction field of a
polynomial ring modulo an *ideal of relations* between the generators,
which is built when the field is created: minimal polynomials, relations
found by lattice reduction and proved, and a Gröbner basis
([](algorithms/exact/field-towers#legacy)). The fields are cached and
identified by their generators. This design, the *legacy design*, is
the one that castle arithmetic succeeds; the main differences are:

| | Calcium (`ca`) | Castle arithmetic (`gr_ca`) |
|---|---|---|
| Fields | formal fields, one per set of generators, with an ideal computed when the field is created | towers of extensions, shared and grown in place, merged when numbers meet |
| Algebraic generators | algebraic numbers in canonical form ([](qqbar)), so that minimal polynomials are needed | steps with any polynomial vanishing at the root; refined when a factorization appears (dynamic evaluation) |
| Relations | searched when a field is created, within limits on the number of generators and degrees | searched by the zero tests that need them, among the generators involved, and recorded in the tower for all later computations |
| Zero test of algebraic numbers | decided within degree limits (falling back on [](qqbar)) | complete |
| Zero test with $\exp$, $\log$ | many cases, with heuristics | Richardson's algorithm: complete under Schanuel's conjecture |
| Special functions | elementary functions, $\Gamma$, $\operatorname{erf}$ and a few others | also polygamma and zeta values, polylogarithms, Dirichlet $L$-functions, hypergeometric, elliptic and modular functions, with their identities |
| Special values | infinities and *Undefined* are elements | fields of numbers; poles and undefined values are errors (`GR_DOMAIN`) |
| Representation | fractions of multivariate polynomials, one representation | rational numbers, flat fractions, dense forms in number fields, factored denominators ([](algorithms/exact/castle-representations)) |
| Interface | `ca_t` functions, `gr` contexts `gr_ctx_init_complex_ca` and others | `gr` contexts `gr_ctx_init_ca` |

Calcium remains available, and its documentation describes the legacy
design in its own terms ([](ca), [](calcium)); new code should use
castle arithmetic. The example programs `machin`, `huge_expr`,
`hilbert_matrix_ca` and `dft` run in both designs (with the option
`-castle` for castle arithmetic), see [](tutorials/exact). On a
benchmark of 744 zero-testing and simplification problems drawn from the
literature and from earlier systems, castle arithmetic decided every
case, Calcium 561, SymPy 522 (one wrongly) and [](qqbar) 313 of the 353
algebraic cases it accepts; on the cases that both decide, the two
designs take about the same time. [](algorithms/exact/castle-design)
describes the design as a whole, with the full comparison
[@Joh2026tower].

## Conventions {#conventions}

- A castle context is thread-safe: operations on it hold a recursive
  mutex of the context, so that elements may be shared between threads;
  the arithmetic of rational numbers and of numbers in one number field
  takes lock-free paths. Distinct contexts are independent.
- Elements of a castle field refer to its towers: they must be used
  with that context (or one of its views) and cleared before it. An
  element of another castle field is converted with `gr_set_other`,
  which transfers the generators it involves by their definitions.
- A castle field never holds infinities or undefined values: a pole or
  an undefined value is an error (`GR_DOMAIN`), as in the other fields
  of [](gr).
- Generic algorithms over castle fields ([](gr_poly), [](gr_mat)) are
  exact: they return `GR_UNABLE` rather than a wrong result when a
  predicate they need is unknown.
