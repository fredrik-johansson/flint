---
title: gr_tower.h
header: gr_tower.h
description: Towers of algebraic and transcendental extensions of the rationals, with dynamic evaluation, complete zero tests, relation searches and a flat representation: the core of castle arithmetic.
category: Exact real and complex numbers
order: 0.2
---

A `gr_tower_t` represents a tower of fields

$$
F_0 \subset F_1 \subset \cdots \subset F_n, \quad F_k = F_{k-1}[x] / (m_k),
$$

where each step adjoins an algebraic element $a_k$ with monic minimal
polynomial $m_k \in F_{k-1}[x]$, and the base field $F_0 = K(t_1, \ldots,
t_r)$ is a rational function field over the constant field $K$
(currently $\mathbb{Q}$) in transcendental generators such as $\pi$,
$\exp(u)$, $\log(u)$ and values of special functions, where $u$ is an
element of the tower defined by the generators preceding $t_j$. Every
generator is embedded in $\mathbb{C}$ by a numerical enclosure, which
fixes the choice of conjugate or branch, so that the tower represents a
definite subfield $F_n \subset \mathbb{C}$.

This module is the core of *castle arithmetic* ([](gr_ca),
[](guide/types/castle-arithmetic)): a castle field keeps a collection of
towers over the common ground field $\mathbb{Q}$ and builds, extends and
merges them as numbers are computed. The tower API can also be used
directly, to build a field step by step and compute in it through
`gr_ctx_init_tower_field`.

The generators are ordered by their *definition order* (the order in
which they were adjoined, compatible with dependencies); the first $p$
generators in this order form a *prefix* of the tower, which is itself a
tower. The algebraic generators have a numbering $a_1, \ldots, a_n$ of
their own and the transcendental generators $t_1, \ldots, t_r$
likewise. Transcendental generators come with a status: $\pi$ is
transcendental, while $\exp(u)$ and $\log(u)$ are transcendental over
the field below only conjecturally (Schanuel's conjecture) in general,
and the values of special functions beyond that.

Elements of $F_n$ have two representations. In the *nested*
representation, an element of $F_k$ is a polynomial over $F_{k-1}$ of
degree less than $\deg m_k$ (see `gr_ctx_init_gr_poly_quotient`); this
is the structural representation, which carries the triangular set, the
refinements and the certificates. In the *flat* representation, an
element is a multivariate rational function in all the generators
([](fmpz_mpoly_q)), reduced modulo the triangular set, which is a
Gröbner basis for the lexicographic order; this is the representation
in which castle fields compute. Both are efficient for towers of many
small extensions, such as $\mathbb{Q}(\sqrt{2}, \sqrt{3}, \ldots,
\sqrt{47})$ of degree $2^{15}$, where a primitive element would be
impractical, as well as for a single extension of high degree.

The minimal polynomials $m_k$ need not be known to be irreducible: a step
may be adjoined with a polynomial that merely vanishes at the root, and
the tower is refined when an operation exposes a factorization (*dynamic
evaluation*, [](algorithms/exact/dynamic-evaluation)). This makes the
zero test of algebraic numbers complete without irreducibility proofs.
With transcendental generators, an element which is nonzero as a
rational function of the generators may be zero as a number; the zero
test then searches for relations between the generators
([](algorithms/exact/zero-testing#castle)), which are proved exactly and
recorded in the tower by eliminating a generator.
## Types and macros {#types-and-macros}


::: type
gr_tower_gen_struct

A generator of a tower: its current kind (algebraic, or one of
`GR_TOWER_EXP`, `GR_TOWER_LOG`, `GR_TOWER_PI`, `GR_TOWER_TAN`,
`GR_TOWER_ATAN`), the kind of its
definition (which differs from the current kind after an elimination;
for algebraic generators also `GR_TOWER_ROOT`,
`GR_TOWER_ROOT_OF_UNITY` and `GR_TOWER_TAN_PI`),
a status flag, its index ($k$ or $j$), the context of $F_k$ if
algebraic, the argument of $\exp$ or $\log$ as a flat element, an
enclosure, a name, a definition id and its definition order.
`GR_TOWER_GEN(T, d)` is the generator with definition order `d`,
`GR_TOWER_STEP(T, k)` the algebraic generator $a_{k+1}$ and
`GR_TOWER_TRANS(T, j)` the transcendental generator $t_{j+1}$.
:::

::: type
gr_tower_struct
gr_tower_t

The tower object.
:::

::: macro
GR_TOWER_TAN_PI

The definition kind of the real algebraic generator
$\tan(\pi / M)$ (`def_param` $M$) of the tangent normal form of the
real trigonometric constants (see the lazy real field).
:::

::: macro
GR_TOWER_STATUS_PROVEN
GR_TOWER_STATUS_DYNAMIC
GR_TOWER_STATUS_INDEPENDENT
GR_TOWER_STATUS_SCHANUEL
GR_TOWER_STATUS_CONJECTURAL

Status of a generator. For a step: the minimal polynomial is known
to be irreducible over the field below (`PROVEN`), or merely
assumed to be so (`DYNAMIC`). For a transcendental generator:
`PROVEN` (its alias `INDEPENDENT`) means that the generator is
known to be algebraically independent of all the generators
preceding it in the definition order, which is what the structural
arguments of the module use (currently only $\pi$ at the head of a
tower has this status: a proof of transcendence alone does not give
it when another transcendental generator precedes); `SCHANUEL`
and `CONJECTURAL` mean that the independence is conjectural (under
Schanuel's conjecture, or beyond it for the values of special
functions).
:::

## Memory management and basic access {#memory-management-and-basic-access}


::: function
void gr_tower_init(gr_tower_t T, gr_ctx_t base)

Initializes `T` to the trivial tower over the field `base`.
The base context must remain valid during the lifetime of `T`.
:::

::: function
void gr_tower_clear(gr_tower_t T)
gr_tower_struct * gr_tower_heap_init(gr_ctx_t base)
void gr_tower_heap_clear(gr_tower_struct * T)
int gr_tower_get_str(char ** s, const gr_tower_t T)
int gr_tower_gen_get(gr_ptr res, const gr_tower_t T, slong d)
gr_ctx_struct * gr_tower_field_ptr(const gr_tower_t T)
slong gr_tower_length_si(const gr_tower_t T)
slong gr_tower_num_gens_si(const gr_tower_t T)
slong gr_tower_num_trans_si(const gr_tower_t T)
ulong gr_tower_version(const gr_tower_t T)
const char * gr_tower_gen_name(const gr_tower_t T, slong d)
int gr_tower_gen_status(const gr_tower_t T, slong d)
int gr_tower_gen_kind(const gr_tower_t T, slong d)

Heap-allocated towers, the string written by `gr_tower_write`,
the generator with definition order `d` as an element of the top
field, and non-inline accessors (for language bindings; the Python
wrapper's `gr_tower` class is built on them). `gr_tower_version` is
incremented by every refinement or rebuild of the tower.
:::

::: function
void gr_tower_set(gr_tower_t res, const gr_tower_t T)

Sets `res` to a copy of `T`. The copy has its own contexts, so
refinements of one tower do not affect the other.
:::

::: function
slong gr_tower_length(const gr_tower_t T)

The number of algebraic steps $n$. The number of transcendental
generators is `T->num_trans` and the total number of generators
is `gr_tower_num_gens`.
:::

::: function
slong gr_tower_num_gens(const gr_tower_t T)
slong gr_tower_prefix_length(const gr_tower_t T, slong p)
slong gr_tower_prefix_num_trans(const gr_tower_t T, slong p)
slong gr_tower_order_code(const gr_tower_t T, slong d)
slong gr_tower_gid_order(const gr_tower_t T, slong gid)
slong gr_tower_find_def_order(const gr_tower_t T, ulong def_id)
slong gr_tower_find_def(const gr_tower_t T, ulong def_id)
slong gr_tower_find_trans_def(const gr_tower_t T, ulong def_id)

The number of generators, the number of algebraic, respectively
transcendental, generators among the first `p` generators in
definition order, the code of the generator with definition order `d`
($k > 0$ for the algebraic generator $a_k$, $-j$ for the
transcendental generator $t_j$), and the definition order of the
generator with the given gid (a tower-local id which is stable under
reorderings) or definition id (-1 if absent); `gr_tower_find_def` and
`gr_tower_find_trans_def` give the index $k$ of the algebraic generator
$a_k$, respectively $j$ of the transcendental generator $t_j$, with the
given definition id (0 if none). The definition order of a generator can
change when the tower is reordered (by the relation search);
`T->structure_version` is incremented then.
:::

::: function
void gr_tower_gen_set_def_id(gr_tower_t T, slong d, ulong def_id)

Sets the definition id of the generator with definition order `d`.
Definition ids identify generators across towers: the merges
(`gr_tower_absorb`) identify generators with equal nonzero
ids. The castle fields assign them to the generators they create; a
user of the tower API may assign its own.
:::

::: function
void gr_tower_set_prefix(gr_tower_t res, const gr_tower_t T, slong p)
void gr_tower_set_subset(gr_tower_t res, const gr_tower_t T, const int * mark)

Sets `res` to a copy of the first `p` generators of `T` in definition
order, respectively of the generators with nonzero `mark[d]`
(indexed by definition order, in that order). The marked set must be
closed under the dependencies of the definitions: the generators in
the moduli and arguments of marked generators must be marked. The
statuses carry over (a modulus irreducible over the field of all
preceding generators is irreducible over a subfield).
:::

::: function
gr_ctx_struct * gr_tower_base(const gr_tower_t T)
gr_ctx_struct * gr_tower_field_at(const gr_tower_t T, slong k)
gr_ctx_struct * gr_tower_field(const gr_tower_t T)

Contexts of $F_0$, $F_k$ and $F_n$. Elements of the tower are
elements of these contexts. Note that the context of $F_k$ is a
generic quotient ring: its `gr_is_zero` and `gr_inv`
are complete only when the tower has been fully refined; use the
`gr_tower` versions of these operations, or the context created by
`gr_ctx_init_tower_field`, for complete versions.
:::

::: function
slong gr_tower_step_degree(const gr_tower_t T, slong k)
slong gr_tower_degree_at(const gr_tower_t T, slong k)
slong gr_tower_degree(const gr_tower_t T)
const gr_poly_struct * gr_tower_step_minpoly(const gr_tower_t T, slong k)

The degree $[F_k : F_{k-1}]$, the degrees $[F_k : F_0]$ and $[F_n :
F_0]$ (saturating at `WORD_MAX`), and the polynomial $m_k$. Here $1 \le
k \le n$.
:::

::: function
int gr_tower_write(gr_stream_t out, const gr_tower_t T)
void gr_tower_print(const gr_tower_t T)

:::

## Adjoining generators {#adjoining-generators}


::: function
int gr_tower_adjoin_algebraic(gr_tower_t T, const gr_poly_t m, const acb_t z, int status, const char * name)

Adjoins the root of the monic polynomial `m` over the current top
field $F_n$ nearest to the midpoint of `z`. The root is certified,
by an interval Newton step, to be the only root of `m` in a disk
around the midpoint of `z` whose radius is that of `z`, inflated
by up to $|z|/4$ when `z` is tighter than the roots can be
separated: `z` may be a rigorous enclosure or an approximation to a
few digits (`1.4142` selects $\sqrt 2$ among the roots of
$x^2 - 2$). Returns `GR_UNABLE` if no such disk is found (several
roots near `z`, or insufficient precision), or `GR_DOMAIN` if `m`
is not monic of positive degree. If `name` is `NULL`, a default
name is generated. `status` should be `GR_TOWER_STATUS_PROVEN`
only if `m` is known to be irreducible over $F_n$.
:::

::: function
int gr_tower_adjoin_qqbar(gr_tower_t T, const qqbar_t x, const char * name)

Adjoins the algebraic number `x` using its minimal polynomial over
$\mathbb{Q}$, which may factor over a nontrivial tower (the step is
then dynamic). The base field must be $\mathbb{Q}$.
:::

::: function
int gr_tower_adjoin_root_ui(gr_tower_t T, gr_srcptr x, ulong n, const char * name)
int gr_tower_adjoin_sqrt(gr_tower_t T, gr_srcptr x, const char * name)

Adjoins the principal `n`-th root (square root) of the element `x` of
the top field. The generator records its definition (`def_kind`
`GR_TOWER_ROOT`, `def_param` $n$, and `x` as its argument), which
the castle field uses to identify repeated roots.
:::

::: function
int gr_tower_adjoin_root_of_unity(gr_tower_t T, ulong n, const char * name)
int gr_tower_adjoin_root_fmpz(gr_tower_t T, const fmpz_t p, ulong n, const char * name)

Adjoins the root of unity $\exp(2 \pi i / n)$ (`def_kind`
`GR_TOWER_ROOT_OF_UNITY`, `def_param` $n$) or the principal
`n`-th root of the positive integer `p` (`GR_TOWER_ROOT` with a
constant argument; `p` may be negative, as for $\sqrt{-163}$, the
root of $x^n - p$), using the minimal polynomial over $\mathbb{Q}$.
Such `structured` generators are identified by their definitions
when towers are merged: a root of prime power order $\ell^e$ is a
power of a generator of the same kind of order $\ell^f$, $f \ge e$,
and roots of composite orders are products of roots of prime power
orders ($\zeta_n = \prod \zeta_{q_i}^{c_i}$ with
$c_i = (n/q_i)^{-1} \bmod q_i$, and likewise for radicals up to a
power of `p`). Roots of prime power orders of different primes
generate linearly disjoint fields, so the representation stays
canonical (irreducible moduli) and sparse. When a merge meets a
generator of order $\ell^f$ with $f < e$, a generator of order
$\ell^e$ is adjoined in front of the tower and the old one is
re-expressed as its power (a linear modulus). Generators which end
up after a generator inserted in front of them are defined over a
larger field, so their proofs of irreducibility are dropped unless
a structural rule covers the new prefix: a root of unity of prime
power order $\ell^e$ stays irreducible over roots of unity of orders
prime to $\ell$, and principal roots of positive integers with
pairwise coprime radicands (or equal radicands with coprime orders),
each irreducible over $\mathbb{Q}$, generate a field of degree the
product of the orders, also in the presence of $i$ (Mordell's theorem
on real radicals, which generalizes Besicovitch's theorem on roots of
distinct primes); an $n$-th root of a prime $p$ stays irreducible
over roots of unity of orders prime to $p$ and roots of integers
$a$ of orders $m$ with $p \nmid a m$ (the field is unramified at $p$,
and $x^n - p$ is Eisenstein there). The last rule also lets merges
skip the search for such a root in the target tower whatever the
statuses of its steps ($\sqrt{11}$ over
$\mathbb{Q}(\zeta_{16}, \sqrt{2}, \ldots, \sqrt{7})$ in the
`examples/dft.c -tower -input 1` benchmark, where this skip made
the computation an order of magnitude faster). A generator whose proof was dropped is then looked
for in the field generated by the new prefix (lattice reduction, for
prefixes of degree at most `GR_TOWER_EXPRESS_DEGREE_LIMIT`): for
instance $\sqrt{2}$ is irreducible over $\mathbb{Q}(i)$, but once
$\zeta_8$ has been moved in front it becomes
$\zeta_8 - \zeta_8^3$ (a linear modulus), keeping the representation
canonical; failing that, dynamic evaluation discovers the
factorization when a zero test needs it.
Square roots of integers are also found through Gauss sums in the
roots of unity of a tower: $\sqrt{p^*} = \sum_a (a/p)\, \zeta_p^a$
for odd primes `p` ($p^* = \pm p \equiv 1 \bmod 4$),
$\sqrt{2} = \zeta_8 + \zeta_8^{-1}$ and $\sqrt{-1} = \zeta_4$. A
merge maps $\sqrt{c}$ to such an expression when the roots of
unity are present (except in real views, which keep real
generators), and the exact zero test turns a dynamic step
$x^2 - c$ which comes after them into a step of degree one (for
instance $\sqrt{3}$ after $\zeta_3$ and $i$), instead of inverting
in a ring which is not a field. Before the modular proofs, the zero
test also applies two structural rules which only look at the
values of the generators before a step, whatever their statuses:
$x^n - p$ (`p` prime, unsplit) over roots of unity and roots of
integers generating a field unramified at `p` (Eisenstein), and
$\Phi_{\ell^e}$ over a field unramified at $\ell$ (linear
disjointness from the totally ramified $\mathbb{Q}(\zeta_{\ell^e})$;
$i$ over $\mathbb{Q}(\zeta_{111}, 3^{1/37})$, which modular proofs
miss since primes $\equiv 1 \bmod 4$ split $x^2 + 1$).
:::

::: function
int gr_tower_adjoin_pi(gr_tower_t T, const char * name)
int gr_tower_adjoin_exp(gr_tower_t T, gr_srcptr u, const char * name)
int gr_tower_adjoin_log(gr_tower_t T, gr_srcptr u, const char * name)
int gr_tower_adjoin_exp_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
int gr_tower_adjoin_log_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)

Adjoins a transcendental generator $\pi$, $\exp(u)$ or $\log(u)$,
where `u` is an element of the top field, given in the nested
representation or as a flat element in the given (current or
superseded) context of the tower's flat machinery. The base field
changes, so the nested contexts of the tower are rebuilt: elements
in the nested representation held by the caller (`u` among them)
are no longer elements of `gr_tower_field`, but they remain
valid elements of the superseded contexts, which the tower keeps
until `gr_tower_clear` (so they can be used with those
contexts and cleared after the call), while flat elements remain
valid elements of the tower. Returns `GR_DOMAIN` for
$\exp(0)$, $\log(0)$ and $\log(1)$, and `GR_UNABLE` if it cannot
be decided whether `u` is such a value. No other relation is searched
for at this point: the generator is adjoined as conjecturally
transcendental, and relations are discovered by zero tests.
:::

::: function
int gr_tower_adjoin_free(gr_tower_t T, const char * name)

Adjoins a free transcendental generator (kind `GR_TOWER_FREE`): a
formal variable, algebraically independent of all the generators by
definition (status `GR_TOWER_STATUS_INDEPENDENT`), without a
numerical value (its enclosure is indeterminate, and
`gr_tower_get_acb` fails on elements involving it) and taking
part in no relation search. The top field is then a rational
function field over the previous one, with the complete arithmetic
and zero test of the flat representation; algebraic steps over it
need an enclosure to select their root, so they can be adjoined
only through generators with values.
:::

::: function
int gr_tower_adjoin_tan(gr_tower_t T, gr_srcptr u, const char * name)
int gr_tower_adjoin_atan(gr_tower_t T, gr_srcptr u, const char * name)
int gr_tower_adjoin_tan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)
int gr_tower_adjoin_atan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)

Adjoins the real transcendental generator $\tan(u)$ or $\arctan(u)$
(kinds `GR_TOWER_TAN` and `GR_TOWER_ATAN`) for a real `u`, which
the caller must ensure ($\tan(u)$ also requires $u \not\equiv \pi/2
\pmod \pi$); $u = 0$ gives `GR_DOMAIN`, as does a `u` which is
numerically not real. Their enclosures are real (with an exactly zero
imaginary part), so that the generators count as real elsewhere (the
real forms of the lazy real field). Their relations with each other,
with $\pi$ and with the exponentials and logarithms are found by
Richardson's algorithm on their `angles` (below).
:::

::: function
int gr_tower_adjoin_special(gr_tower_t T, int kind, slong param, gr_srcptr u, const char * name)
int gr_tower_adjoin_special_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name)

Adjoins the special function value of the given kind at `u` as a
transcendental generator: `GR_TOWER_GAMMA` ($\Gamma(u)$),
`GR_TOWER_ERF` ($\operatorname{erf}(u)$ for `param` = 0,
$\operatorname{erfi}(u)$ for `param` = 1), `GR_TOWER_LAMBERTW`
($W_k(u)$, $k$ = `param`), `GR_TOWER_POLYGAMMA` ($\psi^{(m)}(u)$,
$m$ = `param`), `GR_TOWER_POLYLOG` ($\operatorname{Li}_s(u)$,
$s$ = `param`), `GR_TOWER_ZETA` ($\zeta(u)$), `GR_TOWER_ELLIPTIC_K`
and `GR_TOWER_ELLIPTIC_E` ($K(u)$, $E(u)$ with the parameter
$m = u$), or `GR_TOWER_CONSTANT` (no argument; `param` =
`GR_TOWER_CONST_EULER` or `GR_TOWER_CONST_CATALAN`). The status is
`GR_TOWER_STATUS_CONJECTURAL` (`GR_TOWER_STATUS_SCHANUEL` for
Lambert W). `GR_TOWER_HYPGEOM` (${}_pF_q$) has several arguments
and is adjoined with `gr_tower_adjoin_special_multi_flat`. The caller is responsible for canonical arguments (see
above): a value which is algebraic or expressible through other
generators must not be adjoined. Returns `GR_DOMAIN` for $u = 0$
where the value is trivial or a pole ($\operatorname{erf}$, $W$,
$\operatorname{Li}_s$, $\Gamma$, $\psi^{(m)}$), and `GR_UNABLE` if
the value cannot be evaluated numerically. The generators are
identified across towers (merges, maps) by their kind, parameter and
argument.
:::

::: function
int gr_tower_adjoin_special_multi_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_struct * u, slong nargs, const fmpz_mpoly_ctx_t mctx, const char * name)
slong _gr_tower_special_num_args(int kind, slong param)

Adjoins the value of a special function of several arguments: `u`
is a vector of `nargs` flat elements, of which `u[0]` becomes the
argument of the generator and the others its additional arguments.
For `GR_TOWER_HYPGEOM`, `param` is `GR_TOWER_HYPGEOM_PARAM(p, q)`
and $u = (z, a_1, \ldots, a_p, b_1, \ldots, b_q)$ for the value
${}_pF_q(a; b; z)$ (not regularized; numerically through
`acb_hypgeom_2f1` and friends, and for $p = q + 1 \ge 3$ only
inside the unit disk). `nargs` must be
`_gr_tower_special_num_args(kind, param)`. Generators are
identified across towers by their kind, parameter and all their
arguments.
:::

## Numerical evaluation {#numerical-evaluation}


::: function
int gr_tower_trans_get_acb(acb_t res, gr_tower_t T, slong j, slong prec)

Enclosure of the transcendental generator $t_j$.
:::

::: function
int gr_tower_step_get_acb(acb_t res, gr_tower_t T, slong k, slong prec)

Sets `res` to an enclosure of the generator $a_k$, refined to
`prec` bits of relative accuracy using interval Newton iteration.
The enclosure is cached in the tower.
:::

::: function
int gr_tower_get_acb_at(acb_t res, gr_srcptr x, slong k, slong prec, gr_tower_t T)
int gr_tower_get_acb(acb_t res, gr_srcptr x, slong prec, gr_tower_t T)

Sets `res` to an enclosure of the element `x` of $F_k$ (of $F_n$).
:::

## Dynamic refinement and complete operations {#dynamic-refinement-and-complete-operatio}


::: function
int gr_tower_refine(gr_tower_t T)

Processes all zero divisors recorded by the quotient ring contexts,
refining the corresponding minimal polynomials. Returns nonzero if the
tower changed.
:::

::: function
int gr_tower_refine_step(gr_tower_t T, slong k, const gr_poly_t g)

Given a monic proper factor `g` of $m_k$, replaces $m_k$ by `g` or
$m_k / g$, whichever vanishes at the enclosure of $a_k$.
:::

::: function
truth_t gr_tower_is_zero(gr_srcptr x, gr_tower_t T)
truth_t gr_tower_equal(gr_srcptr x, gr_srcptr y, gr_tower_t T)
int gr_tower_inv(gr_ptr res, gr_srcptr x, gr_tower_t T)
int gr_tower_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_tower_t T)
truth_t gr_tower_is_zero_at(gr_srcptr x, slong k, gr_tower_t T)
truth_t gr_tower_equal_at(gr_srcptr x, gr_srcptr y, slong k, gr_tower_t T)
int gr_tower_inv_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
int gr_tower_div_at(gr_ptr res, gr_srcptr x, gr_srcptr y, slong k, gr_tower_t T)

Zero test, equality test, inversion and division of elements of the
top field (of $F_k$ for the `_at` versions). For algebraic towers these are complete: they refine the
tower as needed and always return a definite answer (`GR_DOMAIN`
for division by zero). An element which is nonzero as a rational
function of conjecturally transcendental generators is decided by
Richardson's algorithm on a copy of the tower (up to the precision
limit, after which the answer is `T_UNKNOWN` / `GR_UNABLE`).
:::

## Modular irreducibility proofs {#modular-irreducibility-proofs}


A step adjoined without a proof of irreducibility (a `dynamic` step,
`GR_TOWER_STATUS_DYNAMIC`) is correct by dynamic evaluation, but
a proof (`GR_TOWER_STATUS_PROVEN`) makes some decisions cheaper or
definite: a nonzero element of a tower whose steps are all proven is
invertible without computing the inverse, a non-rational element of such
a tower is definitely non-rational, and structural rules apply to proven
generators. For towers over $\mathbb{Q}$, proofs are attempted at
`places` of the tower: a chain of simple roots $\alpha_1, \ldots,
\alpha_k$ of the moduli in a finite field $\mathbb{F}_{\ell^m}$
(each $\alpha_j$ a root of $m_j(\alpha_1, \ldots, \alpha_{j-1}, x)$)
defines a homomorphism from the ring $\mathbb{Z}[a_1, \ldots, a_k]$
(with the leading coefficients inverted) to $\mathbb{F}_{\ell^m}$,
whose local ring at the kernel is a discrete valuation ring since the
Jacobian of the triangular system is nonsingular there. A monic
factorization of $m_{k+1}$ over the tower would therefore reduce to a
factorization over the residue field, so $m_{k+1}$ is irreducible over
the tower if its image at the place is irreducible; and a $p$-th root
in the tower of an element $a$ would reduce to a $p$-th root of the
image of $a$, so $a$ is not a $p$-th power in the tower if its image is
not one in the residue field. By Capelli's theorem, $X^n - a$ is
irreducible over the tower if $a$ is not a $p$-th power for any prime
$p \mid n$ and, when $4 \mid n$, $-a/4$ is not a fourth power; this
proves the steps adjoined by `gr_tower_adjoin_root_ui` (the
principal $n$-th roots of the castle field) in almost all cases where
they are irreducible. The proofs are relative to the moduli of the
lower steps, as the status of a step is.

With transcendental generators, the place also gives them values in
$\mathbb{F}_{\ell}$. A modulus then has a root at the place only for
some of the values: a radical $X^n - a(t)$ for about one value of $t$
in $n$, so that over the radicals of $t$ and $1 - t$ of orders 4 (the
theta constants $\lambda^{1/4}$, $(1-\lambda)^{1/4}$, say) one random
choice gives a place of degree one about once in sixteen. When a
modulus has no simple root, the values of the generators which first
occur from the earliest step sharing a generator with it are drawn
again (the chain of roots being taken back to that step), up to 48
times per place. In the merges of the towers of the modular functions,
whose lattice searches for roots of radicals not in the target tower
otherwise took a third of the time, the absence of a root is then
usually shown at a place.

::: function
int gr_tower_prove_step_modular(gr_tower_t T, slong k, slong tries)

Tries to prove that $m_k$ is irreducible over $F_{k-1}$ at places of
the tower, with an effort proportional to `tries` (a few times
`tries` primes with places over the prime fields, then a few primes
with places over small extension fields). Sets the status of the
step to proven and returns 1 on success. Reducible moduli are never
proven; irreducible ones may fail to be (a modulus which is
reducible modulo every prime, like $x^4 + 1$, or a chain of roots
which exists only at primes of low density). The steps below $k$
must be proven (except for a modulus which is Eisenstein in a
transcendental generator): over a quotient ring which is not a field,
a place sees only one of its factors, and the irreducibility over
that factor says nothing about the one the tower is refined to later.
:::

::: function
int gr_tower_prove_modular(gr_tower_t T, slong tries)

Attempts a proof for every dynamic step of the tower which has not
been tried with the current moduli (recorded in `proof_version` of
the generator; the attempt is repeated after a refinement), in
order, the steps above a step which remains dynamic being left for
a later attempt. Returns 1
if all generators are proven afterwards. The zero test calls this
before resorting to an inversion, and the castle field before leaving
a rationality test undecided; the option `GR_TOWER_OPT_MODULAR_TRIES`
is the number of primes. After a failed modular attempt on a step
of a tower over $\mathbb{Q}$, Trager's method (below) is tried,
which decides the question for small degrees: an irreducible $x^4 + 1$ is proven, and a modulus
which is reducible over the tower is refined.
:::

::: function
int gr_tower_prove_step_trager(gr_tower_t T, slong k, slong degree_limit)

Decides the irreducibility of $m_k$ over $F_{k-1}$ by Trager's
method (see the factorization below), for a tower whose steps below
$k$ are proven, over $\mathbb{Q}$ or over a rational function field
$\mathbb{Q}(t_1, \ldots, t_r)$, provided that $\deg(m_k) [F_{k-1} :
F_0] \le$ `degree_limit`. If $m_k$ is irreducible, the step is marked
proven; otherwise the tower is refined with its irreducible factor
vanishing at the generator (the factors which do not vanish being
divided out one at a time), and the refined step is marked proven.
Linear steps are proven directly. Returns 1 if the step is proven
(after a refinement, if any), 0 if the method does not apply, the
degree limit is exceeded, the modulus is not squarefree, or no
shift making the norm squarefree was found. Unlike the modular
proofs this is a decision procedure, but the norm has degree
$\deg(m_k) [F_{k-1} : F_0]$ and is factored over $\mathbb{Z}$ (or
as a multivariate polynomial over $\mathbb{Z}$), so it is limited to
small absolute degrees: the automatic proofs of
`gr_tower_prove_modular` use it for towers over $\mathbb{Q}$
only, with the limit of the option
`GR_TOWER_OPT_TRAGER_DEGREE_LIMIT` (48 by default).
:::

::: function
int gr_tower_binomial_irreducible_modular(gr_tower_t T, const fmpz_mpoly_q_t a, const fmpz_mpoly_ctx_t actx, ulong n, slong tries)

Proves $X^n - a$ irreducible over the field of the tower by the
criterion above, for the element `a` in the flat representation
(in the context `actx` of the tower's flat machinery). Returns 1 on
success; the steps of the tower must be proven.

::: details
(Internally, the same test reports a prime $p \mid n$ for which $a$
was a $p$-th power at many informative places, the evidence used
by `GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT`.)
:::
:::

::: function
int gr_tower_poly_no_roots_modular(const gr_poly_t g, gr_tower_t T, slong tries)

Returns 1 if the polynomial `g` over the top field of `T` is shown to
have no root in that field, and 0 if this is not known. About `tries`
places over prime fields are tried (the transcendental generators
being given random values): if `g` had a root $r$, then at a place
where the coefficients of `g` are integral and its leading
coefficient does not vanish, $\phi(r)$ would be a root of the
reduction of `g` (the local ring at a place where the moduli have
simple roots is regular, hence integrally closed). The statement is
relative to the moduli of the steps being irreducible; the function
serves to skip lattice searches for roots (when towers are merged and
in root finding), where a wrong answer only costs a dynamic step.
:::

The effort of the automatic uses is set by the options
`GR_TOWER_OPT_NO_ROOTS_TRIES` (places tried by
`gr_tower_poly_no_roots_modular`), `GR_TOWER_OPT_MODULAR_TRIES`
(the budget of a modular proof: a schedule of places over prime fields,
then over extension fields, of a size proportional to it) and
`GR_TOWER_OPT_TRAGER_DEGREE_LIMIT` (the largest absolute degree
$\deg(m_k) [F_{k-1} : \mathbb{Q}]$ for which the automatic proofs of
towers over $\mathbb{Q}$ fall back to Trager's method).

::: macro
GR_TOWER_DEFERRED_IDEAL_SIZE

Number of base coefficients above which the flat conversion of a
modulus is deferred (see the flat representation below).
:::

## Factorization of polynomials {#factorization-of-polynomials}


Polynomials over a field $F_k$ of a tower are factored by Trager's
method. The `norm` $N(h) \in F_0[x]$ of $h \in F_k[x]$ is the product of
the conjugates of $h$ over $F_0$, computed as the iterated resultant of
$h$ with the moduli $m_k, \ldots, m_1$ (the coefficients of $h$ viewed
as polynomials in $a_k$ over $F_{k-1}[x]$, and so on down); it has
degree $\deg(h) [F_k : F_0]$. For squarefree $h$, the norm of
$h(x - \theta)$ is squarefree for all but finitely many shifts
$\theta = \sum_j s_j a_j$; for such a shift, the irreducible factors of
$h$ are $\gcd(h(x - \theta), N_i)(x + \theta)$ for the irreducible
factors $N_i$ of the norm over $F_0$. Over $F_0 = \mathbb{Q}$ the norm
is factored by `fmpz_poly_factor`; over $F_0 = \mathbb{Q}(t_1,
\ldots, t_r)$ its denominators are cleared and it is factored as a
polynomial in $t_1, \ldots, t_r, x$ by `fmpz_mpoly_factor`.
Squarefree decomposition is by Yun's algorithm, and the gcds are
computed by the subresultant PRS. The factorization is exact: it
requires the steps up to $k$ to be proven (which is attempted first:
modular proofs, then Trager's method for the steps themselves), so
that $F_k$ is a field.

The option `GR_TOWER_OPT_FACTOR_DEGREE_LIMIT` bounds $\deg(h)
[F_k : F_0]$ for the squarefree parts $h$ factored by
`gr_tower_poly_factor`.

::: function
int gr_tower_poly_factor(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T)
int gr_tower_poly_factor_limit(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, slong degree_limit, gr_tower_t T)

Factors the nonzero polynomial `f` over $F_k$ as $c \prod_i
f_i^{e_i}$ with $c \in F_k$ (the leading coefficient) and distinct
monic irreducible $f_i$, stored in `fac` (a vector of polynomials
over $F_k$) and `mult`. Returns `GR_UNABLE` if a step up to $k$
cannot be proven irreducible, if a squarefree part exceeds the
degree limit, if no suitable shift is found, or if the base field is
not $\mathbb{Q}$ or $\mathbb{Q}(t_1, \ldots, t_r)$; `GR_DOMAIN`
for the zero polynomial. The field contexts
(`gr_ctx_init_tower_field`) implement `gr_factor` for
polynomials (`GR_METHOD_POLY_FACTOR`) by this function.
:::

::: function
int gr_tower_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T)

The roots of `f` in $F_k$ with their multiplicities (the linear
factors of the factorization). The field contexts implement
`gr_poly_roots` by this function.
:::

## Absolute representation {#absolute-representation}


::: function
int gr_tower_get_coeffs_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
int gr_tower_set_coeffs_at(gr_ptr res, gr_srcptr coeffs, slong k, gr_tower_t T)

Conversion between an element of $F_k$ and its vector of $[F_k : F_0]$
coordinates in the monomial basis
$a_1^{i_1} \cdots a_k^{i_k}$, ordered with $i_1$ varying fastest.
:::

::: function
int gr_tower_multiplication_matrix(gr_mat_t res, gr_srcptr x, gr_tower_t T)
int gr_tower_charpoly(gr_poly_t res, gr_srcptr x, gr_tower_t T)

The matrix of multiplication by `x` with respect to the monomial
basis of $F_n$ over $F_0$, and its characteristic polynomial.
:::

::: function
int gr_tower_get_fmpz_poly_minpoly(fmpz_poly_t res, gr_srcptr x, gr_tower_t T)
int gr_tower_get_qqbar(qqbar_t res, gr_srcptr x, gr_tower_t T)
int gr_tower_promote(gr_ptr res, gr_srcptr x, slong j, slong k, gr_tower_t T)

The minimal polynomial of the number represented by `x` over
$\mathbb{Q}$, respectively the number itself as a `qqbar_t`.
The base field must be $\mathbb{Q}$. The annihilating polynomial of
least degree is computed by incremental elimination on the powers of
`x`; if the tower is not yet fully refined this polynomial may be
reducible, in which case the correct factor is selected numerically.
`gr_tower_promote` maps an element of $F_j$ to $F_k$, $j \le k$.
:::

## Expressing numbers in a tower {#expressing-numbers-in-a-tower}


::: function
int gr_tower_sqrt_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T)
int gr_tower_sqrt(gr_ptr res, gr_srcptr x, gr_tower_t T)

Sets `res` to the principal square root of `x` if it lies in $F_k$
($F_n$): with $F_k = F_{k-1}(a)$, $a^2 = r$ after completing the
square, $x = u + v a$ is a square iff $u^2 - r v^2$ is a square $s$ in
$F_{k-1}$ and $(u \pm s)/2$ is a square $p^2$ there, and then
$\sqrt{x} = p + v a / (2p)$; this recurses to the base field, where
rational numbers and rational functions are tested directly. No
lattice reduction is involved. Returns `GR_DOMAIN` if `x` is not a
square in $F_k$, and `GR_UNABLE` if a step below $k$ has degree
other than 2.
:::

Above the degree `GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT` (an option),
lattice-based expression searches are skipped when adjoining roots and
merging towers; roots are then adjoined without searching, dynamic
evaluation keeping zero tests complete.

::: function
int gr_tower_express(gr_ptr res, const gr_poly_t q, const acb_t z, gr_tower_t T)
int gr_tower_express_limit(gr_ptr res, const gr_poly_t q, const acb_t z, slong prec_limit, gr_tower_t T)

Given a polynomial `q` over the top field and an enclosure `z`
isolating a single root of `q`, attempts to find the element of the
top field which is that root. A candidate is found by an integer
relation search (LLL) over the monomial basis, then verified exactly
($q(x) = 0$, using the complete zero test) and identified numerically
with the root in `z`. Returns `GR_SUCCESS` if found and `GR_UNABLE`
if no representation was found within the precision limit (which
does not prove that the number is not in the field).
:::

::: function
int gr_tower_express_qqbar(gr_ptr res, const qqbar_t x, gr_tower_t T)

Expresses the algebraic number `x` in the tower (base field $\mathbb{Q}$).
Returns `GR_DOMAIN` if the degree of `x` does not divide the degree
of the tower, which proves that `x` is not in the field.
:::

## Maps between towers {#maps-between-towers}


::: type
gr_tower_map_struct
gr_tower_map_t

A homomorphism from a source tower to a target tower over the same
constant field, given by the images of the first `length` generators
of the source (in definition order) as flat elements of the target.
Flat elements survive extensions of the target, so a map remains valid
when the target grows (its images may need conversion to the current
polynomial context, see `gr_tower_map_sync`).
:::

::: function
void gr_tower_map_init(gr_tower_map_t map, gr_tower_t source, gr_tower_t target)
void gr_tower_map_clear(gr_tower_map_t map)
void gr_tower_map_set(gr_tower_map_t res, const gr_tower_map_t map)
int gr_tower_map_fit_length(gr_tower_map_t map, slong len)
void gr_tower_map_sync(gr_tower_map_t map)
fmpz_mpoly_q_struct * gr_tower_map_image(gr_tower_map_t map, slong d)
int gr_tower_map_set_image(gr_tower_map_t map, slong d, gr_srcptr x, slong k)
int gr_tower_map_set_image_flat(gr_tower_map_t map, slong d, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t mctx)

A newly initialized map is defined on no generators; images are set
by the functions below or by `set_image` (from an element of $F_k$
of the target) and `set_image_flat`. The image of the generator with
definition order `d` is an element of the context `map->mctx`,
which `sync` brings to the current context of the target.
:::

::: function
int gr_tower_map_set_inclusion(gr_tower_map_t map)

Sets the map to the inclusion of the source into a target whose first
generators (in definition order) are copies of the source's generators.
:::

::: function
int gr_tower_map_apply_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_map_t map)
int gr_tower_map_apply(gr_ptr res, gr_srcptr x, gr_tower_map_t map)
int gr_tower_map_apply_flat(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_tower_map_t map)
int gr_tower_map_apply_poly(gr_poly_t res, const gr_poly_t f, gr_tower_map_t map)

Applies the map to an element of $F_k$ of the source, of its top
field, to a flat element of the source, or to a polynomial over the
top field of the source. The result is an element of (a polynomial
over) the top field of the target, respectively a flat element in
the context `map->mctx` (after the synchronization performed by
this function). Returns `GR_DOMAIN` if `x` involves generators on
which the map is not defined.
:::

::: function
int gr_tower_absorb(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, int flags)
int gr_tower_absorb_prefix(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, slong p, int flags)

Extends `U` (the target of `map`) by the generators of `B` with
definition order from `length` (of the map) to `p` - 1 (or all), and
extends the map accordingly. A generator with the same definition id
as a generator of `U` is identified with it. Otherwise, an algebraic
generator whose minimal polynomial has become linear is substituted;
a minimal polynomial whose image has coefficients in the base field
of `U` is factored there (the factor vanishing at the generator
replaces it); if `flags` contains `GR_TOWER_MERGE_EXPRESS`, a
low-effort expression search is tried, unless a place of `U` shows
that the polynomial has no root in `U`
(`gr_tower_poly_no_roots_modular`); otherwise the generator is
adjoined (as a dynamic step). A transcendental generator is adjoined
with the image of its argument.
:::

::: function
int gr_tower_merge(gr_tower_t U, gr_tower_map_t mapA, gr_tower_map_t mapB, gr_tower_t A, gr_tower_t B, int flags)

Sets `U` (which must be initialized) to a tower containing `A` and `B`,
with maps (initialized by this function) from `A` and `B` into `U`.
:::

::: function
int gr_tower_eliminate(gr_tower_t U, gr_tower_map_t map, gr_tower_t T)

Sets `U` to a tower isomorphic to `T` with the steps of degree 1
removed, with the isomorphism as `map`.
:::

## Tower field context {#tower-field-context}


A fixed tower is built step by step and its top field used as a `gr`
context; a typical use is a field given in advance, in which many
computations take place (linear algebra over $\mathbb{Q}(\sqrt 2,
\sqrt 3)$, say), where the castle field's bookkeeping of towers is not
wanted:

```c
gr_ctx_t QQ, K;
gr_tower_t T;
gr_ptr x;

gr_ctx_init_fmpq(QQ);
gr_tower_init(T, QQ);

/* adjoin sqrt(2) and sqrt(3): the arguments are elements of the
   current top field */
{
    gr_ctx_struct * top = gr_tower_field(T);
    GR_TMP_INIT(x, top);
    GR_MUST_SUCCEED(gr_set_ui(x, 2, top));
    GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
    GR_TMP_CLEAR(x, top);
}
{
    gr_ctx_struct * top = gr_tower_field(T);
    GR_TMP_INIT(x, top);
    GR_MUST_SUCCEED(gr_set_ui(x, 3, top));
    GR_MUST_SUCCEED(gr_tower_adjoin_sqrt(T, x, NULL));
    GR_TMP_CLEAR(x, top);
}

/* the field Q(sqrt 2, sqrt 3) as a gr context: complete zero tests,
   division, generators a1, a2 */
gr_ctx_init_tower_field(K, T);
{
    gr_vec_t gens;
    gr_vec_init(gens, 0, K);
    GR_MUST_SUCCEED(gr_gens(gens, K));
    GR_TMP_INIT(x, K);
    GR_MUST_SUCCEED(gr_add(x, gr_vec_entry_ptr(gens, 0, K), gr_vec_entry_ptr(gens, 1, K), K));
    GR_MUST_SUCCEED(gr_inv(x, x, K));     /* 1/(sqrt 2 + sqrt 3) = sqrt 3 - sqrt 2 */
    gr_println(x, K);
    GR_TMP_CLEAR(x, K);
    gr_vec_clear(gens, K);
}
gr_ctx_clear(K);
gr_tower_clear(T);
```

Adjoining a step whose polynomial is only known to vanish at the root
(`gr_tower_adjoin_algebraic` with `GR_TOWER_STATUS_DYNAMIC`) is
allowed: the field context stays a field, and the tower is refined when
an operation exposes a factorization. Transcendental generators
(`gr_tower_adjoin_pi`, `gr_tower_adjoin_exp`,
`gr_tower_adjoin_log`) rebuild the tower over a rational function
field, which invalidates the contexts created before (see below). The
Python wrapper's `gr_tower` class follows this pattern with
`adjoin_*` methods and `field()`.

::: function
void gr_ctx_init_tower_field(gr_ctx_t ctx, gr_tower_t T)
int gr_tower_field_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx)

Initializes `ctx` to a `gr` context for the field $F_n$, $n$ being
the number of algebraic steps of `T` at that moment. Elements
have the same representation as elements of `gr_tower_field`,
and arithmetic is delegated to that context, but zero tests, equality
tests and inversions are complete (dynamically refining `T`).
The context refers to `T`, which must outlive it. It remains valid
when further algebraic steps are appended to `T` afterwards: its
elements then lie in the subfield $F_n$ of the new top field, and
`gr_set_other` converts them into a context made for the
larger field (and, generally, between contexts of the same tower
from a lower to a higher level). Adjoining a transcendental
generator rebuilds the tower over a new base field, which makes the
contexts created before `stale`: their operations fail with
`GR_UNABLE`, while clearing their elements remains valid (the
tower keeps the superseded nested contexts until it is cleared). Conversion to an `acb`
context (`gr_set_other`) is supported.
:::

## Flat representation {#flat-representation}


The nested representation is the structural one (it carries the
triangular set, dynamic refinement and the certificates), but for
arithmetic it is often better to represent elements as multivariate
rational functions in the generators, reduced modulo the triangular set
with denominators cleared. Since the leading monomials $c_k x_k^{d_k}$ of
the triangular set are pairwise coprime, it is a Gröbner basis for the
lex order with $x_n > \cdots > x_1$, so reduced numerators are canonical
up to scaling. In this representation inversion is free (swap numerator
and denominator), sparse expressions stay sparse and multiplication uses
`fmpz_mpoly` arithmetic, at the cost of non-canonical denominators
(equality tests cross-multiply). Both representations coexist; conversion
in either direction goes through the coordinate vector.

The flat `ideal` (the moduli converted to polynomials in all the
generators) is built from the nested moduli. Over $\mathbb{Q}$, the
conversion of a modulus with non-constant coefficients and more than
`GR_TOWER_DEFERRED_IDEAL_SIZE` base coefficients is deferred until
a reduction actually needs it: the reduction of an element then proceeds
by one division per step from the top down (a reduction by the modulus
of step $k$ raises only degrees in the variables of lower steps, so the
steps needed are known as the reduction goes), converting the moduli on
first use. The moduli of a splitting tower have binomial size in the
middle, and symmetric functions of the roots -- which reduce through the
last, linear steps -- never touch them. A polynomial whose degree in
every variable is already below that of the corresponding modulus is
reduced (the reduction by a modulus of step $k$ raises only degrees in
the variables of lower steps, never that of $a_k$ above $\deg m_k$),
so `gr_tower_flat_reduce` returns at once for such input; this
is the common case after a multiplication of elements which are linear
in the generators, and avoids a division per step.

::: type
gr_tower_flat_struct
gr_tower_flat_t

The flat machinery for one tower: a polynomial context in which the
variable `GR_TOWER_FLAT_VAR_D(F, d)` = `cap - 1 - d` represents the
generator with definition order `d` (`GR_TOWER_FLAT_VAR(F, k)` and
`GR_TOWER_FLAT_TVAR(F, j)` are the variables of $a_k$ and $t_j$),
so that later generators are more significant in lex order and the
triangular set is a Groebner basis whatever the kinds of the
generators are; the tower can grow up to the capacity without
renumbering. The reduced triangular set, and superseded polynomial
contexts (kept alive, each with the list of generators its variables
stood for, so that elements created with them remain convertible even
after reorderings). Every tower has an embedded flat machinery,
`gr_tower_flat`. Reduction is by pseudo-division from the top
generator down; the leading coefficients of the moduli are polynomials
in the transcendental variables only, so the fraction is multiplied
through by them.
:::

::: function
void gr_tower_flat_init(gr_tower_flat_t F, gr_tower_t T, slong cap)
void gr_tower_flat_clear(gr_tower_flat_t F)
int gr_tower_flat_ensure(gr_tower_flat_t F)
gr_tower_flat_struct * gr_tower_flat(gr_tower_t T)

Initialization and maintenance: `ensure` replaces the polynomial
context if the tower has grown beyond the capacity or has been
reordered (returning 1) and rebuilds the triangular set if the tower
has been refined or restructured.
:::

::: function
void gr_tower_flat_convert(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t old_mctx, gr_tower_flat_t F)
int gr_tower_flat_reduce(fmpz_mpoly_q_t x, gr_tower_flat_t F)
int gr_tower_flat_set_nested_at(fmpz_mpoly_q_t res, gr_srcptr x, slong k, gr_tower_flat_t F)
int gr_tower_flat_poly_get_nested_at(gr_ptr res, const fmpz_mpoly_t f, slong k, gr_tower_flat_t F)
int gr_tower_flat_get_nested_at(gr_ptr res, const fmpz_mpoly_q_t x, slong k, gr_tower_flat_t F)
int gr_tower_flat_get_acb(acb_t res, const fmpz_mpoly_q_t x, slong prec, gr_tower_flat_t F)
slong gr_tower_flat_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
slong gr_tower_flat_alg_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
slong gr_tower_flat_def_order(const fmpz_mpoly_q_t x, gr_tower_flat_t F)
int gr_tower_flat_compose(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, fmpz_mpoly_q_struct ** imgs, gr_tower_flat_t F)
truth_t gr_tower_flat_num_is_zero(fmpz_mpoly_q_t x, gr_tower_flat_t F)
int gr_tower_flat_rationalize(fmpz_mpoly_q_t x, gr_tower_flat_t F)
int gr_tower_flat_has_alg_var(const fmpz_mpoly_t p, gr_tower_flat_t F)
slong gr_tower_flat_max_step(const fmpz_mpoly_t p, gr_tower_flat_t F)

Conversion from a superseded context, reduction, conversions to and
from the nested representation of $F_k$, numerical evaluation, the
length of the shortest prefix containing all generators occurring in
`x` (respectively the highest algebraic generator, and the highest
definition order), substitution of elements of `F`'s tower for the
variables of another context, and the complete zero test of the
numerator (numerics, then the nested test, which may refine the tower).
`gr_tower_flat_rationalize` makes the denominator of `x` free of
algebraic generators (`GR_UNABLE` beyond the size limit
`GR_TOWER_OPT_RATIONALIZE_LIMIT`), `gr_tower_flat_has_alg_var` tells
whether the polynomial `p` involves an algebraic generator, and
`gr_tower_flat_max_step` gives the largest index $k$ of an algebraic
generator $a_k$ occurring in it (0 if none).
:::

::: function
void gr_ctx_init_tower_field_flat(gr_ctx_t ctx, gr_tower_t T)
int gr_tower_flat_set_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_tower_flat_get_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_tower_field_flat_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx)

A `gr` context for the top field of a fixed tower in the flat
representation, presented as an extension of the constant field with
generators $t_1, \ldots, t_r, a_1, \ldots, a_n$, conversions to and
from `gr_ctx_init_tower_field`, and numerical evaluation.
:::

## Tuning options {#tuning-options}


The limits and policies of the algorithms are read from a table of
options (indexed by the constants below) that a tower refers to through
its member `options`: `gr_tower_default_options` for a tower made
by `gr_tower_init`, the table of the castle field for its towers
(see `gr_ca_ctx_set_option`). Towers derived from a tower
(prefixes, copies, merges) share its table. Changing an option affects
the computations that follow; the elements created before remain valid.

::: variable
const slong gr_tower_default_options[GR_TOWER_OPT_NUM_OPTIONS]

:::

::: function
void gr_tower_set_options(gr_tower_t T, const slong * options)

Makes `T` read its options from the table `options` (not copied), or
from the defaults if `options` is `NULL`.
:::

::: macro
GR_TOWER_OPTION(T, k)

The value of the option `k` for the tower `T`.
:::

::: function
const char * gr_tower_option_name(slong option)
slong gr_tower_option_find(const char * name)
slong gr_tower_option_default(slong option)
int gr_tower_option_valid(slong option, slong value)
void gr_tower_options_init(slong * options)
int gr_tower_options_set(slong * options, slong option, slong value)

The name of an option (`"prec_limit"` for
`GR_TOWER_OPT_PREC_LIMIT`; `NULL` for an invalid index), the
index of the option with the given name (-1 if none), the default
value, and whether a value is in the range of the option. A table
for `gr_tower_set_options` (`GR_TOWER_OPT_NUM_OPTIONS`
entries) is filled with the defaults by `gr_tower_options_init`
and changed by `gr_tower_options_set`, which returns
`GR_DOMAIN` for a value out of range. (The Python interface takes
the options by their names.)
:::

::: macro
GR_TOWER_OPT_VERBOSE
GR_TOWER_OPT_PRINT_FLAGS
GR_TOWER_OPT_PRINT_DIGITS
GR_TOWER_OPT_PREC_LIMIT
GR_TOWER_OPT_CERTIFY_PREC_LIMIT
GR_TOWER_OPT_NUMERIC_PREC_LIMIT
GR_TOWER_OPT_SMOOTH_LIMIT
GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT
GR_TOWER_OPT_EXPRESS_PREC
GR_TOWER_OPT_EXPRESS_ROUNDS
GR_TOWER_OPT_TRAGER_DEGREE_LIMIT
GR_TOWER_OPT_FACTOR_DEGREE_LIMIT
GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT
GR_TOWER_OPT_MODULAR_TRIES
GR_TOWER_OPT_MODULAR_TERMS_LIMIT
GR_TOWER_OPT_NO_ROOTS_TRIES
GR_TOWER_OPT_RATIONALIZE_LIMIT
GR_TOWER_OPT_RELATION_COST_LIMIT
GR_TOWER_OPT_SQRT_BUDGET
GR_TOWER_OPT_GAUSS_SUM_LIMIT
GR_TOWER_OPT_DENSE_LIMIT
GR_TOWER_OPT_INV_DENSE_DEGREE_LIMIT
GR_TOWER_OPT_INV_LINEAR_DEGREE_LIMIT
GR_TOWER_OPT_DEFERRED_IDEAL_SIZE
GR_TOWER_OPT_GROW_THRESHOLD
GR_TOWER_OPT_ROOT_OF_UNITY_ORDER_LIMIT
GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT
GR_TOWER_OPT_REALIFY_ORDER_LIMIT
GR_TOWER_OPT_TRIG_PI_LIMIT
GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT
GR_TOWER_OPT_TRIG_FORM
GR_TOWER_OPT_SPLIT_IMAGINARY
GR_TOWER_OPT_COMPOSITE_RADICALS
GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT
GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT
GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT
GR_TOWER_OPT_GAMMA_LATTICE_LIMIT
GR_TOWER_OPT_GAMMA_LINE_LIMIT
GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT
GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT
GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT
GR_TOWER_OPT_GAUSS_DIGAMMA_LIMIT
GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT
GR_TOWER_OPT_DENSE_FORM_SPARSITY
GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT
GR_TOWER_OPT_SPLIT_DEGREE_LIMIT
GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT
GR_TOWER_OPT_INV_DENSE_ALG
GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT
GR_TOWER_OPT_FACTORED_DENOMINATORS
GR_TOWER_OPT_UNLIMITED
GR_TOWER_OPT_NUM_OPTIONS

The indices of the options; `GR_TOWER_OPT_NUM_OPTIONS` is their number.
:::

The options (default values in parentheses; precisions are in bits):

* `GR_TOWER_OPT_VERBOSE` (0): with 1, the relations found
  between generators (generators becoming of degree one) are printed;
  with 2, also statistics.
* `GR_TOWER_OPT_PRINT_FLAGS` (`GR_TOWER_PRINT_SYMBOLIC | GR_TOWER_PRINT_DEFS`), `GR_TOWER_OPT_PRINT_DIGITS` (6): printing
  of castle field elements (`gr_ca_ctx_set_print`).
* `GR_TOWER_OPT_PREC_LIMIT` (256): precision up to which a cheap
  numerical separation from zero is attempted before an exact test
  (zero tests, comparisons, signs).
* `GR_TOWER_OPT_CERTIFY_PREC_LIMIT` (4096): precision limit of
  the certifications (roots, signs and real parts of generators, the
  selection among candidates by their enclosures) and of the relation
  searches of the zero test, after which the answer is unknown.
* `GR_TOWER_OPT_NUMERIC_PREC_LIMIT` (65536): when the exact
  methods leave a zero test or a sign undecided (no structure theorem
  covers the element), the numerical evaluation is pushed to this
  precision before giving up: $1 - \operatorname{erf}(100)$ is
  separated from zero at 14500 bits.
* `GR_TOWER_OPT_SMOOTH_LIMIT` (3512): number of primes for the
  trial division of the radicands of roots in castle fields, and of the
  Gaussian integers whose logarithms are decomposed over the prime
  logarithms; a cofactor beyond it is kept as one radicand (or one
  opaque logarithm).
* `GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT` (64),
  `GR_TOWER_OPT_EXPRESS_PREC` (256): lattice-based expression
  searches in fields up to this degree, at a precision of this many
  bits plus four per degree of the field
  (`GR_TOWER_MERGE_EXPRESS_PREC`).
* `GR_TOWER_OPT_EXPRESS_ROUNDS` (2): the lattice-based
  expression searches try at most this many precisions, doubling from
  about twice the degree of the field (0: up to the precision limit).
  A relation is nearly always found at the first or second precision,
  while a search which finds none runs up to the limit, its last
  rounds being the costliest: in the test suite, of 309 searches which
  succeeded, 277 did at the first precision, 31 at the second and 1 at
  the third, and the 319 which failed (92% of the time of the searches)
  never met a spurious relation. With 2, the test suite is about 12%
  faster than with 0 and the profile catalog 7%; with 1, 35% and 12%
  (all results unchanged), an expression missed being left to the zero
  tests, as in the fields where the search does not apply.
* `GR_TOWER_OPT_TRAGER_DEGREE_LIMIT` (48),
  `GR_TOWER_OPT_MODULAR_TRIES` (6),
  `GR_TOWER_OPT_MODULAR_TERMS_LIMIT` (50000),
  `GR_TOWER_OPT_NO_ROOTS_TRIES` (24): effort of the automatic
  irreducibility proofs (Trager's method up to this degree; the number
  of places tried by the modular proofs, which are not attempted for
  moduli of more than this many terms; the places tried by the modular
  proofs of the absence of roots).
* `GR_TOWER_OPT_FACTOR_DEGREE_LIMIT` (128): factorization over
  towers.
* `GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT` (64): in castle fields,
  the roots of a polynomial are found by factoring it over the tower
  of its coefficients when the degree of the tower times that of the
  polynomial is at most this (3/8 of it with transcendental generators,
  over which factoring is costlier).
* `GR_TOWER_OPT_RATIONALIZE_LIMIT` (2000): size of the
  denominators rationalized.
* `GR_TOWER_OPT_RELATION_COST_LIMIT` (200000): multiplicative
  relations costlier to verify are not verified.
* `GR_TOWER_OPT_SQRT_BUDGET` (2000): recursive calls of the
  square root search in towers of quadratic steps (the search branches,
  up to three calls one level down; 200 was too few for the ten
  quadratic steps behind cos(pi/257) as nested square roots).
* `GR_TOWER_OPT_GAUSS_SUM_LIMIT` (10000): square roots of primes
  up to this bound expressed as Gauss sums in the cyclotomic fields
  present.
* `GR_TOWER_OPT_DENSE_LIMIT` (100000): size of the Kronecker arrays
  of dense products.
* `GR_TOWER_OPT_INV_DENSE_DEGREE_LIMIT` (512),
  `GR_TOWER_OPT_INV_LINEAR_DEGREE_LIMIT` (64): inverses by the
  dense multivariate method in fields up to the first degree, and (in
  castle fields) by a linear system up to the second.
  `GR_TOWER_OPT_INV_DENSE_ALG` (0) selects the dense inverse
  algorithm in fields of several generators for testing: 0 the modular
  one, then linear algebra; 1 the modular one only; 2 linear algebra
  only.
* `GR_TOWER_OPT_DEFERRED_IDEAL_SIZE` (12): moduli with more
  coefficients than this over the base field are converted to the flat
  representation only when needed.
* `GR_TOWER_OPT_GROW_THRESHOLD` (6): number of generators beyond
  which a castle field extends towers in place instead of merging them.
* `GR_TOWER_OPT_ROOT_OF_UNITY_ORDER_LIMIT` (360): logarithms of
  roots of unity of orders up to this are rational multiples of
  $\pi i$.
* `GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT` (100000): roots of unity
  of orders up to this may become generators (one generator for the
  orders present, see `GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT`).
* `GR_TOWER_OPT_REALIFY_ORDER_LIMIT` (480): real elements written
  with roots of unity of orders up to this (and whose common level is
  at most this) are rewritten in real generators in real fields.
* `GR_TOWER_OPT_TRIG_PI_LIMIT` (240): sine, cosine and tangent of
  $r \pi$ for denominators of $r$ up to this in the tangent normal form.
* `GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT` (1000): beyond the
  tangent normal form, the trigonometric values at $r \pi$ are
  algebraic numbers (roots of unity, or real algebraic numbers) for
  denominators of $r$ up to this; beyond it, they are built from
  exponentials.
* `GR_TOWER_OPT_TRIG_FORM` (`GR_TOWER_TRIG_EXPONENTIAL`):
  with `GR_TOWER_TRIG_TANGENT`, sine, cosine and tangent of real
  arguments in complex castle fields take the real forms (through
  $\tan(x/2)$, and the tangent normal form at rational multiples of
  $\pi$), as in real fields, instead of exponentials.
* `GR_TOWER_OPT_SPLIT_IMAGINARY` (0),
  `GR_TOWER_OPT_COMPOSITE_RADICALS` (0): the choice of generators
  for square roots (`gr_ca_ctx_set_gen_flags`).
* `GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT` (0): the largest degree
  $\varphi(N)$ of a root of unity $\zeta_N$ of composite order used as
  one generator for the roots of unity whose orders divide $N$; beyond
  it (and with 0), roots of unity of prime power orders are used. One
  field is at best about twice as fast for dense arithmetic within it
  (sums of products of random elements of $\mathbb{Q}(\zeta_{60})$ and
  $\mathbb{Q}(\zeta_{120})$), and slower for orders with several odd
  prime factors ($\mathbb{Q}(\zeta_{105})$, $\mathbb{Q}(\zeta_{255})$:
  twice as slow), and computations involving many orders build large
  fields (the DFT benchmark with $x_k = (k+2)^{1000}$, $N = 105$ is
  an order of magnitude slower with a cap of 64 or 128), whence the
  default.
* `GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT` (120): when towers
  are merged, conjugate roots of one integer polynomial of degree up
  to this are recognized by it (and build a splitting tower with the
  Cauchy moduli).
* `GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT` (1024): annihilating
  polynomials of elements of towers of degree up to this are found by
  linear algebra over the flat representation.
* `GR_TOWER_OPT_GAMMA_LATTICE_LIMIT` (36),
  `GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT` (240): normal forms of
  $\Gamma$ and Hurwitz zeta values at rationals up to these denominators.
* `GR_TOWER_OPT_GAMMA_LINE_LIMIT` (12): the relations between
  $\Gamma$ values on a rational line $a w + b$ are searched for
  multipliers $|a|$ up to this.
* `GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT` (64): Hurwitz zeta values
  and polygamma functions up to this weight take part in the relation
  searches and the normal forms.
* `GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT` (240): levels of the
  relations between $\Gamma$ and Hurwitz zeta generators searched by the
  zero test.
* `GR_TOWER_OPT_GAUSS_DIGAMMA_LIMIT` (30): Gauss's digamma theorem
  for denominators up to this.
* `GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT` (1048576): the dense
  forms of elements of number fields of one generator up to this degree
  (0: none), see the castle field below.
* `GR_TOWER_OPT_DENSE_FORM_SPARSITY` (8): an element given by a
  polynomial in the generator longer than 256 takes the dense form only
  if at least one coefficient in this many is nonzero (0: whatever its
  sparsity), see the castle field below.
* `GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT` (0): in castle fields, a
  merged tower of several proven steps over $\mathbb{Q}$ with monic
  integral univariate moduli (radicals, generic algebraic numbers; no
  roots of unity, tangents or transcendental generators) of degree at
  most this is replaced by $\mathbb{Q}(\theta)$ for a primitive element
  $\theta = \ell \sum c_i a_i$ (small integers $c_i$ making
  $1, \theta, \ldots, \theta^{D-1}$ independent, $\ell$ the least
  integer making $\theta$ integral), the definitions of the tower being
  recorded as polynomials in $\theta$; its elements then have the dense
  form. It loses the structure of the steps (sparse elements become
  dense with large coefficients), so it pays only for small degrees:
  sums of products of random elements are 1.7 times faster in
  $\mathbb{Q}(\sqrt{2}, \sqrt{3})$, but 7 times slower in
  $\mathbb{Q}(\sqrt{2}, \sqrt{3}, \sqrt{5})$, 30 times slower with
  $\sqrt{7}$ added, and fields of degree 20 later merged with
  transcendental generators make zero tests by inversion very slow.
* `GR_TOWER_OPT_SPLIT_DEGREE_LIMIT` (4096),
  `GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT` (512): in castle fields, the
  zero test of an algebraic number $f + g$ in a tower of nominal degree
  beyond the first limit, not all proven, with $f$ and $g$ involving
  disjoint sets of generators (the difference of two numbers built
  independently, whose towers a merge has stacked), decides $f = -g$
  in towers of their own generators: by $P(-g) = 0$ for the minimal
  polynomial $P$ of $f$ (or the other way round; the smaller tower, of
  degree up to the second limit), a zero test in the tower of $-g$
  alone, and the numerical isolation of the roots of $P$. With 0, never.
* `GR_TOWER_OPT_POWER_CHECK_DEGREE_LIMIT` (64): in castle fields,
  a root $\sqrt[n]{x}$ which none of the cheaper rules finds in the
  field is looked for there when $x$ is likely a $p$-th power for a
  prime $p \mid n$ and $p$ times the degree of the tower is at most
  this (with transcendental generators, 3/8 of it, halved for each
  generator beyond the first: the norms of the factorization are
  multivariate): the steps of the tower are proven first if needed (the
  modular proofs, and Trager's method within its limit), and the
  places of the modular proof of irreducibility of $X^n - x$ at which
  $x$ is a unit and $p \mid q - 1$ must all have shown $x$ to be a
  $p$-th power, at least eight of them (more places are tried in that
  case only; an element which is not a $p$-th power is one at about
  one such place in $p$). The roots of $X^p - x$ are then found by
  factoring over the tower, the principal one is identified by its
  enclosure, and $\sqrt[n]{x} = \sqrt[n/p]{\sqrt[p]{x}}$ (principal
  roots). Otherwise, as with 0, a new generator is adjoined whose
  modulus may be reducible, which the zero tests reduce when they meet
  it ($\sqrt{1/(\pi^{1/6} + 2\pi + 1)^2}$ is the inverse of $\pi^{1/6}
  + 2\pi + 1$ at once, rather than a generator of degree 2 over
  $\mathbb{Q}(\pi^{1/6})$). Over the profile catalog, the check costs
  about 0.2 seconds in all (most of it the modular test, which the
  adjunction makes anyway).
* `GR_TOWER_OPT_FACTORED_DENOMINATORS` (0): with 1, the arithmetic of
  castle field elements in towers without algebraic generators (rational
  functions of $\pi$, $e$, logarithms, ...) keeps the denominators as
  products of irreducible factors ([](fmpz_mpoly_q#factored)), and the
  polynomial gcds become trial divisions. It is off by default: it wins
  in some linear algebra over these fields (determinants and systems
  over $\mathbb{Q}(\pi, e, \log 2)$ two to two and a half times
  faster) but loses in sums of many terms with distinct denominators
  (up to fifteen times slower) and in iterated rational maps
  ([](algorithms/exact/castle-representations#factored)).
* `GR_TOWER_OPT_UNLIMITED` (0): with 1, the effort limits below grow
  until an operation succeeds (the *unlimited mode*, see
  [](#effort)); set by the flag `GR_CA_UNLIMITED` of the castle
  fields.

## Effort limits and the unlimited mode {#effort}

Some options bound algorithms whose exhaustion makes an operation fail
(`GR_UNABLE`, `T_UNKNOWN`): the precision limits, the degrees up to
which factorizations and proofs are attempted, the numbers of tries, the
sizes of the special values computed exactly. They are read as *effort
limits*: the value of the option (or of a built-in bound of the same
nature) times $2^e$, where $e$ is the *effort round* in progress in the
thread, and an exhausted effort limit is recorded. Outside the unlimited
mode the effort round is always 0, so that the limits are the options
themselves. In the unlimited mode (`GR_TOWER_OPT_UNLIMITED`), an
outermost operation of a castle field which fails after hitting an
effort limit is repeated in the next effort round, until it succeeds or
fails without hitting one: every heuristic of a sequence "try A up to
its bound, then B up to its bound" is retried with larger bounds, rather
than one of them running forever. The options which only select between
algorithms or representations are read directly, and are not scaled.
See [](algorithms/exact/effort-limits).

::: macro
GR_TOWER_LIMIT(T, k)
GR_TOWER_WITHIN(T, v, k)
GR_TOWER_WITHIN_N(v, n)
GR_TOWER_UNLIMITED(T)

The effort limit of the option `k` of the tower `T`; whether `v` is
within it (otherwise the exhaustion is recorded and the result is 0);
the same for a built-in bound `n`; and whether the unlimited mode is on.
:::

::: function
slong gr_tower_effort_scale(slong v)

The bound `v` scaled by the effort round of the thread: $v \cdot 2^e$,
saturating at `WORD_MAX / 64`; `v` itself when $e = 0$ or $v \le 0$.
:::
