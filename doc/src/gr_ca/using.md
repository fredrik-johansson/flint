---
title: Using castle fields
summary: How to create castle fields and numbers, read and print them, what happens to the towers as numbers are computed, how equality is decided and what it costs, number fields and their fast representations, polynomial roots and splitting towers, real fields, special functions, memory, threads and the Python interface.
order: 1
---

This page describes how to compute with the castle fields of
[](gr_ca). The concepts are introduced in
[](guide/types/castle-arithmetic), and the algorithms are described in
[](algorithms/exact).

## Contexts and numbers {#contexts}

A castle field is created with `gr_ctx_init_ca` over the rational
numbers and cleared with `gr_ctx_clear` after all its elements:

```c
gr_ctx_t QQ, C, R;
gr_ctx_init_fmpq(QQ);
gr_ctx_init_ca(C, QQ, 0);                 /* the complex field */
gr_ctx_init_ca_view(R, C, GR_CA_REAL);    /* its real subfield */
...
gr_ctx_clear(R);
gr_ctx_clear(C);
gr_ctx_clear(QQ);
```

The flags `GR_CA_REAL` and `GR_CA_ALGEBRAIC` select the real and the
algebraic subfields, and `GR_CA_UNLIMITED` the unlimited mode
([](gr_ca#unlimited)). A *view* shares the towers, the definitions and
the lock of its parent, so that elements move between the two by a copy
(`gr_set_other` checks that the element lies in the target field).
Between castle fields with different states, `gr_set_other` transfers an
element structurally: the generators it involves are recreated in the
target from their definitions, roots being matched by their enclosures.

Numbers are created with the usual functions of [](gr): `gr_set_si`,
`gr_set_fmpq`, `gr_pi`, `gr_i`, `gr_sqrt`, `gr_exp`, `gr_gamma`, and so
on, from other rings with `gr_set_other` (rational numbers, [](qqbar)
numbers, elements of other castle fields), and from strings with
`gr_set_str`, which parses expressions with the elementary and special
functions (`sqrt(2) + exp(pi*i/3)`, `gamma(1/3)`, `root(x^5-x-1,
1.1673)`).

## Printing and reading back {#printing}

A number prints as an expression in the generators of its tower,
followed by the definitions of the generators it involves, in order of
creation:

    3*t2+t1+2 {t1 = exp(2); t2 = exp(t1)}

The algebraic generators are named `a1`, `a2`, ... and the
transcendental ones `t1`, `t2`, ... (with `pi` and `i`), uniquely in
the context, so that numbers print consistently whatever tower they
live in. The definitions are `sqrt(x)`, `root(x, n)`, `exp(x)`,
`log(x)`, `exp(2*pi*i/n)`, the special functions (`gamma(x)`,
`zeta(x)`, `polylog(s, x)`, `modular_j(tau)`, ...) and, for an
algebraic generator without a closed form, `root(m(a), approximation)`
with its minimal polynomial over the field below and an approximation
which separates it from the other roots. `gr_ca_ctx_set_print` selects
the symbolic form, the definitions, a numerical value, or a combination:

    4863.92 {3*t2+t1+2 where t1 = exp(2); t2 = exp(t1)}

The printed form with the definitions is read back by `gr_set_str` in
any castle field: the names of a string are local to it. A string
without definitions is read in the generators of the context by their
names. `gr_get_fexpr` gives the number as a symbolic expression
([](fexpr)), with `Where` and `Def` for the definitions, from which
LaTeX is produced by `fexpr_get_str_latex`.

## How towers are built {#towers}

Every number lives in a tower, and the shortest prefix of the tower in
which it can be represented is recorded with it. The operations build
towers as follows.

- **Rational numbers** live in the trivial tower $\mathbb{Q}$, and
  arithmetic between them costs about as much as with `fmpq_t`.
- **Radicals and roots of unity** of rational numbers are products of
  canonical generators of prime power orders, one tower each:
  $\sqrt{12} = 2 \sqrt 2 \sqrt 3$, $\sqrt{-5} = i \sqrt 5$ unless
  $\sqrt{-5}$ is one generator (the default for negative numbers:
  $\mathbb{Q}(\sqrt{-163})$ has degree 2), $\zeta_{15} = \zeta_3^a
  \zeta_5^b$. The choice is configurable (`gr_ca_ctx_set_gen_flags`).
- **Square roots and $n$-th roots** of other numbers are first sought in
  the tower of the argument (exactly, through the quadratic steps of
  the tower, then by a lattice search, and by factoring when the
  argument is likely a power), and adjoined to it otherwise; a rational
  content is taken out ($\sqrt{4\pi} = 2 \sqrt \pi$).
- **Exponentials and logarithms** are transcendental generators,
  unless a relation with the generators present expresses them: a
  relation search runs when they are adjoined, so that $\exp(\log x) =
  x$ at once, and $\exp(\sum c_j \log u_j + r \pi i)$ with rational
  $c_j, r$ and algebraic $u_j$ is a product of principal powers and a
  root of unity, without a new generator. Logarithms of rational and
  Gaussian rational numbers are sums of logarithms of primes and of
  Gaussian primes, plus a multiple of $\pi i$.
- **Trigonometric functions** go through exponentials in the complex
  field, and through tangents in the real field
  ([](algorithms/exact/castle-real-fields)); their values at rational
  multiples of $\pi$ are algebraic numbers.
- **Special functions** are brought to canonical arguments before a
  generator is created, so that related values share generators
  ([](algorithms/exact/castle-special-functions)).
- **Binary operations** on numbers of different towers take place in a
  tower containing the definitions of both: one of the two towers if it
  contains the other's, otherwise a merge, which identifies the
  generators by their definitions and maps the operands by
  substitution ([](algorithms/exact/field-towers#merging)).

Towers are shared by all the numbers built in them and can grow long
in a long computation. To keep history from slowing down new
computations, a transcendental generator for a number stored with a
short prefix of a long tower is adjoined to a copy of the prefix rather
than to the long tower, and the tower of an operand is used for a
result only if no generators of unrelated computations lie between the
definitions involved. Relations between generators of different copies
are found when the copies meet.

## Equality and its cost {#equality}

`gr_is_zero`, `gr_equal` and every predicate built on them (`gr_cmp`,
`gr_sgn`, `gr_is_real`, `gr_floor`, `gr_get_fmpq`, ...) are exact, and
proceed from cheap to expensive:

1. Rational numbers and numbers in canonical form (a proven number
   field of one generator, [](#number-fields)) compare directly.
2. Numerical enclosures at increasing precision separate numbers which
   are not equal.
3. The *zero test* of the flat representation decides the rest: the
   reduced numerator of an algebraic number is zero exactly when the
   number is zero, after the tower is refined where a factorization
   appears ([](algorithms/exact/dynamic-evaluation)); with
   transcendental generators, a nonzero element is a nonzero number when
   the generators it involves are independent, which is decided by
   searching for relations between them
   ([](algorithms/exact/zero-testing#castle)).

A relation found by a zero test is recorded in the tower: one generator
becomes algebraic over the others (of degree one, or a radical), every
number built on the tower is reduced accordingly, and later zero tests
need no search. In the default mode, a test which exceeds its effort
limits returns `T_UNKNOWN`; the options are listed in
[](gr_tower#tuning-options), and the unlimited mode removes the limits
([](gr_ca#unlimited)).

## Number fields {#number-fields}

Computations in one number field are much faster than the generality of
the towers suggests. A number of a tower whose only generator is a
proven algebraic number with a monic integral minimal polynomial
($\mathbb{Q}(\sqrt 2)$, $\mathbb{Q}(\zeta_5)$, $\mathbb{Q}(2^{1/3})$,
$\mathbb{Q}(\zeta_{9973})$) is stored in a *dense form*, a polynomial in
the generator, on which arithmetic needs neither the towers nor the
lock: the cost of a multiplication followed by an addition is
comparable to that of [](nf_elem) for small degrees and lower for large
degrees. Towers of several algebraic generators with univariate integral
moduli (cyclotomic fields as products of prime power cyclotomic fields,
radicals, quadratic fields) multiply densely by Kronecker substitution,
and invert by modular methods. Polynomial and matrix products over such
fields pack all the coefficients into integer polynomials. See
[](algorithms/exact/castle-representations).

## Polynomials, roots and factors {#roots}

The generic polynomial functions of [](gr_poly) work over castle
fields: the gcd and the resultant use subresultant sequences, and
`gr_poly_roots` and `gr_factor` find exact roots:

- a quadratic is solved by the formula, a binomial $a x^n + b$ by a
  principal root and the roots of unity;
- roots which are algebraic numbers of small degree and height are
  recognized numerically and verified exactly;
- otherwise the polynomial is factored exactly over the tower of its
  coefficients (reduced to the generators they involve), and a root of
  an irreducible factor of degree at least 3 is adjoined as a proven
  step, over which the cofactor is factored again: a *splitting tower*
  emerges, in which symmetric functions of the roots reduce to the
  coefficients;
- beyond the degree limits, roots are adjoined as dynamic steps.

In a real field only the real roots are returned, and `gr_factor` gives
quadratic factors for the pairs of nonreal roots. See
[](algorithms/exact/castle-roots).

## Linear algebra {#linear-algebra}

The generic algorithms of [](gr_mat) work over castle fields, and the
castle fields choose among them: determinants use `fmpq_mat` for
rational matrices, the cofactor expansion up to $4 \times 4$, LU
decomposition over number fields of degree smaller than the dimension,
and the division-free Berkowitz algorithm otherwise; linear systems use
fraction-free LU over number fields of large degree and LU otherwise
([](algorithms/exact/linear-algebra)).

## Real fields and trigonometric functions {#real}

A real castle field rewrites real results in real generators: the
trigonometric functions of real arguments are rational functions of
tangents of half angles, their values at rational multiples of $\pi$
are polynomials in one tangent $\tan(\pi/M)$ per level, and real numbers
written with roots of unity are converted to this form. Comparisons,
signs and the rounding functions are decided exactly. See
[](algorithms/exact/castle-real-fields).

## Special, hypergeometric and modular functions {#special}

The special functions of [](gr_special) are available: the gamma
function and its relatives, the error functions, the Lambert $W$
function, zeta, Hurwitz zeta, Dirichlet $L$-functions and the Lerch
transcendent, polylogarithms, complete elliptic integrals,
hypergeometric functions ${}_pF_q$, and the modular and elliptic
functions (the $j$-invariant, $\lambda$, $\eta$, Eisenstein series,
Jacobi theta functions and Weierstrass functions). Their values are
generators at canonical arguments, related by the identities that the
castle fields know; algebraic values (special values at rational
arguments, $j$ at CM points) are returned as algebraic numbers. See
[](algorithms/exact/castle-special-functions),
[](algorithms/exact/castle-hypergeometric) and
[](algorithms/exact/castle-modular).

## Memory {#memory}

Every number counts as a reference to its tower. Towers whose count has
dropped to zero are freed at the end of the outermost operation (in
batches), except those holding the canonical generators of the context
($\pi$, roots of unity, radicals of integers); the definitions of a
collected tower which no other tower contains are forgotten. Memory
therefore stays bounded when a long session keeps few numbers alive:
solving hundreds of polynomials over $\mathbb{Q}(\sqrt 2)$ leaves a
handful of towers. `gr_ca_ctx_num_towers` and `gr_ca_ctx_stats` report
the towers of a context.

## Threads {#threads}

A castle context is thread-safe when FLINT is built with threads: every
operation holds a recursive mutex of the context, since operations may
refine or merge the shared towers and update the cached representations
of elements. Operations on one context are therefore serialized, except
the arithmetic of rational numbers and of dense numbers of one number
field, which takes lock-free paths. Distinct contexts are independent;
for parallel work on unrelated numbers, use one context per thread.

## Python {#python}

The Python interface `src/python/flint_ctypes.py` offers the four castle
fields as `ComplexField_castle` (`CC_castle`), `RealField_castle`
(`RR_castle`), `ComplexAlgebraicField_castle` (`QQbar_castle`) and
`RealAlgebraicField_castle` (`AA_castle`), with the tuning options as
keyword arguments of the constructors
(`ComplexField_castle(unlimited=1)`), `set_print` for the display,
`gens()`, `fexpr()` and `latex()` on elements, and `tower()` to inspect
the tower of an element; the fixed towers of [](gr_tower) are the
classes `gr_tower` and `TowerField`.

```python
>>> from flint_ctypes import *
>>> C = ComplexField_castle()
>>> x = C(2).sqrt() + C(3).sqrt()
>>> x**2
2*a2*a1+5 {a1 = sqrt(2); a2 = sqrt(3)}
>>> (x**2 - 5 - 2*C(6).sqrt()) == 0
True
```

The interface runs the test cases of Calcium on the castle fields, the
open issues of the Calcium issue tracker which castle arithmetic
settles, and examples from the documentation of SageMath's `QQbar`
(`test_castle`, `test_castle_calcium_issues`,
`test_castle_sage_examples`).
