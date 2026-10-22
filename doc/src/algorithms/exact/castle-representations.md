---
title: Representations in castle arithmetic
summary: How the elements of the castle fields are stored and computed with — rational numbers, flat fractions over the polynomial ring of a tower, dense forms for simple number fields (used by density) and factored denominators for rational function fields; canonical fractions by rationalization with norms, compared with Calcium and with rational simplification modulo an ideal; fast arithmetic in products of number fields by Kronecker and Good–Thomas packing, monomial shifts and multimodular inverses at split primes; and the sharing, copying, naming and collection of towers.
complexity: A product in a simple number field of degree $d$ costs one polynomial product of length $d$ and a reduction; in a product of number fields of total degree $D$, one integer polynomial product of length about $\prod (2d_k - 1)$ (or a cyclic convolution of length $n$ for $\mathbb{Q}(\zeta_n)$); an inverse there costs $O(D \sum d_k)$ per prime over a number of primes proportional to the size of the result.
order: 36
---

The elements of a castle field ([](gr_ca)) are elements of towers
([](algorithms/exact/field-towers)) shared between them: each element
carries a pointer to a tower and the length of the prefix of that tower
in which it lives, and is stored in one of four representations. This
page describes the representations and the arithmetic in them, and how
the towers are shared, named and collected.

## Elements and towers {#elements}

Towers are created and extended in place as needed: `gr_sqrt` first
tries to express a root in the tower of its argument
([](algorithms/exact/field-towers#express)), and `gr_exp`, `gr_log` and
`gr_pi` adjoin transcendental generators, with generators of equal
definitions shared. Arithmetic between elements of different towers goes
through a tower containing both, found by the definition ids of the
generators involved, or constructed by a merge
([](algorithms/exact/field-towers#merging)); the operands are mapped by
substitution. Results are stored with the shortest prefix of their tower
in which they can be represented, and rational results move to the
*trivial tower*, which has no generators, never changes and is not
reference counted.

## The four representations {#repr}

An element holds exactly one of four
representations (a tagged union): a rational number (`fmpq`), a
flat element (`fmpz_mpoly_q` over the flat context of its tower,
the general case), a dense form, or a flat element with a factored
denominator ([](#factored)). Rational numbers live in the
trivial tower, which never changes and is not reference counted, so
that arithmetic, comparisons and conversions between rational
elements are done directly on their numerators and denominators,
without the lock (rational results of operations in other towers are
moved there).

Each element also points to its tower and records its *level* (the
length of the prefix of the tower in which it lies), the structure
version of the tower for which that level is valid, the ideal version
for which it is known to be reduced, and a flag marking shallow copies
whose data is shared and must not be modified in place. Levels are
recomputed lazily from the leading monomial (the lexicographic order
puts the most significant variable in the leading term) when the
structure of the tower changes. The element type is public
(`gr_ca_struct`), documented as read-only, with locked accessors for the
representation tag, the flat data, and the polynomial in the first
generator with its modulus, so that the coordinates of number field
elements can be read without copies through the generic interface
([](gr_ca#elements)). An earlier design kept the flat form and the dense
form simultaneously valid, which made shared operands safe to read
without locks; the union is smaller, and the lock-free paths were
preserved by never rewriting an operand in place: the locked
implementations read every operand through a temporary flat view kept
in a per-thread scratch pool, and only an owned result is converted in
place.

*Example (element forms).* In a fresh castle field, $3/7$ is rational.
The root of unity $\zeta_5 = e^{2\pi i/5}$ is flat as created; after
$(1 + \zeta_5)^3$, both the result and the operand $\zeta_5$ are dense,
and `gr_ca_get_fmpq_poly` returns $x^3 + 3x^2 + 3x + 1$ with the modulus
$x^4 + x^3 + x^2 + x + 1$. The number $\pi + \sqrt 2$ involves two
generators and stays flat, and the accessor returns `GR_DOMAIN` for it.
The tag records the history of an element, not only its value: the
same $\zeta_5$ is flat before the operation and dense after it.

## Dense forms {#dense}

An element of a tower whose first generator $a$ is
algebraic, proven, of degree $2 \le d \le$
`GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT` ($2^{20}$ by default;
0 disables dense forms) with a monic integral modulus, and which
involves only $a$ with an integer denominator ($\mathbb{Q}(\sqrt{2})$,
$\mathbb{Q}(\zeta_5)$, $\mathbb{Q}(2^{1/3})$, $\mathbb{Q}(\zeta_{9973})$),
can be stored in dense form: an `fmpq_poly` in $a$ of length less
than $d$ (normalized, without padding, with a positive denominator
coprime to the content of the numerator), referring to an immutable
descriptor of $a$ and its modulus kept by the tower. Products are
computed with FLINT's polynomial multiplication followed by a
reduction by the modulus (folding for prime cyclotomic moduli,
division with remainder for dense moduli of large degree), inverses
by the resultant, sums and differences by `fmpq_poly` additions.
Sums, differences, products, inverses and quotients, operations with
rational scalars, comparisons and the predicates of elements with
the same descriptor and tower, or with rational elements, are
computed in this form without the lock when the result already lives
in that tower (a fresh result goes through the locked operation
once). A polynomial longer than 256 takes the dense form only when
at least one coefficient in
`GR_TOWER_OPT_DENSE_FORM_SPARSITY` (8 by default) is nonzero,
so that sparse elements of large cyclotomic fields (a cosine
$(\zeta^k + \zeta^{-k})/2$ in $\mathbb{Q}(\zeta_{9973})$) stay flat.
The dense form is made at the end of the locked operations on the
outermost level, for the operands and the result, and only when
every flat operand can take it: if one stays flat, nothing is
converted (a running sum of sparse terms would otherwise change form
at every step). Operations needing the flat form of a rational or
dense operand work on a temporary flat copy (released when the
outermost locked operation returns), so that an operand is never
rewritten while another thread may be reading it without the lock;
a flat operand converted to the dense form is read by other threads
only under the lock.

The dense form is the analogue of the `nf_elem` of [](nf) for a field
$\mathbb{Q}(a)$ with one proven algebraic generator; Calcium uses
`nf_elem` for its single-generator fields. We report the cost of element
arithmetic in *Horner steps*: a Horner step is the update $x \gets xy +
c$, with $y$ a fixed element (usually the generator of the field) and $c$
a small rational constant, the inner loop of the evaluation of a
polynomial at $y$ by Horner's rule. It measures one multiplication and
one addition of a constant, including all per-operation overhead
(dispatch, locking, levels, allocation, normalization), which dominates
for small fields. Rational elements are processed without the context
lock or any tower machinery: a Horner step on rationals costs about as
much as with raw `fmpq`. Dense elements of the same descriptor add,
multiply (a shift for monomials, a scalar path for constants,
`fmpz_poly` multiplication otherwise, followed by an all-ones fold for
prime cyclotomic moduli, `fmpz_poly` division with remainder for dense
moduli of degree at least 64, or a sequential fold) and invert (a closed
form for degree two, a resultant-based extended gcd otherwise) without
the lock as well, using FLINT's polynomial kernels at any degree. Results
reuse the storage of the destination, or a per-thread scratch polynomial
swapped with it, so that no allocation occurs in steady state. The
following times (microseconds, measured during development with
`examples/number_field_bench.c`) compare the flat form, the dense form
and `nf_elem` (through `gr_ctx_init_nf`) in prime cyclotomic fields: a
Horner step by the generator, a product of two dense elements with small
coefficients, a sum, and the inverse of $x^2 + k$.

| Field | Horner step | Product | Sum | Inverse |
|---|---|---|---|---|
| $\mathbb{Q}(\zeta_5)$ | 0.57 / 0.13 / 0.17 | 0.50 / 0.18 / 0.22 | 0.38 / 0.06 / 0.03 | 8.0 / 3.1 / 3.5 |
| $\mathbb{Q}(\zeta_{67})$ | 2.2 / 0.53 / 1.5 | 11.6 / 1.6 / 10.0 | 0.50 / 0.21 / 0.19 | 938 / 700 / 742 |
| $\mathbb{Q}(\zeta_{131})$ | 3.0 / 0.72 / 3.8 | 25.7 / 4.4 / 36.7 | 1.15 / 0.33 / 0.37 | 2359 / 2292 / 2578 |
| $\mathbb{Q}(\zeta_{257})$ | 8.2 / 2.4 / 4.7 | 94 / 9.9 / 68 | 3.2 / 1.2 / 1.0 | 11161 / 9253 / 8344 |
| $\mathbb{Q}(\zeta_{521})$ | 10.4 / 2.8 / 14.9 | 139 / 15.5 / 232 | 3.8 / 1.5 / 2.1 | 32201 / 33273 / 33105 |

In large degree the dense form is much faster than `nf_elem`, which
reduces by precomputed powers of the generator: products in
$\mathbb{Q}(\zeta_{521})$ are 9 times faster than in the flat form and 15
times faster than with `nf_elem`, and inverses cost the same in the
three forms.

### When the dense form is used {#density}

The dense form was first an inline coefficient array of at most 64
entries, so that large single-generator fields used the flat form; it is
now an `fmpq_poly`, and the degree limit is an option. The degree is not
the right criterion; the density of the element is. An element of
$\mathbb{Q}(\zeta_{9973})$ such as $\cos(2\pi k/9973) = (\zeta^k +
\zeta^{-k})/2$ with $2 \le k \le 9971$ has two terms in the flat form and
nearly $10^4$ coefficients in the dense one. The rule is therefore: a
polynomial of length at most 256 always takes the dense form (operations
on it are cheap in either form, and conversions are not; the DFT
benchmark in $\mathbb{Q}(\zeta_{256})$ is fastest with every element
dense), and a longer one only when at least one coefficient in
`GR_TOWER_OPT_DENSE_FORM_SPARSITY` (8) is nonzero. The second part of the
rule concerns operations rather than elements: the dense form is made at
the end of a locked operation, for the operands and the result, and if
one operand cannot take it, nothing is converted. Without this, a running
sum of sparse terms in $\mathbb{Q}(\zeta_{9973})$ became dense as soon as
it was dense enough and was then read through a temporary flat view at
every following addition, since the next term stayed flat: the sum of
Gauss's cubic period of 9973 (1662 cosines) took 0.87 s against 0.33 s in
the flat form, and 0.31 s with the rule. Neither form was slow; the
conversions were. The tests of the rule read only the leading term and
the number of terms of each operand before any coefficient is
extracted.

*Example (the density rule in $\mathbb{Q}(\zeta_{9973})$).* Let $z =
e^{2\pi i/9973}$, of degree 9972. Each element below is made flat by the
parser in a fresh field (the operations of the parser are inner
operations, which do not make dense forms), and then multiplied by 3;
the table gives the length of its polynomial in $z$, the number of
nonzero coefficients, and the form of the product and of the operand
afterwards.

| Element | Length | Nonzero | Form after $3x$ |
|---|---|---|---|
| $z^{200} + 1$ | 201 | 2 | dense (length at most 256) |
| $z^{300} + 1$ | 301 | 2 | flat |
| $(z^{800} - 1)/(z^8 - 1) = \sum_{j<100} z^{8j}$ | 793 | 100 | dense ($8 \cdot 100 \ge 793$) |
| $(z^{900} - 1)/(z^9 - 1) = \sum_{j<100} z^{9j}$ | 892 | 100 | flat ($8 \cdot 100 < 892$) |
| $z^2 + z^{-2} = z^2 + z^{9971}$ | 9972 | 2 | flat |
| $z + z^{-1}$ | 9972 | 9971 | dense |

The last line is the exception to $\cos(2\pi k/9973)$ being sparse:
$z^{-1} = z^{9972} = -(1 + z + \cdots + z^{9971})$ modulo $\Phi_{9973}$,
so that $\cos(2\pi/9973)$ has 9971 nonzero coefficients in either form.
The rule applies when a dense form is made from a flat one; arithmetic
between dense elements keeps the dense form whatever the density. If $z$
itself has been made dense (as an operand of $z + 1$), then $z^{300} +
1$ computed from it is dense.

## Products of number fields {#number-fields}

When the tower of the operands has only algebraic generators whose
moduli are monic univariate integer polynomials, each in its own
variable (cyclotomic polynomials of prime power orders, quadratic and
binomial moduli: $\mathbb{Q}(\zeta_N)$ with the roots of unity of prime
power orders, $\mathbb{Q}(\sqrt{-163})$, radicals of integers, and
products of such fields), its elements are polynomials with a tensor
structure, and dense kernels replace the generic multivariate
arithmetic for elements with integer denominators. The generators of
these fields are structured
([](algorithms/exact/field-towers#structured)), so that the moduli are
univariate and the kernels apply. They were developed against a
benchmark comparing Calcium, castle arithmetic and Antic's `nf_elem` on
element steps, polynomial and matrix products, determinants and linear
systems over $\mathbb{Q}(\zeta_N)$ and $\mathbb{Q}(\sqrt D)$
(`examples/number_field_bench.c`); at the start, the towers were 5–10
times slower than Calcium on elements and 20–40 times on matrix
products, and a $6 \times 6$ solve over $\mathbb{Q}(\zeta_{128})$ took
127 s.

### Kronecker substitution {#kronecker}

For generators of degrees $d_1, \ldots, d_m$, an element with exponent
$e_k < d_k$ in generator $k$ is packed into a univariate polynomial with
strides

$$
S_1 = 1, \qquad S_{k+1} = S_k (2d_k - 1),
$$

so that the exponents of a product, which are $< 2d_k - 1$, do not
overlap [@Kro1882; @Har2009a]. One `fmpz_poly` product is followed by a
fold by the moduli, digit by digit from the most significant (vector
`submul` over the contiguous lower digits, only nonzero modulus terms),
and an unpacking in decreasing lexicographic order. The array has length
$L = \prod_k (2d_k - 1) \approx 2^m D$ for $D = \prod_k d_k$.

*Example (Kronecker substitution in $\mathbb{Q}(\sqrt 2, \sqrt 3)$).*
Here $d_1 = d_2 = 2$, so $S_1 = 1$, $S_2 = 3$, $L = 9$, and the monomial
$\sqrt 2^{e_1} \sqrt 3^{e_2}$ is packed as $x^{e_1 + 3e_2}$. The elements
$a = 1 + 2\sqrt 2 + 3\sqrt 3 + \sqrt 6$ and $b = \sqrt 2 + \sqrt 3$
become $A = 1 + 2x + 3x^3 + x^4$ and $B = x + x^3$, and

$$
AB = x + 2x^2 + x^3 + 5x^4 + x^5 + 3x^6 + x^7,
$$

whose exponents have digits $e \bmod 3$ and $\lfloor e/3 \rfloor$ at most
$2 = 2d_k - 2$. The fold by $x_1^2 = 2$ moves digit 2 of each block of
three entries to digit 0 with the factor 2 ($2x^2 \to 4$, $x^5 \to
2x^3$), giving $4 + x + 3x^3 + 5x^4 + 3x^6 + x^7$; the fold by $x_2^2 =
3$ moves the top block $3x^6 + x^7$ down with the factor 3, giving $13 +
4x + 3x^3 + 5x^4$. Unpacking gives $ab = 13 + 4\sqrt 2 + 3\sqrt 3 + 5\sqrt
6$.

### The cyclic layout {#cyclic}

For roots of unity of pairwise coprime orders $q_k$, the map $x_k \mapsto
x^{\epsilon_k}$ with $\epsilon_k \equiv \delta_{jk} \pmod{q_j}$ identifies
$\bigotimes_k \mathbb{Z}[x_k]/(x_k^{q_k} - 1)$ with $\mathbb{Z}[x]/(x^n -
1)$, $n = \prod_k q_k$, by the Chinese remainder theorem (the Good–Thomas
prime factor mapping [@Goo1958]). A product becomes one convolution of
length $n$ with a wrap, followed by the folds by the $\Phi_{q_k}$; this
is used when $2n - 1 < L$, and polynomial and matrix products use the
same layout. Horner steps in $\mathbb{Q}(\zeta_{1155})$ went from 33 to
7.8 µs with it during development.

*Example (the cyclic layout for $\mathbb{Q}(\zeta_3, \zeta_5)$).* For
$q_1 = 3$, $q_2 = 5$ and $n = 15$, $\epsilon_1 = 10$ ($\equiv 1 \bmod 3$,
$\equiv 0 \bmod 5$) and $\epsilon_2 = 6$, so $\zeta_3 \mapsto x^{10}$,
$\zeta_5 \mapsto x^6$, and $x$ itself corresponds to $\zeta_3\zeta_5 =
\zeta_{15}^8$. The elements $u = \zeta_3 + \zeta_5^3$ and $v = \zeta_3 +
\zeta_5^2$ become $x^{10} + x^3$ (as $18 \equiv 3$) and $x^{10} + x^{12}$,
and their product modulo $x^{15} - 1$ is $1 + x^5 + x^7 + x^{13}$. The
exponents map back to the digits $(e \bmod 3, e \bmod 5) = (0, 0), (2,
0), (1, 2), (1, 3)$, that is, to $1 + \zeta_3^2 + \zeta_3\zeta_5^2 +
\zeta_3\zeta_5^3$, and the fold by $\Phi_3$ gives $uv = -\zeta_3 +
\zeta_3\zeta_5^2 + \zeta_3\zeta_5^3$. For two factors the Kronecker array
is shorter ($L = 3 \cdot 7 = 21 < 2n - 1 = 29$) and is the one used for
this field; the cyclic layout takes over from three factors: $209 < 231$
for $\mathbb{Q}(\zeta_3, \zeta_5, \zeta_7)$ ($\epsilon_k = 70, 21, 15$)
and $2309 < 4389$ for $\mathbb{Q}(\zeta_{1155})$.

### Sparse and monomial products {#monomial}

When $\operatorname{len}(x) \operatorname{len}(y) \cdot 8 < L$, the
sparse product is used instead (a 12-term product in a tower of eight
square roots: 4.8 ms to 0.05 ms). Products by monomials are shifts, and
powers of roots of unity are recognized from the exponents of their
terms and applied as shifts: $\zeta_n^j$ is a tensor product of
$\zeta_{q}^{e}$, each a monomial or, for $e \ge d$ ($q = p^f$ an odd
prime power), $-a^{e-d}(1 + a^{q/p} + \cdots + a^{(p-2)q/p})$. The DFT
benchmark, dominated by products by powers of roots of unity, depends on
this.

*Example (powers of roots of unity as monomials).* The canonical
decomposition writes $\zeta_{15} = \zeta_3^2\zeta_5^2$ (as $2/3 + 2/5
\equiv 1/15 \pmod 1$), and since $\zeta_3^2 = -1 - \zeta_3$ is not
reduced, $\zeta_{15}$ has two terms, $-a_1a_2^2 - a_2^2$ for $a_1 =
\zeta_3$, $a_2 = \zeta_5$. The two terms have equal coefficients and take
all $p - 1 = 2$ exponents of the expansion $-a_1^{e-d}(1 + a_1)$ ($q = p =
3$, $e = d = 2$) in $a_1$, so the product recognizes the monomial
$a_1^2a_2^2$ and multiplies by it as a shift by 2 along each axis
followed by the folds: $\zeta_{15}(1 + \zeta_3\zeta_5^3) =
\zeta_3^2\zeta_5^2 + \zeta_3^3\zeta_5^5 = 1 - \zeta_5^2 -
\zeta_3\zeta_5^2$. For $q = 9$ ($d = 6$), $\zeta_9^7$ is stored as
$-\zeta_9(1 + \zeta_9^3) = -\zeta_9 - \zeta_9^4$.

*Reduction.* Measured on monic dense moduli with small coefficients, the
sequential fold beats reduction by precomputed powers of the generator
(Antic's approach) almost everywhere, because the powers have larger
coefficients, and `fmpz_poly_rem` wins from degree 64; the
implementation uses the fold below 64 and the division above, with
special folds for sparse moduli (only the nonzero coefficients of a
modulus such as $x^{64} + 1$ are used) and prime cyclotomic moduli ($1 +
x + \cdots + x^{q-1}$ applied as $x^q = 1$ followed by one subtraction).
Sums with integer denominators only cancel the gcd of the denominators,
and arithmetic with rational scalars (`gr_add_si`, `gr_mul_fmpq`,
`gr_div_ui`, ...) works on the element directly.

### Products of polynomials and matrices {#packed-products}

Polynomial products (`gr_poly_mullow`, from length 4) and matrix
products (`gr_mat_mul`, from $3 \times 3$) over such fields find a
common tower of all entries, clear common denominators (by rows of $A$
and columns of $B$ for matrices), perform one `fmpz_poly` (or
`fmpz_poly_mat`) product of the packed entries, and reduce the outputs.
This is Calcium's `ca_mat_mul_same_nf` approach, extended to
multi-generator fields. Products with rational numbers are scalar
multiplications.

*Example (a matrix product over $\mathbb{Q}(\zeta_5)$).* For $\zeta =
\zeta_5$ and

$$
A = \begin{pmatrix} \frac12 & \frac{\zeta}{2} & 0 \\ 1 & 1 & \zeta^2 \\ \frac{\zeta^3}{3} & 0 & 1 \end{pmatrix}, \qquad
B = \begin{pmatrix} \zeta & 0 & 1 \\ 1 & \frac12 & 0 \\ 0 & \zeta & \frac13 \end{pmatrix},
$$

the row denominators of $A$ are $r = (2, 1, 3)$ and the column
denominators of $B$ are $c = (1, 2, 3)$. The scaled entries, packed
($\zeta \mapsto x$), form two matrices over $\mathbb{Z}[x]$, and one
`fmpz_poly_mat` product gives

$$
\begin{pmatrix} 1 & x & 0 \\ 1 & 1 & x^2 \\ x^3 & 0 & 3 \end{pmatrix}
\begin{pmatrix} x & 0 & 3 \\ 1 & 1 & 0 \\ 0 & 2x & 1 \end{pmatrix} =
\begin{pmatrix} 2x & x & 3 \\ x + 1 & 2x^3 + 1 & x^2 + 3 \\ x^4 & 6x & 3x^3 + 3 \end{pmatrix}.
$$

Each entry is reduced modulo $\Phi_5$ (here only $x^4 \to -1 - x - x^2 -
x^3$) and divided by $r_ic_j$: for instance $(AB)_{22} = (1 +
2\zeta^3)/2$ and $(AB)_{31} = -(1 + \zeta + \zeta^2 + \zeta^3)/3$.

### Inverses {#inverses}

A single generator uses the extended gcd of `fmpq_poly`, or a closed
form in degree two; the dense element form uses the resultant-based
extended gcd of `fmpz_poly`, which yields the inverse with the norm as
its denominator ([](#canonical)). For several generators, the
implementation inverts modulo primes at which every modulus splits into
distinct linear factors.

::: proposition Multimodular inverses {#prop-modular-inverse}
Let $A$ be an element of a tower of univariate monic integer moduli
$M_1, \ldots, M_m$, let $N \in \mathbb{Z}$ and let $X$ be an integer
coefficient vector with $XA \equiv N$ modulo the moduli and modulo
primes $p_1, \ldots, p_s$. Let $G = \prod_k R_k$ with $R_k =
\max_i \sum_{e < 2d_k - 1} |[x^i](x^e \bmod M_k)|$. If $P = \prod_j p_j
> 2(\|X\|_1 \|A\|_\infty G + |N|)$, then $XA = N$ exactly, so $X/N$ is the
inverse of $A$.
:::

::: proof
The coefficients of the reduced product $XA$ are bounded by $\|X\|_1
\|A\|_\infty G$: the product before reduction is bounded by $\|X\|_1
\|A\|_\infty$ coefficientwise, and the reduction modulo $M_k$ multiplies
the bound by at most $R_k$. Each coefficient of $XA - N$ is divisible by
$P$ and smaller than $P/2$ in absolute value, hence zero.
:::

Modulo a split prime $p$ (primes $p \equiv 1 \pmod n$ for roots of unity
of order $n$, and $p \equiv 1 \pmod{8|c|}$ for $\sqrt c$, by quadratic
reciprocity), the algebra is $\mathbb{F}_p^D$, reached by one Vandermonde
transform per axis (cyclotomic roots as powers of a primitive root,
quadratic roots by modular square roots, others by `nmod_poly_roots`), so
the inverse is pointwise. The adjugate $X$ and the norm $N$ are
reconstructed by batched Chinese remaindering with a doubling number of
primes, and [](#prop-modular-inverse) replaces an exact product check
that had cost 70% of the time. Primes are drawn from progressions $1
\bmod \text{step}$ just below $2^{60}$ (so that a first batch of four
primes covers bounds of about 240 bits) and cached with their roots per
ideal version. The cost is $O(D(d_1 + \cdots + d_m))$ per prime, against
$O(D^3)$ for a linear system with the multiplication matrix; a failed
search falls back to Dixon's $p$-adic solver [@Dix1982] on the
multiplication matrix (up to degree 512), and to norms with nested
Euclid beyond. An inverse in $\mathbb{Q}(\zeta_{13}, \zeta_{11},
\zeta_7)$ (degree 720) went from more than 300 s (nested Euclid) to
14 ms. Delayed-reduction dot products and Shoup's precomputed
multiplication [@Har2014] are used in the transforms.

*Example (a modular inverse in $\mathbb{Q}(\sqrt 2, \sqrt 3)$).* Let $a =
1 + 2\sqrt 2 + 3\sqrt 3 + \sqrt 6$ as in [](#kronecker). Its norm is $N =
184 = 2^3 \cdot 23$ and its adjugate is $X = N/a = 32 - 34\sqrt 2 +
44\sqrt 3 - 18\sqrt 6$. Both moduli split modulo the primes $p \equiv 1
\bmod 48$ (the step is $\operatorname{lcm}(2, 8 \cdot 2, 8 \cdot 3)$).
For $p = 97$, $\sqrt 2 \equiv \pm 14$ and $\sqrt 3 \equiv \pm 10$, and the
transform evaluates $a$ at $(14, 10), (14, -10), (-14, 10), (-14,
-10)$, with the values $5, 53, 57, 83$. Their product is $N \bmod 97 =
87$, the pointwise quotients $N/a$ are $95, 84, 73, 70$, and the inverse
transform gives $X \bmod 97 = (32, 63, 44, 79)$ on the basis $1, \sqrt 2,
\sqrt 3, \sqrt 6$. With $R_1 = 3$, $R_2 = 4$ (from $x^2 \equiv 2$ and $x^2
\equiv 3$), $G = 12$, $\|X\|_1 = 128$ and $\|a\|_\infty = 3$, the bound
of [](#prop-modular-inverse) is $2(128 \cdot 3 \cdot 12 + 184) = 9584$,
so the primes 97 and 193 suffice: Chinese remaindering in the symmetric
range recovers $X$ and $N$, and

$$
\frac1a = \frac{X}{N} = \frac{16 - 17\sqrt 2 + 22\sqrt 3 - 9\sqrt 6}{92},
$$

which is what the implementation returns. (The implementation itself
draws its primes from the same progression just below $2^{60}$.)

Towers whose moduli are not univariate (triangular steps such as
$\sqrt{1 + \sqrt 2}$, or roots of polynomials with algebraic
coefficients) of degree at most 64 over $\mathbb{Q}$ invert through the
multiplication matrix on the monomial basis and a rational linear solve,
its columns $x b_j$ obtained from one another by multiplications by
single generators, which is 2–3 times faster than rationalization
through the nested representation, which remains the fallback. The last
dense inverse is cached in the tower, since fraction-free elimination
divides every entry of a step by the same pivot.

## Canonical fractions: rationalization by norms {#canonical}

Let the tower have transcendental generators $t = (t_1, \ldots, t_m)$
and algebraic generators $a = (a_1, \ldots, a_n)$, so that $F = K[a]/I$
with $K = \mathbb{Q}(t)$ and $I \subset K[a]$ generated by the
triangular set $m_1, \ldots, m_n$, $m_k$ of degree $d_k$ in $a_k$ with
a leading coefficient in $K$ (the flat moduli are these polynomials
with their denominators cleared). The flat form of an element is a
fraction $N/D$ with

- $D \in \mathbb{Z}[t]$: the denominator involves no algebraic
  generator;
- $N \in \mathbb{Z}[t][a]$ reduced: $\deg_{a_k} N < d_k$ for every $k$;
- $N$ and $D$ without common factor in $\mathbb{Z}[t]$, and $D$ with a
  positive leading coefficient.

::: proposition Canonical fractions {#prop-canonical}
If every algebraic step is proven irreducible, every element of $F$ has
exactly one such representation.
:::

::: proof
Uniqueness: if $N/D = N'/D'$ in $F$, then $ND' - N'D \in I$; this
polynomial is reduced (the reduced polynomials form a $K$-subspace, and
$D, D' \in K$), and the only reduced element of $I$ is 0 because the
triangular set is a Gröbner basis
([](algorithms/exact/field-towers#prop-groebner)). Hence $ND' = N'D$ in
$\mathbb{Z}[t][a]$, and since $\mathbb{Z}[t]$ is a unique factorization
domain and $N, D$ (respectively $N', D'$) have no common factor in
$\mathbb{Z}[t]$, $D$ and $D'$ divide each other; with the sign
normalization $D = D'$ and $N = N'$. Existence: the steps being
irreducible, $F$ is a field and each $F_k = F_{k-1}[x_k]/(m_k)$ is a
vector space over $F_{k-1}$ with basis $1, a_k, \ldots, a_k^{d_k - 1}$;
so the reduced monomials $a^e$ ($e_k < d_k$) form a basis of $F$ over
$K$, and the coordinates of an element are elements of $K$. Clearing
their denominators gives the form above.
:::

The computation of the coordinates of a quotient $A/B$ is where norms
come in. For $B \in F_k$, multiplication by $B$ is an $F_{k-1}$-linear
map of $F_k$ with matrix $M_B$ (of size $d_k$), and

$$
B \cdot \bigl(\operatorname{adj}(M_B) e_0\bigr) = \det M_B = N_{F_k/F_{k-1}}(B) \in F_{k-1},
$$

so that $1/B = \operatorname{adj}(M_B) e_0 / N_{F_k/F_{k-1}}(B)$: the
algebraic denominator $B$ is replaced by its relative norm, an element
of the field below, and the procedure recurses until the denominator
lies in $K$. For a single generator the same is done in one step with
the resultant: $sm + uB = \operatorname{Res}(m, B) = N_{F/\mathbb{Q}}(B)$
(up to sign) gives $1/B = u/\operatorname{Res}(m, B)$, which is how
inverses of dense elements are computed ([](#dense)). A monomial in
algebraic generators is inverted through the moduli ($1/a = -(c_1 + c_2a
+ \cdots + a^{d-1})/c_0$ for $m = c_0 + c_1x + \cdots + x^d$,
recursively), and a denominator in $K$ costs nothing. In the castle
fields, every flat operation produces this form directly: sums and
products of fractions with denominators in $\mathbb{Z}[t]$ need a
reduction of the numerator and a gcd with a polynomial of $\mathbb{Z}[t]$
only (an integer gcd when there are no transcendental generators), and a
quotient is a product by an inverse computed as above. If a norm turns
out to be zero, $B$ is a zero divisor of a pretend field, and the nested
division triggers a refinement
([](algorithms/exact/dynamic-evaluation#alg-d5)).

*Example (rationalizing $1/(1 + \pi + \sqrt 2 - \pi\sqrt 3)$).* In a
fresh complex castle field, $B = 1 + \pi + \sqrt 2 - \pi\sqrt 3$ lies in
the tower $\mathbb{Q}(\pi)(a_1)(a_2)$ with $a_1 = \sqrt 2$, $a_2 = \sqrt
3$, and $1/B$ prints as

```text
(-2*a2*a1*pi^2-2*a2*a1*pi-2*a2*pi^3+2*a2*pi^2+3*a2*pi-4*a1*pi^2-2*a1*pi+a1-2*pi^3+pi-1)/(4*pi^4-8*pi^3-16*pi^2-4*pi+1)
```

with the definitions `a1 = sqrt(2); a2 = sqrt(3)`. The denominator is
obtained by two relative norms. Over $F_1 = \mathbb{Q}(\pi)(\sqrt 2)$,
the conjugate of $B$ is $B' = 1 + \pi + \sqrt 2 + \pi\sqrt 3$, and

$$
N_{F_2/F_1}(B) = BB' = (1 + \pi + \sqrt 2)^2 - 3\pi^2 = 3 + 2\pi - 2\pi^2 + 2(1 + \pi)\sqrt 2 \in F_1
$$

no longer involves $\sqrt 3$. Over $F_0 = \mathbb{Q}(\pi)$, the norm of
this element eliminates $\sqrt 2$:

$$
N_{F_1/F_0}\bigl(N_{F_2/F_1}(B)\bigr) = (3 + 2\pi - 2\pi^2)^2 - 8(1 + \pi)^2 = 4\pi^4 - 8\pi^3 - 16\pi^2 - 4\pi + 1,
$$

the product of the four conjugates of $B$ under the sign changes of
$\sqrt 2$ and $\sqrt 3$, and irreducible over $\mathbb{Q}$. The
numerator is $B'$ times the conjugate $3 + 2\pi - 2\pi^2 - 2(1 +
\pi)\sqrt 2$ of the intermediate norm. The two forms of $1/B$ have the
same value, but only the second satisfies the conditions above.

*Costs and the size limit.* Rationalization trades size for canonicity:
$\operatorname{adj}(M_B) e_0$ is dense in the basis, and the norm has
degree up to $d_k \deg_t B$ in $t$, so that the representation of $1/B$
can be much larger than $B$ (the inverse of a sparse element of a large
cyclotomic tower is dense). The adjugate is therefore computed only
within a size limit (`GR_TOWER_OPT_RATIONALIZE_LIMIT`); beyond it, the
denominator remains an element of $\mathbb{Z}[t][a]$ with numerator and
denominator reduced separately, which is Calcium's representation. Such
fractions are valid but not canonical: equality with them goes through
the zero test of the difference. The benefit is not only canonicity:
without rationalization, sums of many fractions with algebraic
denominators accumulate products that no polynomial gcd can cancel,
because the cancellation happens modulo $I$ and not in $\mathbb{Z}[t,
a]$. One benchmark (input 5 of the DFT benchmark at $N = 16$) went from
6.9 s to 0.28 s when rationalization was introduced during development.

*Example (the size of a rationalized inverse).* Let $\zeta = e^{2\pi
i/17}$, of degree $d = 16$, and $B = \pi + \zeta$, with two terms. The
norm of $B$ is $\prod_{k=1}^{16} (\pi + \zeta^k) = \Phi_{17}(-\pi) =
\pi^{16} - \pi^{15} + \cdots - \pi + 1$, of degree $d \deg_\pi B = 16$,
and the rationalized inverse $1/B$ has a numerator with $136 = 16 \cdot
17/2$ terms, of degree 15 in both $\pi$ and $\zeta$. In a number field
alone the same growth affects the coefficients: $1/(2 + e^{2\pi i/101})$
has all 100 coefficients nonzero, with the common denominator
$\Phi_{101}(-2) = (2^{101} + 1)/3$ and integer numerators of up to 99
bits.

*Comparison with Calcium.* Calcium represents an element of a
multivariate field as an `fmpz_mpoly_q` whose numerator and denominator
are reduced separately modulo the reduction ideal, with their
polynomial gcd cancelled ([](algorithms/exact/field-towers#reduction)).
The representation is not canonical even when the ideal contains all
relations: with $a = \sqrt 2$ and $t = \pi$, both $1/(a + t)$ and $(t -
a)/(t^2 - 2)$ are reduced and gcd-free, and either may occur, so
equality requires a cross multiplication and a reduction, and sums
accumulate algebraic denominators. In Calcium, $1/(\pi + \sqrt 2)$
prints as `(1)/(a+b)` and $(\pi - \sqrt 2)/(\pi^2 - 2)$ as
`(a-b)/(a^2-2)`, with $a = \pi$ and $b = \sqrt 2$ in both: the two
numbers are equal, which Calcium proves, but their representatives
differ. In a castle field both print as `(-a1+pi)/(pi^2-2)` with $a_1 =
\sqrt 2$, and the printed strings are identical. Together with the
triangular set, which makes reduction canonical without Gröbner
completion, rationalization is what makes equality of flat elements
syntactic in proven towers; the zero test
([](algorithms/exact/zero-testing#castle)) is then needed only to
discover relations, not to compare representatives.

*Comparison with rational simplification modulo an ideal.* Monagan and
Pearce [@MP2006] simplify a fraction $a/b$ modulo a prime ideal $I
\subset k[x_1, \ldots, x_n]$ of any dimension, given by a Gröbner basis.
Their first algorithm computes a canonical form for a monomial order:
the generators of the ideal quotient $(\langle b \rangle + I) : a$ give
pairs $[c, d]$ with $ad \equiv bc \bmod I$, and the reduced Gröbner basis
of the module they generate with $I$ (in a term-over-position order)
contains the pair with the least leading monomial. Their second
algorithm finds a fraction of minimal total degree $\deg c + \deg d$ by
solving homogeneous linear systems for increasing degree pairs, with
unknowns the coefficients of $c$ and $d$; its result need not be unique.
Neither algorithm cancels all common factors in general. The two
approaches answer different questions on different data. In a tower,
the transcendence basis $t$ and a triangular presentation of $F$ over
$K = \mathbb{Q}(t)$ are known by construction, so the fraction field of
$\mathbb{Q}[t, a]/I$ is the finite extension $K[a]/I$ and every
denominator can be taken in $K$ ([](#prop-canonical)); the canonical
form is obtained with norms in the tower (determinants or resultants),
and no Gröbner basis of a module or linear system over monomials of
growing degree is needed. For a general prime ideal, Monagan and Pearce's
setting, such a presentation would first have to be computed (a Gröbner
basis in an elimination order), which is the cost that the tower avoids.
Conversely, their forms can be much smaller: on the circle $s^2 + c^2 =
1$, with $c$ transcendental and $s$ algebraic over $\mathbb{Q}(c)$, the
rationalized forms of $(1 - c)/s$ and $1/s$ are $s/(1 + c)$ (also of
minimal degree) and $s/(1 - c^2)$ (degree 3 instead of 1).
Minimal-degree fractions are a natural output form, and a canonical form
of the first kind would be a better fallback than unsimplified fractions
for elements beyond the size limit above; neither is implemented. The
rationalized form is preferred as the working representation because it
is maintained by every arithmetic operation at the cost of one norm per
division, and because equality and hashing become syntactic.

Linear algebra over the castle fields uses division-free algorithms by
default, since inverses of tower elements are large (see
[](algorithms/exact/linear-algebra)).

## Factored denominators {#factored}

In a tower without algebraic generators, a field $\mathbb{Q}(t_1,
\ldots, t_r)$ of rational functions of transcendental generators ($\pi$,
$e$, logarithms, values of special functions), the flat elements need no
reduction, and their arithmetic is that of [](fmpz_mpoly_q), dominated by
the polynomial gcds which keep the fractions canonical. With the option
`GR_TOWER_OPT_FACTORED_DENOMINATORS`, such elements use the factored
variant ([](fmpz_mpoly_q#factored)) instead: the denominator is kept as a
product of irreducible factors, and the gcds become trial divisions by
the known factors. A flat operand enters the factored form at its first
operation (its denominator is factored once), and the results stay
factored; every other routine reads elements through their flat view,
which multiplies the denominator out. An operation which would merge
towers goes through the flat representation.

The option is off by default, because neither form dominates. Measured
against the flat arithmetic over $\mathbb{Q}(\pi)$, $\mathbb{Q}(\pi,
e)$ and $\mathbb{Q}(\pi, e, \log 2)$:

| Workload | Factored time / flat time |
|---|---|
| Berkowitz determinant, $5 \times 5$ reciprocals of linear forms in three generators | 0.54 |
| LU solving, same matrices | 0.38 |
| Determinant, $8 \times 8$ reciprocals in two generators | 1.26 |
| $\sum_{k \le 200} 1/(\pi + k)$ | 15 |
| iterating $x \mapsto (x^2 + \pi)/(x + e)$ | 14 |

Sums of many terms with distinct denominators lose because the factored
denominator accumulates factors; iterated rational maps lose because
their divisions factor the growing numerators of the divisors. Neither
cost guards (on the sizes of the polynomials to factor, or the number of
factors of a sum) nor restricting the option to the linear algebra
routines gave a combination that avoids the losses while keeping the
gains, and the catalog of `src/python/gr_tower_profile.py`, about half
of whose cases do some arithmetic of this kind (mostly multiples of
$\pi$ before an exponential), takes the same time either way. The
option is therefore a choice for applications dominated by linear
algebra over function fields.

## Sharing, copying and collecting towers {#sharing}

Towers are shared by all elements built in them, so a tower can grow
long (many generators accumulated by earlier, unrelated computations).
To keep such history from slowing down new computations, a
transcendental generator ($\exp$, $\log$, $\tan$, $\arctan$) for an
argument stored with a short prefix of a long tower (at least four
generators beyond the prefix) is not appended to the long tower but
adjoined to a copy of the prefix, keeping the definition id and name
of an equal generator already defined elsewhere; the relation search
for it then runs over the generators of the prefix only, and its
values are expressed without the generators of the rest of the tower
(for instance, after sixteen identities for $\sum_{n=1}^N \cos(n a)$
have filled a tower with the tangents of multiples of $a/2$, the
identity $\tan a \tan b = (\cos(a-b) - \cos(a+b)) / (\cos(a-b) +
\cos(a+b))$ for other angles is then found without a relation
search). Radicals (`gr_sqrt`, $n$-th roots) are still
adjoined to the tower of their argument, since the expressions and
factorizations found for them there are too costly to rediscover
(polynomial roots are sought over the generators their coefficients
involve, see
[](algorithms/exact/castle-roots#strategies)), unless the tower continues with at least
twelve generators beyond the prefix (Watson's sum after Dixon's in
one context took four minutes, the roots of its sines appended
after the gamma values of Dixon's sum, which the relation searches
of its merges then visited; a second). For the same reason the tower of
one operand of an operation, containing the definitions of the
other, is used for the result only if no generators of other
computations lie between them ($i \pi$ with $i$ in a tower of its own grown in
place by an earlier computation, say: the result would carry that
computation's generators in its prefix), and the cached generators
(the roots of integers and of unity, algebraic numbers) whose
tower is collected move to a tower where they come first.
Relations between generators of different copies are found when
the copies meet: a merge which adds
transcendental generators to a tower which already has some is
followed by a relation search at the default precision. Likewise,
$\exp(r \pi i + y)$ with rational $r$ is split as $\exp(r \pi i)
\exp(y)$ when the
root of unity $\exp(r \pi i)$ is already in the tower of the
argument, rather than defining a new exponential. When a relation
makes an exponential a radical $x^c = u$ with $u = r w^c$
($r$ rational, $w$ a monomial in the generators, or a square up to
the content for $c = 2$), its modulus is the factor of $x^c - u$
vanishing at it, $w^e f(x/w)$ with $f$ the factor of $y^c - r$ over
$\mathbb{Q}$ (degree 1 when $x/w = \pm 1$).

Every definition (a generator adjoined to some tower, or found
expressible in one) gets a name unique in the context, `a1`,
`a2`, ... for algebraic definitions and `t1`, `t2`, ... for
exponentials and logarithms in order of creation ($\pi$ and $i$
keep their names), which follows the generator into every tower
containing it, so that elements print consistently across towers
and the symbols can be read back; `gr_gens` lists the
definitions made so far which are still present in some tower, in
order. Generators created inside a tower by the relation search (a
primitive exponential $\exp(u/m)$ inserted before its powers, say)
are named when they are printed or listed at the latest.

Elements print with the definitions of the generators they involve, in
dependency order (each definition after those it refers to, so that the
string reads back left to right), for example
`3*t6+t5+2 {t5 = exp(2); t6 = exp(t5)}`, optionally with a numerical
value. Definitions use `sqrt(x)`, `root(x, n)`, `exp(x)`, `log(x)`,
`exp(2*pi*i/n)`, `tan(pi/M)`, special function names, and `root(m(a),
approx)` for a root of a polynomial nearest to a decimal approximation
([](algorithms/exact/castle-roots#printing)). The printed form reads back
through `gr_set_str`, with names local to the string
(`gr_generic_set_str_expr_terminals`), and the symbolic expression output
(`gr_get_fexpr`) produces `Where(expr, Def(a_1, Sqrt(2)), ...)`, from
which LaTeX is obtained. Elements are converted between castle contexts
structurally: generators are rebuilt from their definitions rather than
re-recognized numerically, which matters for roots that agree to many
digits. This printing style follows Calcium's, which prints extension
numbers with their definitions; the round trip through strings is new.

Towers in which no element lives any more are collected: every
element counts as a reference to its tower, and the towers whose
count has dropped to zero are freed at the end of the outermost
operation on the context (in batches). Towers holding the canonical
definitions of the registries ($\pi$ and other constant
transcendentals, roots of unity, radicals of integers) are kept; the
cache of algebraic numbers and the recorded expressions of
definitions are moved to other towers containing them, or dropped.
The definitions of a collected tower which no other tower contains
are forgotten: they no longer appear in `gr_gens`, and their
names can no longer be parsed (printed forms with definitions remain
readable). Memory therefore stays bounded when a long session keeps
few elements alive (solving hundreds of polynomials over
$\mathbb{Q}(\sqrt{2})$ leaves a handful of towers). Towers whose count
reaches zero are queued and collected at the end of the outermost
operation, so that no internal pointer to them is live; the trivial
tower and the towers of the registries of canonical constants are
pinned, since dropping them could duplicate definitions.

*Example (levels and collection).* The element $x = \sqrt 2 + \sqrt 3 +
\pi + e^{\sqrt 2}$ lives in a tower with generators $\sqrt 2$, $\sqrt 3$,
$\pi$, $e^{\sqrt 2}$ (in definition order) and has level 4. The
differences $x - \sqrt 3 - \pi - e^{\sqrt 2}$ and $x - \sqrt 2 - \pi -
e^{\sqrt 2}$ stay in the same tower, with levels 1 and 2, and $\sqrt 2 x -
\sqrt 6 - \sqrt 2\pi - \sqrt 2 e^{\sqrt 2} = 2$ has level 0. In a fresh
context, which holds one tower (the trivial one), creating $\sqrt 2$,
$\sqrt 3$, $\sqrt 2 + \sqrt 3$, $\pi$ and $\sqrt 2 + \sqrt 3 + \pi$ gives 6
towers; deleting the two sums collects their towers (4 remain), and
deleting $\sqrt 2$, $\sqrt 3$ and $\pi$ leaves the count at 4, since
$\mathbb{Q}(\sqrt 2)$, $\mathbb{Q}(\sqrt 3)$ and $\mathbb{Q}(\pi)$ hold
canonical constants and are pinned. Computing the roots of $x^3 - \sqrt 2
x - k$ for $k = 0, \ldots, 199$ and discarding them then leaves 5 towers.

## Threads {#threads}

The context is thread-safe when FLINT is built with pthreads: every
operation on it (including the `gr_ca_*` functions below)
holds a recursive mutex of the context, since operations may refine
or merge the towers of the shared registry and update the cached
representation of elements. Operations on one context are thus
serialized, except for the arithmetic of the two kinds of elements
below; distinct contexts are independent. The lower-level
`gr_ctx_init_tower_field` context is not thread-safe.

*Example.* In a fresh context, $1/3$ is in the rational representation,
and $\sqrt 2$ is created flat, as a generator of its tower. Once it is an
operand of an arithmetic operation it is converted to the dense form (a
polynomial in $\sqrt 2$ over $\mathbb{Q}$), and so is the result $(1 +
\sqrt 2)^{10} = 3363 + 2378\sqrt 2$; $\sqrt 2 + \sqrt 3$ and $\sqrt 2 + \pi$
are flat. Further arithmetic on $1/3$, $\sqrt 2$ and $(1 + \sqrt 2)^{10}$
needs no lock: a Horner step in $\mathbb{Q}(\sqrt{-163})$ costs about
0.2 µs, against 4 µs for a product of two flat elements of
$\mathbb{Q}(\zeta_{60})$ (which has the generators $\zeta_3$, $\zeta_4$,
$\zeta_5$; see [](algorithms/exact/castle-design#open) for the gap between
single- and multi-generator fields).

## Performance {#performance}

`examples/dft.c` runs Calcium's exact DFT benchmark: the DFT of length
$N$ of $x_0, \ldots, x_{N-1}$ followed by the inverse DFT and an exact
comparison with the input. The table gives CPU seconds of castle
arithmetic (`-castle`) and of Calcium for a few lengths (one machine,
indicative only).

| Input $x_n$ | $N$ | Castle | Calcium |
|---|---|---|---|
| 0: $n + 2$ | 16 / 64 | 0.003 / 0.006 | 0.004 / 0.015 |
| 1: $\sqrt{n + 2}$ | 16 / 64 | 0.015 / 1.0 | 0.023 / 1.4 |
| 2: $\log(n + 2)$ | 16 / 64 | 0.006 / 0.14 | 0.016 / 1.4 |
| 3: $e^{2\pi i/(n+2)}$ | 16 / 64 / 128 | 0.011 / 0.68 / 8.2 | 0.073 / fails / — |
| 4: $1/(1 + (n + 2)\pi)$ | 16 / 64 | 0.026 / 16.6 | 0.038 / 18.0 |
| 5: $1/(1 + \sqrt{n + 2}\,\pi)$ | 16 / 32 | 0.13 / 6.1 | 5.9 / > 300 |
| 6: $(n + 2)^{1000}$ | 16 / 64 | 0.003 / 0.067 | 0.002 / 0.20 |
| 7: $(n + 2)^{10000}$ | 16 / 64 | 0.005 / 0.25 | 0.005 / 1.0 |

Input 3 shows the prime-power decomposition of roots of unity: at $N =
16$, the inputs and the powers of $\omega = e^{-2\pi i/16}$ involve the
roots of unity of all orders up to 17, and their sum lies in a tower with
one generator for each of the prime powers 17, 16, 13, 11, 9, 7, 5
($\zeta_8$, $i$ and $\zeta_3$ are powers of $\zeta_{16}$ and $\zeta_9$), of
degree $16 \cdot 8 \cdot 12 \cdot 10 \cdot 6 \cdot 6 \cdot 4 = 2211840$,
which is the degree of $\mathbb{Q}(\zeta_{12252240})$, $12252240 =
\operatorname{lcm}(2, \ldots, 17)$. Elements stay sparse in this tower:
$\omega = -\zeta_{16}^7$, $x_{10} = \zeta_{12} = -i\zeta_3$ and $\omega^5
x_{10} = \zeta_9^3\zeta_{16}^7$ are monomials, and the products by powers
of $\omega$ are shifts ([](#monomial)); Calcium fails on a single
cyclotomic field of huge degree. Input 4 is a pure rational function
computation in $\pi$, with the same `fmpz_mpoly_q` arithmetic as Calcium,
and the times agree. Inputs 6 and 7 have huge integer coefficients, where
Calcium's time goes into full products of number field elements by
powers of $\zeta$ that castle arithmetic performs as shifts. Input 5 at
$N = 64$ is out of reach of both (multivariate gcds in fraction
arithmetic; [](algorithms/exact/castle-design#open)).

`examples/machin.c -castle` proves all the Machin-like formulas
tabulated in [](mp_real) (up to 48 terms) in about three seconds, where
Calcium proves only those with few terms
([](algorithms/exact/field-towers#logarithms)).

## Use in FLINT {#use}

- `gr_ca_repr`, `gr_ca_get_data` and `gr_ca_get_fmpq_poly` give access
  to the representations ([](gr_ca#elements)).
- The options `GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT`,
  `GR_TOWER_OPT_DENSE_FORM_SPARSITY`, `GR_TOWER_OPT_RATIONALIZE_LIMIT`
  and `GR_TOWER_OPT_FACTORED_DENOMINATORS` select the representations
  ([](gr_tower#tuning-options)).
- `src/gr_ca/dense.c`, `src/gr_ca/factored.c` and `src/gr_ca/multi.c`
  implement the dense forms, the factored denominators and the
  arithmetic in products of number fields.
