---
title: Linear algebra over exact numbers
summary: Gaussian elimination, fraction-free elimination, division-free determinants and characteristic polynomials, eigenvalues, Jordan forms and matrix functions over exact real and complex numbers, where every pivot decision is a zero test; the choice of algorithms and the packed products of castle fields, and how the legacy ca_mat builds on the generic gr_mat algorithms with its own pivoting and specialisations.
complexity: Elimination costs $O(n^3)$ field operations and up to $O(n^2)$ zero tests; the Berkowitz algorithm costs $O(n^4)$ ring operations and no zero tests; eigenvalue and Jordan form computations add root finding for the characteristic polynomial and $O(n)$ rank computations per eigenvalue.
order: 60
---

Linear algebra over a field needs to decide, at each step of Gaussian
elimination, whether a candidate pivot is zero. Over exact real and
complex numbers this decision is a zero test (see
[](algorithms/exact/zero-testing)), and arithmetic is expensive and
grows the representations. FLINT relies on the generic algorithms of
[](gr_mat), which report a failed zero test through their status codes,
combined with algorithms that avoid zero tests, and with specialised
arithmetic when the entries permit: in the castle fields ([](gr_ca)),
through methods of the context ([](#castle)), and in the legacy
[](ca_mat) module through its own wrappers ([](#legacy)).

## Castle fields {#castle}

*Example (one determinant over two exact fields).* The generic
`gr_mat_det`, applied to

$$
M = \begin{pmatrix} 1 & \sqrt 2 & \pi \\ \sqrt 2 & \pi & 1 \\ \pi & 1 & \sqrt 2 \end{pmatrix}
$$

over the complex castle field, returns $\det M = -\pi^3 + 3\sqrt 2\pi -
2\sqrt 2 - 1$, printed as `-pi^3+3*pi*a1-2*a1-1 {a1 = sqrt(2)}`: an
element of a tower of degree 2 over $\mathbb{Q}(\pi)$. The same call over
Calcium's `gr` wrapper returns the same polynomial in $\pi$ and $\sqrt 2$.
The determinant code only calls the arithmetic and the predicates of the
context.

The generic functions of [](gr_mat) (`gr_mat_det`, `gr_mat_rank`,
`gr_mat_rref`, `gr_mat_nonsingular_solve`, `gr_mat_charpoly`,
`gr_mat_eigenvalues`, `gr_mat_jordan_blocks`, ...) work over a castle
field directly. Since the zero test of algebraic numbers is complete
([](algorithms/exact/dynamic-evaluation#thm-complete)), every pivot
decision over algebraic entries succeeds, and the rank, the reduced row
echelon form and the solvability of a system of algebraic numbers are
always decided; with transcendental entries, a pivot search fails
(`GR_UNABLE`) only when every candidate is undecided. The context
overrides three methods of the matrix interface:

- *Multiplication* (`gr_mat_mul`, from $3 \times 3$): over the products
  of number fields of [](algorithms/exact/castle-representations#number-fields),
  the entries are packed into integer polynomials with common
  denominators by rows and columns and multiplied by one
  `fmpz_poly_mat_mul`; products with rational matrices are scalar
  multiplications. Polynomial products (`gr_poly_mullow`, from length
  4) use the same packing.
- *Determinants*: rational matrices go to `fmpq_mat_det`
  ([](algorithms/linear-algebra/multimodular)); matrices up to $4
  \times 4$ use the cofactor expansion; over a number field of degree
  $D$ with dense arithmetic, LU decomposition (recursive, with the fast
  products) when $D < n$, and the division-free Berkowitz algorithm
  ([](#division-free)) otherwise, also for entries with transcendental
  generators, where elimination accumulates denominators.
- *Linear systems* (`gr_mat_nonsingular_solve`): fraction-free LU
  ([](#fflu)) over number fields of degree $D \ge n/2$, LU otherwise.

The thresholds come from measurements on random matrices with small
entries (`examples/number_field_bench.c`): over $\mathbb{Q}(\zeta_{64})$,
a $16 \times 16$ determinant is more than an order of magnitude faster
by Berkowitz than by LU, since the fractions of elimination have large
heights when the degree is large compared to the size of the matrix,
while over $\mathbb{Q}(\sqrt{-163})$ LU is fastest. Over fields of
rational functions in transcendental generators, the factored
denominators of [](algorithms/exact/castle-representations#factored)
make Berkowitz determinants about twice as fast and LU solving about
two and a half times as fast, at the expense of other workloads, and are
therefore optional.

Over general towers (with transcendental generators, or through the
nested representation), division-free algorithms are mandatory: in
early experiments, Gaussian elimination with divisions was catastrophic
in every representation, each division by a pivot producing an inverse
with large coefficients that the subsequent operations never cancel.
Over number fields with the dense kernels, inverses are cheap enough
that the choice depends on $n$ and $D$. The thresholds were measured
during development: for a $16 \times 16$ determinant over
$\mathbb{Q}(\zeta_{64})$, elimination took 1.5 s, fraction-free
elimination 0.46 s and Berkowitz 0.07 s; over $\mathbb{Q}(\sqrt{-163})$
at $40 \times 40$, elimination took 0.16 s and Berkowitz 1.2 s. Calcium
uses elimination for number field entries throughout; its $6 \times 6$
determinant over $\mathbb{Q}(\zeta_{128})$ took 1.6 s against 2.7 ms in
a castle field.

*Example (heights in elimination).* Let $M$ be a $10 \times 10$ matrix
over $\mathbb{Q}(\zeta_{32})$ ($D = 16$) whose entries have all 16
coefficients drawn at random from $\{-3, \ldots, 3\}$. Its determinant
has integer coefficients of up to 13 digits. Gaussian elimination
produces the pivots $\Delta_k/\Delta_{k-1}$ (quotients of leading
principal minors), whose denominators divide the norms of the
$\Delta_{k-1}$, of degree 16 in their coefficients: they have 14, 29, 43,
67, 86, 106, 127, 148 and 171 digits for $k = 2, \ldots, 10$. Since $D
\ge n$, castle arithmetic uses Berkowitz's algorithm, which stays in
$\mathbb{Z}[\zeta_{32}]$; on one machine, elimination took 26 ms,
fraction-free elimination 9 ms and Berkowitz's algorithm 2.6 ms. In the
other direction, the $7 \times 7$ Fourier matrix $(\zeta_7^{jk})$ has $D =
6 < n = 7$ and goes to elimination; its determinant is returned as
$343(1 + 2\zeta_7 + 2\zeta_7^2 + 2\zeta_7^4) = 343\sqrt{-7}$.

## Fraction-free elimination {#fflu}

Gaussian elimination with divisions produces fractions whose numerators
and denominators grow. Bareiss's algorithm [@Bar1968] avoids this: with
$a^{(0)} = A$ and $a_{-1,-1}^{(-1)} = 1$ it updates, for $i, j > k$,

$$
a_{ij}^{(k+1)} = \frac{a_{kk}^{(k)} a_{ij}^{(k)} - a_{ik}^{(k)} a_{kj}^{(k)}}{a_{k-1,k-1}^{(k-1)}},
$$

where by Sylvester's identity every $a_{ij}^{(k)}$ is a $(k+1) \times
(k+1)$ minor of the input, so the division is exact (see
[](algorithms/linear-algebra/fraction-free) for the integer case). `gr_mat_fflu` implements
this for any field, and `ca_mat_fflu` with `ca_mat_find_pivot` for
pivoting; the last pivot is the determinant (up to sign) and is returned
as the denominator. Over a
field such as $\mathbb{Q}(\pi)$ the entries are rational functions in
the generators, the exact divisions keep them polynomial in the input
entries, and the representations stay smaller than with ordinary
elimination. In castle fields the exact division is a division of
polynomials in the flat representation, which is exact as well; in
[](ca) it is a field division, so fractions can appear when the formal
representation does not see that it is exact.

## Algorithms without zero tests {#division-free}

Some quantities can be computed without any zero test, which makes them
always succeed:

- The Berkowitz algorithm [@Ber1984] computes the characteristic
  polynomial with $O(n^4)$ ring operations and no divisions, from
  products of Toeplitz matrices built from the leading principal
  submatrices (`gr_mat_charpoly_berkowitz`); determinants over castle
  fields with transcendental generators or of large degree use it
  ([](#castle)), `ca_mat_charpoly` uses it by default and
  `ca_mat_det_berkowitz` takes the constant coefficient.
- Cofactor expansion (`ca_mat_det_cofactor`, `ca_mat_adjugate_cofactor`)
  is used for matrices of size up to 4 or 5.
- The adjugate is obtained from the characteristic polynomial $p(x) =
  \det(xI - A) = x q(x) + p(0)$ as $\operatorname{adj}(A) = (-1)^{n+1}
  q(A)$, by the Cayley–Hamilton theorem (`ca_mat_adjugate_charpoly`).
  Then $A^{-1} = \operatorname{adj}(A)/\det A$ and $A^{-1} B$ need a
  single zero test, of the determinant, which is how `ca_mat_inv` and
  `ca_mat_nonsingular_solve_adjugate` proceed for small matrices or
  matrices with entries in several fields.

The Danilevsky algorithm [@Dan1937], which reduces $A$ to companion form
by similarity transformations with $O(n^3)$ operations, is faster but
needs pivots and therefore zero tests (`ca_mat_charpoly_danilevsky`).
See [](algorithms/linear-algebra/charpoly) for these algorithms over
general rings.

## Eigenvalues and the Jordan form {#jordan}

`gr_mat_eigenvalues` computes the characteristic polynomial and its
roots with multiplicities. Over a castle field the roots come from
`gr_poly_roots` ([](algorithms/exact/castle-roots)), which adjoins the
roots of factors of any degree as algebraic generators, so that the
eigenvalues of any matrix of algebraic numbers are found (the
symmetric functions of the eigenvalues of a rational matrix then reduce
to its coefficients in a splitting tower). `ca_mat_eigenvalues` uses
`ca_poly_roots` ([](ca_poly/roots)) and fails if a factor of degree 4 or
more with non-rational coefficients occurs; for an integer or rational
matrix it always succeeds, with algebraic eigenvalues.

Given the eigenvalues $\lambda$ with algebraic multiplicities $m$,
`ca_mat_jordan_blocks` (through `gr_mat_jordan_blocks`) computes the
ranks $r_k = \operatorname{rank}(A - \lambda I)^k$ for $k = 0, 1, \ldots$
until they stabilise at $n - m$. The number of Jordan blocks of size at
least $k$ for $\lambda$ is $r_{k-1} - r_k$, so the sizes of the blocks
form the conjugate of the partition $(r_0 - r_1, r_1 - r_2, \ldots)$
(the transpose of its Ferrers diagram). `ca_mat_jordan_transformation`
then constructs a basis of generalised eigenvectors from nullspaces, and
`ca_mat_diagonalization` checks whether every geometric multiplicity
$\dim \ker(A - \lambda I)$ equals the algebraic multiplicity.

Matrix functions are computed from the Jordan decomposition $A = P J
P^{-1}$: for a Jordan block $J_\lambda$ of size $s$,

$$
f(J_\lambda) = \sum_{k=0}^{s-1} \frac{f^{(k)}(\lambda)}{k!} N^k,
$$

where $N$ is the nilpotent part, so that $f(A) = P f(J) P^{-1}$ only
needs $f(\lambda)$ and its derivatives at the eigenvalues
(`gr_mat_func_jordan`, used by `ca_mat_exp` with $f = \exp$ and by
`ca_mat_log` with the principal logarithm, which fails for singular
matrices). Polynomials in a matrix
are evaluated with the Paterson–Stockmeyer algorithm [@Paterson1973]
(`ca_mat_ca_poly_evaluate`), with about $2\sqrt{\deg f}$ matrix
multiplications.

## The legacy ca_mat (Calcium) {#legacy}

### Generic algorithms {#generic}

Most functions of [](ca_mat) are wrappers that call the corresponding
`gr_mat` function with a `ca` context: `ca_mat_lu` calls
`gr_mat_lu_generic` (classical Gaussian elimination below dimension 5,
block recursive elimination above), `ca_mat_rank`, `ca_mat_rref`,
`ca_mat_right_kernel`, `ca_mat_nonsingular_solve_lu`, the Berkowitz and
Danilevsky characteristic polynomials, the Jordan form and the matrix
exponential and logarithm all come from [](gr_mat). The `ca` context
overrides three methods of the generic matrix interface:
multiplication (`ca_mat_mul`), determinants (`ca_mat_det`) and the pivot
search (`ca_mat_find_pivot`), so that the generic algorithms profit from
the specialisations below. A failed zero test makes a generic function
return `GR_UNABLE`, which the `ca_mat` wrappers translate into a return
value 0 or `T_UNKNOWN`.

### Choosing pivots {#pivots}

`ca_mat_find_pivot` searches column $j$ for a nonzero entry in two
passes. The first pass ignores entries that are rational zeros and takes
the remaining entry of least complexity of representation (according to
`ca_cmp_repr`): a simple entry, such as a small integer or an element of
a small number field, is cheap to test and likely to be nonzero, so a
single zero test usually suffices. If that entry is not proved nonzero,
the second pass tests every entry and takes the simplest one that is
proved nonzero. Entries that are proved zero along the way are replaced
by an exact zero, which simplifies later steps. If no entry is proved
nonzero and not all are proved zero, the search returns `T_UNKNOWN` and
the elimination fails: the rank is then unknown.

### Rational and number field matrices {#specialised}

Matrices whose entries are all represented as rational numbers are
recognised with `_ca_mat_same_field` and the computation is delegated
to [](fmpz_mat) or [](fmpq_mat): `ca_mat_mul` (for inner dimension at
least 3) and `ca_mat_det` (for size at least 3), with their multimodular
algorithms (see [](algorithms/linear-algebra/multimodular)). Matrices
whose entries lie in one number field $\mathbb{Q}(a)$ are multiplied by
`ca_mat_mul_same_nf` (when the dimensions are at least 3 and the inner
dimension at least 4), which writes the entries as
polynomials in $a$ over a common denominator, multiplies the resulting
integer polynomial matrices with `fmpz_poly_mat_mul` and reduces modulo
the minimal polynomial; for such matrices, determinants and inverses use
LU decomposition, since zero tests in number fields never fail.

### Use in FLINT (legacy) {#legacy-use}

- `ca_mat_find_pivot`, `ca_mat_fflu`, `ca_mat_rref_fflu` and
  `ca_mat_nonsingular_solve_fflu` are specific to [](ca_mat); `ca_mat_lu`,
  `ca_mat_rank`, `ca_mat_rref`, `ca_mat_right_kernel` and
  `ca_mat_nonsingular_solve_lu` wrap the generic [](gr_mat) algorithms.
- `ca_mat_charpoly`, `ca_mat_det`, `ca_mat_adjugate` and `ca_mat_inv`
  choose between the methods above; `ca_mat_mul` and
  `ca_mat_mul_same_nf` specialise multiplication.
- `ca_mat_eigenvalues`, `ca_mat_diagonalization`, the Jordan form
  functions, `ca_mat_exp` and `ca_mat_log` build on root finding and
  nullspaces.

## Use in FLINT {#use}

- `gr_mat_det`, `gr_mat_nonsingular_solve` and `gr_mat_mul` over the
  castle fields dispatch to the methods of [](#castle)
  (`src/gr_ca/dense.c`, `src/gr_ca/multi.c`); the other functions of
  [](gr_mat) are the generic ones.
- `gr_mat_fflu`, `gr_mat_nonsingular_solve_fflu`,
  `gr_mat_charpoly_berkowitz` and `gr_mat_det_berkowitz` are the
  fraction-free and division-free algorithms.
- `qqbar_eigenvalues_fmpz_mat` computes eigenvalues of integer matrices
  directly as algebraic numbers.
