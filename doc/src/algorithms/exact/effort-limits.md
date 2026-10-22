---
title: Effort limits and the unlimited mode
summary: How castle arithmetic bounds the work of its semi-decision procedures — effort limits read from the options and scaled by an effort round per thread, the recording of exhausted limits, and the unlimited mode, which repeats a failed outermost operation with doubled limits until it succeeds or fails for a reason other than effort, so that no fast heuristic is skipped and no operation gives up prematurely.
complexity: An operation which succeeds in effort round $e$ has repeated the work of rounds $0, \ldots, e-1$, whose limits are smaller by factors of 2; for limits on which the cost depends at least linearly, this at most doubles the cost of the last round.
order: 41
---

The zero test of castle arithmetic ([](algorithms/exact/zero-testing))
is a semi-decision procedure: proving $x \ne 0$ needs a precision of
about $-\log_2 |x|$ bits, and proving $x = 0$ needs relations whose size
is not known in advance. Many other operations are semi-decision
procedures or contain heuristics with bounds: the irreducibility proofs
at places ([](algorithms/exact/dynamic-evaluation#places)), the
factorizations by Trager's method within degree limits, the expression
searches of merges and root finding, the normal forms of special values
at rational arguments up to a level. Each is bounded, and an operation
whose bound is exhausted fails with `GR_UNABLE` or `T_UNKNOWN`. This page
describes how these bounds are organized, and how the *unlimited mode*
removes them.

## Effort limits {#limits}

The bounds which can make an operation fail are *effort limits*. They
are of two kinds:

- the tuning options of the tower ([](gr_tower#tuning-options)) which
  bound an algorithm: the precision limits
  (`GR_TOWER_OPT_PREC_LIMIT`, `GR_TOWER_OPT_CERTIFY_PREC_LIMIT`,
  `GR_TOWER_OPT_NUMERIC_PREC_LIMIT`), the degrees up to which proofs,
  factorizations and minimal polynomials are attempted
  (`GR_TOWER_OPT_TRAGER_DEGREE_LIMIT`, `GR_TOWER_OPT_FACTOR_DEGREE_LIMIT`,
  `GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT`,
  `GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT`), the numbers of tries
  (`GR_TOWER_OPT_MODULAR_TRIES`), the costs of relation searches and
  rationalizations (`GR_TOWER_OPT_RELATION_COST_LIMIT`,
  `GR_TOWER_OPT_RATIONALIZE_LIMIT`, `GR_TOWER_OPT_MODULAR_TERMS_LIMIT`,
  `GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT`) and the levels of the special
  relations (`GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT`,
  `GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT`, `GR_TOWER_OPT_GAMMA_LINE_LIMIT`);
- built-in bounds of the same nature in the implementations: working
  precisions of internal evaluations, the numbers of terms and shifts of
  the hypergeometric normal form, the orders of the polylogarithm
  distribution relations, the conductors of the Dirichlet characters
  expanded exactly.

The options which only choose between algorithms or representations
(the dense form and factored denominator options, the thresholds of
linear algebra) are not effort limits: exceeding them changes how a
result is computed, not whether it is.

Every effort limit is read through one of the macros `GR_TOWER_LIMIT`,
`GR_TOWER_WITHIN` and `GR_TOWER_WITHIN_N` ([](gr_tower#effort)), which
return the value $v$ of the bound times $2^e$, where $e$ is the *effort
round* in progress in the thread (`gr_tower_effort_scale`), saturating at
`WORD_MAX / 64` so that products with small factors do not overflow. A
test that finds a quantity beyond its limit records the exhaustion in a
thread-local flag before giving up. Outside the unlimited mode, the
effort round is always 0: the limits are the values of the options, and
the flag is ignored.

## The unlimited mode {#unlimited}

In the unlimited mode (`GR_TOWER_OPT_UNLIMITED`, set by the flag
`GR_CA_UNLIMITED` of [](gr_ca#unlimited)), every outermost operation on
a castle field runs as follows.

::: algorithm Retry with growing effort {#alg-retry}
**Input:** an operation $f$ of a castle field in the unlimited mode,
called from outside the field.

**Output:** the result of $f$.

1. Set the effort round $e \gets 0$ and clear the exhaustion flag.
2. Run $f$. If it succeeds (a status other than `GR_UNABLE`, a truth
   value other than `T_UNKNOWN`), return its result.
3. If no effort limit was exhausted during the run, or $e$ has reached
   `FLINT_BITS` $- 7$ (where all scaled limits saturate), return the
   failure.
4. Clear the flag, set $e \gets e + 1$ and go to step 2.
:::

The retry happens only at the outermost level: the operations which
an operation calls internally on the same context run within its effort
round, so that one failure does not start nested retries. An operation
which fails for a reason that no effort can remove (a division by zero,
`GR_DOMAIN` outside the domain of a function, an operation not
implemented) records no exhaustion and fails at once, as in the limited
mode.

::: proposition {#prop-unlimited}
In the unlimited mode, an operation fails because of the effort limits
only if it fails in every effort round up to the saturation of the
limits. In particular, a zero test of an element $x \ne 0$ succeeds
whenever the precision at which the enclosure of $x$ excludes zero is
below the saturated precision limit, and a zero test of an algebraic
number always succeeds.
:::

::: proof
A failure is returned in round $e$ only if no limit was exhausted in
that round, which means that the failure does not depend on the limits,
or if $e$ is the last round. A nonzero $x$ is separated from zero at a
finite precision ([](algorithms/exact/zero-testing#numerical)), which a
round reaches since the precision limits double with each round; the
zero test of algebraic numbers by dynamic evaluation needs no limit
other than the precisions of the refinements and the proofs
([](algorithms/exact/dynamic-evaluation#thm-complete)).
:::

In practice, the saturation is never reached: the limits grow
exponentially with the rounds, and an operation which still fails after
a few rounds runs out of time or memory instead. The unlimited mode thus
turns the zero test into a procedure that terminates exactly when a
proof exists within the methods of castle arithmetic (and, for exp-log
constants, under Schanuel's conjecture), as for [](qqbar), which also
never aborts prematurely. It does not decide more identities than the
limited mode would with large enough options: an identity outside the
known relations ([](algorithms/exact/zero-testing#decided)) is not
decided, and its zero test runs forever.

## Fast heuristics stay fast {#heuristics}

A simpler way to remove the limits would be to set every option to
infinity. This would be much slower, because the operations are
sequences of heuristics of increasing cost — try A up to its bound, then
B up to its bound, then the general method — and infinite bounds would
make the first heuristic run forever on inputs that the second decides
at once. With growing effort rounds, every heuristic of a sequence is
tried in each round with bounds that grow together, so that the cheap
methods keep deciding what they decide in round 0, at the cost of round
0, and an operation that needs a larger bound for one method pays for
the earlier rounds only. For the same reason, bounded loops over
increasing precisions or degrees inside a heuristic are written as loops
which continue past their nominal bound only when an exhaustion was
recorded, rather than with a fixed number of steps.

The options remain in force in the unlimited mode, as the starting
values of the limits: lowering one makes the first rounds cheaper and
the retries more frequent; raising one avoids rounds for inputs known to
need it.

## Measurements {#measurements}

With all precision limits set to 64 bits, $1 - \operatorname{erf}(30)
\approx 2.6 \cdot 10^{-393}$ is not separated from zero in the limited
mode (`T_UNKNOWN`), and is proved nonzero in the unlimited mode after
five rounds. On the catalog of cases of `src/python/gr_tower_profile.py`
with all limits lowered by the same factor, the limited mode leaves 12
predicates undecided and 3 operations unable; the unlimited mode decides
or computes all of them, while 12 cases run past the time limit instead of
returning `T_UNKNOWN`. With the default options, the unlimited mode
costs nothing for operations which exhaust no limit.

## Use in FLINT {#use}

- `GR_TOWER_OPT_UNLIMITED` and the flag `GR_CA_UNLIMITED` select the
  unlimited mode ([](gr_tower#tuning-options), [](gr_ca#unlimited)).
- `GR_TOWER_LIMIT`, `GR_TOWER_WITHIN`, `GR_TOWER_WITHIN_N`,
  `GR_TOWER_UNLIMITED` and `gr_tower_effort_scale` read the effort
  limits ([](gr_tower#effort)).
- The retry of [](#alg-retry) is performed by the locking wrappers of
  the castle fields (`src/gr_ca/impl.h`); the verbosity option
  `GR_TOWER_OPT_VERBOSE` at level 2 prints the effort rounds.
