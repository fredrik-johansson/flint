---
title: Exact real and complex numbers
short: Exact numbers
summary: Algebraic numbers and exact real and complex numbers: castle arithmetic in towers of fields, the legacy formal fields of Calcium, and individual algebraic numbers; representations, arithmetic and decision procedures.
order: 14
---

FLINT computes exactly with real and complex numbers in three ways.
*Castle arithmetic* ([](gr_ca), [](guide/types/castle-arithmetic))
represents numbers built from rational numbers by arithmetic, roots,
exponentials, logarithms and special functions as elements of towers of
fields ([](gr_tower)) built on demand; its zero test is complete for
algebraic numbers and decides exp-log numbers by Richardson's algorithm.
The legacy design of Calcium ([](ca)) represents such numbers as
elements of formal fields with ideals of relations [@Joh2021calcium].
And [](qqbar) represents individual algebraic numbers canonically, by
minimal polynomials and isolating balls.

The sections below describe how towers are built and merged
([](algorithms/exact/field-towers)), how algebraic steps are computed
with before they are known to be irreducible
([](algorithms/exact/dynamic-evaluation)), how elements are stored
([](algorithms/exact/castle-representations)), how equality is decided
([](algorithms/exact/zero-testing),
[](algorithms/exact/effort-limits)), how roots, real fields and special
functions are handled, and how linear algebra and simplification work
on top of them. Each page describes the legacy Calcium design in a
section of its own. The [overview](algorithms/exact/overview) compares
the representations and says which algorithm each operation uses, and
[](algorithms/exact/castle-design) gathers the design of castle
arithmetic as a whole: what it inherited from Calcium and what changed,
the lessons of its development, measurements and open problems
[@Joh2026tower].
