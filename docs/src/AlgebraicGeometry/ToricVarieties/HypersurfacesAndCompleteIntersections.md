# Hypersurfaces and Complete Intersections

## Introduction

A hypersurface or complete intersection in a normal toric variety is described
by homogeneous equations in its Cox ring. The types
`HypersurfaceInToricVariety` and `CompleteIntersectionInToricVariety` retain
these equations and their ambient toric variety. Different collections of
equations may define the same closed subscheme; the chosen presentation
determines the stored equations and their divisor classes.

The ambient fan must have only full-dimensional maximal cones, as required by
the current covered-scheme implementation. The constructors check this even
when `check_codimension=false`; varieties with lower-dimensional maximal
cones are not yet supported.

By default, the constructors check that the associated scheme is either
empty or has codimension equal to the number of equations. Empty zero loci
are valid and accepted. The check uses OSCAR's covered-scheme functionality
and can be expensive. Set `check_codimension=false` to skip it when the
equations are already known to define a complete intersection.
**Without this verification, results relying on the complete-intersection
assumption may be incorrect.**

Neither smoothness nor reducedness is required: the defining equations may
describe a singular or nonreduced scheme. A nonzero constant equation is also
allowed and defines an empty zero locus.

### Complete intersections in singular ambient varieties

On a smooth ambient variety, the defining equations are locally equations of
Cartier divisors. Their expected codimension then certifies that the local
equations form a regular sequence.

On a singular toric variety, a homogeneous Cox-ring equation can instead
describe a Weil divisor that is not Cartier. Such a divisor need not be
locally cut out by one equation on the toric variety. Consequently, the
expected codimension alone does not certify a *local complete intersection*.
For this reason, with `check_codimension=true`, the constructors also require
each defining divisor class to be Cartier unless the ambient variety is
smooth or the zero locus is empty. Setting `check_codimension=false` skips
this certification as well.

## Constructors

```@docs
hypersurface_in_toric_variety(X::NormalToricVariety, equation::MPolyDecRingElem)
complete_intersection_in_toric_variety(X::NormalToricVariety, equations::AbstractVector)
```

## Attributes

```@docs
toric_variety(Y::AbsCompleteIntersectionInToricVariety)
defining_equations(Y::AbsCompleteIntersectionInToricVariety)
defining_equation(Y::HypersurfaceInToricVariety)
defining_divisor_classes(Y::AbsCompleteIntersectionInToricVariety)
defining_divisor_class(Y::HypersurfaceInToricVariety)
defining_ideal(Y::CompleteIntersectionInToricVariety)
covered_scheme(Y::AbsCompleteIntersectionInToricVariety)
closed_embedding(Y::CompleteIntersectionInToricVariety)
```

## Scheme functionality

Standard scheme operations work directly on hypersurfaces and complete
intersections. In particular, the following are available:

- `dim`, `base_ring`, `default_covering`, `affine_charts`
- `is_empty`, `is_smooth`, `is_reduced`, `is_irreducible`, `is_integral`

Use `covered_scheme(Y)` to obtain the corresponding covered scheme and
`Oscar.closed_embedding(Y)` to obtain its embedding into the ambient toric
variety.

## Properties

The following property tests whether the chosen equations have the
expected codimension:

```@docs
has_expected_codimension(Y::CompleteIntersectionInToricVariety)
```

## Fundamental class

For a simplicial ambient toric variety, the product of the defining
divisor classes represents the fundamental class of a complete intersection,
subject to the expected-codimension and Cartier conditions described above.
The fundamental class of an empty zero locus is zero.

```@docs
fundamental_class(Y::AbsCompleteIntersectionInToricVariety)
```
