# Tropical linear spaces

## Introduction
A tropical linear space is a balanced polyhedral complex supported on a finite intersection of linear tropical hypersurfaces with all multiplicities one.  It is dual to a matroid subdivision of a hypersimplex, and may arise as tropicalizations of linear ideals. For more on tropical linear spaces, see
- Chapter 4.4 in [MS15](@cite)
- Chapter 10 in [Jos21](@cite)

Objects of type `TropicalLinearSpace` need to be embedded, abstract tropical linear spaces are currently not supported.


## Constructors
In addition to converting from `TropicalVariety`, objects of type `TropicalLinearSpace` can be constructed from:
1. Pluecker vectors over a tropical semiring: uses a low-level implementation in `polymake`
2. Pluecker vectors over a field and a tropical semiring map: computes coordinatewise valuation and uses constructor (1.)
3. matrices over a tropical semiring: computes tropical minors and uses constructor (1.)
4. matrices over a field and a tropical semiring map.
    - if matrix over `QQ` and tropical semiring map is trivial, uses an implementation of Rincon's algorithm [Rin13](@cite) in `polymake`
    - for general input, computes minors and uses constructor (2.)
```@docs
tropical_linear_space
```

## Properties
In addition to the properties inherited from `TropicalVariety`, objects of type `TropicalLinearSpace` have the following exclusive properties:
```@docs
pluecker_indices(TropL::TropicalLinearSpace)
tropical_pluecker_vector(TropL::TropicalLinearSpace)
algebraic_pluecker_vector(TropL::TropicalLinearSpace)
tropical_semiring_map(TropL::TropicalLinearSpace)
tropical_matrix(TropL::TropicalLinearSpace)
algebraic_matrix(TropL::TropicalLinearSpace)
algebraic_ideal(TropL::TropicalLinearSpace)
```
