# Tropical hypersurfaces

## Introduction
A tropical hypersurface is a balanced polyhedral complex of codimension one.  It is dual to a regular subdivision of a Newton polytope. For more on tropical hypersurfaces, see
- Chapter 3.1 in [MS15](@cite)
- Chapter 1 in [Jos21](@cite)

#### Note:
- Objects of type `TropicalHypersurface` need to be embedded, abstract tropical hypersurfaces are currently not supported.
- The type `TropicalHypersurface` can be thought of as subtype of `TropicalVariety` in the sense that it should have all properties and features of the latter.


## Constructors
In addition to converting from `TropicalVariety`, objects of type `TropicalHypersurface` can be constructed from:
- polynomials over a tropical semiring,
- polynomials over a field and a tropical semiring map,
- subdivision of points and a choice of min- or max-convention.
```@docs
tropical_hypersurface
```

## Properties
In addition to the properties inherited from `TropicalVariety`, objects of type `TropicalHypersurface` have the following exclusive properties:
```@docs
algebraic_polynomial(TropH::TropicalHypersurface)
tropical_polynomial(TropH::TropicalHypersurface)
dual_subdivision(TropH::TropicalHypersurface)
```
