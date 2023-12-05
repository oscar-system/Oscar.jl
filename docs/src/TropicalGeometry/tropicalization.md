# Tropicalization of polynomial ideals

## Introduction
Tropical varieties can arise as tropicalizations of polynomial ideals. For a general introduction, see
- Chapter 3 in [MS15](@cite)

For algorithmic details, see
- [BJSST07](@cite)
- [MR20](@cite)

## Main function
```@docs
tropical_variety(I::MPolyIdeal, nu::Union{TropicalSemiringMap,Nothing}=nothing; weighted_polyhedral_complex_only::Bool=false, skip_saturation::Bool=false, skip_primary_decomposition::Bool=false)
```
