# Groebner theory

## Introduction
Groebner bases and related notions in tropical geometry are a generalization of its classical counterparts to fields with valuation.  The Groebner bases take the valuation into account, the initial ideals live over the residue field, and consequently the Groebner complex is not necessarily a polyhedral fan.  For details, see
- Chapter 2.4 and 2.5 in [MS15](@cite)

## Groebner bases
```@docs
groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}})
```

## Initial forms and initial ideals
```@docs
initial(f::MPolyRingElem, nu::TropicalSemiringMap, w::Vector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}})
initial(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}}; skip_groebner_basis_computation::Bool=false)
```
