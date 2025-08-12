```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Tropical Groebner theory

## Introduction
Tropical algebraic geometry incorporates the valuation of the underlying ground field, and therefore so does its Groebner theory. Tropical Groebner theory is a generalization of its classical counterpart to fields with valuation, and the classical Groebner theory is a specialization of tropical Groebner theory to fields with trivial valuation.  Instead of monomial orderings there are term orderings which take coefficient valuations into account, and initial forms and ideals live over the residue field.  For details, see Chapter 2.4 in [MS15](@cite).

## Groebner bases
In [MS15](@cite), (tropical) Groebner bases are only defined for homogeneous ideals and they are finite sets whose initial forms generate the initial ideal. In OSCAR, Groebner bases are defined for all polynomial ideals and they are finite generating sets whose initial forms generate the initial ideal. For homogeneous ideals generating the initial ideal implies generating the original ideal, so both notions coincide. Currently, the computation of Groebner bases is only supported for homogeneous and principal ideals.

```@docs
groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}})
```

## Initial forms and initial ideals
```@docs
initial(f::MPolyRingElem, nu::TropicalSemiringMap, w::AbstractVector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}})
initial(I::MPolyIdeal, nu::TropicalSemiringMap, w::AbstractVector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}}; skip_groebner_basis_computation::Bool=false)
```
