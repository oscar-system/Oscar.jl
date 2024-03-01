# Groebner theory

## Introduction
Tropical algebraic geometry incorporates the valuation of the underlying ground field, and therefore so does its Groebner theory. Tropical Groebner theory is a generalization of its classical counterpart to fields with valuation, and the classical Groebner theory is a specialization of tropical Groebner theory to fields with trivial valuation.  Instead of monomial orderings there are term orderings which take the valuation into account, and initial forms and ideals live over the residue field.  For details, see Chapter 2.4 in [MS15](@cite).

## Groebner bases
Groebner bases in [MS15](@cite) are only defined for homogeneous ideals and they are finite sets whose initial forms generate the initial ideal.  Groebner bases in OSCAR are defined for all ideals and they are finite generating sets whose initial forms generate the initial ideal. For homogeneous ideals generating the initial ideal implies generating the original ideal, so both notions coincide.

No general algorithm exists for computing Groebner bases for inhomogeneous ideals.  The command `groebner_basis` computes a Groebner basis if it knows how, and raises an error if it does not.  The cases in which a Groebner basis of an ideal `I` with respect to the tropical semiring map `nu` and weight vector `w` is returned are:
* `I` looks principal
* `I` looks binomial and `w` lies in `tropical_variety(I)`
* `I` looks affine linear and `nu` is trivial
* `nu` is trivial and (`w` is negative if `convention(nu)==min` or `w` is positive if `convention(nu)==max`)
* `I` is homogeneous
Here, "looks" means `gens(I)` having the desired property.  We are not computing a classical reduced Groebner basis to verify that `I` actually has the desired property.
```@docs
groebner_basis(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}})
```

## Initial forms and initial ideals
```@docs
initial(f::MPolyRingElem, nu::TropicalSemiringMap, w::Vector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}})
initial(I::MPolyIdeal, nu::TropicalSemiringMap, w::Vector{<:Union{QQFieldElem,ZZRingElem,Rational,Integer}}; skip_groebner_basis_computation::Bool=false)
```
