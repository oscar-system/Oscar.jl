```@meta
CurrentModule = Oscar
```

# Cycles and divisors 

## Algebraic Cycles
```@docs 
AbsAlgebraicCycle{CoveredSchemeType<:AbsCoveredScheme, CoefficientRingType<:AbstractAlgebra.Ring}
```
### Constructors
```@docs
algebraic_cycle(X::AbsCoveredScheme, R::Ring)
algebraic_cycle(I::AbsIdealSheaf, R::Ring)
algebraic_cycle(I::AbsIdealSheaf)
```
### Properties
```@docs
ambient_scheme(D::AbsAlgebraicCycle)
components(D::AbsAlgebraicCycle)
dim(D::AbsAlgebraicCycle)
irreducible_decomposition(D::AbsAlgebraicCycle)
integral(W::AbsAlgebraicCycle; check::Bool=true)
```
### Attributes
```@docs 
is_effective(A::AbsAlgebraicCycle)
is_prime(D::AbsAlgebraicCycle)
```
### Methods 
```@docs 
Base.:<=(A::AbsAlgebraicCycle,B::AbsAlgebraicCycle)
```

## Weil Divisors
```@docs
AbsWeilDivisor{CoveredSchemeType, CoefficientRingType}
```
### Constructors
```@docs
weil_divisor(X::AbsCoveredScheme, R::Ring)
weil_divisor(I::AbsIdealSheaf; check::Bool=true)
weil_divisor(I::AbsIdealSheaf, R::Ring; check::Bool=true)
```
### Methods 
Besides the methods for [`AbsAlgebraicCycle`](@ref)
the following are available.
```@docs
is_in_linear_system(f::VarietyFunctionFieldElem, D::AbsWeilDivisor; regular_on_complement::Bool=false, check::Bool=true)
order_of_vanishing(f::VarietyFunctionFieldElem, D::AbsWeilDivisor; check::Bool=true)
intersect(D::AbsWeilDivisor, E::AbsWeilDivisor; covering::Covering=default_covering(scheme(D)))
```

## Linear Systems
```@docs
LinearSystem{DivisorType<:AbsWeilDivisor}
weil_divisor(L::LinearSystem)
variety(L::LinearSystem)
subsystem(L::LinearSystem, D::AbsWeilDivisor)
```

## Cartier Divisors 
```@docs 
CartierDivisor{CoveredSchemeType<:AbsCoveredScheme, CoeffType<:RingElem}
EffectiveCartierDivisor{CoveredSchemeType<:AbsCoveredScheme}
```
Cartier divisors support elementary arithmetic.
### Constructors 
```@docs 
effective_cartier_divisor(I::AbsIdealSheaf; trivializing_covering::Covering = default_covering(scheme(I)), check::Bool = true)
effective_cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
cartier_divisor(E::EffectiveCartierDivisor)
cartier_divisor(IP::AbsProjectiveScheme, f::Union{MPolyDecRingElem, MPolyQuoRingElem})
```
### Attributes
```@docs
ideal_sheaf(C::EffectiveCartierDivisor)
ambient_scheme(C::EffectiveCartierDivisor)
ambient_scheme(C::CartierDivisor)
coefficient_ring(C::CartierDivisor)
components(C::CartierDivisor)
trivializing_covering(C::EffectiveCartierDivisor)
```
