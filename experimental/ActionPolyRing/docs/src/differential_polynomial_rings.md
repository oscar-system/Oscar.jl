```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Differential polynomial rings](@id differentialpolyring)

## Construction

```@docs
differential_polynomial_ring
```

!!! warning
    After calling one of the two following methods, all jet variables that arise within their computation will
    be tracked afterwards.

## Action maps
The action maps of a differential polynomial ring over the commutative ring `R` are `R`-linear derivations.

```@docs
diff_action(p::DifferentialPolyRingElem, i::Int)
diff_action(dpre::DifferencePolyRingElem{T}, d::Vector{Int}) where {T}
```
