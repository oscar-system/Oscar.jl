```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Difference polynomial rings](@id differencepolyring)

A difference polynomial ring over the commutative ring ``R`` is an action polynomial ring ``A`` whose action maps are (injective) endomorphisms of ``A``, i.e. ``R``-linear maps that are also multiplicative.

## Construction

```@docs
difference_polynomial_ring(R::Ring, n_elementary_symbols::Int, n_action_maps::Int; kwargs...)
difference_polynomial_ring(R::Ring, x::Symbol, n_action_maps::Int; kwargs...)
```

## Action maps
The action maps of a difference polynomial ring over the commutative ring `R` are injective `R`-algebra homomorphism.

!!! warning
    After calling one of the two following methods, all jet variables that arise within their computation will
    be tracked afterwards.

```@docs
diff_action(p::DifferencePolyRingElem, i::Int)
diff_action(p::DifferencePolyRingElem, d::Vector{Int})
```

