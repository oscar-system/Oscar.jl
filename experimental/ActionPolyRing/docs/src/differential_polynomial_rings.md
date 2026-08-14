# [Differential polynomial rings](@id differentialpolyring)

A differential polynomial ring over the commutative ring ``R`` is an action polynomial ring ``A`` whose action maps are
derivations of ``A``, i.e. ``R``-linear maps that also satisfy the Leibniz-rule.

## [Construction](@id differentialpolyring_construction)
We provide the following methods to construct derivations on the coefficient ring `R`. Using these is necessary, if
one want to use differential polynomial rings with nontrivial derivations.

```@docs
action_derivation(R::Ring)
action_derivation(m::Map{D, D}; check::Bool=true) where {D <: Ring}
```

We provide the following constructors to create differential polynomial rings where all derivations are trivial,
i.e. the zero map on the coefficient ring.

```@docs
differential_polynomial_ring(R::Ring, n_action_indeterminates::Int, n_action_maps::Int; kwargs...)
differential_polynomial_ring(R::Ring, x::Symbol, n_action_maps::Int; kwargs...)
```

We provide the following constructors to create differential polynomial rings with arbitrary commuting derivations. We verify that
the derivations commute on generators of the coefficient ring `R`. However, if `gens` is not applicable to it then no check is performed at all.

```@docs
differential_polynomial_ring(R::D, n_action_indeterminates::Int, action_maps::Vector{<:ActionDerivation{D}}; kwargs...) where {D <: Ring}
differential_polynomial_ring(R::D, action_indeterminate::Symbol, action_maps::Vector{<:ActionDerivation{D}}; kwargs...) where {D <: Ring}
```

## Applying derivations
You can use the following methods to apply the derivations of a differential polynomial ring to its elements:

```@docs
apply_action(p::DifferentialPolyRingElem, i::Int)
apply_action(p::DifferentialPolyRingElem, d::Vector{Int})
```

!!! warning
    After calling one of these methods, all jet variables that arise within their computation will
    be tracked afterwards.

