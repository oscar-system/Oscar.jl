# [Difference polynomial rings](@id differencepolyring)

A difference polynomial ring over the commutative ring ``R`` is an action polynomial ring ``A`` whose action maps
are (injective) endomorphisms of ``A``, i.e. ``R``-linear maps that are also multiplicative. We also call them
shift operators.

## [Construction](@id differencepolyring_construction)
We provide the following methods to construct shift operators on the coefficient ring `R`. Using these is necessary, if
one want to use difference polynomial rings with nontrivial shift operators.

```@docs
action_shift(R::Ring)
action_shift(m::Map{D, D}; check::Bool=true) where {D <: Ring}
```

We provide the following methods to create difference polynomial rings where all shift operators are trivial,
i.e. the identity map on the coefficient ring.

```@docs
difference_polynomial_ring(R::Ring, n_action_indeterminates::Int, n_action_maps::Int; kwargs...)
difference_polynomial_ring(R::Ring, x::Symbol, n_action_maps::Int; kwargs...)
```

We provide the following constructors to create difference polynomial rings with arbitrary commuting shift operators. We verify that
the shift operators commute on generators of the coefficient ring `R`. However, if `gens` is not applicable to it then no check is performed at all.

```@docs
difference_polynomial_ring(R::D, n_action_indeterminates::Int, action_maps::Vector{<:ActionShift{D}}; kwargs...) where {D <: Ring}
difference_polynomial_ring(R::D, action_indeterminate::Symbol, action_maps::Vector{<:ActionShift{D}}; kwargs...) where {D <: Ring}
```

## Applying shift operators
You can use the following methods to apply the shift operators of a difference polynomial ring to its elements:

```@docs
apply_action(p::DifferencePolyRingElem, i::Int)
apply_action(p::DifferencePolyRingElem, d::Vector{Int})
```

!!! warning
    After calling one of these methods, all jet variables that arise within their computation will
    be tracked afterwards.
