### Constructors
@doc raw"""
    action_derivation(R::Ring)

Construct the zero derivation on the ring `R`.
"""
action_derivation(R::Ring) = TrivialActionDerivation{typeof(R)}(R)

@doc raw"""
    action_derivation(m::Map{D, D}) where {D <: Ring}

Wrap the map `m` into an `ActionDerivation`. This does not check whether `m` is actually a derivation.
"""
action_derivation(m::Map{D, D}) where {D <: Ring} = NontrivialActionDerivation{D}(m)

@doc raw"""
    action_shift(R::Ring)

Construct the trivial shift, i.e. the identity map on the ring `R`.
"""
action_shift(R::Ring) = TrivialActionShift{typeof(R)}(R)

@doc raw"""
    action_shift(m::Map{D, D}) where {D <: Ring}

Wrap the map `m` into an `ActionShift`. This does not check whether `m` is actually a shift operator.
"""
action_shift(m::Map{D, D}) where {D <: Ring} = NontrivialActionShift{D}(m)

### Getters
domain(m::Union{TrivialActionDerivation, TrivialActionShift}) = m.domain
codomain(m::Union{TrivialActionDerivation, TrivialActionShift}) = domain(m)

__underlying_map(m::Union{NontrivialActionDerivation, NontrivialActionShift}) = m.underlying_map

### Standard functionality
domain(m::Union{NontrivialActionDerivation, NontrivialActionShift}) = domain(__underlying_map(m))
codomain(m::Union{NontrivialActionDerivation, NontrivialActionShift}) = codomain(__underlying_map(m))

function (m::TrivialActionDerivation)(x)
  R = domain(m)
  _ = R(x) # Check coercivity
  return zero(R)
end
(m::TrivialActionShift)(x) = domain(m)(x)
(m::Union{NontrivialActionDerivation, NontrivialActionShift})(x) = __underlying_map(m)(domain(m)(x))

