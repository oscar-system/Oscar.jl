### Constructors
@doc raw"""
    action_derivation(R::Ring)

Construct the zero derivation on the ring `R`.
"""
action_derivation(R::Ring) = TrivialActionDerivation{typeof(R)}(R)

@doc raw"""
    action_derivation(m::Map{D, D}; check::Bool=true) where {D <: Ring}

Wrap the map `m` into an `ActionDerivation`. If `check` is true, a heuristic validation on the generators is performed.
"""
function action_derivation(m::Map{D, D}; check::Bool=true) where {D <: Ring}
  check && @req __is_probably_valid_derivation(m) "The provided map fails the Leibniz rule or additivity on generators; it is not a valid derivation"
  
  __is_trivial_derivation(m) && return TrivialActionDerivation{D}(domain(m))
  return NontrivialActionDerivation{D}(m)
end

@doc raw"""
    action_shift(R::Ring)

Construct the trivial shift, i.e. the identity map on the ring `R`.
"""
action_shift(R::Ring) = TrivialActionShift{typeof(R)}(R)

@doc raw"""
    action_shift(m::Map{D, D}; check::Bool=true) where {D <: Ring}

Wrap the map `m` into an `ActionShift`. If `check` is true, a heuristic validation on the generators is performed.
"""
function action_shift(m::Map{D, D}; check::Bool=true) where {D <: Ring}
  check && @req __is_probably_valid_shift(m) "The provided map fails multiplicativity or additivity on generators; it is not a valid shift"
  
  __is_trivial_shift(m) && return TrivialActionShift{D}(domain(m))
  return NontrivialActionShift{D}(m)
end

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

### Triviality Helpers
function __is_trivial_shift(m::Map{D, D}) where {D <: Ring}
  m isa AbstractAlgebra.Generic.IdentityMap && return true
  
  R = domain(m)
  !applicable(gens, R) && return false 
  
  return all(g -> m(g) == g, gens(R))
end

function __is_trivial_derivation(m::Map{D, D}) where {D <: Ring}
  R = domain(m)
  !applicable(gens, R) && return false 
  return all(g -> iszero(m(g)), gens(R))
end

### Verifiers
function __is_probably_valid_shift(m::Map{D, D}) where {D <: Ring}
  m isa AbstractAlgebra.Generic.IdentityMap && return true
  R = domain(m)
    
  (!is_one(m(one(R))) || !is_zero(m(zero(R)))) && return false
  
  !applicable(gens, R) && return true # Impossible to check so trust the user
  G = gens(R)
  mG = [m(g) for g in G]
  n = length(G)
  @inbounds for i in 1:n
    @inbounds for j in i:n
      (m(G[i] + G[j]) != mG[i] + mG[j]) && return false
      (m(G[i] * G[j]) != mG[i] * mG[j]) && return false
    end
  end
  return true
end

function __is_probably_valid_derivation(m::Map{D, D}) where {D <: Ring}
  R = domain(m)
    
  !is_zero(m(one(R))) && return false
  
  !applicable(gens, R) && return true # Impossible to check so trust the user
  G = gens(R)
  mG = [m(g) for g in G]
  n = length(G)
  @inbounds for i in 1:n
    @inbounds for j in i:n
      (m(G[i] + G[j]) != mG[i] + mG[j]) && return false
      (m(G[i] * G[j]) != mG[i] * G[j] + G[i] * mG[j]) && return false
    end
  end
  
  return true
end

