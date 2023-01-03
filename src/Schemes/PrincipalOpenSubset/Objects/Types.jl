export PrincipalOpenSubset

########################################################################
# Principal open subsets of affine schemes                             #
########################################################################
@attributes mutable struct PrincipalOpenSubset{BRT, RT, AmbientType} <: AbsSpec{BRT, RT}
  X::AmbientType
  U::Spec{BRT, RT}
  f::RingElem
  inc::OpenInclusion

  function PrincipalOpenSubset(X::AbsSpec, f::RingElem)
    parent(f) == OO(X) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, [f])
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, f)
  end

  function PrincipalOpenSubset(X::AbsSpec, f::Vector{<:RingElem})
    all(x->(parent(x) == OO(X)), f) || return PrincipalOpenSubset(X, OO(X).(f))
    U = hypersurface_complement(X, f)
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, (length(f)>0 ? prod(f) : one(OO(X))))
  end
end

