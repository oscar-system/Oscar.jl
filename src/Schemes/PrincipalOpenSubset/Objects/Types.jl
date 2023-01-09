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

  # The following constructors allow to pass a ring as an extra argument.
  # It is assumed that the complement of f in X is isomorphic to Spec(R).
  function PrincipalOpenSubset(X::AbsSpec, R::Ring, f::RingElem;
      check::Bool=true
    )
    U = Spec(R)
    check && (U == hypersurface_complement(X, f) || error("scheme is not isomorphic to the anticipated open subset"))
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, f)
  end
  function PrincipalOpenSubset(X::AbsSpec, R::Ring, f::Vector{T};
      check::Bool=true
    ) where {T<:RingElem}
    d = prod(length(f) > 0 ? f : one(OO(X)))
    return PrincipalOpenSubset(X, R, d)
  end
end

