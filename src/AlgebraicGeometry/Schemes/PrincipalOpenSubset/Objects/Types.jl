########################################################################
# Principal open subsets of affine schemes                             #
########################################################################
@attributes mutable struct PrincipalOpenSubset{BRT, RT<:Ring, AmbientType} <: AbsAffineScheme{BRT, RT}
  X::AmbientType
  U::AffineScheme{BRT, RT}
  f::RingElem
  inc::AbsAffineSchemeMor # A PrincipalOpenInclusion, really; but this can not be
                  # said here because of the inclusion order.

  function PrincipalOpenSubset(X::AbsAffineScheme, f::RingElem)
    parent(f) == OO(X) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, [f])
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, f)
  end

  function PrincipalOpenSubset(X::AbsAffineScheme, f::Vector{<:RingElem})
    all(x->(parent(x) == OO(X)), f) || return PrincipalOpenSubset(X, OO(X).(f))
    U = hypersurface_complement(X, f)
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, (length(f)>0 ? prod(f) : one(OO(X))))
  end

  # The following constructors allow to pass a ring as an extra argument.
  # It is assumed that the complement of f in X is isomorphic to Spec(R).
  function PrincipalOpenSubset(X::AbsAffineScheme, R::Ring, f::RingElem;
      check::Bool=true
    )
    U = spec(R)
    @check U == hypersurface_complement(X, f) "scheme is not isomorphic to the anticipated open subset"
    return new{base_ring_type(X), ring_type(U), typeof(X)}(X, U, f)
  end
  function PrincipalOpenSubset(X::AbsAffineScheme, R::Ring, f::Vector{T};
      check::Bool=true
    ) where {T<:RingElem}
    d = prod(length(f) > 0 ? f : one(OO(X)))
    return PrincipalOpenSubset(X, R, d, check=check)
  end
end

