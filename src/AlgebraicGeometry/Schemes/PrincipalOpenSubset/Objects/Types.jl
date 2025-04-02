########################################################################
# Principal open subsets of affine schemes                             #
########################################################################
@attributes mutable struct PrincipalOpenSubset{BRT, RT<:Ring} <: AbsAffineScheme{BRT, RT}
  X::AbsAffineScheme
  U::AffineScheme{BRT, RT}
  f::Vector{<:MPolyRingElem} # factors of a representative of the complement equation
  h::RingElem # The equation in the ambient scheme defining the complement of this
  inc::AbsAffineSchemeMor # A PrincipalOpenInclusion, really; but this can not be
                  # said here because of the inclusion order.

  function PrincipalOpenSubset(X::AbsAffineScheme, f::RingElem)
    parent(f) == OO(X) || error("element does not belong to the correct ring")
    U = hypersurface_complement(X, [f])
    return new{base_ring_type(X), ring_type(U)}(X, U, [lifted_numerator(f)], f)
  end

  function PrincipalOpenSubset(X::AbsAffineScheme, f::Vector{<:RingElem})
    all(x->(parent(x) == OO(X)), f) || return PrincipalOpenSubset(X, OO(X).(f))
    U = hypersurface_complement(X, f)
    return new{base_ring_type(X), ring_type(U)}(X, U, lifted_numerator.(f), (length(f)>0 ? OO(X)(prod(lifted_numerator.(f))) : one(OO(X))))
  end

  # The following constructors allow to pass a ring as an extra argument.
  # It is assumed that the complement of f in X is isomorphic to Spec(R).
  function PrincipalOpenSubset(X::AbsAffineScheme, R::Ring, f::RingElem;
      check::Bool=true
    )
    U = spec(R)
    @check U == hypersurface_complement(X, f) "scheme is not isomorphic to the anticipated open subset"
    return new{base_ring_type(X), ring_type(U)}(X, U, lifted_numerator.(f))
  end
  function PrincipalOpenSubset(X::AbsAffineScheme, R::Ring, f::Vector{T};
      check::Bool=true
    ) where {T<:RingElem}
    U = spec(R)
    @check begin
        d = prod(length(f) > 0 ? f : one(OO(X)))
        U == hypersurface_complement(X, d)
    end
    return new{base_ring_type(X), ring_type(U)}(X, U, lifted_numerator.(f))
  end
  function PrincipalOpenSubset(X::AbsAffineScheme, U::AbsAffineScheme, f::Vector{T};
      check::Bool=true
    ) where {T<:RingElem}
    @check begin
        d = prod(length(f) > 0 ? f : one(OO(X)))
        U == hypersurface_complement(X, d)
    end
    return new{base_ring_type(X), ring_type(U)}(X, U, lifted_numerator.(f))
  end
end

