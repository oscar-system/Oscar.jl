@attributes mutable struct AffineAlgebraicSet{BaseRing<:Field, RingType<:MPolyAnyRing} <: AbsAffineAlgebraicSet{BaseRing, RingType}
  X::Spec
  function AffineAlgebraicSet(X::Spec; check::Bool=true)
    if check
      is_geometrically_reduced(X) || error("algebraic sets must be geometrically reduced")
    else
      # unlock the scheme methods for geometrically reduced schemes.
      set_attribute!(X, :is_geometrically_reduced, true)
      set_attribute!(X, :is_reduced, true)
    end
    new{typeof(base_ring(X)), typeof(OO(X))}(X)
  end
end


