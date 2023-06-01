@attributes mutable struct ProjectiveAlgebraicSet{BaseRing<:Field, GradedRingType<:Ring} <: AbsProjectiveAlgebraicSet{BaseRing, GradedRingType}
  X::ProjectiveScheme{BaseRing,GradedRingType}
  function ProjectiveAlgebraicSet(X::ProjectiveScheme; check::Bool=true)
    @check
    if check
      is_geometrically_reduced(X) || error("algebraic sets must be geometrically reduced")
    else
      # unlock the scheme methods for geometrically reduced schemes.
      set_attribute!(X, :is_geometrically_reduced, true)
      set_attribute!(X, :is_reduced, true)
    end
    new{typeof(base_ring(X)), typeof(homogeneous_coordinate_ring(X))}(X)
  end
end


