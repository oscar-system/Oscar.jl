@attributes mutable struct ProjectiveAlgebraicSet{BaseRing<:Field, GradedRingType<:Ring} <: AbsProjectiveAlgebraicSet{BaseRing, GradedRingType}
  X::ProjectiveScheme{BaseRing, GradedRingType}
  Xred::ProjectiveScheme{BaseRing, GradedRingType}
  function ProjectiveAlgebraicSet(X::ProjectiveScheme{S, T}; is_reduced::Bool=false,
                                  check::Bool=true) where {S, T}
    Y = new{S, T}()
    Y.X = X
    if is_reduced
      @check is_geometrically_reduced(X) || error("algebraic sets must be geometrically reduced")
      Y.Xred = X
    else
      # unlock the scheme methods for geometrically reduced schemes.
      set_attribute!(Y, :is_geometrically_reduced, true)
      set_attribute!(Y, :is_reduced, true)
    end
    return Y
  end
end


