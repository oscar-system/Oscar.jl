@attributes mutable struct ProjectiveAlgebraicSet{BaseRing<:Field, GradedRingType<:Ring} <: AbsProjectiveAlgebraicSet{BaseRing, GradedRingType}
  X::ProjectiveScheme{BaseRing, GradedRingType}
  Xred::ProjectiveScheme{BaseRing, GradedRingType}
  function ProjectiveAlgebraicSet(X::ProjectiveScheme{S, T}; is_reduced::Bool=false,
                                  check::Bool=true) where {S, T}
    Y = new{S, T}()
    Y.X = X
    if is_reduced
      Y.Xred = X
      set_attribute!(Y.Xred, :is_reduced, true)
      set_attribute!(Y.Xred, :is_geometrically_reduced, true)
    end
    set_attribute!(Y, :is_reduced, true)
    @check is_geometrically_reduced(Y) "Algebraic sets must be geometrically reduced"
    set_attribute!(Y, :is_geometrically_reduced, true)
    return Y
  end
end

