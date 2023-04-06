@attributes mutable struct ProjectiveVariety{BaseRing<:Field, GradedRingType} <: AbsProjectiveVariety{BaseRing, GradedRingType}
  X::ProjectiveAlgebraicSet{BaseRing,GradedRingType}

  function ProjectiveVariety(X::ProjectiveScheme{S, T}; check::Bool=true) where {S, T}
    if check
      S <:QQField || error("varieties must be geometrically integral, but we test this only over QQ at the moment; disable this check if you know the variety is geometrically integral or proceed at your own risk")
      is_geometrically_integral(X) || error("varieties must be geometrically integral")
    else
      set_attribute!(X, :is_geometrically_integral => true)
      set_attribute!(X, :is_integral => true)
      set_attribute!(X, :is_geometrically_reduced => true)
      set_attribute!(X, :is_reduced => true)
    end
    Y = ProjectiveAlgebraicSet(X, check=false)
    new{typeof(base_ring(Y)), typeof(homogeneous_coordinate_ring(Y))}(Y)
  end

  function ProjectiveVariety(X::ProjectiveAlgebraicSet{S, T}; check::Bool=true) where {S, T}
    if check
      S <:QQField || error("varieties must be geometrically integral, but we test this only over QQ at the moment; disable this check if you know the variety is geometrically integral or proceed at your own risk")
      is_geometrically_integral(X) || error("varieties must be geometrically integral")
    else
      set_attribute!(X, :is_geometrically_integral => true)
      set_attribute!(X, :is_integral => true)
      set_attribute!(X, :is_geometrically_reduced => true)
      set_attribute!(X, :is_reduced => true)
    end
    new{typeof(base_ring(X)), typeof(homogeneous_coordinate_ring(X))}(X)
  end
end
