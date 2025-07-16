@attributes mutable struct ProjectiveVariety{BaseRing<:Field, GradedRingType} <: AbsProjectiveVariety{BaseRing, GradedRingType}
  X::ProjectiveAlgebraicSet{BaseRing,GradedRingType}

  function ProjectiveVariety(X::ProjectiveScheme{S, T}; check::Bool=true) where {S, T}
    @check begin
      S <:QQField || error("varieties must be geometrically integral, but we test this only over QQ at the moment; disable this check if you know the variety is geometrically integral or proceed at your own risk")
      is_geometrically_integral(X) || error("varieties must be geometrically integral")
    end
    Y = ProjectiveAlgebraicSet(X; check=false)
    return ProjectiveVariety(Y; check=false)
  end

  function ProjectiveVariety(X::ProjectiveAlgebraicSet{S, T}; check::Bool=true) where {S, T}
    @check begin
      S <:QQField || error("varieties must be geometrically integral, but we test this only over QQ at the moment; disable this check if you know the variety is geometrically integral or proceed at your own risk")
      is_geometrically_integral(X) || error("varieties must be geometrically integral")
    end
    Z = new{S, T}(X)
    Z.__attrs = copy(X.__attrs)
    set_attribute!(Z, :is_geometrically_integral => true)
    set_attribute!(Z, :is_integral => true)
    set_attribute!(Z, :is_geometrically_reduced => true)
    set_attribute!(Z, :is_reduced => true)
  end
end
