@attributes mutable struct AffineVariety{BaseRing<:Field, RingType<:MPolyAnyRing} <: AbsAffineVariety{BaseRing, RingType}
  X::Spec

  function AffineVariety(X::Spec; check::Bool)
    if check
      error("varieties must be geometrically integral, but we cannot test this at the moment, disable this check if you know X is geometrically integral or proceed at your own risk")
      is_geometrically_integral(X) || error("varieties must be geometrically integral")
    else
      set_attribute!(X, :is_geometrically_integral => true)
      set_attribute!(X, :is_integral => true)
      set_attribute!(X, :is_geometrically_reduced => true)
      set_attribute!(X, :is_reduced => true)
    end
    new{typeof(base_ring(X)), typeof(OO(X))}(X)
  end

  function AffineVariety(X::Spec{<:QQField,T}; check::Bool) where T
    if check
      is_geometrically_integral(X) || error("varieties must be geometrically integral")
    else
      set_attribute!(X, :is_geometrically_integral => true)
      set_attribute!(X, :is_integral => true)
      set_attribute!(X, :is_geometrically_reduced => true)
      set_attribute!(X, :is_reduced => true)
    end
    new{typeof(base_ring(X)), typeof(OO(X))}(X)
  end
end




