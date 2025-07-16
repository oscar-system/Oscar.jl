@attributes mutable struct AffineVariety{BaseRing<:Field, RingType<:MPolyAnyRing} <: AbsAffineVariety{BaseRing, RingType}
  X::AbsAffineAlgebraicSet

  function AffineVariety(X::AbsAffineAlgebraicSet{S, T}; check::Bool) where {S, T}
    if check
      S <:QQField || error("varieties must be geometrically integral, but we test this only over QQ at the moment, disable this check if you know the variety is geometrically integral or proceed at your own risk")
      is_geometrically_integral(X) || error("varieties must be geometrically integral")
    end
    Y = new{S,T}(X)
    Y.__attrs = copy(X.__attrs)
    set_attribute!(Y, :is_geometrically_integral => true)
    set_attribute!(Y, :is_integral => true)
    set_attribute!(Y, :is_geometrically_reduced => true)
    set_attribute!(Y, :is_reduced => true)
    return Y
  end

end




