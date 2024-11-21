function components(::Type{T}, D::AbsAlgebraicCycle) where {T <: AbsAlgebraicCycle}
  X = scheme(D)
  R = coefficient_ring(D)
  return [AlgebraicCycle(X, R, IdDict{AbsIdealSheaf, elem_type(R)}([I=>one(R)]); check=false) for I in components(D)]
end

function components(::Type{T}, D::AbsWeilDivisor) where {T <: AbsWeilDivisor}
  X = scheme(D)
  R = coefficient_ring(D)
  return [WeilDivisor(X, R, IdDict{AbsIdealSheaf, elem_type(R)}([I=>one(R)]); check=false) for I in components(D)]
end

function getindex(D::AbsAlgebraicCycle, C::AbsAlgebraicCycle)
  comps = components(C)
  @assert isone(length(comps))
  return D[first(comps)]
end

function canonical_divisor(X:: AbsAffineScheme{<:Field, <:MPolyRing}) 
  return CartierDivisor(covered_scheme(X), ZZ)
end

function self_intersection_via_adjunction(K_X::CartierDivisor, C::EffectiveCartierDivisor, g::Int)
  return 2g-2-intersect(K_X, C)
end
function self_intersection_via_adjunction(K_X::CartierDivisor, C::CartierDivisor, g::Int)
    return 2g-2-sum(a*intersect(K_X, comp) for (a, comp) in coefficient_dict(C); init=0)
end




