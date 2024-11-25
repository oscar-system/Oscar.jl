function canonical_divisor(X:: AbsAffineScheme{<:Field, <:MPolyRing}) 
  return CartierDivisor(covered_scheme(X), ZZ)
end

function self_intersection_via_adjunction(K_X::CartierDivisor, C::EffectiveCartierDivisor, g::Int)
  return 2g-2-intersect(K_X, C)
end
function self_intersection_via_adjunction(K_X::CartierDivisor, C::CartierDivisor, g::Int)
    return 2g-2-sum(a*intersect(K_X, comp) for (a, comp) in coefficient_dict(C); init=0)
end

