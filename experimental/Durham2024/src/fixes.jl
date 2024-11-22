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




scheme(f::VarietyFunctionFieldElem) = scheme(parent(f))
variety(f::VarietyFunctionFieldElem) = variety(parent(f))

function _colength_in_localization(I::Ideal, P::Ideal)
  return _colength_in_localization(saturated_ideal(I), saturated_ideal(P))
end

function _colength_in_localization(I::MPolyIdeal, P::MPolyIdeal)
  R = base_ring(I)
  @assert R === base_ring(P)
  U = MPolyComplementOfPrimeIdeal(P)
  L, loc = localization(R, U)
  I_loc = loc(I)
  @assert !isone(I_loc)
  F = free_module(L, 1)
  IF, inc_IF = I_loc*F
  M = cokernel(inc_IF)
  return length(M)
end


function weil_divisor(
    f::VarietyFunctionFieldElem; 
    ring::Ring=ZZ, covering::Covering=default_covering(scheme(f))
  )
  X = scheme(f)
  ideal_dict = IdDict{AbsIdealSheaf, elem_type(ring)}()
  for U in patches(covering)
    inc_dict = IdDict{AbsIdealSheaf, elem_type(ring)}()
    f_loc = f[U]
    num = numerator(f_loc)
    den = denominator(f_loc)
    num_ideal = ideal(OO(U), num)
    den_ideal = ideal(OO(U), den)
    num_dec = primary_decomposition(num_ideal)
    den_dec = primary_decomposition(den_ideal)
    for (_, P) in num_dec
      # If this component was already seen in another patch, skip it.
      new_comp = PrimeIdealSheafFromChart(X, U, P)
      any(new_comp == PP for PP in keys(ideal_dict)) && continue 
      c = _colength_in_localization(num_ideal, P)
      inc_dict[new_comp] = c
    end
    for (_, P) in den_dec
      # If this component was already seen in another patch, skip it.
      new_comp = PrimeIdealSheafFromChart(X, U, P)
      any(new_comp == PP for PP in keys(ideal_dict)) && continue 
      c = _colength_in_localization(den_ideal, P)
      key_list = collect(keys(inc_dict))
      k = findfirst(==(new_comp), key_list)
      if k === nothing
        is_zero(c) && continue
        inc_dict[new_comp] = -c
      else
        d = inc_dict[key_list[k]]
        if c == d
          delete!(inc_dict, key_list[k])
          continue
        end
        inc_dict[key_list[k]] = d - c
      end
    end
    for (pp, c) in inc_dict
      ideal_dict[pp] = c
    end
  end
  return WeilDivisor(X, ring, ideal_dict; check=false)
end

function move_divisor(D::AbsWeilDivisor; check::Bool=false)
  X = scheme(D)
  is_zero(D) && return D

  if !is_prime(D)
    R = coefficient_ring(D)
    return sum(a*move_divisor(WeilDivisor(D, R; check=false)) for (D, a) in coefficient_dict(irreducible_decomposition(D)); init=WeilDivisor(X, R))
  end

  P = first(components(D))
  i = findfirst(U->!isone(P(U)), affine_charts(X))
  i === nothing && error("divisor is trivial")
  U = affine_charts(X)[i]
  I = P(U)
  L, loc = localization(OO(U), complement_of_prime_ideal(I))
  LP = loc(I)
  g = gens(saturated_ideal(I))
  g = sort!(g, by=total_degree)
  i = findfirst(f->(ideal(L, f) == LP), g)
  x = g[i]
  f = function_field(X; check)(x)
  return irreducible_decomposition(D - weil_divisor(f))
end

function is_zero(D::AbsAlgebraicCycle)
  return all(is_zero(c) || is_one(I) for (I, c) in coefficient_dict(D))
end

function intersect(C::CartierDivisor, D::AbsWeilDivisor)
  X = scheme(C)
  @assert X === scheme(D)
  R = coefficient_ring(C)
  @assert R === coefficient_ring(D)
  result = AlgebraicCycle(X, R)
  for (E, c) in coefficient_dict(C)
    result = result + c*intersect(E, D)
  end
  return result
end

function intersect(E::EffectiveCartierDivisor, D::AbsWeilDivisor; check::Bool=true)
  X = scheme(E)
  @assert X === scheme(D)
  R = coefficient_ring(D)
  result = AlgebraicCycle(X, R)
  DD = irreducible_decomposition(D)
  for (P, a) in coefficient_dict(DD)
    _, inc_P = sub(P)
    if is_zero(pullback(inc_P, ideal_sheaf(E)))
      P_moved = move_divisor(WeilDivisor(X, R, 
        IdDict{AbsIdealSheaf, elem_type(R)}(
          [P=>one(R)]); check=false); check=false)
      result = result + a*intersect(E, P_moved)
    else
      result = result + a*AlgebraicCycle(X, R, 
        IdDict{AbsIdealSheaf, elem_type(R)}(
          [P + ideal_sheaf(E)=>one(R)]); check=false)
    end
  end
  return result
end

function is_zero(II::AbsIdealSheaf)
  return all(iszero(II(U)) for U in affine_charts(scheme(II)))
end

function +(C::CartierDivisor, E::EffectiveCartierDivisor)
  return C + 1*E
end

function -(C::CartierDivisor, E::EffectiveCartierDivisor)
  return C - 1*E
end
