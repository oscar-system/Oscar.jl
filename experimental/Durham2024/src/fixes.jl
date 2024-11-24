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

# Compute the colength of I_P in the localization R_P.
# This assumes that R itself is a domain and that I is of height 1. 
# Then there exists a regular 
# point on Spec(R) and hence R_P is also regular and a UFD. 
# Being of height 1, I_P must then be principal.
function _colength_in_localization(I::Ideal, P::Ideal)
  R = base_ring(I)
  @assert R === base_ring(P)
  U = MPolyComplementOfPrimeIdeal(saturated_ideal(P))
  L, loc = localization(R, U)
  I_loc = loc(I)
  @assert base_ring(I_loc) === L
  P_loc = loc(P)
  x = _find_principal_generator(P_loc)
  y = one(L)
  k = 0
  while true
    (y in I_loc) && return k
    y = y*x
    k += 1
  end
    
  return _colength_in_localization(saturated_ideal(I), saturated_ideal(P))
end

# same assumptions as above apply.
function _find_principal_generator(I::Union{<:MPolyLocalizedIdeal, <:MPolyQuoLocalizedIdeal})
  L = base_ring(I)
  g = gens(I)
  g = sort!(g, by=x->total_degree(lifted_numerator(x)))
  for x in g
    is_zero(x) && continue
    ideal(L, x) == I && return x
  end
  error("no principal generator found")
end

function weil_divisor(
    f::VarietyFunctionFieldElem; 
    ring::Ring=ZZ, covering::Covering=default_covering(scheme(f))
  )
  @vprint :Divisors 4 "calculating principal divisor for $f\n"
  X = scheme(f)
  ideal_dict = IdDict{AbsIdealSheaf, elem_type(ring)}()
  for U in patches(covering)
    @vprint :Divisors 4 "doing patch $U\n"
    inc_dict = IdDict{AbsIdealSheaf, elem_type(ring)}()
    f_loc = f[U]
    num = numerator(f_loc)
    den = denominator(f_loc)
    num_ideal = ideal(OO(U), num)
    den_ideal = ideal(OO(U), den)
    num_dec = primary_decomposition(num_ideal)
    den_dec = primary_decomposition(den_ideal)
    @vprint :Divisors 4 "  numerator:\n"
    for (_, P) in num_dec
      @vprint :Divisors 4 "    $P\n"
      # If this component was already seen in another patch, skip it.
      new_comp = PrimeIdealSheafFromChart(X, U, P)
      @vprint :Divisors 4 "    $(any(new_comp == PP for PP in keys(ideal_dict)) ? "already found" : "new component")\n"
      any(new_comp == PP for PP in keys(ideal_dict)) && continue 
      c = _colength_in_localization(num_ideal, P)
      @vprint :Divisors 4 "    multiplicity $c\n"
      inc_dict[new_comp] = c
    end
    @vprint :Divisors 4 "  denominator:\n"
    for (_, P) in den_dec
      # If this component was already seen in another patch, skip it.
      new_comp = PrimeIdealSheafFromChart(X, U, P)
      @vprint :Divisors 4 "    $(any(new_comp == PP for PP in keys(ideal_dict)) ? "already found" : "new component")\n"
      any(new_comp == PP for PP in keys(ideal_dict)) && continue 
      c = _colength_in_localization(den_ideal, P)
      @vprint :Divisors 4 "    multiplicity $c\n"
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
  cpcd = copy(coefficient_dict(D))
  DD = irreducible_decomposition(D)
  cpcd = copy(coefficient_dict(DD))
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

function -(C::CartierDivisor, E::EffectiveCartierDivisor)
  return C - 1*E
end

