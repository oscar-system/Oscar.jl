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
      any(P == PP for PP in keys(ideal_dict)) && continue 
      c = _colength_in_localization(num_ideal, P)
      new_comp = PrimeIdealSheafFromChart(X, U, P)
      if !any(new_comp == P for P in keys(inc_dict))
        is_zero(c) && continue
        inc_dict[new_comp] = -c
      else
        d = inc_dict[new_comp]
        if c == d
          delete!(inc_dict, new_comp)
          continue
        end
        inc_dict[new_comp] = d - c
      end
    end
    for (pp, c) in inc_dict
      ideal_dict[pp] = c
    end
  end
  return WeilDivisor(X, ring, ideal_dict; check=false)
end

