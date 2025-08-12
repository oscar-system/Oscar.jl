function base_change(
    phi::Any, C::AbsAlgebraicCycle; 
    scheme_base_change::Map=base_change(phi, scheme(C)[2])
  )
  X = scheme(C)
  @assert codomain(scheme_base_change) === X
  Y = domain(scheme_base_change)
  kk = coefficient_ring(C)
  ideal_dict = IdDict{AbsIdealSheaf, elem_type(kk)}()
  for (I, c) in coefficient_dict(C)
    I_red = pullback(scheme_base_change, I)
    has_attribute(I, :_self_intersection) && set_attribute!(I_red, :_self_intersection=>(get_attribute(I, :_self_intersection)::Int))
    ideal_dict[I_red] = c
  end
  return AlgebraicCycle(Y, kk, ideal_dict; check=false)
end

function base_change(
    phi::Any, D::WeilDivisor;
    scheme_base_change::Map=base_change(phi, scheme(D))[2]
  )
  und = base_change(phi, underlying_cycle(D); scheme_base_change)
  return WeilDivisor(und; check=false)
end

