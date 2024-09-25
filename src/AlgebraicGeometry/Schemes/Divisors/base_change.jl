function base_change(
    phi::Any, C::AbsAlgebraicCycle; 
    scheme_base_change::Map=base_change(phi, scheme(C)[2])
  )
  X = scheme(C)
  @assert codomain(scheme_base_change) === X
  Y = domain(scheme_base_change)
  kk = coefficient_ring(C)
  return AlgebraicCycle(Y, kk, IdDict{AbsIdealSheaf, elem_type(kk)}(pullback(scheme_base_change, I)=>c for (I, c) in coefficient_dict(C)); check=false)
end

function base_change(
    phi::Any, D::WeilDivisor;
    scheme_base_change::Map=base_change(phi, scheme(D))[2]
  )
  und = base_change(phi, underlying_cycle(D); scheme_base_change)
  return WeilDivisor(und; check=false)
end

