export length
export decomposition_series

function length(M::ModuleFP{RingElemType}) where {RingElemType<:AbsLocalizedRingElem{<:Any, <:Any, <:MPolyComplementOfPrimeIdeal}}
  if iszero(M) 
    set_attribute!(M, :decomposition_series, Vector{elem_type(M)}())
    return 0
  end
  R = base_ring(M)
  P = prime_ideal(inverted_set(R))
  F = FreeMod(R, 1)
  N, _ = quo(F, (R(P)*F)[1])
  k = 0
  MM = M
  p = identity_map(M)
  NtoM, interp = hom(N, MM)
  decomp = Vector{elem_type(M)}()
  while !iszero(NtoM)
    g = gens(NtoM)
    i = findfirst(x->!(iszero(x)), g)
    Msub, inc = sub(MM, [interp(g[i])(N[1])])
    push!(decomp, preimage(p, inc(Msub[1])))
    MM, pp = quo(MM, Msub)
    p = compose(p, pp)
    NtoM, interp = hom(N, MM)
    k = k + 1
  end
  iszero(MM) || error("module does not have finite length")
  set_attribute!(M, :decomposition_series, decomp)
  return k
end

function decomposition_series(
    M::ModuleFP{RingElemType}
  ) where {RingElemType<:AbsLocalizedRingElem{<:Any, <:Any, <:MPolyComplementOfPrimeIdeal}}
  if !has_attribute(M, :decomposition_series)
    length(M)
  end
  return get_attribute(M, :decomposition_series)::Vector{elem_type(M)}
end

