### Lifting maps through projective resolutions
is_chain_complex(c::AbsHyperComplex) = (dim(c) == 1 ? direction(c, 1) == (:chain) : error("complex is not one-dimensional"))
is_cochain_complex(c::AbsHyperComplex) = (dim(c) == 1 ? direction(c, 1) == (:cochain) : error("complex is not one-dimensional"))


function (l::MapLifter{T})(Phi::AbsHyperComplexMorphism, I::Tuple) where {T<:ModuleFPHom}
  i = first(I)
  i == l.start_index && return l.orig_map

  inc = 1 
  if is_chain_complex(domain(Phi))
    inc = -1
  end
  psi = Phi[Tuple(i+inc)]
  dom_map = map(domain(Phi), 1, I)
  cod_map = map(codomain(Phi), 1, Tuple(i + l.offset))
  @assert codomain(dom_map) === domain(psi)
  @assert codomain(cod_map) === codomain(psi)

  x = gens(domain(dom_map))
  y = dom_map.(x)
  z = psi.(y)
  w = [preimage(cod_map, a) for a in z]
  return hom(domain(dom_map), domain(cod_map), w, check=l.check)
end

function can_compute(l::MapLifter, Phi::AbsHyperComplexMorphism, I::Tuple)
  i = first(I)
  direction(domain(Phi), 1) == :chain  && return i >= l.start_index
  direction(domain(Phi), 1) == :cochain && return i <= l.start_index
end


function lift_map(dom::AbsHyperComplex{DomChainType}, cod::AbsHyperComplex{CodChainType}, phi::Map;
    start_index::Int=0, offset::Int=0, check::Bool=true
  ) where {DomChainType, CodChainType}
  l = MapLifter(morphism_type(DomChainType, CodChainType), dom, cod, phi; start_index, offset, check)
  return HyperComplexMorphism(dom, cod, l, offset=[offset])
end
  
morphism_type(::Type{T1}, ::Type{T2}) where {T1<:ModuleFP, T2<:ModuleFP} = ModuleFPHom{<:T1, <:T2}

### Cartan Eilenberg resolutions

