struct TensorProductChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  list::Vector{<:ComplexOfMorphisms}

  function TensorProductChainFactory(list::Vector{T}) where {ChainType, T<:ComplexOfMorphisms{ChainType}}
    return new{ChainType}(list)
  end
end

function (fac::TensorProductChainFactory{T})(HC::AbsHyperComplex, i::Tuple) where {T<:ModuleFP}
  d = length(i)
  if all(k->i[k] in range(fac.list[k]), 1:d)
    return tensor_product([fac.list[k][i[k]] for k in 1:d])
  end
  C = first(list)
  k = first(range(C))
  R = base_ring(C[k])
  return FreeMod(R, 0)
end

function can_compute(fac::TensorProductChainFactory, HC::AbsHyperComplex, i::Tuple)
  return true
end

struct TensorProductMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  list::Vector{<:ComplexOfMorphisms}

  function TensorProductMapFactory(list::Vector{T}) where {ChainType, T<:ComplexOfMorphisms{ChainType}}
    MorphismType = morphism_type(ChainType)
    return new{MorphismType}(list)
  end
end

function (fac::TensorProductMapFactory{T})(HC::AbsHyperComplex, p::Int, i::Tuple) where {T<:ModuleFPHom}
  d = length(i)
  v = collect(i)
  v[p] = (direction(HC, p) == :chain ? v[p] - 1 : v[p] + 1)
  j = Tuple(v)
  if all(k->i[k] in range(fac.list[k]), 1:d) && all(k->j[k] in range(fac.list[k]), 1:d)
    map_vec = T[identity_map(fac.list[k][i[k]]) for k in 1:p-1]
    push!(map_vec, map(fac.list[p], i[p]))
    map_vec = vcat(map_vec, T[identity_map(fac.list[k][i[k]]) for k in p+1:d])
    return hom_tensor(HC[i], HC[j], map_vec)
  end
  dom = HC[i]
  cod = HC[j]
  return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
end

function can_compute(fac::TensorProductMapFactory, HC::AbsHyperComplex, i::Tuple)
  return true
end

function tensor_product(list::Vector{T}) where {T<:ComplexOfMorphisms}
  d = length(list)
  upper_bounds = [(typ(C) == :chain ? last(range(C)) : first(range(C))) for C in list]
  lower_bounds = [(typ(C) == :cochain ? last(range(C)) : first(range(C))) for C in list]
  typs = [typ(C) for C in list]

  chain_factory = TensorProductChainFactory(list)
  map_factory = TensorProductMapFactory(list)
  return HyperComplex(d, chain_factory, map_factory, typs, upper_bounds=upper_bounds, 
                      lower_bounds=lower_bounds
                     )
end

