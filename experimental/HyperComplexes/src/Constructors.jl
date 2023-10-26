struct ChainFactoryFromComplex{ChainType} <: HyperComplexChainFactory{ChainType}
  C::ComplexOfMorphisms{ChainType}

  function ChainFactoryFromComplex(C::ComplexOfMorphisms{ChainType}) where ChainType
    return new{ChainType}(C)
  end
end

function (fac::ChainFactoryFromComplex{T})(HC::AbsHyperComplex, i::Tuple) where {T<:ModuleFP}
  @assert length(i) == 1 "wrong type of index"
  k = i[1]
  if k in range(fac.C)
    return fac.C[k]
  end
  R = base_ring(fac.C[first(range(fac.C))])
  return FreeMod(R, 0)
end

function can_compute(fac::ChainFactoryFromComplex, HC::AbsHyperComplex, i::Tuple)
  @assert length(i) == 1 "wrong type of index"
  return true
end

struct MapFactoryFromComplex{MorphismType} <: HyperComplexMapFactory{MorphismType}
  C::ComplexOfMorphisms

  function MapFactoryFromComplex(C::ComplexOfMorphisms{ChainType}) where {ChainType}
    MorphismType = morphism_type(ChainType)
    return new{MorphismType}(C)
  end
end

function (fac::MapFactoryFromComplex{T})(HC::AbsHyperComplex, p::Int, i::Tuple) where {T<:ModuleFPHom}
  @assert length(i) == 1 "wrong type of index"
  @assert p == 1 "complex is one-dimensional"
  k = i[1]
  if k in map_range(fac.C)
    return map(fac.C, k)
  end
  dom = HC[(k,)]
  cod = HC[(k-1,)]
  return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
end


function hyper_complex(C::ComplexOfMorphisms)
  chain_factory = ChainFactoryFromComplex(C)
  map_factory = MapFactoryFromComplex(C)
  upper_bound = (typ(C) == :chain ? last(range(C)) : first(range(C)))
  lower_bound = (typ(C) == :cochain ? last(range(C)) : first(range(C)))
  result = HyperComplex(1, chain_factory, map_factory, [typ(C)],
                        upper_bounds = [upper_bound], lower_bounds = [lower_bound]
                       )
  return result
end


