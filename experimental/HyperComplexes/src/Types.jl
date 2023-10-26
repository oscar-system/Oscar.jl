### Abstract type and interface for double complexes
abstract type AbsHyperComplex{ChainType, MapType} end

### asking for the entries of the complex
getindex(HC::AbsHyperComplex, i::Tuple) = underlying_complex(HC)[i]
has_index(HC::AbsHyperComplex, i::Tuple) = has_index(underlying_complex(HC), i)
can_compute_index(HC::AbsHyperComplex, i::Tuple) = can_compute_index(underlying_complex(HC), i)

### accessing the maps
map(HC::AbsHyperComplex, p::Int, i::Tuple) = map(underlying_complex(HC), p, i)
has_map(HC::AbsHyperComplex, p::Int, i::Tuple) = has_map(underlying_complex(HC), p, i)
can_compute_map(HC::AbsHyperComplex, p::Int, i::Tuple) = can_compute_map(underlying_complex(HC), p, i)

### properties
direction(HC::AbsHyperComplex, p::Int) = direction(underlying_complex(HC), p)
dim(HC::AbsHyperComplex) = dim(underlying_complex(HC))
is_complete(HC::AbsHyperComplex) = is_complete(underlying_complex(HC))
has_upper_bound(HC::AbsHyperComplex, p::Int) = has_upper_bound(underlying_complex(HC), p)
has_lower_bound(HC::AbsHyperComplex, p::Int) = has_lower_bound(underlying_complex(HC), p)
upper_bound(HC::AbsHyperComplex, p::Int) = upper_bound(underlying_complex(HC), p)
lower_bound(HC::AbsHyperComplex, p::Int) = lower_bound(underlying_complex(HC), p)

### factories for the chains
abstract type HyperComplexChainFactory{ChainType} end

function (fac::HyperComplexChainFactory)(HC::AbsHyperComplex, i::Tuple)
  error("production of the $i-th entry not implemented for factory of type $(typeof(fac))")
end

function can_compute(fac::HyperComplexChainFactory, HC::AbsHyperComplex, i::Tuple)
  error("testing whether the $i-th entry can be computed using $fac is not implemented; please overwrite this method")
end

### factories for the maps
abstract type HyperComplexMapFactory{MorphismType} end

function (fac::HyperComplexMapFactory)(HC::AbsHyperComplex, p::Int, i::Tuple)
  error("production of the $i-th map in the $p-th direction not implemented for factory of type $(typeof(fac))")
end

function can_compute(fac::HyperComplexMapFactory, HC::AbsHyperComplex, p::Int, i::Tuple)
  error("testing whether the $i-th map in the $p-th direction can be computed using $fac is not implemented; please overwrite this method")
end

mutable struct HyperComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType}
  d::Int 
  chains::Dict{Tuple, <:ChainType}
  morphisms::Dict{Tuple, Dict{Int, <:MorphismType}}

  chain_factory::HyperComplexChainFactory{ChainType}
  map_factory::HyperComplexMapFactory{MorphismType}

  directions::Vector{Symbol}

  upper_bounds::Vector{Union{Int, Nothing}}
  lower_bounds::Vector{Union{Int, Nothing}}

  # fields for caching
  is_complete::Bool

  function HyperComplex(
      d::Int,
      chain_factory::HyperComplexChainFactory{ChainType},
      map_factory::HyperComplexMapFactory{MorphismType},
      directions::Vector{Symbol};
      upper_bounds::Vector=[nothing for i in 1:d],
      lower_bounds::Vector=[nothing for i in 1:d]
    ) where {ChainType, MorphismType}
    @assert d > 0 "can not create zero or negative dimensional hypercomplex"
    chains = Dict{Tuple, ChainType}()
    morphisms = Dict{Tuple, Dict{Int, <:MorphismType}}()
    return new{ChainType, MorphismType}(d, chains, morphisms, 
                                        chain_factory, map_factory, directions, 
                                        Vector{Union{Int, Nothing}}(upper_bounds), 
                                        Vector{Union{Int, Nothing}}(lower_bounds)
                                       )
  end
end


