struct ChainFactoryFromComplex{ChainType} <: HyperComplexChainFactory{ChainType}
  C::ComplexOfMorphisms{ChainType}
  auto_extend::Bool

  function ChainFactoryFromComplex(
      C::ComplexOfMorphisms{ChainType}; 
      auto_extend::Bool=false
    ) where ChainType
    return new{ChainType}(C, auto_extend)
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
  fac.auto_extend && return true
  return first(i) in range(fac.C)
end

struct MapFactoryFromComplex{MorphismType} <: HyperComplexMapFactory{MorphismType}
  C::ComplexOfMorphisms
  auto_extend::Bool

  function MapFactoryFromComplex(
      C::ComplexOfMorphisms{ChainType};
      auto_extend::Bool=false
    ) where {ChainType}
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

function can_compute(fac::MapFactoryFromComplex, HC::AbsHyperComplex, p::Int, i::Tuple)
  @assert isone(p) "request out of bounds"
  @assert isone(length(i)) "index out of bounds"
  fac.auto_extend && return true
  return first(i) in map_range(fac.C)
end

# Wrap a conventional complex into a hypercomplex
function hyper_complex(C::ComplexOfMorphisms; auto_extend::Bool=false)
  chain_factory = ChainFactoryFromComplex(C; auto_extend)
  map_factory = MapFactoryFromComplex(C; auto_extend)
  upper_bound = (typ(C) == :cochain ? last(range(C)) : first(range(C)))
  lower_bound = (typ(C) == :chain ? last(range(C)) : first(range(C)))
  result = HyperComplex(1, chain_factory, map_factory, [typ(C)],
                        upper_bounds = [upper_bound], lower_bounds = [lower_bound]
                       )
  return result
end

# Conversion of a conventional complex to an AbsSimpleComplex compatible 
# with the hypercomplex framework
function SimpleComplexWrapper(C::ComplexOfMorphisms; auto_extend::Bool=false)
  hc = hyper_complex(C; auto_extend)
  return SimpleComplexWrapper(hc)
end

### Shifts
shift(c::AbsHyperComplex, i::Tuple) = ShiftedHyperComplex(c, i)
shift(c::AbsHyperComplex, i::Int...) = shift(c, i)

shift(d::AbsDoubleComplexOfMorphisms, i::Tuple) = DoubleComplexWrapper(ShiftedHyperComplex(d, i))

shift(d::AbsSimpleComplex, i::Tuple) = SimpleComplexWrapper(ShiftedHyperComplex(d, i))


########################################################################
# Constructors for views (slices) of hypercomplexes                    #
########################################################################

# Suppose C is a double complex. Then C[5, 2:4] should return a 
# simple complex concentrated in degrees 2:4 given by the fifth
# column of C in this range.
function getindex(c::AbsHyperComplex, u::Union{Int, UnitRange}...)
  return c[Vector{Union{Int, UnitRange{Int}}}(collect(u))]
end

function getindex(c::AbsHyperComplex, u::Vector{Union{Int, UnitRange{Int}}})
  new_dims = 0
  v = Vector{Tuple{Int, Union{Int, UnitRange{Int}}}}()
  for r in u
    if r isa UnitRange
      new_dims = new_dims + 1
      push!(v, (new_dims, r))
    else
      push!(v, (0, r))
    end
  end
  return HyperComplexView(c, v)
end

function extend_dimension(c::AbsHyperComplex, d::Int)
  d >= dim(c) || error("new dimension must be at least the given one")
  all(k->has_upper_bound(c, k), 1:dim(c)) || error("extension only implemented for bounded complexes")
  all(k->has_lower_bound(c, k), 1:dim(c)) || error("extension only implemented for bounded complexes")
  u = [(k, lower_bound(c, k):upper_bound(c, k)) for k in 1:dim(c)]
  u = vcat(u, [(i, 0:0) for i in dim(c)+1:d])
  return HyperComplexView(c, u)
end

