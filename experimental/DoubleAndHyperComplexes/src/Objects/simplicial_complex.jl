#=
# Since some code coverage tool blocks the tests from passing if this file is 
# present with proper julia code inside, the code below is commented out. 
#
# If you want to implement your own hyper complex class, you may start with 
# the template below and replace every occurrence of `SimplicialCo` with your favourite
# name for your new class. Then you have to fill in the gaps according to your 
# needs. We provide a sample implementation below.
=#


### Production of the chains
struct SimplicialCoChainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  R::Ring
  K::SimplicialComplex

  function SimplicialCoChainFactory(R::Ring, K::SimplicialComplex)
    return new{FreeMod{elem_type(R)}}(R, K)
  end
end

function (fac::SimplicialCoChainFactory)(self::AbsHyperComplex, i::Tuple)
  return FreeMod(fac.R, length(faces(fac.K, only(i))))
end

function can_compute(fac::SimplicialCoChainFactory, self::AbsHyperComplex, i::Tuple)
  return 0 <= only(i) <= dim(fac.K)
end

### Production of the morphisms 
struct SimplicialCoMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  function SimplicialCoMapFactory()
    return new{FreeModuleHom}()
  end
end

function (::SimplicialCoMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  fac = chain_factory(self)
  i = only(I)
  dom = self[i]
  cod = self[i+1]
  return hom(dom, cod, transpose(matrix(fac.R, Polymake.topaz.boundary_matrix(fac.K.pm_simplicialcomplex, i+1))))
end

function can_compute(fac::SimplicialCoMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  fac = chain_factory(self)
  p == 1 || return false
  return 0 <= only(i) < dim(fac.K)
end

### The concrete struct
@attributes mutable struct SimplicialCoComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}
  face_to_index_map::Dict{Int, Dict{Set{Int}, Int}}

  function SimplicialCoComplex(R::Ring, K::SimplicialComplex)
    chain_fac = SimplicialCoChainFactory(R::Ring, K::SimplicialComplex)
    map_fac = SimplicialCoMapFactory()

    internal_complex = HyperComplex(1, chain_fac, map_fac, [:cochain])
    ChainType = FreeMod{elem_type(R)}
    MorphismType = FreeModuleHom
    return new{ChainType, MorphismType}(internal_complex)
  end
end

### Show the simplices belonging to the cochain a
function show_elem(io::IO, C::SimplicialCoComplex, a::FreeModElem, p::Int)
  @req parent(a) === C[p] "parent mismatch"
  fst = true
  for (c, g) in coordinates(a)
    !fst && print(io, " + ")
    fst = false
    print(io, g, "*", sort(collect(faces(simplicial_complex(C), p)[c])))
  end
end

### Generic version (TODO: move elsewhere)
function show_elem(io::IO, C::AbsHyperComplex, a::FreeModElem, p::Int)
  @req parent(a) === C[p] "parent mismatch"
  print(io, a)
end

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::SimplicialCoComplex) = c.internal_complex

function face_to_index_map(c::SimplicialCoComplex, p::Int)
  if !isdefined(c, :face_to_index_map)
    c.face_to_index_map = Dict{Int, Dict{Set{Int}, Int}}()
  end
  return get!(c.face_to_index_map, p) do
    Dict{Set{Int}, Int}(s => i for (i,s) in enumerate(faces(simplicial_complex(c), p)))
  end
end

### additional getters
simplicial_complex(C::SimplicialCoComplex) = chain_factory(C).K
base_ring(C::SimplicialCoComplex) = chain_factory(C).R

### cup product for (composable) dual simplices
function mul_cochains(C::SimplicialCoComplex, a::FreeModElem, p::Int, b::FreeModElem, q::Int)
  @req parent(a) === C[p] && parent(b) === C[q] "parent mismatch"
  K = simplicial_complex(C)
  cochain = zero(C[p+q])
  f = face_to_index_map(C, p+q)
  for (ga, ca) in coordinates(a), (gb, cb) in coordinates(b)
    sa = faces(K, p)[ga]
    sb = faces(K, q)[gb]
    s = union(sa, sb)
    (maximum(sa) == minimum(sb) && s in keys(f)) || continue
    cochain += ca * cb * gens(C[p+q])[f[s]]
  end
  return cochain
end
