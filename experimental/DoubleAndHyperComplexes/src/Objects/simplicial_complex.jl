#= Simplicial cochain complexes from simplicial complexes
# This is to extract from a simplicial complex (i.e. the topological space given 
# by means of combinatorial data) the corresponding cochain complex for singular 
# cohomology.
=#
### Production of the chains
struct SimplicialCochainFactory{ChainType} <: HyperComplexChainFactory{ChainType}
  R::Ring
  K::SimplicialComplex

  function SimplicialCochainFactory(R::Ring, K::SimplicialComplex)
    return new{FreeMod{elem_type(R)}}(R, K)
  end
end

### Initializes the free modules over the ring of appropriate rank in the relevant degrees. 
function (fac::SimplicialCochainFactory)(self::AbsHyperComplex, i::Tuple)
  return FreeMod(fac.R, length(faces(fac.K, only(i))))
end

### The complex is only computable in at most the dimension of the simplicial complex
function can_compute(fac::SimplicialCochainFactory, self::AbsHyperComplex, i::Tuple)
  return 0 <= only(i) <= dim(fac.K)
end

### Production of the morphisms 
struct SimplicialCochainMapFactory{MorphismType} <: HyperComplexMapFactory{MorphismType}
  function SimplicialCochainMapFactory()
    return new{FreeModuleHom}()
  end
end

### Populates morphisms with matrices defined by the boundary matrices of the simplicial complexes. 
function (::SimplicialCochainMapFactory)(self::AbsHyperComplex, p::Int, I::Tuple)
  fac = chain_factory(self)
  i = only(I)
  dom = self[i]
  cod = self[i+1]
  return hom(dom, cod, transpose(matrix(fac.R, Polymake.topaz.boundary_matrix(fac.K.pm_simplicialcomplex, i+1))))
end

### 
function can_compute(fac::SimplicialCochainMapFactory, self::AbsHyperComplex, p::Int, i::Tuple)
  fac = chain_factory(self)
  p == 1 || return false
  return 0 <= only(i) < dim(fac.K)
end

### The concrete struct
@doc raw"""
    SimplicialCochainComplex(R::Ring, K::SimplicialComplex)

Given a ring and a simplicial complex (in OSCAR), produces the cochain complex which computes simplicial cohomology of the complex. 

This relies on the abstract type for hypercomplexes. 
"""
@attributes mutable struct SimplicialCochainComplex{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType} 
  internal_complex::HyperComplex{ChainType, MorphismType}
  face_to_index_map::Dict{Int, Dict{Set{Int}, Int}}
  multiplication_dict::Dict{Tuple{Int, Int, Int, Int}, <:FreeModElem}

  function SimplicialCochainComplex(R::Ring, K::SimplicialComplex)
    chain_fac = SimplicialCochainFactory(R::Ring, K::SimplicialComplex)
    map_fac = SimplicialCochainMapFactory()

    internal_complex = HyperComplex(1, chain_fac, map_fac, [:cochain]; upper_bounds=Union{Int, Nothing}[dim(K)], lower_bounds=Union{Int, Nothing}[0])
    ChainType = FreeMod{elem_type(R)}
    MorphismType = FreeModuleHom
    return new{ChainType, MorphismType}(internal_complex)
  end
end

### Show the simplices belonging to the cochain a
function show_elem(io::IO, C::SimplicialCochainComplex, a::FreeModElem, p::Int)
  @req parent(a) === C[p] "parent mismatch"
  fst = true
  for (c, g) in coordinates(a)
    !fst && print(io, " + ")
    fst = false
    print(io, g, "*", sort(collect(faces(simplicial_complex(C), p)[c])))
  end
end

### Generic version (TODO: move elsewhere)
#=
function show_elem(io::IO, C::AbsHyperComplex, a::FreeModElem, p::Int)
  @req parent(a) === C[p] "parent mismatch"
  print(io, a)
end
=#

### Implementing the AbsHyperComplex interface via `underlying_complex`
underlying_complex(c::SimplicialCochainComplex) = c.internal_complex

### 
function face_to_index_map(c::SimplicialCochainComplex, p::Int)
  if !isdefined(c, :face_to_index_map)
    c.face_to_index_map = Dict{Int, Dict{Set{Int}, Int}}()
  end
  return get!(c.face_to_index_map, p) do
    Dict{Set{Int}, Int}(s => i for (i,s) in enumerate(faces(simplicial_complex(c), p)))
  end
end

# a dynamic table for faster multiplication
function multiplication_dict(c::SimplicialCochainComplex)
  if !isdefined(c, :multiplication_dict)
    c.multiplication_dict = Dict{Tuple{Int, Int, Int, Int}, FreeModElem{elem_type(base_ring(c))}}()
  end
  return c.multiplication_dict::Dict{Tuple{Int, Int, Int, Int}, FreeModElem{elem_type(base_ring(c))}}
end

### additional getters
@doc raw"""
    simplicial_complex(C::SimplicialCochainComplex)

Given a simplicial cochain complex, recovers the underlying simplicial complex. 
"""
simplicial_complex(C::SimplicialCochainComplex) = chain_factory(C).K

@doc raw"""
    base_ring(C::SimplicialCochainComplex)

Given a simplicial cochain complex, recovers the ring over which the cochain complex is defined. 
"""
base_ring(C::SimplicialCochainComplex) = chain_factory(C).R

### cup product for (composable) dual simplices
@doc raw"""
    mul_cochains(C::SimplicialCochainComplex, a::FreeModElem, p::Int, b::FreeModElem, q::Int)

Computes the cup product of two classes in simplicial cohomology follwing the formula in Chapter 3.2 of Hatcher's text on algebraic topology (https://pi.math.cornell.edu/~hatcher/AT/ATch3.pdf). 
"""
function mul_cochains(C::SimplicialCochainComplex, a::FreeModElem, p::Int, b::FreeModElem, q::Int)
  @req parent(a) === C[p] && parent(b) === C[q] "parent mismatch"
  K = simplicial_complex(C)
  cochain = zero(C[p+q])
  f = face_to_index_map(C, p+q)
  mdict = multiplication_dict(C)
  for (ga, ca) in coordinates(a), (gb, cb) in coordinates(b)
    res = get!(mdict, (ga, p, gb, q)) do
      sa = faces(K, p)[ga]
      sb = faces(K, q)[gb]
      s = union(sa, sb)
      maximum(sa) == minimum(sb) && s in keys(f) && return gen(C[p+q], f[s])
      zero(C[p+q])
    end
    @assert parent(res) === C[p+q]
    if !iszero(res)
      cochain += ca * cb * res
    end
  end
  return cochain
end

