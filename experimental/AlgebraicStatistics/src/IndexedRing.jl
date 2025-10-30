# this should find a home outside of AlgStats at some point
struct IndexedRing{T, U} <: NCRing
  R::Ring
  index_to_gen::Dict{T, U}
  gen_to_index::Dict{U, T}

  function IndexedRing(S::Ring, varnames::Vector{<:VarName}; kw...)
    R, x =  polynomial_ring(S, varnames; kw...)
    index_to_gen = Dict{Int, elem_type(R)}(i => x[i] for i in 1:ngens(R))
    gen_to_index = Dict{elem_type(R), Int}(v => k for (k, v) in index_to_gen)
    return new{Int, elem_type(R)}(R, index_to_gen, gen_to_index)
  end

  function IndexedRing(S::Ring, varnames::Pair{String, <: Array{T}}; kw...) where T
    gen_names = ["$(varnames.first)[$(join(x, ","))]" for x in varnames.second]
    R, r =  polynomial_ring(S, gen_names; kw...)
    U = elem_type(R)
    index_to_gen = Dict{T, U}(t => r[i] for (i, t) in enumerate(varnames.second))
    gen_to_index = Dict{U, T}(v => k for (k, v) in index_to_gen)
    return new{T, U}(R, index_to_gen, gen_to_index)
  end
end

_ring(R::IndexedRing) = R.R
symbols(R::IndexedRing) = symbols(_ring(R))
base_ring(R::IndexedRing) = base_ring(_ring(R))
elem_type(R::IndexedRing) = elem_type(_ring(R))
gens(R::IndexedRing) = R.index_to_gen
ngens(R::IndexedRing) = ngens(_ring(R))

function Base.getindex(R::IndexedRing, r::RingElem)
  @req r in keys(R.gen_to_index) "$r is not a generator of $R"
  return R.gen_to_index[r]
end

@doc raw"""
     indexed_ring(R::Ring, varnames::Vector{<:VarName}; kw...)
     indexed_ring(R::Ring, varnames::Pair{String, <:Array{T}}; kw...) where T

Returns a polynomial ring with coefficient ring `R` and the generators as a `Dict`.
The keyword arguments `kw` will be passed to the underlying [`polynomial_ring`](@ref) function.
Additional to the usual polynomial ring functionality one can also ask for the index of a given generator.

#Examples
```jldoctest
julia> R, x = indexed_ring(QQ, "x" => collect(Iterators.product(1:5, 1:3)));

julia> x[1, 2]
x[1,2]

julia> x[(1, 2)]
x[1,2]

julia> R[x[1, 2]]
(1, 2)

julia> R, (y, z) = indexed_ring(QQ, [:y, :z]);

julia> R, d = indexed_ring(QQ, [:y, :z]);

julia> d[1]
y

julia> R[d[1]]
1
```
"""
function indexed_ring(R::Ring, varnames; kw...)
  MR = IndexedRing(R, varnames; kw...)
  return MR, gens(MR)
end

hom(R::IndexedRing, S::NCRing, images::Vector; check::Bool = true) =  hom(_ring(R), S, images)
hom(R::IndexedRing, S::IndexedRing, images::Vector; check::Bool = true) = hom(R, _ring(S), images)

(MR::IndexedRing)(r) = _ring(MR)(r)

