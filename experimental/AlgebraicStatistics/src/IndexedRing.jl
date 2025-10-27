# this should find a home outside of AlgStats at some point
struct IndexedRing{T, U} <: NCRing
  R::Ring
  index_to_gen::Dict{T, U}
  gen_to_index::Dict{U, T}

  function IndexedRing(S::Ring, varnames::Vector{<:VarName}; cached=cached)
    R, x =  polynomial_ring(S, varnames; cached=cached)
    index_to_gen = Dict{Int, elem_type(R)}(i => x[i] for i in 1:ngens(R))
    gen_to_index = Dict{elem_type(R), Int}(v => k for (k, v) in index_to_gen)
    return new{Int, elem_type(R)}(R, index_to_gen, gen_to_index)
  end

  function IndexedRing(S::Ring, varnames::Pair{String, <: Array{T}}; cached=cached) where T
    gen_names = ["$(varnames.first)[$(join(x, ","))]" for x in varnames.second]
    R, r =  polynomial_ring(S, gen_names; cached=cached)
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

function indexed_ring(R::Ring, varnames; cached=false)
  MR = IndexedRing(R, varnames; cached=cached)
  return MR, gens(MR)
end

hom(R::IndexedRing, S::NCRing, images::Vector; check::Bool = true) =  hom(_ring(R), S, images)
hom(R::IndexedRing, S::IndexedRing, images::Vector; check::Bool = true) = hom(R, _ring(S), images)

(MR::IndexedRing)(r) = _ring(MR)(r)

