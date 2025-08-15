struct ModelRing{T, U}
  R::Ring
  index_to_gen::Dict{T, U}
  gen_to_index::Dict{U, T}

  function ModelRing(S::Ring, varnames::Vector{VarName}; cached=cached)
    R, x =  polynomial_ring(S, varnames; cached=cached)
    index_to_gen = Dict{Int, MPolyRingElem}(i => x[i] for i in 1:ngens(R))
    gen_to_index = Dict{MPolyRingElem, Int}(v => k for (k, v) in index_to_gen)
    return new{Int, elem_type(R)}(R, index_to_gen, gen_to_index)
  end

  function ModelRing(S::Ring, varnames::Pair{String, <: Array{T}}; cached=cached) where T
    gen_names = ["$(varnames.first)[$(join(x, ","))]" for x in varnames.second]
    R, r =  polynomial_ring(S, gen_names; cached=cached)
    U = elem_type(R)
    index_to_gen = Dict{T, U}(t => r[i] for (i, t) in enumerate(varnames.second))
    gen_to_index = Dict{U, T}(v => k for (k, v) in index_to_gen)
    return new{T, U}(R, index_to_gen, gen_to_index)
  end
end

gens(R::ModelRing) = R.index_to_gen
ngens(R::ModelRing) = ngens(R.R)

function model_ring(R::Ring, varnames; cached=false)
  MR = ModelRing(R, varnames; cached=cached)
  return MR, gens(MR)
end
