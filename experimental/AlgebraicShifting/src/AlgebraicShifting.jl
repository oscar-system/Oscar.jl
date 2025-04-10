include("SimplicialComplexesExtended.jl")
include("UniformHypergraph.jl")
include("PartialShift.jl")
include("SymmetricShift.jl")
include("PartialShiftGraph.jl")

export UniformHypergraph

export compound_matrix
export contracted_partial_shift_graph
export exterior_shift
export face_size
export generic_unipotent_matrix
export partial_shift_graph
export partial_shift_graph_vertices
export rothe_matrix
export uniform_hypergraph
export eargmin
export efindmin

function independent_columns(A::MatElem)
  col_indices = Int[]
  row_index = 1
  row_bound = nrows(A)
  for col_index in 1:ncols(A)
    if !is_zero(A[row_index, col_index])
      push!(col_indices, col_index)
      row_index += 1
    end
    row_index > row_bound && break
  end
  return col_indices
end

""" eargmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)

  Extended argmin function. Allows custum filter, default value and comparator.
"""
function eargmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)
  best = nothing
  for x in xs
    filter(x) || continue
    if isnothing(best) || lt(f(x), best)
      best = f(x)
    end
  end
  return i_best
end

""" efindmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)

  Extended findmin function. Allows custum filter, default value and comparator.
  #Examples
  ```jldoctest
    julia> efindmin(length, ["", "a", "ab", "abc"]; filter=(!=)(""))
    (1, 2)
  ```
  because `1 == length("a")`, which is the shortest non-empty element of the array, and occurs at index 2.
"""
function efindmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)
  best = nothing
  i_best = default
  for (i, x) in enumerate(xs)
    filter(x) || continue
    if isnothing(best) || lt(f(x), best)
      best = f(x)
      i_best = i
    end
  end
  return (best, i_best)
end
efindmin(xs; filter=_->true, default=nothing, lt=Base.isless) = efindmin(identity, xs; filter=filter, default=default, lt=lt)

################################################################################
# linear orderings
function isless_hasone_revlex(S1::Vector{Int}, S2::Vector{Int})
  @req length(S1) == length(S2) "Vectors should be of same size for comparison"
  1 in S1 && 1 in S2 && return S1 >= S2
  1 in S1 && return true
  1 in S2 && return false
  return S1 >= S2
end

function isless_hasone_revlex(S1::Set{Int}, S2::Set{Int})
  return isless_hasone_revlex(sort(collect(S1)), sort(collect(S2)))
end

function isless_lex(S1::Set{Set{Int}}, S2::Set{Set{Int}})
  S_diff = collect(symdiff(S1, S2))
  isempty(S_diff) && return false
  set_cmp(a, b) = min(symdiff(a, b)...) in a
  return sort(S_diff;lt=set_cmp)[1] in S1
end

function isless_lex(K1::SimplicialComplex, K2::SimplicialComplex)
  return isless_lex(Set(facets(K1)), Set(facets(K2)))
end

function isless_lex(K1::UniformHypergraph, K2::UniformHypergraph)
  return faces(K1) < faces(K2)
end

"""
  Given two simplicial complexes `K1`, `K2` return true if
  `K1` is lexicographically less than `K2`
"""
isless_lex(K1::ComplexOrHypergraph, K2::ComplexOrHypergraph) = isless_lex(Set(facets(K1)), Set(facets(K2)))

function isless_induced(lt::Function, ::Type{Set{Set{Int}}})
  function induced(S1::Set{Set{Int}}, S2::Set{Set{Int}})
    S_diff = collect(symdiff(S1, S2))
    isempty(S_diff) && return false
    return sort(S_diff;lt=lt)[1] in S1
  end
  return induced
end

function isless_induced(lt::Function, ::Type{T}) where T <: SimplicialComplex
  f = isless_induced(lt, Set{Set{Int}})
  function induced(K1::T, K2::T)
    return f(Set(facets(K1)), Set(facets(K2)))
  end
  return induced
end
