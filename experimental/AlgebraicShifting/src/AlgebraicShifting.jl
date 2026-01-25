include("UniformHypergraph.jl")
include("PartialShift.jl")
include("SymmetricShift.jl")
include("PartialShiftGraph.jl")
include("serialization.jl")

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


function pivot_columns(A::MatElem)
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

"""
    eargmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)

Extended argmin function. Allows custom filter, default value and comparator.
"""
function eargmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)
  best = nothing
  for x in xs
    filter(x) || continue
    if isnothing(best) || lt(f(x), best)
      best = f(x)
    end
  end
  return best
end

"""
    efindmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)

Extended findmin function. Allows custum filter, default value and comparator.
#Examples
```jldoctest
julia> Oscar.efindmin(length, ["", "a", "ab", "abc"]; filter=(!=)(""))
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

function _domination(s1::Combination{Int}, s2::Combination{Int})
  v1 = sort(s1)
  v2 = sort(s2)
  return all([v1[i] <= v2[i] for i in 1:length(v1)])
end
