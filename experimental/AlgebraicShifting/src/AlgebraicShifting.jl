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
  The collection xs need not have linear nor 1-based indexing.
"""
function eargmin(f, xs; filter=_->true, default=nothing, lt=Base.isless)
  best = nothing
  i_best = default
  for i in eachindex(xs)
    x = xs[i]
    filter(x) || continue
    if isnothing(best) || lt(f(x), best)
      best = f(x)
      i_best = i
    end
  end
  return i_best
end
