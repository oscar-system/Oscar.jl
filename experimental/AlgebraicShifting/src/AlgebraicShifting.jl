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

# uses permutation lemma to find other dependent columns
function permuted_dependent_cols(omega::GSet,
                                 cols::Vector{Vector{Int}},
                                 lex_min_indep_sets::Vector{Vector{Int}},
                                 max_col::Vector{Int})
  dependent_cols = Set{Vector{Int}}()
  I = first(cols)
  S, _ = stabilizer(omega, I)
  
  for (j, J) in enumerate(cols[2:end])
    exists, g = is_conjugate_with_data(omega, I, J)
    for p in right_coset(S, g)
      if exists && all(<(sort(p.(J))), sort.(on_sets(lex_min_indep_sets, p)))
        p.(J) <= max_col && push!(dependent_cols, p.(J))
      end
    end
  end
  println(dependent_cols)
  return dependent_cols
end

function lex_min_col_basis(m::AbstractAlgebra.Generic.MatSpaceElem{T},
                           n::Int, k::Int;
                           n_dependent_columns::Int=-1) where T <: MPolyRingElem
  v = identity_matrix(base_ring(m), size(m, 1))
  r = 0
  I = Int[]
  n_current_dep_cols = 0
  col_sets = sort(subsets(n, k))[1:ncols(m)]
  G = symmetric_group(n)
  omega = gset(G, on_sets, col_sets)
  max_col = col_sets[end]
  println(max_col)
  for j = 1:size(m, 2)
    # Evaluate j-th column of m * v
    c = v[r + 1:end, :] * m[:, j:j]
    m[:, j:j] = zero(m[:, j:j])
    if iszero(c)
      n_current_dep_cols += 1
      if n_dependent_columns == n_current_dep_cols
        append!(I, j + 1: ncols(m))
        break
      end
      p_dependent_cols = permuted_dependent_cols(omega, col_sets[j:end], col_sets[I], max_col)
      p_dep_ind = [findfirst(==(p_J), col_sets) for p_J in p_dependent_cols]
      println(col_sets)
      println(p_dep_ind, j)
      # m[:, p_dependent_cols] = zero(m[:, p_dependent_cols])
      continue
    end
    m[r + 1, j] = one(base_ring(m))
    push!(I, j)
    # Break if this is the last necessary column
    if r + 1 == size(m, 1)
      r += 1
      break
    end
    # Find the shortest non-zero polynomial in c
    _, i = efindmin(length, c[:, 1]; filter=!iszero)
    # Use that as pivot; move corresponding row into row r+1
    if i > 1
      swap_rows!(v, r + i, r + 1)
      swap_rows!(c,     i,     1)
    end
    # Eliminate other entries of c
    for i in 2:size(m, 1) - r
      if iszero(c[i])
        continue
      end
      _, a, b = gcd_with_cofactors(c[i], c[1])
      v[r + i, :] = b * v[r + i:r + i,:] - a * v[r + 1:r + 1, :]
      v[r + i, :] = divexact(v[r + i:r + i,:], content(v[r + i:r + i, :]))
    end
    r += 1
  end
  println(I)
  return I
end

