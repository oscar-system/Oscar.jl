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

function _domination(s1::Vector{Int}, s2::Vector{Int})
  v1 = sort(s1)
  v2 = sort(s2)
  return all([v1[i] <= v2[i] for i in 1:length(v1)])
end

function larger_cols(c::Vector{Int}, col_sets::Vector{Vector{Int}})
  col_sets[findall(x -> _domination(c, x), col_sets)]
end

# uses permutation lemma to find other dependent columns
function permuted_dependent_cols(omega::GSet,
                                 src_col::Vector{Int},
                                 target_cols::Vector{Vector{Int}},
                                 lex_min_indep_sets::Vector{Vector{Int}})

  dependent_cols = Set{Vector{Int}}()
  S, _ = stabilizer(omega, src_col)
  for J in target_cols
    exists, g = is_conjugate_with_data(omega, src_col, J)
    !exists && continue
    for p in right_coset(S, g)
      if all(<(J), sort.(on_sets(lex_min_indep_sets, p)))
        push!(dependent_cols, J)
        break
      end
    end
  end
  return dependent_cols
end

function lex_min_col_basis_cf(m::AbstractAlgebra.Generic.MatSpaceElem{T},
                              n::Int, k::Int;
                              n_dependent_columns::Int=-1) where T <: MPolyRingElem
  for i=1:nrows(m)
    c = content(m[i:i, :])
    if !isone(c)
      m[i, :] = divexact(m[i:i, :], c)
    end
  end

  I = Int[]
  dep_col_ind = Set{Int}([])
  n_current_dep_cols = 0
  col_sets = sort(subsets(n, k))[1:ncols(m)]
  G = symmetric_group(n)
  omega = gset(G, on_sets, col_sets)
  max_col = col_sets[end]
  j = 1
  for i=1:nrows(m)
    best_j = 0
    best_t = typemax(Int)
    while j <= ncols(m) 
      best_i = 0
      best_t = 0
      for ii = i:nrows(m)
        if is_zero_entry(m, ii, j)
          continue
        end
        if best_i == 0
          best_i = ii
          best_t = length(m[ii,j])
        elseif best_t > length(m[ii, j])
          best_t = length(m[ii, j])
          best_i = ii
        end
      end
      if best_i == 0
        # dependent col
        push!(dep_col_ind, j)
        n_current_dep_cols = length(dep_col_ind)
        if n_dependent_columns == n_current_dep_cols
          append!(I, [i for i in j + 1:ncols(m) if !(i in dep_col_ind)])
          return I
        end

        possible_col_ind = [i for i in j + 1:ncols(m) if !(i in dep_col_ind)]
        l_cols = larger_cols(col_sets[j], col_sets[possible_col_ind])
        
        if !isempty(l_cols)
          l_dep_ind = [findfirst(==(J), col_sets) for J in l_cols]
          dep_col_ind = union(dep_col_ind, l_dep_ind)
          m[:, l_dep_ind] = zero(m[:, l_dep_ind])
          
          next_j_index = findfirst(!in(dep_col_ind), j + 1:ncols(m))
          if isnothing(next_j_index)
            append!(I, [i for i in j + 1:ncols(m) if !(i in dep_col_ind)])
            return I
          end
          j = collect(j + 1:ncols(m))[next_j_index]
          continue
        end

        # # use permutation to find other dependent cols
        target_cols = col_sets[[i for i in j:ncols(m) if !(i in dep_col_ind)]]
        src_col = col_sets[j]
        p_dependent_cols = permuted_dependent_cols(omega, src_col, target_cols, col_sets[I])
        isempty(p_dependent_cols) && continue

        # # get the columnn indices
        p_dep_ind = [findfirst(==(J), col_sets) for J in p_dependent_cols]

        # # update set of dependent columns and the number of them
        dep_col_ind = union(dep_col_ind, p_dep_ind)
        m[:, p_dep_ind] = zero(m[:, p_dep_ind])

        next_j_index = findfirst(!in(dep_col_ind), j + 1:ncols(m))
        if isnothing(next_j_index)
          append!(I, [i for i in j + 1:ncols(m) if !(i in dep_col_ind)])
          return I
        end
        j = collect(j + 1:ncols(m))[next_j_index]
        continue
      end
      if best_i > i
        m = swap_rows!(m, i, best_i)
      end
      break
    end
    push!(I, j)

    for k=i+1:nrows(m)
      if iszero(m[k, j])
        continue
      end
      g, a, b = gcd_with_cofactors(m[k, j], m[i, j])
      m[k, :] = b*m[k:k, :] - a * m[i:i, :]
      m[k, :] = divexact(m[k:k, :], content(m[k:k, :]))
    end
    j += 1
  end

  append!(I, [i for i in j + 1:ncols(m) if !(i in dep_col_ind)])
  return I
end
