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
  println("permuted")
  dependent_cols = Set{Vector{Int}}()
  I = first(cols)
  S, _ = stabilizer(omega, I)
  println("# of possible dependent cols", length(cols[2:end]))
  for (index, J) in enumerate(cols[2:end])
    println(index)
    exists, g = is_conjugate_with_data(omega, I, J)
     for p in right_coset(S, g)
       if exists && all(<(J), sort.(on_sets(lex_min_indep_sets, p)))
         J <= max_col && push!(dependent_cols, J)
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
    println(j)
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
        println("n_dependent", n_dependent_columns)
        n_current_dep_cols += 1
        if j in dep_col_ind
          n_current_dep_cols += 1
        end
        if n_dependent_columns == n_current_dep_cols
          append!(I, j + 1:ncols(m))
          return I
        end

        if j in dep_col_ind
          continue
        end
        # 
        # # use permutation to find other dependent cols
        target_cols = col_sets[[i for i in j:ncols(m) if !(i in dep_col_ind)]]
        p_dependent_cols = permuted_dependent_cols(omega, target_cols, col_sets[I], max_col)
        isempty(p_dependent_cols) && continue
        # 
        # # get the columnn indices
        p_dep_ind = [findfirst(==(J), col_sets) for J in p_dependent_cols]
        # 
        # # update set of dependent columns and the number of them
        dep_col_ind = union(dep_col_ind, p_dep_ind)
        println(dep_col_ind, " ", j, " ", n_dependent_columns)
        # n_current_dep_cols = length(dep_col_ind)
        # if n_current_dep_cols == n_dependent_columns
          # println(dep_col_ind)
          # break
        # end
        # println()

        m[:, p_dep_ind] = zero(m[:, p_dep_ind])
        j += 1
        continue
      end
      if best_i > i
        m = swap_rows!(m, i, best_i)
      end
      break
    end
    if j > ncols(m)
      return I
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
  return I
end

function lex_min_col_basis_fl(m::AbstractAlgebra.Generic.MatSpaceElem{T},
                              n::Int, k::Int;
                              n_dependent_columns::Int=-1) where T <: MPolyRingElem
  v = identity_matrix(base_ring(m), size(m, 1))
  r = 0
  I = Int[]
  dep_col_ind = Set{Int}([])
  n_current_dep_cols = 0
  col_sets = sort(subsets(n, k))[1:ncols(m)]
  G = symmetric_group(n)
  omega = gset(G, on_sets, col_sets)
  max_col = col_sets[end]
  for j = 1:size(m, 2)
    # Evaluate j-th column of m * v
    j in dep_col_ind && continue
    c = v[r + 1:end, :] * m[:, j:j]
    m[:, j:j] = zero(m[:, j:j])
    if iszero(c)
      n_current_dep_cols += 1
      push!(dep_col_ind, j)
      if n_dependent_columns == n_current_dep_cols
        append!(I, j + 1:ncols(m))
        return I
      end
      
      # use permutation to find other dependent cols
      p_dependent_cols = permuted_dependent_cols(omega, col_sets[j:end], col_sets[I], max_col)#
      isempty(p_dependent_cols) && continue
      # 
      # # get the columnn indices
      p_dep_ind = [findfirst(==(J), col_sets) for J in p_dependent_cols]
      # 
      # # update set of dependent columns and the number of them
      # dep_col_ind = union(dep_col_ind, p_dep_ind)
      # n_current_dep_cols = length(dep_col_ind)

      # if n_current_dep_cols == n_dependent_columns
      #   append!(I, j + 1:ncols(m))
      #   return collect(symdiff(Set(I), dep_col_ind))
      # end
      m[:, p_dep_ind] = zero(m[:, p_dep_ind])
      continue
    end
    m[r + 1, j] = one(base_ring(m))
    # Break if this is the last necessary column
    push!(I, j)
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
  return I
end

