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
  return best
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

function _domination(s1::Combination{Int}, s2::Combination{Int})
  v1 = sort(s1)
  v2 = sort(s2)
  return all([v1[i] <= v2[i] for i in 1:length(v1)])
end

function larger_cols(c::Combination{Int}, col_sets::Vector{Combination{Int}})
  col_sets[findall(x -> _domination(c, x), col_sets)]
end

function _find_evaluations(dep_col_index::Int,
                           cols::Vector{Combination{Int}},
                           generic_rothe::AbstractAlgebra.Generic.MatSpaceElem{<:MPolyRingElem{T}},
                           random_rothe::MatElem{T}) where T <: FieldElem
  evaluations = Set{Tuple{MPolyRingElem{T}, T}}([])
  generic_cols = cols[1:dep_col_index - 1]
  for g_col in generic_cols
    for i in 1:nrows(generic_rothe), j in g_col
      (is_unit(generic_rothe[i, j]) || iszero(generic_rothe[i, j])) && continue
      push!(evaluations, (generic_rothe[i, j], random_rothe[i, j]))
    end
  end
  return collect.(zip(evaluations...))
end

function _check_schur_complements(M::MatElem)
  n = ncols(M)
  P = M[1:n - 1, 1:n - 1]
  Q = M[1:n - 1, n]

  F = fraction_field(base_ring(M))
  y = solve(F.(P), F.(Q); side=:right)

  non_singular = false
  best_i = -1
  sum_terms = typemax(Int)
  for i in n:nrows(M)
    r = M[i, 1:n - 1]
    s = M[i, n]
    if !iszero(dot(r, y) - s)
      non_singular = true
      s = sum(length.(M[i, :]))
      if sum_terms > s
        best_i, sum_terms = i, s
      end
    end
  end
  return non_singular, best_i
end

function check_dep_cols(m::AbstractAlgebra.Generic.MatSpaceElem{T},
                        cols_to_check::Vector{Int},
                        dependent_columns::Vector{Int},
                        col_sets::Vector{Combination{Int}}) where T <: MPolyRingElem

  for col in cols_to_check
    if !iszero(m[:, col])
      col_indices = [col_index for col_index in 1:col if !(col_index in dependent_columns)]
      push!(col_indices, col)
      m_t = transpose(m[:, col_indices])
      rank_m_t = rank_dropped(m_t)
      !rank_m_t && return false
      m[:, col_indices] = transpose(m_t)

      for i=1:nrows(m)
        c = content(m[i:i, :])
        if !isone(c)
          m[i, :] = divexact(m[i:i, :], c)
        end
      end
    end
  end
  return true
end

function lex_min_col_basis(m::AbstractAlgebra.Generic.MatSpaceElem{T},
                           uhg::UniformHypergraph,
                           cols_to_check::Vector{Int},
                           dependent_columns::Vector{Int};
                           full_shift=false) where T <: MPolyRingElem
  for i=1:nrows(m)
    c = content(m[i:i, :])
    if !isone(c)
      m[i, :] = divexact(m[i:i, :], c)
    end
  end
  for j=1:ncols(m)
    c = content(m[:, j:j])
    if !isone(c)
      println("dividing")
      m[:, j] = divexact(m[:, j:j], c)
    end
  end

  n = n_vertices(uhg)
  k = face_size(uhg)
  col_sets = collect(combinations(n, k))[1:ncols(m)]

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
        elseif is_unit(m[ii, j])
          best_i = ii
          break
        elseif best_t > length(m[ii, j]) 
          best_t = length(m[ii, j])
          best_i = ii
        end
      end
      if best_i == 0
        # j is dependent on the columns before it
        m[:, j] .= zero(base_ring(m))
        if full_shift
          for col in cols_to_check
            if _domination(col_sets[j], col_sets[col])
              m[:, col] .= zero(base_ring(m))
            end
          end
        end
        iszero(m[:, cols_to_check]) && return true
        j += 1
        continue
      end
      if best_i > i
        m = swap_rows!(m, i, best_i)
      end
      break
    end

    for k=i+1:nrows(m)
      if iszero(m[k, j])
        continue
      end
      g, a, b = gcd_with_cofactors(m[k, j], m[i, j])
      m[k, :] = b*m[k:k, :] - a * m[i:i, :]
      m[k, :] = divexact(m[k:k, :], content(m[k:k, :]))
    end

    iszero(m[:, cols_to_check]) && return true

    for j=1:ncols(m)
      c = content(m[:, j:j])
      if !isone(c)
        m[:, j] = divexact(m[:, j:j], c)
      end
    end

    j += 1
  end

  return false
end

function rank_dropped(M::MatElem{<:MPolyRingElem})
  rk = 0
  #for i=1:nrows(M)
  #  c = content(M[i:i, :])
  #  if !isone(c)
  #    M[i, :] = divexact(M[i:i, :], c)
  #  end
  #end
  j = 1
  for i=1:nrows(M)
    best_j = 0
    best_t = typemax(Int)
    while j <= ncols(M) 
      best_i = [0]
      best_t = 0
      for ii = i:nrows(M)
        if is_zero_entry(M, ii, j)
          continue
        end
        if best_i == [0]
          best_i = ii
          best_t = length(M[ii,j])
        elseif is_unit(M[ii, j])
          best_i = ii
          break
        elseif best_t > length(M[ii, j])
          best_t = length(M[ii, j])
          best_i = ii
        # elseif best_t == length(M[ii, j])
          # push!(best_i, ii)
        end
      end

      if first(best_i) == 0
        j += 1
        continue
      end
      best_row_index = best_i
      # best_row_index = eargmin(ii -> -sum(iszero.(M[ii, j:end])), best_i;)
      if best_row_index > i
        M = swap_rows!(M, i, best_row_index)
      end
      break
    end
    if j > ncols(M)
      return iszero(M[end, :])
    end
    rk += 1

    for k=i+1:nrows(M)
      if iszero(M[k, j])
        continue
      end
      g, a, b = gcd_with_cofactors(M[k, j], M[i, j])
      M[k, :] = b*M[k:k, :] - a * M[i:i, :]
      M[k, :] = divexact(M[k:k, :], content(M[k:k, :]))
      if iszero(M[k, :])
        M = swap_rows!(M, k, nrows(M))
        return true
      end
    end
    
    j += 1
  end
  return false
end
