###############################################################################
#
#   Cartan Matrix Helpers
#   
###############################################################################

@doc raw"""
    cartan_matrix(fam::Symbol, rk::Int) -> ZZMatrix

Returns the Cartan matrix of finite type, where `fam` is the family ($A$, $B$, $C$, $D$, $E$, $F$ $G$)
and `rk` is the rank of the associated the root system.
The convention is $(a_ij) = (\langle \alpha_i^\vee, \alpha_j \rangle)$ for simple roots $\alpha_i$.
```jldoctest
julia> cartan_matrix(:B, 2)
[ 2   -1]
[-2    2]
```
"""
function cartan_matrix(fam::Symbol, rk::Int)
  if fam == :A
    @req rk >= 1 "type An requires rank rk to be at least 1"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
  elseif fam == :B
    @req rk >= 2 "type Bn requires rank rk to be at least  2"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk, rk - 1] = -2
  elseif fam == :C
    @req rk >= 2 "type Cn requires rank rk to be at least 2"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk - 1, rk] = -2
  elseif fam == :D
    @req rk >= 4 "type Dn requires rank to be at least 4"

    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 2)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk - 2, rk] = -1
    mat[rk, rk - 2] = -1
  elseif fam == :E
    @req rk in 6:8 "type En requires rank to be one of 6, 7, 8"
    mat = matrix(
      ZZ,
      [
        2 0 -1 0 0 0 0 0
        0 2 0 -1 0 0 0 0
        -1 0 2 -1 0 0 0 0
        0 -1 -1 2 -1 0 0 0
        0 0 0 -1 2 -1 0 0
        0 0 0 0 -1 2 -1 0
        0 0 0 0 0 -1 2 -1
        0 0 0 0 0 0 -1 2
      ],
    )

    if rk == 6
      mat = mat[1:6, 1:6]
    elseif rk == 7
      mat = mat[1:7, 1:7]
    end
  elseif fam == :F
    @req rk == 4 "type Fn requires rank to be 4"
    mat = matrix(ZZ, [2 -1 0 0; -1 2 -1 0; 0 -2 2 -1; 0 0 -1 2])
  elseif fam == :G
    @req rk == 2 "type Gn requires rank to be 2"
    mat = matrix(ZZ, [2 -3; -1 2])
  else
    throw(ArgumentError("unknown family"))
  end

  return mat
end

@doc raw"""
    cartan_matrix(types::Tuple{Symbol,Int}...) -> ZZMatrix
"""
function cartan_matrix(types::Tuple{Symbol,Int}...)
  blocks = ZZMatrix[cartan_matrix(type...) for type in types]

  return block_diagonal_matrix(blocks)
end

@doc raw"""
    cartan_to_coxeter_matrix(mat::ZZMatrix; check::Bool=true) -> Bool


"""
function cartan_to_coxeter_matrix(gcm; check::Bool=true)
  @req !check || is_cartan_matrix(gcm) "requires a generalized Cartan matrix"

  cm = identity_matrix(gcm)
  rk, _ = size(gcm)
  for i in 1:rk, j in 1:rk
    if i == j
      continue
    end

    if gcm[i, j] == 0
      cm[i, j] = 2
    elseif gcm[i, j] == -1
      cm[i, j] = 3
    elseif gcm[i, j] == -2
      cm[i, j] = 4
    elseif gcm[i, j] == -3
      cm[i, j] = 6
    end
  end

  return cm
end

@doc raw"""
    is_cartan_matrix(mat::ZZMatrix; generalized::Bool=true) -> Bool

Test if `mat` is a generalized Cartan matrix. The keyword argument `generalized`
can be set to `false` to restrict this to Cartan matrices of finite type.
```jldoctest
julia> is_cartan_matrix(ZZ[2 -2; -2 2])
true

julia> is_cartan_matrix(ZZ[2 -2; -2 2]; generalized=false)
false
"""
function is_cartan_matrix(mat::ZZMatrix; generalized::Bool=true)
  n, m = size(mat)
  if n != m
    return false
  end

  if !all(mat[i, i] == 2 for i in 1:n)
    return false
  end

  for i in 1:n
    for j in (i + 1):n
      if is_zero_entry(mat, i, j) && is_zero_entry(mat, j, i)
        continue
      end

      # mat[i,j] != 0 or mat[j,i] != 0, so both entries must be < 0
      if mat[i, j] >= 0 || mat[j, i] >= 0
        return false
      end

      # if we only want a generalized Cartan matrix, we are done with this pair of entries
      if generalized
        continue
      end

      if !(mat[i, j] * mat[j, i] in 1:3)
        return false
      end
    end
  end

  return true
end

@doc raw"""
    cartan_matrix_symmetrizer(gcm::ZZMatrix; check::Bool=true) -> Vector{ZZRingElem}

Return a vector $d$ such that $d_i a_{ij}$ is a symmtric matrix, where $a_{ij}$ are the entries of `gcm`.
`check` optionally verifies that `gcm` is indeed a generalized Cartan matrix.
```jldoctest
julia> cartan_symmetrizer(cartan_matrix(:B, 2))
2-element Vector{ZZRingElem}:
2
1
"""
function cartan_symmetrizer(gcm::ZZMatrix; check::Bool=true)
  @req !check || is_cartan_matrix(gcm) "requires a generalized Cartan matrix"
  rk = nrows(gcm)
  diag = ones(ZZRingElem, rk)
  for i in 1:rk, j in (i + 1):rk
    if !is_zero_entry(gcm, i, j)
      diag[i] = lcm(diag[i], gcm[i, j], gcm[j, i])
    end
  end

  return diag
end

@doc raw"""
    cartan_bilinear_form(gcm::ZZMatrix; check::Bool=true) -> ZZMatrix
    
```jldoctest
julia> cartan_bilinear_form(cartan_matrix(:B, 2))
[ 4   -2]
[-2    2]
"""
function cartan_bilinear_form(gcm::ZZMatrix; check::Bool=true)
  @req !check || is_cartan_matrix(gcm) "requires a generalized Cartan matrix"

  bil = deepcopy(gcm)
  sym = cartan_symmetrizer(gcm; check=false)
  for i in 1:length(sym)
    mul!(bil[i, :], sym[i])
  end
  return bil
end

@doc raw"""
    cartan_type(gcm::ZZMatrix; check::Bool=true) -> Vector{Tuple{Symbol, Int}}

Currently only works for Cartan matrices of finite type.
"""
function cartan_type(gcm::ZZMatrix; check::Bool=true)
  type, _ = cartan_type_with_ordering(gcm; check=check)
  return type
end

@doc raw"""
    cartan_type(gcm::ZZMatrix; check::Bool=true) -> Tuple{Vector{Tuple{Symbol, Int}}, Vector{Int}}

Currently only works for Cartan matrices of finite type.
"""
function cartan_type_with_ordering(gcm::ZZMatrix; check::Bool=true)
  @req !check || is_cartan_matrix(gcm; generalized=false) "requires Cartan matrix of finite type"

  rk = nrows(gcm)
  type = Tuple{Symbol,Int}[]

  # global information
  ord = sizehint!(Int[], rk) # ordering of the roots
  adj = [sizehint!(Int[], 3) for _ in 1:rk] # store adjacent roots

  # information about current type
  num = 0 # number of roots
  roots = sizehint!(Int[], rk)

  # used for traversal
  undone = trues(rk)
  plan = zeros(Int, 3) # a root is connected to up to 3 others
  head = 0 # index up to which we planned (cyclic)
  tail = 0 # index of plan which will be done next

  i, j = 1, 2
  while true
    while i == j || (j <= rk && is_zero_entry(gcm, i, j))
      j += 1
    end
    if j == rk + 1
      num += 1
      undone[i] = false
      push!(roots, i)

      if tail != head
        tail += 1
        if tail == 4
          tail = 1
        end
        i = plan[tail]
      else # nothing further is planned
        offset = length(ord) + 1
        if num == 1 # rank 1
          push!(type, (:A, 1))
          push!(ord, i)
        elseif num == 2 # rank 2
          i, j = roots[1], roots[2]
          if gcm[i, j] * gcm[j, i] == 1
            push!(type, (:A, 2))
            push!(ord, i, j)
          elseif gcm[i, j] == -2
            push!(type, (:C, 2))
            push!(ord, i, j)
          elseif gcm[j, i] == -2
            push!(type, (:B, 2))
            push!(ord, i, j)
          elseif gcm[i, j] == -3
            push!(type, (:G, 2))
            push!(ord, j, i)
          else
            push!(type, (:G, 2))
            push!(ord, i, j)
          end
        else # rank > 2
          # find start of the Dynkin graph
          v = roots[1] # default to the first root (this only matters for D4)
          for i in roots
            j = adj[i][1]
            if length(adj[i]) == 1 && length(adj[j]) < 3 && gcm[i, j] * gcm[j, i] == 1
              v = i
              break
            end
          end
          push!(ord, v)

          n = 1
          adj3 = false # true if found a root with 3 adjacents
          while n < num
            nv = v
            for vv in adj[v]
              filter!(x -> x != v, adj[vv])
              push!(ord, vv)
              n += 1
              if length(adj[vv]) > 0
                nv = vv
                if length(adj[vv]) == 2 # +1 for the predecessor
                  adj3 = true
                end
              end
            end
            v = nv
          end

          if adj3
            if isempty(adj[ord[end]]) && isempty(adj[ord[end - 1]])
              push!(type, (:D, num))
            else
              push!(type, (:E, num))
            end
          elseif num == 4 &&
            gcm[ord[end - 2], ord[end - 1]] * gcm[ord[end - 1], ord[end - 2]] > 1
            push!(type, (:F, 4))
          elseif gcm[ord[end - 1], ord[end]] * gcm[ord[end], ord[end - 1]] == 1
            push!(type, (:A, num))
          elseif gcm[ord[end - 1], ord[end]] == -1
            push!(type, (:B, num))
          else
            push!(type, (:C, num))
          end
        end

        # find next component
        i = findfirst(undone)
        if isnothing(i)
          break
        end

        # reset number of roots
        num = 0
        empty!(roots)
      end

      j = 1
      continue
    end

    # plan to visit j if undone
    if undone[j]
      head += 1
      if head == 4
        head = 1
      end
      plan[head] = j
    end

    push!(adj[i], j)
    j += 1
  end

  return type, ord
end
