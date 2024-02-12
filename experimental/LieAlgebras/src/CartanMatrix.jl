###############################################################################
#
#   Cartan Matrix Functionality
#   - Cartan matrices (and root orderings) are taken from Bourbaki
#
###############################################################################

@doc raw"""
    cartan_matrix(fam::Symbol, rk::Int) -> ZZMatrix

Returns the Cartan matrix of finite type, where `fam` is the family ($A$, $B$, $C$, $D$, $E$, $F$ $G$)
and `rk` is the rank of the associated the root system; for $B$ and $C$ the rank has to be at least 2, for $D$ at least 4.
The convention is $(a_{ij}) = (\langle \alpha_i^\vee, \alpha_j \rangle)$ for simple roots $\alpha_i$.

# Example
```jldoctest
julia> cartan_matrix(:B, 2)
[ 2   -1]
[-2    2]

julia> cartan_matrix(:C, 2)
[ 2   -2]
[-1    2]
```
"""
function cartan_matrix(fam::Symbol, rk::Int)
  @req is_cartan_type(fam, rk) "Invalid cartan type"
  if fam == :A
    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
  elseif fam == :B
    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk, rk - 1] = -2
  elseif fam == :C
    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 1)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk - 1, rk] = -2
  elseif fam == :D
    mat = diagonal_matrix(ZZ(2), rk)
    for i in 1:(rk - 2)
      mat[i + 1, i], mat[i, i + 1] = -1, -1
    end
    mat[rk - 2, rk] = -1
    mat[rk, rk - 2] = -1
  elseif fam == :E
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
    mat = matrix(ZZ, [2 -1 0 0; -1 2 -1 0; 0 -2 2 -1; 0 0 -1 2])
  elseif fam == :G
    mat = matrix(ZZ, [2 -3; -1 2])
  else
    error("Unreachable")
  end

  return mat
end

@doc raw"""
    cartan_matrix(type::Tuple{Symbol,Int}...) -> ZZMatrix

Returns a block diagonal matrix of indecomposable Cartan matrices as defined by type.
For allowed values see `cartan_matrix(fam::Symbol, rk::Int)`.

# Example
```jldoctest
julia> cartan_matrix((:A, 2), (:B, 2))
[ 2   -1    0    0]
[-1    2    0    0]
[ 0    0    2   -1]
[ 0    0   -2    2]
```
"""
function cartan_matrix(type::Tuple{Symbol,Int}...)
  @req length(type) > 0 "At least one type is required"

  blocks = [cartan_matrix(t...) for t in type]
  return block_diagonal_matrix(blocks)
end

@doc raw"""
    is_cartan_matrix(mat::ZZMatrix; generalized::Bool=true) -> Bool

Checks if `mat` is a generalized Cartan matrix. The keyword argument `generalized`
can be set to `false` to restrict this to Cartan matrices of finite type.

# Example
```jldoctest
julia> is_cartan_matrix(ZZ[2 -2; -2 2])
true

julia> is_cartan_matrix(ZZ[2 -2; -2 2]; generalized=false)
false
```
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
    end
  end

  if generalized
    return true
  end

  return is_positive_definite(mat)
end

@doc raw"""
    cartan_symmetrizer(gcm::ZZMatrix; check::Bool=true) -> Vector{ZZRingElem}

Return a vector $d$ of coprime integers such that $(d_i a_{ij})_{ij}$ is a symmetric matrix,
where $a_{ij}$ are the entries of the Cartan matrix `gcm`.
The keyword argument `check` can be set to `false` to skip verification whether `gcm` is indeed a generalized Cartan matrix.

# Example
```jldoctest
julia> cartan_symmetrizer(cartan_matrix(:B, 2))
2-element Vector{ZZRingElem}:
 2
 1
```
"""
function cartan_symmetrizer(gcm::ZZMatrix; check::Bool=true)
  @req !check || is_cartan_matrix(gcm) "Requires a generalized Cartan matrix"
  rk = nrows(gcm)
  diag = ones(ZZRingElem, rk)

  # used for traversal
  undone = trues(rk)
  plan = zeros(Int, rk) # roots planned sorted asc grouped by component
  head = 0
  tail = 0

  # we collect roots of the same length
  # once we know if they are short or long we scale approriately
  while any(undone)
    if head == tail
      head += 1
      plan[head] = findfirst(undone)
      undone[plan[head]] = false
    end

    prev = head
    i = plan[head]
    for j in 1:rk
      if i == j
        continue
      end

      if !undone[j] || is_zero_entry(gcm, i, j)
        continue
      end

      head += 1
      plan[head] = j
      undone[j] = false

      if diag[i] * gcm[i, j] == diag[j] * gcm[j, i]
        continue
      elseif gcm[i, j] == gcm[j, i]
        diag[i] = lcm(diag[i], diag[j])
        diag[j] = diag[i]
        continue
      end

      if gcm[j, i] < -1
        tail += 1
        v = -gcm[j, i]
        while tail < head
          diag[plan[tail]] *= v
          tail += 1
        end
      end
      if gcm[i, j] < -1
        diag[j] *= -gcm[i, j]
        tail = head - 1
      end
    end

    # we found new roots, meaning we are done with this component of the root system
    if prev == head
      tail = head
    end
  end

  return diag
end

@doc raw"""
    cartan_bilinear_form(gcm::ZZMatrix; check::Bool=true) -> ZZMatrix

Returns the matrix of the symmetric bilinear form associated to the Cartan matrix from `cartan_symmetrizer`.
The keyword argument `check` can be set to `false` to skip verification whether `gcm` is indeed a generalized Cartan matrix.

# Example
```jldoctest
julia> cartan_bilinear_form(cartan_matrix(:B, 2))
[ 4   -2]
[-2    2]
```
"""
function cartan_bilinear_form(gcm::ZZMatrix; check::Bool=true)
  sym = cartan_symmetrizer(gcm; check)
  bil = deepcopy(gcm)
  for i in 1:length(sym)
    mul!(view(bil, i:i, :), sym[i])
  end
  return bil
end

@doc raw"""
    cartan_type(gcm::ZZMatrix; check::Bool=true) -> Vector{Tuple{Symbol, Int}}

Returns the Cartan type of a Cartan matrix `gcm` (currently only Cartan matrices of finite type are supported).
This function is left inverse to `cartan_matrix`, i.e. in the case of isomorphic types (e.g. $B_2$ and $C_2$)
the ordering of the roots does matter (see the example below).
The keyword argument `check` can be set to `false` to skip verification whether `gcm` is indeed a Cartan matrix of finite type.

# Example
```jldoctest
julia> cartan_type(ZZ[2 -1; -2 2])
1-element Vector{Tuple{Symbol, Int64}}:
 (:B, 2)

julia> cartan_type(ZZ[2 -2; -1 2])
1-element Vector{Tuple{Symbol, Int64}}:
 (:C, 2)
```
"""
function cartan_type(gcm::ZZMatrix; check::Bool=true)
  type, _ = cartan_type_with_ordering(gcm; check=check)
  return type
end

@doc raw"""
    cartan_type_with_ordering(gcm::ZZMatrix; check::Bool=true) -> Vector{Tuple{Symbol, Int}}, Vector{Int}

Returns the Cartan type of a Cartan matrix `gcm` together with a vector indicating a canonical ordering
of the roots in the Dynkin diagram (currently only Cartan matrices of finite type are supported).
The keyword argument `check` can be set to `false` to skip verification whether `gcm` is indeed a Cartan matrix of finite type.
# Example
```jldoctest
julia> cartan_matrix(:E, 6)
[ 2    0   -1    0    0    0]
[ 0    2    0   -1    0    0]
[-1    0    2   -1    0    0]
[ 0   -1   -1    2   -1    0]
[ 0    0    0   -1    2   -1]
[ 0    0    0    0   -1    2]

julia> cartan_type_with_ordering(cartan_matrix(:E, 6))
([(:E, 6)], [1, 3, 4, 2, 5, 6])

julia> cartan_type_with_ordering(ZZ[2 0 -1 0; 0 2 0 -2; -2 0 2 0; 0 -1 0 2])
([(:B, 2), (:C, 2)], [1, 3, 2, 4])
```
"""
function cartan_type_with_ordering(gcm::ZZMatrix; check::Bool=true)
  @req !check || is_cartan_matrix(gcm; generalized=false) "requires Cartan matrix of finite type"

  rk = nrows(gcm)
  type = Tuple{Symbol,Int}[]

  # global information
  ord = sizehint!(Int[], rk) # ordering of the roots
  adj = [sizehint!(Int[], 3) for _ in 1:rk] # store adjacent roots, adj[i] is ordered asc

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
            push!(ord, i, j)
          else
            push!(type, (:G, 2))
            push!(ord, j, i)
          end
        else # rank > 2
          # find start of the Dynkin graph
          v = 0
          for i in roots
            j = adj[i][1]
            if length(adj[i]) == 1 && gcm[i, j] * gcm[j, i] == 1
              if length(adj[j]) == 1 || (length(adj[j]) == 2 && gcm[j, adj[j][2]] == -1)
                v = i
                break
              elseif v == 0
                v = i
              end
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

@doc raw"""
    is_cartan_type(fam::Symbol, rk::Int) -> Bool

Checks if the pair (`fam`, `rk`) is a valid Cartan type,
i.e. one of `A_l` (l >= 1), `B_l` (l >= 2), `C_l` (l >= 2), `D_l` (l >= 4), `E_6`, `E_7`, `E_8`, `F_4`, `G_2`.
"""
function is_cartan_type(fam::Symbol, rk::Int)
  fam in [:A, :B, :C, :D, :E, :F, :G] || return false
  fam == :A && rk >= 1 && return true
  fam == :B && rk >= 2 && return true
  fam == :C && rk >= 2 && return true
  fam == :D && rk >= 4 && return true
  fam == :E && rk in [6, 7, 8] && return true
  fam == :F && rk == 4 && return true
  fam == :G && rk == 2 && return true
  return false
end
