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
    cartan_matrix(type::Vector{Tuple{Symbol,Int}}) -> ZZMatrix

Returns a block diagonal matrix of indecomposable Cartan matrices as defined by type.
For allowed values see `cartan_matrix(fam::Symbol, rk::Int)`.

# Example
```jldoctest
julia> cartan_matrix([(:A, 2), (:B, 2)])
[ 2   -1    0    0]
[-1    2    0    0]
[ 0    0    2   -1]
[ 0    0   -2    2]
```
"""
function cartan_matrix(type::Vector{Tuple{Symbol,Int}})
  @req length(type) > 0 "At least one type is required"

  blocks = [cartan_matrix(t...) for t in type]
  return block_diagonal_matrix(blocks)
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

  return cartan_matrix(collect(type))
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

The order of returned components is, in general, not unique and might change between versions.
If this function is called with the output of `cartan_matrix(type)`, it will keep the order of `type`.

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
The keyword argument `check` can be set to `false` to skip verification whether `gcm` is indeed a
Cartan matrix of finite type.

The order of returned components and the ordering is, in general, not unique and might change between versions.
If this function is called with the output of `cartan_matrix(type)`, it will keep the order of `type` and the
returned ordering will be the identity.

# Example
```jldoctest
julia> cartan_type_with_ordering(cartan_matrix(:E, 6))
([(:E, 6)], [1, 2, 3, 4, 5, 6])

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
  adj = [[j for j in 1:rk if i != j && !is_zero_entry(gcm, i, j)] for i in 1:rk] # adjacency list
  done = falses(rk) # whether a root is already in a component

  for v0 in 1:rk
    done[v0] && continue

    # rank 1
    if length(adj[v0]) == 0
      push!(type, (:A, 1))
      push!(ord, v0)
      done[v0] = true
      continue
    end

    # rank 2
    if length(adj[v0]) == 1 && length(adj[only(adj[v0])]) == 1
      v1 = only(adj[v0])
      if gcm[v0, v1] * gcm[v1, v0] == 1
        push!(type, (:A, 2))
        push!(ord, v0, v1)
      elseif gcm[v0, v1] == -2
        push!(type, (:C, 2))
        push!(ord, v0, v1)
      elseif gcm[v1, v0] == -2
        push!(type, (:B, 2))
        push!(ord, v0, v1)
      elseif gcm[v0, v1] == -3
        push!(type, (:G, 2))
        push!(ord, v0, v1)
      elseif gcm[v1, v0] == -3
        push!(type, (:G, 2))
        push!(ord, v1, v0)
      else
        error("unreachable")
      end
      done[v0] = true
      done[v1] = true
      continue
    end

    # rank > 2
    # do a DFS to find the whole component
    comp = [v0]
    todo = [v0]
    done[v0] = true
    while !isempty(todo)
      v = pop!(todo)
      for w in adj[v]
        if !done[w]
          push!(comp, w)
          push!(todo, w)
          done[w] = true
        end
      end
    end
    sort!(comp)
    len_comp = length(comp)

    deg3 = findfirst(v -> length(adj[v]) == 3, comp)
    if isnothing(deg3)
      # case A, B, C, F

      # find start of the Dynkin graph
      start = 0
      for v1 in filter(v -> length(adj[v]) == 1, comp)
        v2 = only(adj[v1])
        gcm[v1, v2] * gcm[v2, v1] == 1 || continue          # discard right end of B and C
        if len_comp == 4
          v3 = only(filter(!=(v1), adj[v2]))
          gcm[v2, v3] == -1 || continue                   # discard right end of F
        end

        # found start
        start = v1
        break
      end
      @assert start != 0

      # find the path
      path = [start, only(adj[start])]
      for _ in 1:(len_comp - 2)
        push!(path, only(filter(!=(path[end - 1]), adj[path[end]])))
      end
      # determine type
      if len_comp == 4 && gcm[path[3], path[2]] == -2
        push!(type, (:F, 4))
      elseif gcm[path[end - 1], path[end]] == -2
        push!(type, (:C, len_comp))
      elseif gcm[path[end], path[end - 1]] == -2
        push!(type, (:B, len_comp))
      else
        push!(type, (:A, len_comp))
      end
      append!(ord, path)
    else
      # case D or E

      # find the three paths
      v_deg3 = comp[deg3]
      paths = [[v_deg3, v_n] for v_n in adj[v_deg3]]
      for path in paths
        while length(adj[path[end]]) == 2
          push!(path, only(filter(!=(path[end - 1]), adj[path[end]])))
        end
        popfirst!(path)
      end
      sort!(paths; by=length)
      @assert sum(length, paths) + 1 == len_comp
      # determine type
      if length(paths[2]) == 1
        # case D
        push!(type, (:D, len_comp))
        if len_comp == 4
          push!(ord, only(paths[1]), v_deg3, only(paths[2]), only(paths[3]))
        else
          append!(ord, reverse!(paths[3]))
          push!(ord, v_deg3, only(paths[1]), only(paths[2]))
        end
      elseif length(paths[2]) == 2
        # case E
        push!(type, (:E, len_comp))
        push!(ord, paths[2][2], only(paths[1]), paths[2][1], v_deg3)
        append!(ord, paths[3])
      else
        error("unreachable")
      end
    end
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
