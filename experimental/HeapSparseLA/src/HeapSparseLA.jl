# A heap based implementation of sparse row vectors for elements
# which are big in the sense that copying them around in memory 
# for sorting is expensive.
include("Types.jl")

### `HeapNode`
# A minimal data structure for nodes in binary trees holding the data
# for a sparse vector. 
#
# The nodes in the tree must be partially ordered according to `_heap_leq` 
# below. Every node contains an `index` and a non-zero `content`, i.e. 
# the value of the vector at this index. 
#
# There is an additional field 
# `cofactor`. If it is not `nothing`, then it indicates that the whole 
# subsequent tree, including this node, must be multiplied with the value 
# in `cofactor` to obtain the actual values.
#
# Note: Any index may occur multiple times among the nodes in a tree!
# In general the design is made to allow for lazy evaluation.

### user facing getters for `HeapNode`s
content(n::HeapNode) = n.content
cofactor(n::HeapNode{T}) where T = n.cofactor === nothing ? one(content(n)) : n.cofactor::T
index(n::HeapNode) = n.index
has_left_child(n::HeapNode) = (n.left !== nothing)
has_right_child(n::HeapNode) = (n.right !== nothing)
is_leaf(n::HeapNode) = !has_left_child(n) && !has_right_child(n)
left_child(n::HeapNode{T}) where {T} = n.left::HeapNode{T}
right_child(n::HeapNode{T}) where {T} = n.right::HeapNode{T}

### basic functionality for `HeapNode`s
tree_size(n::HeapNode) = (is_leaf(n) ? 1 : 1 + tree_size(n.left) + tree_size(n.right))
tree_size(n::Nothing) = 0

function ancestry_depth(n::HeapNode)
  !has_left_child(n) && !has_right_child(n) && return 0
  n1 = n2 = 0
  has_left_child(n) && (n1 = ancestry_depth(left_child(n)))
  has_right_child(n) && (n2 = ancestry_depth(right_child(n)))
  return maximum([n1, n2]) + 1
end

# TODO: Overwrite hashing!
function ==(a::HeapNode, b::HeapNode)
  return index(a) == index(b) && content(a) == content(b)
end

function show(io::IO, n::HeapNode; cofactor=1)
  show_recursive(io, n)
end

function show_recursive(io::IO, n::HeapNode; indent::Int=0, cofactor=1)
  new_cof = n.cofactor === nothing ? cofactor : cofactor*n.cofactor
  offset = prod(" " for i in 1:indent; init = "")
  println(offset*"$(index(n)) -> $(new_cof*content(n))")
  has_left_child(n) && show_recursive(io, left_child(n), indent=indent+1, cofactor=new_cof)
  has_right_child(n) && show_recursive(io, right_child(n), indent=indent+1, cofactor=new_cof)
end

function copy_tree(n::HeapNode)
  result = HeapNode(index(n), content(n), cofactor=n.cofactor)
  result.left = copy_tree(n.left)
  result.right = copy_tree(n.right)
  return result
end

copy_tree(n::Nothing) = nothing

###
# Comparison of nodes.
#
# We put an ordering on the nodes which first compares 
# the indices. This can be refined by an arbitrary 
# ordering on the elements of the ring by overwriting 
# the method `_heap_leq` for this type. For instance, 
# one could sort polynomials by degree.
function _heap_leq(a::HeapNode, b::HeapNode)
  index(a) < index(b) && return true
  index(a) > index(b) && return false
  return _heap_leq(content(a), content(b))
end

function _heap_leq(a::Any, b::Any)
  return true
end

function _heap_leq(a::HeapNode, b::Nothing)
  return true
end

function _heap_leq(a::Nothing, b::HeapNode)
  return false
end

function _heap_leq(a::Nothing, b::Nothing)
  return true
end

### Merging trees (identified by the root nodes) and preserving 
# the partial ordering by `_heap_leq`.
function Base.merge!(a::HeapNode{T}, b::HeapNode{T}) where {T}
  if _heap_leq(a, b)
    # a is the new top
    a.left = merge!(a.left, a.right)
    scale_cofactor!(a.left, a.cofactor)
    if a.cofactor!==nothing
      a.content *= a.cofactor
      a.cofactor = nothing
    end
    a.right = b
    return a
  else
    b.left = merge!(b.left, b.right)
    scale_cofactor!(b.left, b.cofactor)
    if b.cofactor!==nothing
      b.content *= b.cofactor
      b.cofactor = nothing
    end
    b.right = a
    return b
  end
end

function Base.merge!(a::HeapNode, b::Nothing)
  return a
end

function Base.merge!(a::Nothing, b::HeapNode)
  return b
end

function Base.merge!(a::Nothing, b::Nothing)
  return nothing
end

function all_indices(n::HeapNode)
  is_leaf(n) && return [index(n)]
  return vcat([index(n)], all_indices(n.left), all_indices(n.right))
end

all_indices(n::Nothing) = Int[]

### Write out the full vector to a dictionary without duplicate entries.
as_dictionary(n::Nothing, result::Dict{Int, T}, c::Any=nothing) where {T} = result

function as_dictionary(n::HeapNode{T}, result::Dict{Int, T}, cofactor::Nothing=nothing) where {T}
  new_cof = n.cofactor === nothing ? one(content(n)) : n.cofactor
  i = index(n)
  if i in keys(result)
    result[i] = result[i] + new_cof * content(n)
  else
    result[i] = new_cof * content(n)
  end
  as_dictionary(n.left, result, n.cofactor)
  as_dictionary(n.right, result, n.cofactor)
  return result
end

function as_dictionary(n::HeapNode{T}, result::Dict{Int, T}, cofactor::T) where {T}
  new_cof = n.cofactor === nothing ? cofactor : n.cofactor*cofactor
  i = index(n)
  if i in keys(result)
    result[i] += new_cof * content(n)
  else
    result[i] = new_cof * content(n)
  end
  as_dictionary(n.left, result, new_cof)
  as_dictionary(n.right, result, new_cof)
  return result
end

### `HeapSRow`
# A data structure to realize sparse vectors with entries in a ring 
# based on binary trees (using `HeapNode` above) and supporting lazy 
# evaluation.

### basic getters for `HeapSRow`
base_ring(v::HeapSRow) = v.entry_parent
top_node(v::HeapSRow{T}) where {T} = v.top::HeapNode{T}

function is_consolidated(v::HeapSRow) 
  v.is_consolidated && return true
  if is_empty(v)
    v.is_consolidated = true
    return true
  end

  n = top_node(v)
  result, _ = _has_duplicate(n, Int[])
  if !result
    v.result = true
    return true
  end
  return false
end

function getindex(v::HeapSRow{T}, i::Int) where {T}
  is_empty(v) && return zero(base_ring(v))
  haskey(index_cache(v), i) && return index_cache(v)[i]
  n = top_node(v)
  v.top_is_consolidated && index(n) == i && return content(n)
  c, rem = pop_index!(n, i)
  index_cache(v)[i] = c
  if iszero(c)
    v.top = rem
  else
    v.top = merge!(rem, HeapNode(i, c))
  end
  return c
end


function index_cache(v::HeapSRow{T}) where {T}
  if !isdefined(v, :index_cache)
    v.index_cache = Dict{Int, T}()
  end
  return v.index_cache
end

function pivot(v::HeapSRow)
  consolidate_top!(v)
  isempty(v) && return 0, zero(base_ring(v))
  return index(top_node(v)), content(top_node(v))
end

### Constructors for `HeapSRow`
function HeapSRow(R::NCRing, n::Nothing)
  return HeapSRow(R)
end

function HeapSRow(R::NCRing, list::Vector{Pair{Int, T}}) where {T}
  return HeapSRow(R, [(a, b) for (a, b) in list])
end

function HeapSRow(v::SRow{T}) where {T}
  result = HeapSRow(base_ring(v))
  result.top_is_consolidated = false
  open = HeapNode{T}[]
  
  side = false
  for (i, c) in v
    n = HeapNode(i, c)
    if result.top === nothing 
      result.top = n
      push!(open, n)
      continue
    end
    if !side
      first(open).left = n
      push!(open, n)
      side = true
    else
      popfirst!(open).right = n
      push!(open, n)
      side = false
    end
  end
  return result
end

function HeapSRow(R::NCRing, list::Vector{Tuple{Int, T}}) where {T}
  result = HeapSRow(R)
  isempty(list) && return result
  if isone(length(list))
    i, c = first(list)
    result.top = HeapNode(i, c)
    return result
  end

  n = length(list)
  h = div(n, 2)
  result = HeapSRow(R, list[1:h]) + HeapSRow(R, list[h+1:end])
  return result
end

### functionality for `HeapSRow`s

# Find the first index `i` for which `v` has a non-zero entry, 
# compute the value `c` for this entry, remove that entry 
# from `v` and return a node holding the information `(i, c)`.
function pop_top_node!(v::HeapSRow)
  isempty(v) && return nothing
  n = top_node(v)
  i = index(n)
  a, left = pop_index!(n.left, i)
  b, right = pop_index!(n.right, i)
  n.content += a
  n.content += b
  n.cofactor !== nothing && (n.content *= n.cofactor)
  v.top = merge!(left, right)
  scale_cofactor!(v.top, n.cofactor)
  n.cofactor = nothing
  n.left = nothing
  n.right = nothing
  return n
end

function pop_top!(v::HeapSRow)
  n = pop_top_node!(v)
  return index(n), content(n)
end

# For every index `i` gather all occurences of `i` in the tree, 
# compute the actual value `c` of `v` for this index and, if it 
# is not zero, store it in a single node. Compile all these 
# polished nodes in a balanced tree.
function consolidate!(v::HeapSRow{T}) where T
  v.is_consolidated && return v
  result = HeapSRow(base_ring(v), pop_top_node!(v))
  side = false
  open = [result.top]
  while !isempty(v)
    n = pop_top_node!(v)
    if !side
      first(open).left = n
      push!(open, n)
      side = true
    else
      popfirst!(open).right = n
      push!(open, n)
      side = false
    end
  end
  result.is_consolidated = true
  return result
end

+(v::HeapSRow, w::HeapSRow) = add!(copy(v), copy(w))

function add!(v::HeapSRow{T}, w::HeapSRow{T}) where {T}
  @assert base_ring(v) === base_ring(w)
  R = base_ring(v)
  isempty(v) && return w
  isempty(w) && return v
  nt = merge!(top_node(v), top_node(w))
  return HeapSRow(R, nt)
end
  
function addmul!(v::HeapSRow{T}, w::HeapSRow{T}, a::T) where {T}
  return add!(v, mul!(a, w))
end

function copy(v::HeapSRow)
  isempty(v) && return HeapSRow(base_ring(v))
  n = top_node(v)
  nn = copy_tree(n)
  return HeapSRow(base_ring(v), nn)
end

function *(a::T, v::HeapSRow{T}) where {T}
  return mul!(a, copy(v))
end

function mul!(a::T, v::HeapSRow{T}) where {T}
  isempty(v) && return v
  iszero(a) && return HeapSRow(base_ring(v))
  isone(a) && return v
  # try to preserve consolidation
  if v.top_is_consolidated
    v.top.content *= a
    scale_cofactor!(v.top.left, a)
    scale_cofactor!(v.top.right, a)
  else
    scale_cofactor!(v.top, a)
  end
  # clear the cache
  v.index_cache = Dict{Int, T}()
  return v
end

*(a, v::HeapSRow) = base_ring(v)(a)*v
mul!(a, v::HeapSRow) = mul!(base_ring(v)(a), v)

# TODO: deprecated?
mult_content_rec_left!(a, n::Nothing) = nothing

function mult_content_rec_left!(a::T, n::HeapNode{T}) where {T}
  iszero(a) && return nothing
  isone(a) && return n
  n.content *= a
  mult_content_rec_left!(a, n.left)
  mult_content_rec_left!(a, n.right)
  is_zero(n.content) && return merge!(n.left, n.right)
  return n
end

-(v::HeapSRow, w::HeapSRow) = addmul!(copy(v), copy(w), -one(base_ring(w)))
-(w::HeapSRow) = -one(base_ring(w))*w

# Remove all nodes with `index` `i` from `v`, 
# sum them up and return a pair `(c, w)` consisting of 
# the value `c` of `v` at `i` and the remainder vector `w`.
function pop_index!(v::HeapSRow, i::Int)
  R = base_ring(v)
  isempty(v) && return zero(R), v
  c, n = pop_index!(top_node(v), i)
  return c, HeapSRow(R, n)
end

function pop_index!(n::HeapNode{T}, i::Int) where {T}
  # early abort: if the index of the node is bigger than i
  # then everything that follows will also have greater index
  index(n) > i && return zero(content(n)), n
  if index(n) == i
    result = content(n)
    a, left = pop_index!(n.left, i)
    b, right = pop_index!(n.right, i)
    rem = merge!(left, right)
    result += a
    result += b
    n.cofactor !== nothing && (result *= n.cofactor)
    return result, scale_cofactor!(rem, n.cofactor)
  end

  a, rem_left = pop_index!(n.left, i)
  b, rem_right = pop_index!(n.right, i)
  result = a + b
  n.cofactor !== nothing && (result *= n.cofactor)
  n.left = rem_left
  n.right = rem_right
  return result, n
end

function pop_index!(n::Nothing, i::Int)
  return 0, nothing
end

function scale_cofactor!(n::HeapNode{T}, a::T) where T
  if n.cofactor === nothing 
    n.cofactor = a
  else
    n.cofactor *= a
  end
  return n
end

scale_cofactor!(n::Nothing, a::Any) = nothing
scale_cofactor!(n::HeapNode, a::Nothing) = n

function consolidate_top!(v::HeapSRow)
  v.top_is_consolidated && return v
  isempty(v) && return v
  n = top_node(v)
  i = index(n)
  cof = n.cofactor
  if n.cofactor !== nothing 
    n.content *= n.cofactor
    n.cofactor = nothing
    scale_cofactor!(n.left, cof)
    scale_cofactor!(n.right, cof)
  end
  a, left = pop_index!(n.left, i)
  b, right = pop_index!(n.right, i)
  n.content += a + b
  n.left = left
  n.right = right

  if iszero(n.content)
    v.top = merge!(n.left, n.right)
    return consolidate_top!(v)
  end

  v.top_is_consolidated = true
  return v
end

isempty(v::HeapSRow) = v.top === nothing

function has_index(v::HeapSRow{T}, i::Int) where {T}
  isempty(v) && return false
  return _has_index(top_node(v), i)
end

function _has_index(n::HeapNode, i::Int)
  index(n) == i && return true
  has_left_child(n) && _has_index(left_child(n), i) && return true
  has_right_child(n) && _has_index(right_child(n), i) && return true
  return false
end

function all_indices(v::HeapSRow)
  isempty(v) && return Int[]
  return all_indices(top_node(v))
end

function as_dictionary(v::HeapSRow)
  result = Dict{Int, elem_type(base_ring(v))}()
  isempty(v) && return result
  result = as_dictionary(top_node(v), result)
  v.index_cache = result
  return result
end

function iszero(v::HeapSRow)
  isempty(v) && return true
  return all(iszero(c) for c in values(as_dictionary(v)))
end

function ==(v::HeapSRow{T}, w::HeapSRow{T}) where {T}
  isempty(v) && return iszero(w)
  isempty(w) && return iszero(v)
  pivot(v) == pivot(w) || return false

  d1 = as_dictionary(v)
  d2 = as_dictionary(w)
  for (k, c) in d1
    if !(k in keys(d2))
      iszero(c) || return false
      continue
    end
    c == d2[k] || return false
    delete!(d2, k)
  end
  return all(iszero(c) for c in values(d2))
end

### Sparse matrices based on `HeapSRow`

### basic getters for `HeapSMat`
base_ring(A::HeapSMat) = A.base_ring
nrows(A::HeapSMat) = A.nrows
ncols(A::HeapSMat) = A.ncols

function getindex(A::HeapSMat, i::Int, j::Int)
  @assert j <= ncols(A)
  return A.rows[i][j]
end

getindex(A::HeapSMat, i::Int) = A.rows[i]

function setindex!(A::HeapSMat, v::HeapSRow, i::Int)
  if A.transpose !== nothing 
    A.transpose.transpose = nothing
    A.transpose = nothing
  end
  A.rows[i] = v
end

function push!(A::HeapSMat, v::HeapSRow)
  if A.transpose !== nothing 
    A.transpose.transpose = nothing
    A.transpose = nothing
  end
  push!(A.rows, v)
  A.nrows += 1
  return A
end

function dense_row(v::HeapSRow, n::Int)
  return [v[i] for i in 1:n]
end

function sparse_row(v::HeapSRow)
  return sparse_row(base_ring(v), collect((i, c) for (i, c) in as_dictionary(v) if !iszero(c)))
end

function consolidate!(A::HeapSMat{T}) where T
  for (i, v) in enumerate(A.rows)
    A[i] = consolidate!(v)
  end
  return A
end

function dense_matrix(A::HeapSMat)
  result = zero(matrix_space(base_ring(A), nrows(A), ncols(A)))
  for i in 1:nrows(A)
    for (j, c) in as_dictionary(A[i])
      result[i, j] = c
    end
  end
  return result
end

function ==(A::HeapSMat, B::HeapSMat)
  base_ring(A) === base_ring(B) || return false
  nrows(A) == nrows(B) || return false
  ncols(A) == ncols(B) || return false
  return all(A[i] == B[i] for i in 1:nrows(A))
end

function *(A::HeapSMat{T}, B::HeapSMat{T}) where {T}
  return mul!(consolidate!(copy(A)), B)
end

function mul!(A::HeapSMat{T}, B::HeapSMat{T}) where {T}
  R = base_ring(A)
  @assert R === base_ring(B)
  for (i, v) in enumerate(A.rows)
    A[i] = mul!(v, B)
  end
  A.ncols = ncols(B)
  return A
end

function *(v::HeapSRow, A::HeapSMat)
  return mul!(copy(v), A)
end

function mul!(v::HeapSRow, A::HeapSMat)
  R = base_ring(A)
  @assert R === base_ring(v)
  result = HeapSRow(R)
  while !isempty(v)
    i, c = pop_top!(v)
    result = add!(result, c*A[i])
  end
  return result
end

function mul!(a, A::HeapSMat{T}) where T
  return mul!(base_ring(A)(a), A)
end

function mul!(a::T, A::HeapSMat{T}) where T
  if A.transpose !== nothing 
    A.transpose.transpose = nothing
    A.transpose = nothing
  end
  for v in A.rows
    mul!(a, v)
  end
  return A
end

function copy(A::HeapSMat)
  R = base_ring(A)
  result = HeapSMat(R, 0, ncols(A))
  for v in A.rows
    push!(result, copy(v))
  end
  return result
end

function *(a, A::HeapSMat)
  return mul!(a, copy(A))
end

function show(io::IO, A::HeapSMat)
  println(io, "sparse ($(nrows(A)) x $(ncols(A)))-matrix over $(base_ring(A))")
end

function add!(A::HeapSMat, B::HeapSMat)
  if A.transpose !== nothing 
    A.transpose.transpose = nothing
    A.transpose = nothing
  end
  base_ring(A) === base_ring(B) || return false
  nrows(A) == nrows(B) || return false
  ncols(A) == ncols(B) || return false
  for (i, v) in enumerate(B.rows)
    A[i] = add!(A[i], v)
  end
  return A
end

+(A::HeapSMat, B::HeapSMat) = add!(copy(A), B)
-(A::HeapSMat, B::HeapSMat) = add!(-B, A)
-(A::HeapSMat) = -one(base_ring(A))*A

function getindex(A::HeapSMat, ind::Vector{Int})
  return HeapSMat(base_ring(A), length(ind), ncols(A), A.rows[ind])
end

function getindex(A::HeapSMat, r::UnitRange)
  return HeapSMat(base_ring(A), length(r), ncols(A), A.rows[r])
end

function transpose(A::HeapSMat{T}) where {T}
  A.transpose !== nothing && return A.transpose::HeapSMat{T}
  result = transpose!(copy(A))
  A.transpose = result
  result.transpose = A
  return result
end

function transpose!(A::HeapSMat{T}) where {T}
  res_list = HeapSRow{T}[]
  m = nrows(A)
  n = ncols(A)
  pivots = Tuple{Int, T}[pivot(v) for v in A.rows]
  for j in 1:n
    next = Tuple{Int, T}[]
    for (i, (k, c)) in enumerate(pivots)
      k == j || continue
      push!(next, (i, c))
      A[i] = pop_index!(A[i], k)[2]
      A[i] = pop_index!(A[i], k)[2]
      pivots[i] = pivot(A[i])
    end
    push!(res_list, HeapSRow(base_ring(A), next))
  end
  return HeapSMat(base_ring(A), n, m, res_list)
end

function _upper_triangular_form!(A::HeapSMat{ZZRingElem})
  # we iterate through the columns of A
  S = HeapSMat(ZZ, 0, nrows(A))
  H = HeapSMat(ZZ, 0, ncols(A))
  T_list = HeapSRow{ZZRingElem}[]
  for j in 1:ncols(A)
    @show dense_matrix(A)
    @show j
    # find the rows whose pivot is in this column
    cand = [i for i in 1:nrows(A) if Oscar.pivot(A[i])[1] == j]
    isempty(cand) && continue # nothing to do here
    @show cand
    B = A[cand]
    # collect the generators
    g = [c for (_, c) in Oscar.pivot.(A.rows[cand])]
    @show g
    # find the principal generator for the ideal and the Bezout coeffs
    res, coeff = Oscar._gcdx(g)
    @show res
    @show coeff
    coeff_vec = Oscar.HeapSRow(ZZ, collect(zip(cand, coeff)))
    # store the transition coefficients in a base change matrix
    push!(S, coeff_vec)
    # compile the new line which will be in place of the previous ones
    new_line = coeff_vec*A
    # store it in the output matrix
    push!(H, new_line)
    # determine the inverse base change and store it
    c = [divexact(g, res) for g in g]
    @show c
    t_vec = Oscar.HeapSRow(ZZ, collect(zip(cand, c)))
    @show t_vec
    op = tensor_product(t_vec, coeff_vec, nrows(A), nrows(A))
    @show dense_matrix(op)
    @show dense_matrix(unit_matrix(HeapSMat, ZZ, nrows(A)))
    op = unit_matrix(HeapSMat, ZZ, nrows(A)) - op
    @show dense_matrix(op)
    push!(T_list, t_vec)
    # perform the reduction in place on A
    # this will create additional lines for the next step
    AA = copy(A)
    for (k, c) in zip(cand, c)
      A[k] = add!(A[k], -c*new_line)
      @assert iszero(pivot(A[k])[1]) || pivot(A[k])[1] > j
    end
    @assert op*AA == A
    @show dense_matrix(A)
  end
  return H, S, transpose!(HeapSMat(ZZ, nrows(H), nrows(A), T_list))
end

function unit_matrix(::Type{HeapSMat}, R::NCRing, n::Int)
  iszero(n) && return HeapSMat(R, 0, 0)
  res_list = [HeapSRow(R, [(i, one(R))]) for i in 1:n]
  return HeapSMat(R, n, n, res_list)
end

function zero_matrix(::Type{HeapSMat}, R::NCRing, m::Int, n::Int)
  return HeapSMat(R, m, n)
end

function HeapSMat(AA::MatElem{T}) where {T}
  R = base_ring(AA)
  m = nrows(AA)
  n = ncols(AA)
  result = HeapSMat(R, 0, n)
  for i in 1:m
    push!(result, HeapSRow(R, [(i, c) for (i, c) in enumerate(AA[i, :]) if !iszero(c)]))
  end
  return result
end

function HeapSMat(AA::SMat{T}) where {T}
  R = base_ring(AA)
  m = nrows(AA)
  n = ncols(AA)
  result = HeapSMat(R, 0, n)
  for i in 1:m
    push!(result, HeapSRow(AA[i]))
  end
  return result
end

function sparse_matrix(A::HeapSMat)
  R = base_ring(A)
  result = sparse_matrix(R, 0, ncols(A))
  for v in A.rows
    push!(result, sparse_row(v))
  end
  return result
end

@doc raw"""
    upper_triangular_form!(A::MatElem)

Given a matrix ``A ∈ Rᵐˣⁿ`` over a Euclidean domain ``R``, 
return a triple `(B, S, T)` consisting of 

  * an upper-triangular matrix ``B ∈ Rᵖˣⁿ``
  * a matrix ``S ∈ Rᵖˣᵐ`` such that ``B = S⋅A``
  * a matrix ``T ∈ Rᵐˣᵖ`` such that ``A = T⋅B``

In other words: We determine an upper-triangular generating 
set for the submodule of ``Rⁿ`` generated by the rows of ``A``,
together with the base change matrices.
"""
function upper_triangular_form!(A::MatElem)
  error("not implemented")
end

function upper_triangular_form!(A::HeapSMat{ZZRingElem})
  consolidate!(A)
  # we iterate through the columns of A
  S = unit_matrix(HeapSMat, ZZ, nrows(A))
  T_trans = unit_matrix(HeapSMat, ZZ, nrows(A))
  m0 = nrows(A)
  for j in 1:ncols(A)
    # find the rows whose pivot is in this column
    cand = [i for i in 1:m0 if Oscar.pivot(A[i])[1] == j]
    isempty(cand) && continue # nothing to do here
    # collect the generators
    g = [c for (_, c) in Oscar.pivot.(A.rows[cand])]
    # find the principal generator for the ideal and the Bezout coeffs
    res, coeff = Oscar._gcdx(g)
    coeff_vec = Oscar.HeapSRow(ZZ, collect(zip(cand, coeff)))
    # compile the new line which will be added to A
    new_line = consolidate!(coeff_vec*A)
    push!(A, new_line)
    # update the base change matrix
    push!(S, coeff_vec*S)

    push!(T_trans, HeapSRow(ZZ))
    # determine the elimination of the pivots in the original A
    c = [divexact(g, res) for g in g]

    w = consolidate!(copy(S[nrows(S)]))
    for (k, mu) in zip(cand, c)
      A[k] = add!(A[k], -mu*new_line)
      S[k] = add!(S[k], -mu*w)
    end
    w = T_trans[nrows(T_trans)]
    for (k, c) in zip(cand, c)
      w = add!(w, c*T_trans[k])
    end
    T_trans[nrows(T_trans)] = w
  end
  return A[m0+1:nrows(A)], S[m0+1:nrows(S)], transpose(T_trans[m0+1:nrows(T_trans)])
end

function upper_triangular_form!(A::SMat{ZZRingElem})
  # we iterate through the columns of A
  S = unit_matrix(SMat, ZZ, nrows(A))
  T_trans = unit_matrix(SMat, ZZ, nrows(A))
  m0 = nrows(A)
  for j in 1:ncols(A)
    # find the rows whose pivot is in this column
    cand = [i for i in 1:m0 if Oscar.pivot(A[i])[1] == j]
    isempty(cand) && continue # nothing to do here
    # collect the generators
    g = [c for (_, c) in pivot.(A.rows[cand])]
    # find the principal generator for the ideal and the Bezout coeffs
    res, coeff = Oscar._gcdx(g)
    coeff_vec = sparse_row(ZZ, collect(zip(cand, coeff)))
    # compile the new line which will be added to A
    new_line = coeff_vec*A
    push!(A, new_line)
    push!(S, coeff_vec*S)

    push!(T_trans, sparse_row(ZZ))
    c = [divexact(g, res) for g in g]

    w = S[nrows(S)]
    for (k, mu) in zip(cand, c)
      A[k] = A[k] - mu*new_line
      S[k] = S[k] - mu*w
    end
    w = T_trans[nrows(T_trans)]
    for (k, c) in zip(cand, c)
      w = w + c*T_trans[k]
    end
    T_trans[nrows(T_trans)] = w
  end
  return A[m0+1:nrows(A)], S[m0+1:nrows(S)], transpose(T_trans[m0+1:nrows(T_trans)])
end

# An extension of `gcdx` to multiple arguments.
function _gcdx(v::Vector)
  @assert all(!iszero(x) for x in v)
  w = sort(v; by=abs) # If we sort, we have to translate the result back
  res, bez = _gcdx_rec(w)
  n = length(v)
  trans = Int[]
  for i in 1:n
    k = findfirst(w[k] == v[i] for k in 1:n)
    k === nothing && error("permutation could not be determined")
    push!(trans, k)
    w[k] = 0
  end
  return res, bez[trans]
end

function _gcdx_rec(v::Vector)
  @assert !isempty(v)
  if isone(length(v))
    return first(v), [one(first(v))]
  elseif length(v) == 2 
    res, a, b = gcdx(v[1], v[2])
    return res, [a, b]
  end

  k = div(length(v), 2)
  u = v[1:k]
  w = v[k+1:end]
  res1, coeff1 = _gcdx_rec(u)
  res2, coeff2 = _gcdx_rec(w)
  res, a, b = gcdx(res1, res2)
  return res, vcat(a*coeff1, b*coeff2)
end


### Some auxiliary additional functionality for `SMat` to make 
# the above code work.
function getindex(A::SMat, r::UnitRange)
  R = base_ring(A)
  result = sparse_matrix(R, 0, ncols(A))
  for i in r
    push!(result, A[i])
  end
  return result
end

function unit_matrix(::Type{SMat}, R::Ring, n::Int)
  result = sparse_matrix(R, 0, n)
  for i in 1:n
    push!(result, sparse_row(R, [(i, one(R))]))
  end
  return result
end

function pivot(v::SRow)
  isempty(v) && return 0, zero(base_ring(v))
  return first(v.pos), first(v.values)
end


