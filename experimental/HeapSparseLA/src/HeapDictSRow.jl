# Binary heaps

# comparison for ordering in the heap; to be overwritten eventually
_binary_heap_leq(a::Any, b::Any) = a <= b

function length(h::BinaryHeap)
  return length(h.content)
end

function Base.push!(h::BinaryHeap{T}, v::T) where {T}
  !h.is_heapified && heapify!(h)
  n = length(h)
  push!(h.content, v) # just to make space
  i = n + 1
  p = i >> 1
  while i > 1 && !_binary_heap_leq(h.content[p], v)
    h.content[i] = h.content[p]
    i = p
    p = i >> 1
  end
  i <= n && (h.content[i] = v)
  return h
end

function Base.pop!(h::BinaryHeap{T}) where T
  !h.is_heapified && heapify!(h)
  # move the top node down
  isempty(h) && return nothing
  result = first(h.content)
  i = 1
  l = i << 1
  r = l + 1
  n = length(h)
  while i <= n
    l > n && break # we're in a leaf node
    if r > n # the left child is the only leaf node
      h.content[i] = h.content[l]
      i = l
      break
    end

    a = h.content[l]
    b = h.content[r]
    if _binary_heap_leq(a, b)
      h.content[i] = a
      i = l
      l = i << 1
      r = l + 1
    else
      h.content[i] = b
      i = r
      l = i << 1
      r = l + 1
    end
  end

  # move the last element to the gap at i and then up
  i0 = i
  v = pop!(h.content) # the last element
  i == n && return result # nothing to do in this case
  p = i >> 1
  while i > 1 && !_binary_heap_leq(h.content[p], v)
    h.content[i] = h.content[p]
    i = p
    p = i >> 1
  end
  !isempty(h.content) && (h.content[i] = v)
  return result
end

function peek(h::BinaryHeap)
  !h.is_heapified && heapify!(h)
  return first(h.content)
end

function heapify!(h::BinaryHeap{T}) where {T}
  for q in length(h.content):-1:1
    i = q
    p = i >> 1
    v = h.content[q]
    while i > 1 && !_binary_heap_leq(h.content[p], v)
      h.content[i] = h.content[p]
      i = p
      p = i >> 1
    end
    i < q && (h.content[i] = v)
  end
  h.is_heapified = true
  return h
end

isempty(h::BinaryHeap) = isempty(h.content)

function copy(h::BinaryHeap)
  return BinaryHeap(h.content; is_heapified = h.is_heapified)
end

### sparse rows based on heaps and dictionaries
function HeapDictSRow(
    R::NCRing, a::Vector{Tuple{IndexType, T}}
  ) where {IndexType, T}
  a = filter!(x->!iszero(x[2]), a)
  return HeapDictSRow(R, [i for (i, _) in a], [v for (_, v) in a])
end

base_ring(v::HeapDictSRow) = v.R

function getindex(v::HeapDictSRow{I, T}, i::I) where {I, T}
  @hassert :SparseLA 3 is_sane(v)
  !haskey(v.coeffs, i) && return zero(base_ring(v))
  return v.coeffs[i]
end

function setindex!(v::HeapDictSRow{I, T}, c::T, i::I) where {I, T}
  @hassert :SparseLA 3 is_sane(v)
  !haskey(v.coeffs, i) && push!(v.indices, i)
  v.coeffs[i] = c
  @hassert :SparseLA 3 is_sane(v)
  return c
end

function Base.pop!(v::HeapDictSRow)
  @hassert :SparseLA 3 is_sane(v)
  i = pop!(v.indices)
  c = pop!(v.coeffs, i)
  @hassert :SparseLA 3 is_sane(v)
  return i, c
end

function isempty(v::HeapDictSRow)
  return isempty(v.indices)
end

function length(v::HeapDictSRow)
  return length(v.indices)
end

function copy(v::HeapDictSRow)
  return HeapDictSRow(v)
end

function Base.delete!(h::BinaryHeap{I}, k::I) where {I}
  isempty(h) && return h
  # try to find k
  i = 1
  i = find_node_rec(h.content, i, k)
  i === nothing && error("key not found")

  # move the element down
  l = i << 1
  r = l + 1
  n = length(h)
  while i <= n
    l > n && break # we're in a leaf node
    if r > n # the left child is the only leaf node
      h.content[i] = h.content[l]
      i = l
      break
    end

    a = h.content[l]
    b = h.content[r]
    if _binary_heap_leq(a, b)
      h.content[i] = a
      i = l
      l = i << 1
      r = l + 1
    else
      h.content[i] = b
      i = r
      l = i << 1
      r = l + 1
    end
  end

  # move the last element to the gap at i and then up
  n = length(h.content)
  v = pop!(h.content) # the last element
  i == n && return h # nothing to do in this case
  p = i >> 1
  while i > 1 && !_binary_heap_leq(h.content[p], v)
    h.content[i] = h.content[p]
    i = p
    p = i >> 1
  end
  !isempty(h.content) && (h.content[i] = v)
  return h
end

function find_node_rec(c::Vector{I}, i::Int, k::I) where {I}
  k == c[i] && return i
  l = i << 1
  r = l + 1
  n = length(c)

  if l <= n && _heap_leq(c[l], k)
    l_res = find_node_rec(c, l, k)
    l_res !== nothing && return l_res
  end

  if r <= n && _heap_leq(c[r], k)
    r_res = find_node_rec(c, r, k)
    r_res !== nothing && return r_res
  end
  return nothing
end



function add!(v::HeapDictSRow{I, T}, w::HeapDictSRow{I, T}) where {I, T}
  @hassert :SparseLA 3 (is_sane(v) && is_sane(w))
  # TODO: use sizehint to avoid incremental allocations
  m = length(w)
  n = length(v)
  new_ind = sizehint!(Vector{I}(), m)
  old_ind = sizehint!(Vector{I}(), m)
  for k in keys(w.coeffs)
    if k in keys(v.coeffs)
      push!(old_ind, k)
    else
      push!(new_ind, k)
    end
  end
  for i in old_ind
    v.coeffs[i] += pop!(w.coeffs, i)
    if iszero(v.coeffs[i]) 
      delete!(v.coeffs, i)
      delete!(v.indices, i)
    end
  end
  # TODO: Can we indicate that there are no duplicate keys and get a speedup?
  v.coeffs = merge!(v.coeffs, w.coeffs)
  # TODO: can we use array-copy functions here?
  sizehint!(v.indices.content, m + length(new_ind))
  for j in new_ind
    push!(v.indices.content, j)
  end
  #v.indices.content = vcat(v.indices.content, new_ind)
  v.indices.is_heapified = false
  heapify!(v.indices)
  @hassert :SparseLA 3 is_sane(v)
  return v
end

function +(v::HeapDictSRow{I, T}, w::HeapDictSRow{I, T}) where {I, T}
  return add!(copy(v), copy(w))
end

*(a::Any, v::HeapDictSRow) = base_ring(a)*v

function *(a::T, v::HeapDictSRow{I, T}) where {I, T}
  return mul!(a, copy(v))
end

function mul!(a::T, v::HeapDictSRow{I, T}) where {I, T}
  @hassert :SparseLA 3 is_sane(v)
  for (i, c) in v.coeffs
    v.coeffs[i] *= a
  end
  @hassert :SparseLA 3 is_sane(v)
  return v
end

-(v::HeapDictSRow) = -one(base_ring(v))*v

function Base.show(io::IO, v::HeapDictSRow)
  h = copy(v.indices)
  while !isempty(h)
    i = pop!(h)
    println(io, "$i => $(v.coeffs[i])")
  end
end

function sparse_row(v::HeapDictSRow)
  return sparse_row(base_ring(v), [(i, c) for (i, c) in v.coeffs])
end

function HeapDictSRow{I, T}(R::NCRing, a::Vector{Tuple{I, T}}) where {I, T}
  return HeapDictSRow(R, a)
end

function pivot(v::HeapDictSRow{I, T}) where {I, T}
  @hassert :SparseLA 3 is_sane(v)
  isempty(v) && return zero(I), zero(base_ring(v))
  i = peek(v.indices)
  return i, v.coeffs[i]
end

function *(v::HeapDictSRow{I, ET}, A::HeapSMat{ET, HeapDictSRow{I, ET}}) where {I, ET}
  return mul!(v, A)
end

function mul!(v::HeapDictSRow{I, ET}, A::HeapSMat{ET, HeapDictSRow{I, ET}}) where {I, ET}
  @hassert :SparseLA 3 is_sane(v)
  @hassert :SparseLA 3 all(is_sane(w) for w in A.rows)
  R = base_ring(v)
  @assert R === base_ring(A)
  m = nrows(A)
  n = ncols(A)
  result = typeof(v)(R)
  return sum(c*A[i] for (i, c) in v.coeffs; init=result)
end

consolidate!(v::HeapDictSRow) = v

function transpose!(A::HeapSMat{ET, T}) where {ET, T<:HeapDictSRow{Int}}
  @hassert :SparseLA 3 all(is_sane(w) for w in A.rows)
  m = nrows(A)
  n = ncols(A)
  R = base_ring(A)
  result = typeof(A)(R, 0, m)

  res_list = T[]
  m = nrows(A)
  n = ncols(A)
  for j in 1:n
    next = Tuple{Int, ET}[]
    for i in 1:m
      if j in keys(A[i].coeffs)
        push!(next, (i, A[i, j]))
      end
    end
    new_line = T(R, next)
    push!(result, new_line)
  end
  @hassert :SparseLA 3 all(is_sane(w) for w in A.rows)
  @hassert :SparseLA 3 all(is_sane(w) for w in result.rows)
  return result
end

function ==(v::T, w::T) where {T <: HeapDictSRow}
  return v.coeffs == w.coeffs
end

as_dictionary(v::HeapDictSRow) = copy(v.coeffs)

function is_sane(v::HeapDictSRow)
  all(x in keys(v.coeffs) for x in v.indices.content) || error("more indices than coefficients")
  all(x in v.indices.content for x in keys(v.coeffs)) || error("more coefficients than indices")
  length(v.indices.content) == length(v.coeffs) || error("double indices")
  any(iszero(x) for x in values(v.coeffs)) && error("zero value found")
  return true
end

