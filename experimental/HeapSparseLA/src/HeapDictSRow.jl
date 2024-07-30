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
  return HeapDictSRow(R, [i for (i, _) in a], [v for (_, v) in a])
end

base_ring(v::HeapDictSRow) = v.R

function getindex(v::HeapDictSRow{I, T}, i::I) where {I, T}
  !haskey(v.coeffs, i) && return zero(base_ring(v))
  return v.coeffs[i]
end

function setindex!(v::HeapDictSRow{I, T}, c::T, i::I) where {I, T}
  !haskey(v.coeffs, i) && push!(v.indices, i)
  v.coeffs[i] = c
  return c
end

function Base.pop!(v::HeapDictSRow)
  i = pop!(v.indices)
  c = pop!(v.coeffs, i)
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

function add!(v::HeapDictSRow{I, T}, w::HeapDictSRow{I, T}) where {I, T}
  # TODO: use sizehint to avoid incremental allocations
  h = v.indices
  sizehint!(h.content, length(h.content) + length(w.indices.content))
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
  end
  # TODO: Can we indicate that there are no duplicate keys and get a speedup?
  v.coeffs = merge!(v.coeffs, w.coeffs)
  v.indices.content = vcat(v.indices.content, new_ind)
  v.indices.is_heapified = false
  heapify!(v.indices)
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
  for (i, c) in v.coeffs
    v.coeffs[i] *= a
  end
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

