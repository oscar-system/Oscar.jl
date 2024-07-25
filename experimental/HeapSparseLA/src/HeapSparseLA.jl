# A heap based implementation of sparse row vectors for elements
# which are big in the sense that copying them around in memory 
# for sorting is expensive.
#=
mutable struct HeapNode{T}
  content::T
  index::Int
  left::Union{HeapNode{T}, Nothing}
  right::Union{HeapNode{T}, Nothing}

  function HeapNode(i::Int, cont::T; 
        left::Union{HeapNode{T}, Nothing}=nothing, 
        right::Union{HeapNode{T}, Nothing}=nothing
      ) where {T}
    return new{T}(cont, i, left, right)
  end
end

content(n::HeapNode) = n.content
index(n::HeapNode) = n.i
has_left_child(n::HeapNode) = (n.left !== nothing)
has_right_child(n::HeapNode) = (n.right !== nothing)
is_leaf(n::HeapNode) = !has_left_child(n) && !has_right_child(n)
left_child(n::HeapNode{T}) where {T} = n.left::T
right_child(n::HeapNode{T}) where {T} = n.right::T


mutable struct HeapSRow{T}
  entry_parent::NCRing
  top::Union{Nothing, HeapNode{T}}
  is_consolidated::Bool

  function HeapSRow(R::NCRing)
    return new{elem_type(R)}(R, nothing, true)
  end

  function HeapSRow(R::NCRing, n::HeapNode{T};
      is_consolidated::Bool=false
    ) where {T}
    return new{T}(R, n, is_consolidated)
  end
end

top(v::HeapSRow{T}) where {T} = v.top::HeapNode{T}
isempty(v::HeapSRow) = v.top === nothing

function is_consolidated(v::HeapSRow) 
  v.is_consolidated && return true
  if is_empty(v)
    v.is_consolidated = true
    return true
  end

  n = top(v)
  result, _ = _has_duplicate(n, Int[])
  if !result
    v.result = true
    return true
  end
  return false
end

function _has_duplicate(n::HeapNode, key_list::Vector{Int})
  i = index(n)
  k = searchsortedfirst(key_list, i, rev=true)
  if key_list[k] != i
    push!(key_list, i)
    if has_left_child(n)
      result, key_list = _has_duplicate(left_child(n), key_list)
      result && return true, key_list
    end
    if has_right_child(n)
      result, key_list = _has_duplicate(right_child(n), key_list)
      result && return true, key_list
    end
    return false, key_list
  end

  push!(key_list, i)
  return true, key_list
end

function getindex(v::HeapSRow{T}, i::Int)
  #TODO: Also clean up the heap if there is more than one value for the key i.
  all_indices = gather_values(top(v), i, is_consolidated(v))
  return sum(all_indices; init=zero(v.entry_parent))
end

function gather_values(n::HeapNode{T}, i::Int, is_consolidated::Bool) where {T}
  result = T[]
  i > index(n) && return result
  if i == index(n) 
    push!(result, content(n))
    is_consolidated && return result # There can only be one entry for i and it is this one
  end
  has_left_child(n) && (result = vcat(result, gather_values(left_child(n), i)))
  has_right_child(n) && (result = vcat(result, gather_values(right_child(n), i)))
  return result
end

function has_index(v::HeapSRow{T}, i::Int)
  isempty(v) && return false
  return _has_index(top(v), i)
end

function _has_index(n::HeapNode, i::Int)
  index(n) == i && return true
  has_left_child(n) && _has_index(left_child(n), i) && return true
  has_right_child(n) && _has_index(right_child(n), i) && return true
  return false
end

=#
