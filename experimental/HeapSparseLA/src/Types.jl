
mutable struct HeapNode{T}
  content::T
  index::Int
  cofactor::Union{T, Nothing} # an element !==nothing here indicates
                              # that the contents here and all of the 
                              # following subtree shall be multiplied 
                              # with `cofactor` to get the actual values.
  left::Union{HeapNode{T}, Nothing}
  right::Union{HeapNode{T}, Nothing}

  function HeapNode(i::Int, cont::T; 
        cofactor::Union{T, Nothing}=nothing,
        left::Union{HeapNode{T}, Nothing}=nothing, 
        right::Union{HeapNode{T}, Nothing}=nothing
      ) where {T}
    return new{T}(cont, i, cofactor, left, right)
  end
end

mutable struct HeapSRow{T}
  entry_parent::NCRing
  top::Union{Nothing, HeapNode{T}}
  top_is_consolidated::Bool
  is_consolidated::Bool
  index_cache::Dict{Int, T}

  function HeapSRow(R::NCRing)
    return new{elem_type(R)}(R, nothing, true, true)
  end

  function HeapSRow(R::NCRing, n::HeapNode{T};
      top_is_consolidated::Bool=false,
      is_consolidated::Bool=false
    ) where {T}
    return new{T}(R, n, top_is_consolidated, is_consolidated)
  end
end

mutable struct HeapSMat{T}
  base_ring::NCRing
  nrows::Int
  ncols::Int
  rows::Vector{HeapSRow{T}}
  transpose::Union{HeapSMat{T}, Nothing}

  function HeapSMat(R::NCRing, nrows::Int, ncols::Int, rows::Vector{HeapSRow{T}}) where {T}
    @assert all(base_ring(v) === R for v in rows)
    @assert nrows >= length(rows)
    crows = copy(rows)
    for i in length(rows)+1:nrows
      push!(crows, HeapSRow(R))
    end
    return new{T}(R, nrows, ncols, crows, nothing)
  end

  function HeapSMat(R::Ring, nrows::Int, ncols::Int)
    rows = HeapSRow{elem_type(R)}[HeapSRow(R) for i in 1:nrows]
    return new{elem_type(R)}(R, nrows, ncols, rows, nothing)
  end
end


