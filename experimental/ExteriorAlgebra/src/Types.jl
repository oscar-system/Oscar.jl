mutable struct ExteriorAlgebra{T} <: NCRing
  base_ring::NCRing
  symbols::Vector{Symbol}
  rank::Int

  print_symbols::Dict{Int, Vector{Symbol}}
  multiplication_hash_tables::Dict{Tuple{Int, Int}, Matrix{Tuple{Int,Int}}}

  function ExteriorAlgebra(R::NCRing, a::Vector{Symbol})
    return new{elem_type(R)}(R, a, length(a))
  end
end

mutable struct ExtAlgElem{T} <: NCRingElem
  parent::ExteriorAlgebra{T}
  components::Dict{Int, SRow{T}}

  function ExtAlgElem(E::ExteriorAlgebra{T}) where {T}
    return new{T}(E)
  end

  function ExtAlgElem(E::ExteriorAlgebra{T}, comp::Dict{Int, SRow{T}}; 
      check::Bool=true
    ) where {T}
    @assert all(base_ring(x) === base_ring(E) for (_, x) in comp)
    @check all(!iszero(x) for (_, x) in comp)
    return new{T}(E, comp)
  end
end


