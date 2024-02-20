################################################################################
#
#  Partition
#
################################################################################

struct Partition{T<:IntegerUnion} <: AbstractVector{T}
  p::Vector{T}
end

struct PartitionSet{T<:IntegerUnion}
  n::T

  function PartitionSet(n::T) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    return new{T}(n)
  end
  
end

################################################################################
#
#  Young Tableaux
#
################################################################################

struct YoungTableau{T} <: AbstractVector{AbstractVector{T}}
  t::Vector{Vector{T}}
end
