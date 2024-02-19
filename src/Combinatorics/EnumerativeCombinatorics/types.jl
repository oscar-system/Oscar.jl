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
end

################################################################################
#
#  Young Tableaux
#
################################################################################

struct YoungTableau{T} <: AbstractVector{AbstractVector{T}}
  t::Vector{Vector{T}}
end
