################################################################################
#
#  Partition
#
################################################################################

struct Partition{T<:IntegerUnion} <: AbstractVector{T}
  p::Vector{T}
end

################################################################################
#
#  WeakComposition
#
################################################################################

struct WeakComposition{T<:IntegerUnion} <: AbstractVector{T}
  c::Vector{T}
end

# Iterator type: all weak compositions of n into k parts
struct WeakCompositions{T<:IntegerUnion}
  n::T
  k::Int # This is the length of the vectors, so cannot be larger than an Int
         # and there is no performance benefit from making it an Int8 (for example)

  function WeakCompositions(n::T, k::Int) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    @req k >= 0 "k >= 0 required"
    return new{T}(n, k)
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
