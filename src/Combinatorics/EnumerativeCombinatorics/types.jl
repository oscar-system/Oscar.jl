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
#  Composition
#
################################################################################

struct Composition{T<:IntegerUnion} <: AbstractVector{T}
  c::Vector{T}
end

# Iterator type: all compositions of n
struct Compositions{T<:IntegerUnion}
  n::T

  function Compositions(n::T) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    return new{T}(n)
  end
end

# Iterator type: all compositions of n into k parts
struct CompositionsFixedNumParts{T<:IntegerUnion}
  n::T
  k::Int # This is the length of the vectors, so cannot be larger than an Int
         # and there is no performance benefit from making it an Int8 (for example)
  weak_comp_iter::WeakCompositions{T}

  function CompositionsFixedNumParts(n::T, k::Int) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    @req k >= 0 "k >= 0 required"

    # We get a composition of n into k parts from a weak composition of n - k
    # into k parts by adding 1 to every entry, so we may reuse the iterator
    # for weak compositions here. If k > n, so n - k < 0, there are no
    # compositions, but we have to cheat a bit to get an empty iterator
    # for the weak compositions.
    if k > n
      # 1 does not have any weak compositions into 0 parts, so this will
      # produce an empty iterator
      nk = 1
      kk = 0
    else
      nk = n - k
      kk = k
    end
    return new{T}(n, k, weak_compositions(nk, kk))
  end
end

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
#  Young Tableaux
#
################################################################################

struct YoungTableau{T} <: AbstractVector{AbstractVector{T}}
  t::Vector{Vector{T}}
end
