################################################################################
#
#  WeakComposition
#
################################################################################

@doc raw"""
    WeakComposition{T<:IntegerUnion} <: AbstractVector{T}

A weak composition consisting of integers of type `T`.
This is a wrapper around `Vector{T}`.

See [`weak_composition`](@ref) for the user-facing constructor and an example.
"""
struct WeakComposition{T<:IntegerUnion} <: AbstractVector{T}
  c::Vector{T}
end

# Iterator type: all weak compositions of n into k parts
struct WeakCompositions{T<:IntegerUnion}
  n::T
  k::Int # This is the length of the vectors, so cannot be larger than an Int
         # and there is no performance benefit from making it an Int8 (for example)
  inplace::Bool # Whether all generated compositions share the same array in
                # memory

  function WeakCompositions(n::T, k::Int, inplace::Bool = false) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    @req k >= 0 "k >= 0 required"
    return new{T}(n, k, inplace)
  end
end

################################################################################
#
#  Composition
#
################################################################################

@doc raw"""
    Composition{T<:IntegerUnion} <: AbstractVector{T}

A composition consisting of integers of type `T`.
This is a wrapper around `Vector{T}`.

See [`composition`](@ref) for the user-facing constructor and an example.
"""
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
#  Ascending composition
#
################################################################################

# Iterator type: all ascending compositions of n
struct AscendingCompositions{T<:IntegerUnion}
  n::T

  function AscendingCompositions(n::T) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    return new{T}(n)
  end
end

# Internal type: state of the iterator
mutable struct AscendingCompositionsState{T<:IntegerUnion}
  a::Vector{T}
  x::T
  y::T
  k::Int

  function AscendingCompositionsState{T}() where {T<:IntegerUnion}
    return new{T}()
  end
end

################################################################################
#
#  Partition
#
################################################################################

@doc raw"""
    Partition{T<:IntegerUnion} <: AbstractVector{T}

A partition consisting of integers of type `T`.
This is a wrapper around `Vector{T}`.

See [`partition`](@ref) for the user-facing constructor and an example.
"""
struct Partition{T<:IntegerUnion} <: AbstractVector{T}
  p::Vector{T}
end

################################################################################
#
#  Young Tableaux
#
################################################################################

@doc raw"""
    YoungTableau{T<:IntegerUnion} <: AbstractVector{AbstractVector{T}}

A Young tableau filled with integers of type `T`.
This is a wrapper around `Vector{Vector{T}}`.

See [`young_tableau`](@ref) for the user-facing constructor and an example.
"""
struct YoungTableau{T<:IntegerUnion} <: AbstractVector{AbstractVector{T}}
  t::Vector{Vector{T}}
end
