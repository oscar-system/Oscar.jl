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
      nk = T(1)
      kk = 0
    else
      nk = n - T(k)
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

# Iterator type: all partitions of an integer n
struct Partitions{T<:IntegerUnion}
  n::T

  function Partitions(n::T) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    return new{T}(n)
  end

end

# Iterator type: partitions of n into k parts, with optional lower/upper bounds on the parts
# If distinct_parts == true, then all parts have distinct values.
struct PartitionsFixedNumParts{T<:IntegerUnion}
  n::T
  k::Int

  lb::T
  ub::T
  distinct_parts::Bool

  function PartitionsFixedNumParts(n::T, k::IntegerUnion, lb::IntegerUnion, ub::IntegerUnion, only_distinct_parts::Bool) where T<:IntegerUnion
    @req n >= 0 "n >= 0 required"
    @req k >= 0 "k >= 0 required"
    @req lb >= 0 "lb >=0 required"
    # If lb == 0 the algorithm will actually create lists containing the
    # entry zero, e.g. partitions(2, 2, 0, 2) will contain [2, 0].
    # This is nonsense, so we set lb = 1 in this case.
    if lb == 0
      lb = 1
    end
    return new{T}(n, Int(k), T(lb), T(ub), only_distinct_parts)
  end

end

function PartitionsFixedNumParts(n::T, k::IntegerUnion; only_distinct_parts::Bool = false) where T<:IntegerUnion
  return PartitionsFixedNumParts(n, k, 1, n, only_distinct_parts)
end

function PartitionsFixedNumParts(n::T, k::IntegerUnion, lb::IntegerUnion, ub::IntegerUnion; only_distinct_parts::Bool = false) where T<:IntegerUnion
  return PartitionsFixedNumParts(n, k, lb, ub, only_distinct_parts)
end

# Iterator type: partitions of n into k parts with values in v and every value
# occurring according to the multiplicities in mu
struct PartitionsFixedNumPartsAndValues{T<:IntegerUnion}
  n::T
  k::Int
  v::Vector{T}
  mu::Vector{Int}

  function PartitionsFixedNumPartsAndValues(n::T, k::Int, v::Vector{T}, mu::Vector{Int}) where {T <: IntegerUnion}
    @req n >= 0 "n >= 0 required"
    @req k >= 0 "k >= 0 required"
    @req length(mu) == length(v) "mu and v should have the same length"

    # Algorithm partb in [RJ76] assumes that v is strictly increasing.
    # Added (and noticed) on Mar 22, 2023.
    @req all([v[i] < v[i + 1] for i in 1:length(v) - 1]) "v must be strictly increasing"

    # Parta allows v[1] = 0 but this is nonsense for entries of a partition
    @req all(>(0), v) "Entries of v must be positive"

    # For safety
    @req all(>(0), mu) "Entries of mu must be positive"
    return new{T}(n, k, v, mu)
  end
end

# Internal type: state of the iterator
mutable struct PartitionsFixedNumPartsAndValuesState{T<:IntegerUnion}
  x::Vector{T}
  y::Vector{T}
  ii::Vector{T}
  N::T
  i::Int
  r::Int
  done::Bool

  function PartitionsFixedNumPartsAndValuesState{T}() where {T<:IntegerUnion}
    return new{T}()
  end
end

# Iterator type: partitions of n with values in v
# Optionally, every value occurs according to the multiplicities in mu
struct PartitionsFixedValues{T<:IntegerUnion}
  n::T
  v::Vector{T}
  mu::Vector{Int}
  kmin::Int # minimum number of parts
  kmax::Int # maximum number of parts

  # Constructor without multiplicities
  function PartitionsFixedValues(n::T, v::Vector{T}) where {T <: IntegerUnion}
    @req n >= 0 "n >= 0 required"
    @req all([v[i] < v[i + 1] for i in 1:length(v) - 1]) "v must be strictly increasing"
    @req all(>(0), v) "Entries of v must be positive"

    if n == 0 || isempty(v)
      # We fill the multiplicities with 1s to get around the @req's.
      # If n == 0, the multiplicities aren't checked and if isempty(v), then
      # the multiplicities are empty (as required).
      return new{T}(n, v, ones(Int, length(v)), 0, 0)
    end

    kmin = div(n, v[end])
    kmax = div(n, v[1])
    return new{T}(n, v, fill(kmax, length(v)), kmin, kmax)
  end

  # Constructor with multiplicities
  function PartitionsFixedValues(n::T, v::Vector{T}, mu::Vector{Int}) where {T <: IntegerUnion}
    @req n >= 0 "n >= 0 required"
    @req length(mu) == length(v) "mu and v should have the same length"
    @req all([v[i] < v[i + 1] for i in 1:length(v) - 1]) "v must be strictly increasing"
    @req all(>(0), v) "Entries of v must be positive"
    @req all(>(0), mu) "Entries of mu must be positive"

    if n == 0 || isempty(v)
      # We fill the multiplicities with 1s to get around the @req's.
      # If n == 0, the multiplicities aren't checked and if isempty(v), then
      # the multiplicities are empty (as required).
      return new{T}(n, v, ones(Int, length(v)), 0, 0)
    end

    kmax = 0
    cursum = 0
    for i in 1:length(v), j in 1:mu[i]
      kmax += 1
      cursum += v[i]
      cursum >= n && break
    end

    kmin = 0
    cursum = 0
    for i in length(v):-1:1, j in 1:mu[i]
      kmin += 1
      cursum += v[i]
      cursum >= n && break
    end

    return new{T}(n, v, mu, kmin, kmax)
  end
end

################################################################################
#
#  Multipartition
#
################################################################################

@doc raw"""
    Multipartition{T<:IntegerUnion} <: AbstractVector{Partition{T}}

Multipartitions are implemented as a subtype of 1-dimensional arrays of partitions. You can use smaller integer types to increase performance.

See [`multipartition`](@ref) for the user-facing constructor and an example.
"""
struct Multipartition{T<:IntegerUnion} <: AbstractVector{Partition{T}}
    mp::Vector{Partition{T}}
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

# Iterator type: all semistandard tableaux of a given shape
struct SemiStandardTableaux{T<:IntegerUnion}
  shape::Partition{T}
  max_val::T

  function SemiStandardTableaux(p::Partition{T}, max_val::T) where {T <: IntegerUnion}
    return new{T}(p, max_val)
  end
end

# Iterator type: all semistandard tableaux with a given number of boxes
struct SemiStandardTableauxFixedBoxNum{T<:IntegerUnion}
  box_num::T
  max_val::T

  function SemiStandardTableauxFixedBoxNum(box_num::T, max_val::T) where {T <: IntegerUnion}
    @req box_num >= 0 "Number of boxes must be non-negative"
    return new{T}(box_num, max_val)
  end
end

# Iterator type: all semistandard tableaux with a given shape and weight
struct SemiStandardTableauxFixedShapeAndWeight{T<:IntegerUnion}
  shape::Partition{T}
  weight::Vector{T}

  function SemiStandardTableauxFixedShapeAndWeight(shape::Partition{T}, weight::Vector{T}) where {T <: IntegerUnion}
    @req sum(shape) == sum(weight) "Sum of shape and weight must agree"
    i = findlast(!iszero, weight) # Trim trailing zeros; they upset the iterator
    if isnothing(i)
      i = 0
    end
    return new{T}(shape, weight[1:i])
  end
end

# Internal type: state of the iterator
mutable struct SemiStandardTableauxFixedShapeAndWeightState{T<:IntegerUnion}
  n::Int
  increaseN::Bool
  tab::YoungTableau{T}
  boxes_filled::Vector{Int}
  n_used_weight::Vector{Int}
  row_pointer::Matrix{Int}

  function SemiStandardTableauxFixedShapeAndWeightState{T}() where {T <: IntegerUnion}
    return new{T}()
  end
end

# Iterator type: all standard tableaux of a given shape
struct StandardTableaux{T<:IntegerUnion}
  shape::Partition{T}

  function StandardTableaux(p::Partition{T}) where {T <: IntegerUnion}
    return new{T}(p)
  end
end

# Internal type: state of the iterator
mutable struct StandardTableauxState{T<:IntegerUnion}
  n::Int
  i::Int
  j::Int
  tab::YoungTableau{T}
  sub_s::Vector{Int}
  tracker_row::Vector{Int}
end

# Iterator type: all standard tableaux with a given number of boxes
struct StandardTableauxFixedBoxNum{T<:IntegerUnion}
  box_num::T

  function StandardTableauxFixedBoxNum(box_num::T) where {T <: IntegerUnion}
    @req box_num >= 0 "Number of boxes must be non-negative"
    return new{T}(box_num)
  end
end

################################################################################
#
#  Combination(s)
#
################################################################################

struct Combination{T} <: AbstractVector{T}
  v::Vector{T}
end

# Iterator type: all combinations of k elements from the vector v
struct Combinations{T, U<:IntegerUnion}
  v::T
  n::U
  k::U

  inplace::Bool # Whether all generated combinations share the same array in
                # memory

  function Combinations(v::T, n::U, k::U, inplace::Bool = false) where {T, U<:IntegerUnion}
    return new{T,U}(v, n, k, inplace)
  end
end


Combinations(v::AbstractArray, k::T) where {T<:IntegerUnion} = Combinations(v, T(length(v)), k)

################################################################################
#
#  Multicombination(s)
#
################################################################################

# Iterator type: all combinations of k elements from the vector v with repetition
struct MultiCombinations{T, U<:IntegerUnion}
  v::T
  n::U
  k::U

  inplace::Bool # Whether all generated combinations share the same array in
                # memory

  function MultiCombinations(v::T, n::U, k::U, inplace::Bool = false) where {T, U<:IntegerUnion}
    return new{T,U}(v, n, k, inplace)
  end
end

MultiCombinations(v::AbstractArray, k::T) where {T<:IntegerUnion} = MultiCombinations(v, T(length(v)), k)
