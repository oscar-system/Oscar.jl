################################################################################
# Partitions of an integer.
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
# by Ulrich Thiel and OSCAR-ified by Claudia He Yun and Matthias Zach.
#
# Mar 2023: Improvements and bug fixes by UT.
################################################################################

################################################################################
#
#  Constructors and basic functionality
#
################################################################################

function _defines_partition(parts::Vector{<: IntegerUnion})
  if isempty(parts)
    return true
  end
  return issorted(parts; rev=true) && is_positive(parts[end])
end

@doc raw"""
    partition([T::Type{<:IntegerUnion}], parts::IntegerUnion...; check::Bool = true)
    partition(parts::Vector{T}; check::Bool = true) where T <: IntegerUnion

Return the partition given by the integer sequence `parts` as an object of type
`Partition{T}`.

The element type `T` may be optionally specified, see also the examples below.

If `check` is `true` (default), it is checked whether the given sequence defines
a partition.

# Examples
```jldoctest
julia> P = partition([6, 4, 4, 2]) # the partition 6 + 4 + 4 + 2 of 16
[6, 4, 4, 2]

julia> P = partition(6, 4, 4, 2) # the same partition
[6, 4, 4, 2]

julia> P = partition(Int8, 6, 4, 4, 2) # save the elements in 8-bit integers
Int8[6, 4, 4, 2]
```
"""
partition

function partition(parts::IntegerUnion...; check::Bool = true)
  return partition(Int, parts..., check = check)
end

function partition(T::Type{<:IntegerUnion}, parts::IntegerUnion...; check::Bool = true)
  return partition(collect(T, parts), check = check)
end

function partition(parts::Vector{T}; check::Bool = true) where {T <: IntegerUnion}
  if check
    @req _defines_partition(parts) "The sequence does not define a partition"
  end
  return Partition{T}(parts)
end

# The empty array is of "Any" type, and this is not what we want.
# We want it to be an array of integers of the default type Int64.
function partition(p::Vector{Any}; check::Bool = true)
  return partition(Vector{Int}(p), check = check)
end

data(P::Partition) = P.p

function Base.show(io::IO, ::MIME"text/plain", P::Partition)
  p = data(P)
  if isempty(p)
    print(io, "Empty partition")
    return
  end
  print(io, p)
end

################################################################################
#
#  Array-like functionality
#
################################################################################

function Base.size(P::Partition)
  return size(data(P))
end

function Base.length(P::Partition)
  return length(data(P))
end

function Base.getindex(P::Partition, i::IntegerUnion)
  return getindex(data(P), Int(i))
end

function Base.setindex!(P::Partition, x::IntegerUnion, i::IntegerUnion)
  return setindex!(data(P), x, Int(i))
end

function Base.copy(P::Partition)
  return partition(copy(data(P)))
end

@doc raw"""
    getindex_safe(P::Partition, i::IntegerUnion)

Return `P[i]` if `i < length(P)` and `0` otherwise.
It is assumed that `i` is positive.

# Examples
```jldoctest
julia> P = partition([3, 2, 1])
[3, 2, 1]

julia> getindex_safe(P, 3)
1

julia> getindex_safe(P, 4)
0
```
"""
function getindex_safe(P::Partition{T}, i::IntegerUnion) where T
  return (i > length(data(P)) ? zero(T) : getindex(data(P), Int(i)))
end

base(P::Partitions) = P.n

Base.eltype(::Type{Partitions{T}}) where T = Partition{T}

function Base.show(io::IO, ::MIME"text/plain", P::Partitions)
  print(pretty(io), "Iterator over the partitions of $(base(P))")
end

Base.length(P::Partitions) = BigInt(number_of_partitions(P.n))



base(P::PartitionsFixedNumParts) = P.n

Base.eltype(::Type{PartitionsFixedNumParts{T}}) where T = Partition{T}

function Base.show(io::IO, ::MIME"text/plain", P::PartitionsFixedNumParts)
  print(pretty(io), "Iterator over the partitions of $(base(P)) into ",
  ItemQuantity(P.k, "part"))
end

# NOTE this will not be accurate in many cases,
# in particular if upper/lower bounds are given,
# or if `only_distinct_parts == true`.
# Base.length(P::PartitionsFixedNumParts) = BigInt(number_of_partitions(P.n, P.k))

Base.IteratorSize(::Type{PartitionsFixedNumParts{T}}) where T = Base.SizeUnknown()
################################################################################
#
# Generating and counting unrestricted partitions
#
################################################################################

@doc raw"""
    number_of_partitions(n::IntegerUnion)

Return the number of integer partitions of `n`. For `n < 0`, return `0`.

# Examples
```jldoctest
julia> number_of_partitions(1000)
24061467864032622473692149727991
```
"""
function number_of_partitions(n::IntegerUnion)
  if n < 0
    return ZZ(0)
  end
  # This function should always return a ZZRingElem, see the discussion here:
  # https://github.com/oscar-system/Oscar.jl/pull/3159 .
  # We hence "overwrite" the Nemo function number_of_partitions(::Int) which
  # returns an Int (or throws an InexactError).
  return Nemo.number_of_partitions(ZZ(n))
end


@doc raw"""
    partitions(n::IntegerUnion)

Return an iterator over all partitions of a non-negative integer `n`, produced
in lexicographically *descending* order.
Using a smaller integer type for `n` (e.g. `Int8`) may increase performance.

The algorithm used is "Algorithm ZS1" by [ZS98](@cite). This algorithm is also
discussed in [Knu11](@cite), Algorithm P (page 392).

# Examples
```jldoctest
julia> p = partitions(4);

julia> first(p)
[4]

julia> collect(p)
5-element Vector{Partition{Int64}}:
 [4]
 [3, 1]
 [2, 2]
 [2, 1, 1]
 [1, 1, 1, 1]

julia> collect(partitions(Int8(4))) # using less memory
5-element Vector{Partition{Int8}}:
 Int8[4]
 Int8[3, 1]
 Int8[2, 2]
 Int8[2, 1, 1]
 Int8[1, 1, 1, 1]
```
"""
function partitions(n::IntegerUnion)
  return Partitions(n)
end


function Base.iterate(P::Partitions{T}) where T
  n = base(P)

  if n == 0
    return partition(T[], check=false), (T[], 0, 0)
  elseif P.n == 1
    return partition(T[1], check=false), (T[1], 1, 0)
  end

  d = fill(T(1), Int(n))
  d[1] = n
  return partition(d[1:1], check=false), (d, 1, 1)
end

@inline function Base.iterate(P::Partitions{T}, state::Tuple{Vector{T}, Int, Int}) where T
  d, k, q = state
  q==0 && return nothing
  if d[q] == 2
    d[q] = 1
    k += 1
    q -= 1
  else
    m = d[q] - 1
    np = k - q +1
    d[q] = m
    while np >= m
      q += 1
      d[q] = m
      np -= m
    end
    if np == 0
      k = q
    else
      k = q + 1
      if np > 1
        q += 1
        d[q] = np
      end
    end
  end
  return partition(d[1:k], check=false), (d, k, q)
end

################################################################################
#
# Generating and counting restricted partitions
#
################################################################################

@doc raw"""
    number_of_partitions(n::IntegerUnion, k::IntegerUnion)

Return the number of integer partitions of the non-negative integer `n` into
`k >= 0` parts.
If `n < 0` or `k < 0`, return `0`.
"""
function number_of_partitions(n::IntegerUnion, k::IntegerUnion)
  if n < 0 || k < 0
    return ZZ(0)
  end

  n = ZZ(n)
  k = ZZ(k)

  # Special cases
  if n == k
    return ZZ(1)
  elseif n < k || k == 0
    return ZZ(0)
  elseif k == 1
    return ZZ(1)

  # See https://oeis.org/A008284
  elseif n < 2*k
    return number_of_partitions(n - k) #n - k >= 0 holds since the case n<k was already handled

  # See https://oeis.org/A008284
  elseif n <= 2 + 3*k
    p = number_of_partitions(n - k) #n - k >= 0 holds since the case n<k was already handled
    for i = 0:Int(n) - 2*Int(k) - 1
      p = p - number_of_partitions(ZZ(i))
    end
    return p

  # Otherwise, use recurrence.
  # The following is taken from the GAP code in lib/combinat.gi
  # It uses the standard recurrence relation but in a more intelligent
  # way without recursion.
  else
    n = Int(n)
    k = Int(k)
    p = fill( ZZ(1), n )
    for l = 2:k
      for k = l + 1:n - l + 1
        p[k] = p[k] + p[k - l]
      end
    end
    return p[n - k + 1]
  end

end

@doc raw"""
    partitions(n::IntegerUnion, k::IntegerUnion; only_distinct_parts::Bool = false)
    partitions(n::IntegerUnion, k::IntegerUnion, lb::IntegerUnion, ub::IntegerUnion; only_distinct_parts::Bool = false)

Return an iterator over all partitions of a non-negative integer `n` into
`k >= 0` parts. Optionally, a lower bound `lb >= 0` and an upper bound `ub` for
the parts can be supplied. In this case, the partitions are produced in
*decreasing* order.

There are two choices for the parameter `only_distinct_parts`:
* `false`: no further restriction (*default*);
* `true`: only distinct parts.

The implemented algorithm is "parta" in [RJ76](@cite).

# Examples
All partitions of 7 into 3 parts:
```jldoctest
julia> collect(partitions(7, 3))
4-element Vector{Partition{Int64}}:
 [5, 1, 1]
 [4, 2, 1]
 [3, 3, 1]
 [3, 2, 2]
```
All partitions of 7 into 3 parts where all parts are between 1 and 4:
```jldoctest
julia> collect(partitions(7, 3, 1, 4))
3-element Vector{Partition{Int64}}:
 [4, 2, 1]
 [3, 3, 1]
 [3, 2, 2]
```
Same as above but requiring all parts to be distinct:
```jldoctest
julia> collect(partitions(7, 3, 1, 4; only_distinct_parts = true))
1-element Vector{Partition{Int64}}:
 [4, 2, 1]
```
"""
function partitions(n::IntegerUnion, k::IntegerUnion, lb::IntegerUnion, ub::IntegerUnion; only_distinct_parts::Bool = false)
  return PartitionsFixedNumParts(n, k, lb, ub; only_distinct_parts = only_distinct_parts)
end

function partitions(n::IntegerUnion, k::IntegerUnion; only_distinct_parts::Bool = false)
  return PartitionsFixedNumParts(n, k; only_distinct_parts = only_distinct_parts)
end

# Algorithm "parta" in [RJ76](@cite), de-gotoed from old ALGOL 60 code by E. Thiel.

# Note that the algorithm is given as partitioning m into n parts,
# but we have refactored to align more closely with standard terminology.
function Base.iterate(P::PartitionsFixedNumParts{T}) where T
  n = P.n
  k = P.k
  lb = P.lb
  ub = P.ub
  only_distinct_parts = P.distinct_parts

  if n == 0 && k == 0
    return partition(T[], check=false), (T[], T[], T(0), T(0), 1, false)
  end

  # This iterator should be empty
  if k == 0 || k > n || ub < lb
    return nothing
  end

  if n == k && lb == 1
    only_distinct_parts && k > 1 && return nothing
    return partition(T[1 for i in 1:n], check=false), (T[], T[], T(0), T(0), 1, false)
  end

  if k == 1 && lb <= n <= ub
    return partition(T[n], check=false), (T[], T[], T(0), T(0), 1, false)
  end

  x = zeros(T,k)
  y = zeros(T,k)
  jj = only_distinct_parts*k*(k-1)
  N = T(n - k*lb - div(jj,2))
  L2 = ub-lb
  0 <= N <= k*L2 - jj || return nothing

  for i in 1:k
    y[i] = x[i] = lb + only_distinct_parts*(k-i)
  end

  i = 1
  L2 = L2 - only_distinct_parts*T(k-1)

  while N > L2
    N -= L2
    x[i] = y[i] + L2
    i += 1
  end
  x[i] = y[i] + N
  return partition(x[1:k], check = false), (x, y, N, L2, i, true)
end

@inline function Base.iterate(P::PartitionsFixedNumParts{T}, state::Tuple{Vector{T}, Vector{T}, T, T, Int, Bool}) where T
  k = P.k
  x, y, N, L2, i, flag = state

  N == 0 && return nothing

  if flag
    if i < k && N > 1
      N = T(1)
      x[i] = x[i] - 1
      i += 1
      x[i] = y[i] + 1
      # We could do
      #   return partition(x[1:k], check = false), (x,y,N,L2,i,false)
      # here, but apparently having only one `return` in combination with
      # `@inline` leads to only half as many allocations. So we have to do
      # `if flag` ... `if !flag`.
    else
      flag = false
    end
  end
  if !flag
    lcycle = false
    for j in i - 1:-1:1
      L2 = x[j] - y[j] - T(1)
      N = N + T(1)
      if N <= (k-j)*L2
        x[j] = y[j] + L2
        lcycle = true
        break
      end
      N = N + L2
      x[i] = y[i]
      i = j
    end
    lcycle || return nothing
    while N > L2
      N -= L2
      x[i] = y[i] + L2
      i += 1
    end
    x[i] = y[i] + N
  end
  return partition(x[1:k], check = false), (x,y,N,L2,i,!flag)
end

@doc raw"""
    partitions(n::T, v::Vector{T}) where T <: IntegerUnion
    partitions(n::T, v::Vector{T}, mu::Vector{<:IntegerUnion}) where T <: IntegerUnion
    partitions(n::T, k::IntegerUnion, v::Vector{T}, mu::Vector{<:IntegerUnion}) where T <: IntegerUnion

Return an iterator over all partitions of a non-negative integer `n` where each
part is an element in the vector `v` of positive integers.
It is assumed that the entries in `v` are strictly increasing.

If the optional vector `mu` is supplied, then each `v[i]` occurs a maximum of
`mu[i] > 0` times per partition.

If the optional integer `k >= 0` is supplied, the partitions will be into `k`
parts. In this case, the partitions are produced in lexicographically *decreasing*
order.

The implemented algorithm is "partb" in [RJ76](@cite).

# Examples
The number of partitions of 100 where the parts are from {1, 2, 5, 10, 20, 50}:
```jldoctest
julia> length(collect(partitions(100, [1, 2, 5, 10, 20, 50])))
4562
```
All partitions of 100 where the parts are from {1, 2, 5, 10, 20, 50} and each
part is allowed to occur at most twice:
```jldoctest
julia> collect(partitions(100, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))
6-element Vector{Partition{Int64}}:
 [50, 50]
 [50, 20, 20, 10]
 [50, 20, 20, 5, 5]
 [50, 20, 10, 10, 5, 5]
 [50, 20, 20, 5, 2, 2, 1]
 [50, 20, 10, 10, 5, 2, 2, 1]
```
The partitions of 100 into seven parts, where the parts are required to be
elements from {1, 2, 5, 10, 20, 50} and each part is allowed to occur at most twice.
```jldoctest
julia> collect(partitions(100, 7, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]))
1-element Vector{Partition{Int64}}:
 [50, 20, 20, 5, 2, 2, 1]
```
"""
function partitions(n::T, v::Vector{T}) where T <: IntegerUnion
  return PartitionsFixedValues(n, v)
end

function partitions(n::T, v::Vector{T}, mu::Vector{S}) where {T <: IntegerUnion, S <: IntegerUnion}
  return PartitionsFixedValues(n, v, convert(Vector{Int}, mu))
end

base(P::PartitionsFixedValues) = P.n
parts_min(P::PartitionsFixedValues) = P.kmin
parts_max(P::PartitionsFixedValues) = P.kmax
values(P::PartitionsFixedValues) = P.v
multiplicities(P::PartitionsFixedValues) = P.mu

Base.eltype(::Type{PartitionsFixedValues{T}}) where T = Partition{T}
Base.IteratorSize(::Type{PartitionsFixedValues{T}}) where T = Base.SizeUnknown()

function Base.show(io::IO, ::MIME"text/plain", P::PartitionsFixedValues)
  print(pretty(io), "Iterator over the partitions of $(base(P)) with fixed values")
end

# We iterate over the different iterators of n into k parts for k in
# parts_min(P), ..., parts_max(P)
@inline function iterate(P::PartitionsFixedValues{T}, state::Union{Nothing, Tuple{PartitionsFixedNumPartsAndValues{T}, PartitionsFixedNumPartsAndValuesState{T}}} = nothing) where T
  n = base(P)

  if isnothing(state)
    Pk = partitions(n, parts_min(P), values(P), multiplicities(P))
    s = nothing
  else
    Pk, s = state
  end
  next = iterate(Pk, s)

  while isnothing(next)
    parts(Pk) == parts_max(P) && return nothing
    Pk = partitions(n, parts(Pk) + 1, values(P), multiplicities(P))
    next = iterate(Pk)
  end

  return next[1], (Pk, next[2])
end

function partitions(n::T, k::IntegerUnion, v::Vector{T}, mu::Vector{S}) where {T <: IntegerUnion, S <: IntegerUnion}
  return PartitionsFixedNumPartsAndValues(n, Int(k), v, convert(Vector{Int}, mu))
end

base(P::PartitionsFixedNumPartsAndValues) = P.n
parts(P::PartitionsFixedNumPartsAndValues) = P.k
values(P::PartitionsFixedNumPartsAndValues) = P.v
multiplicities(P::PartitionsFixedNumPartsAndValues) = P.mu

Base.eltype(::Type{PartitionsFixedNumPartsAndValues{T}}) where T = Partition{T}

function Base.show(io::IO, ::MIME"text/plain", P::PartitionsFixedNumPartsAndValues)
  print(pretty(io), "Iterator over the partitions of $(base(P)) into ",
  ItemQuantity(P.k, "part"), " with fixed values")
end

Base.IteratorSize(::Type{PartitionsFixedNumPartsAndValues{T}}) where T = Base.SizeUnknown()

# Algorithm "partb" in [RJ76](@cite), de-gotoed from old ALGOL 60 code by E. Thiel.
# The algorithm as published in the paper has several issues and we hope to have fixed
# them all, see below for details. Some initial fixing was done by T. Schmit.
# Initialize a `state` object for the iterator
function _initial_state(P::PartitionsFixedNumPartsAndValues{T}) where T
  state = PartitionsFixedNumPartsAndValuesState{T}()
  state.done = false

  if is_empty(multiplicities(P))
    # There are no partitions
    state.done = true
    return state
  end

  n = base(P)
  k = parts(P)
  v = values(P)
  mu = multiplicities(P)

  # Initialize variables
  r = length(v)
  j = 1
  m = mu[1]
  x = zeros(T, k)
  y = zeros(T, k)
  ii = zeros(T, k)
  N = n

  # The first step in the algorithm is to initialize the arrays x and y
  # for the backtrack search.
  # This fills the array x from right with the values from v from the left
  # (the smallest) up to their specified multiplicity.
  # If this is larger than n, there cannot be a partition.
  for i = k:-1:1
    x[i] = v[j]
    y[i] = v[j]
    m = m - 1
    N = N - v[j]
    if m == 0
      if j == r
        # In the original algorithm there's a goto b3 here which means
        # the program will terminate, and thus return an empty list.
        # But this is wrong, an example is partitions(1, 1, [1], [1]).
        # What we need to do instead here is to break the loop.
        break
      end
      j = j + 1
      m = mu[j]
    end
  end

  state.x = x
  state.y = y
  state.ii = ii
  state.N = N
  state.i = 1
  state.r = r
  return state
end

# Check whether there is obviously at most one partition.
# Should be called before the first (actual) iteration and modifies `state` in place.
# Return `flag, p` where `flag` is true if there is at most one partition.
# `p` is that partition or `nothing`
function _is_trivial!(P::PartitionsFixedNumPartsAndValues{T}, state::PartitionsFixedNumPartsAndValuesState{T}) where T
  if parts(P) == 0
    if base(P) == 0
      state.done = true
      return true, partition(T[], check = false)
    else
      return true, nothing
    end
  end

  if isempty(multiplicities(P))
    state.done = true
    return true, nothing
  end

  k = parts(P)
  v = values(P)
  x = state.x
  y = state.y
  N = state.N
  r = state.r

  # This is a necessary condition for existence of a partition
  if N < 0 || N > k * (v[r] - v[1])
    return true, nothing
  end

  # The following is a condition for when only a single partition
  # exists. We added the || i == 1 condition because without it the example
  # partitions(17, 7, [1, 4], [1, 4]) returns [0, 0, 4, 4, 4, 4, 1], which is
  # nonsense. So, we have to make sure that all entries of x were modified in
  # the initial backtracking, which means i was counted down to 1.
  # Noticed on Mar 23, 2023.
  if N == 0 && x[1] != 0
    state.done = true
    return true, partition(copy(x), check = false)
  end

  # The iterator is not obviously trivial, so we have to set up N for the first
  # iteration.
  # Logically, this should be done by `_initial_state`, but then we would have
  # to work with N - y[1] in this function which sounds equally horrible.
  state.N = N + y[1]

  return false, nothing
end

@inline function iterate(P::PartitionsFixedNumPartsAndValues{T}, s::Union{Nothing, PartitionsFixedNumPartsAndValuesState{T}} = nothing) where {T <: IntegerUnion}
  first_round = isnothing(s)
  state = first_round ? _initial_state(P) : s

  if first_round
    fl, p = _is_trivial!(P, state)
    if fl
      if !isnothing(p)
        return p, state
      else
        return nothing
      end
    end
  end

  if state.done
    return nothing
  end

  k = parts(P)
  v = values(P)
  mu = multiplicities(P)

  # Initialize variables
  x = state.x
  y = state.y
  ii = state.ii
  N = state.N
  r = state.r
  i = state.i

  p = nothing
  while isnothing(p) && !state.done
    state.done = true

    j = 1
    while j <= mu[r] && N > v[r]
      x[i] = v[r]
      ii[i] = r - 1

      # Added, otherwise get out of bounds
      if i == k
        r = r - 1
        break # inner while loop
      end
      i += 1
      N = N - v[r] + y[i]
      j += 1
    end

    if j > mu[r]
      r = r - 1
    end

    while r > 0 && v[r] > N # Added additional r > 0, otherwise get out of bounds
      r = r - 1
    end
    r == 0 && return nothing

    if N == v[r]
      x[i] = v[r]
      p = partition(copy(x), check = false)
      r = r - 1
      # Added, otherwise get out of bounds
      r == 0 && break
    end

    m = y[i]
    while true
      # Here comes the most intricate mistake.
      # On "Knuth's" problem
      # 100, 7, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]
      # the published algorithm does not find the valid partition
      # 50 + 20 + 20 + 5 + 2 + 2 + 1.
      # But when we replace v[r] > k by v[r] >= k, everything works.
      # Finding this was a wild guess!
      # Found on Mar 27, 2023.
      if v[r] >= m && N - v[r] <= (k - i)*(v[r] - v[1])
        state.done = false
        break # inner while loop
      end
      x[i] = m
      i -= 1
      i == 0 && break # inner while loop
      r = Int(ii[i])
      N = N + x[i] - m
      m = y[i]
    end
  end
  if isnothing(p)
    # We didn't find a partition
    return nothing
  end

  state.x = x
  state.y = y
  state.ii = ii
  state.N = N
  state.i = i
  state.r = r
  return p, state
end

################################################################################
#
# Relations
#
################################################################################

@doc raw"""
    dominates(lambda::Partition, mu::Partition)

Return `true` if `lambda` dominates `mu`, `false` otherwise.

# Examples
```jldoctest
julia> dominates(partition(3, 1), partition(2, 2))
true

julia> dominates(partition(4, 1), partition(3, 3))
false
```
"""
function dominates(lambda::Partition, mu::Partition)
  dif = 0
  i = 1
  while i <= min(length(lambda), length(mu))
    dif += lambda[i] - mu[i]
    i += 1
    if dif < 0
      return false
    end
  end
  if length(lambda) < length(mu)
    while i <= length(mu)
      dif -= mu[i]
      i += 1
    end
    if dif < 0
      return false
    end
  end
  return true
end

################################################################################
#
# Operations
#
################################################################################

@doc raw"""
    conjugate(lambda::Partition)

Return the conjugate of the partition `lambda`.

# Examples
```jldoctest
julia> conjugate(partition(8, 8, 8, 7, 2, 1, 1))
[7, 5, 4, 4, 4, 4, 4, 3]
```
"""
function conjugate(lambda::Partition{T}) where T <: IntegerUnion
  if isempty(lambda)
    return copy(lambda)
  end

  mu = zeros(T, lambda[1])

  for i = 1:length(lambda)
    for j = 1:lambda[i]
      mu[j] += 1
    end
  end

  return partition(mu, check = false)
end
