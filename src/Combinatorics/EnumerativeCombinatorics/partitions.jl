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
# Constructors and basic functionality
#
################################################################################

function _defines_partition(parts::Vector{<: IntegerUnion})
  if isempty(parts)
    return true
  end
  return all(i -> parts[i] >= parts[i + 1], 1:length(parts) - 1) && is_positive(parts[end])
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
  print(io, data(P))
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

function Base.copy(P::Partition{T}) where T <: IntegerUnion
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

################################################################################
#
# Generating and counting unrestricted partitions
#
################################################################################

@doc raw"""
    number_of_partitions(n::IntegerUnion)

Return the number of integer partitions of `n`.

# Examples
```jldoctest
julia> number_of_partitions(1000)
24061467864032622473692149727991
```
"""
function number_of_partitions(n::IntegerUnion)
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
function partitions(n::T) where T <: IntegerUnion

  #Argument checking
  @req n >= 0 "n >= 0 required"

  # Some trivial cases
  if n == 0
    return (p for p in Partition{T}[ partition(T[], check = false) ])
  elseif n == 1
    return (p for p in Partition{T}[ partition(T[1], check = false) ])
  end

  # Now, the algorithm starts
  P = Partition{T}[]    #this will be the array of all partitions
  k = 1
  q = 1
  d = fill( T(1), n )
  d[1] = n
  push!(P, partition(d[1:1], check = false))
  while q != 0
    if d[q] == 2
      k += 1
      d[q] = 1
      q -= 1
    else
      m = d[q] - 1
      np = k - q + 1
      d[q] = m
      while np >= m
        q += 1
        d[q] = m
        np = np - m
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
    push!(P, partition(d[1:k], check = false))
  end
  return (p for p in P)
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
"""
function number_of_partitions(n::IntegerUnion, k::IntegerUnion)

  @req n >= 0 "n >= 0 required"
  @req k >= 0 "k >= 0 required"

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
      for m = l + 1:n - l + 1
        p[m] = p[m] + p[m - l]
      end
    end
    return p[n - k + 1]
  end

end

@doc raw"""
    partitions(m::IntegerUnion, n::IntegerUnion; only_distinct_parts::Bool = false)
    partitions(m::IntegerUnion, n::IntegerUnion, l1::IntegerUnion, l2::IntegerUnion; only_distinct_parts::Bool = false)

Return an iterator over all partitions of a non-negative integer `m` into
`n >= 0` parts. Optionally, a lower bound `l1 >= 0` and an upper bound `l2` for
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
function partitions(m::T, n::IntegerUnion, l1::IntegerUnion, l2::IntegerUnion; only_distinct_parts::Bool = false) where T <: IntegerUnion
  # Algorithm "parta" in [RJ76](@cite), de-gotoed from old ALGOL 60 code by E. Thiel.

  # Note that we are considering partitions of m here. I would switch m and n
  # but the algorithm was given like that and I would otherwise confuse myself
  # implementing it.

  #Argument checking
  @req m >= 0 "m >= 0 required"
  @req n >= 0 "n >= 0 required"
  @req l1 >= 0 "l1 >= 0 required"

  # Use type of m
  n = convert(T, n)

  # Some trivial cases
  if m == 0 && n == 0
    return (p for p in Partition{T}[ partition(T[], check = false) ])
  end

  if n == 0 || n > m
    return (p for p in Partition{T}[])
  end

  if l2 < l1
    return (p for p in Partition{T}[])
  end

  # If l1 == 0 the algorithm parta will actually create lists containing the
  # entry zero, e.g. partitions(2, 2, 0, 2) will contain [2, 0].
  # This is nonsense, so we set l1 = 1 in this case.
  if l1 == 0
    l1 = 1
  end

  #Algorithm starts here
  P = Partition{T}[]    #this will be the array of all partitions
  x = zeros(T, n)
  y = zeros(T, n)
  j = only_distinct_parts*n*(n - 1)
  m = m - n*l1 - div(j, 2)
  l2 = l2 - l1
  if 0 <= m <= n*l2 - j

    for i = 1:n
      y[i] = x[i] = l1 + only_distinct_parts*(n - i)
    end

    i = 1
    l2 = l2 - only_distinct_parts*(n - 1)

    while true
      while m > l2
        m -= l2
        x[i] = y[i] + l2
        i += 1
      end

      x[i] = y[i] + m
      push!(P, partition(x[1:n], check = false))

      if i < n && m > 1
        m = 1
        x[i] = x[i] - 1
        i += 1
        x[i] = y[i] + 1
        push!(P, partition(x[1:n], check = false))
      end

      lcycle = false
      for j = i - 1:-1:1
        l2 = x[j] - y[j] - 1
        m = m + 1
        if m <= (n - j)*l2
          x[j] = y[j] + l2
          lcycle = true
          break
        end
        m = m + l2
        x[i] = y[i]
        i = j
      end

      if !lcycle
        break
      end
    end
  end

  return (p for p in P)
end

function partitions(m::T, n::IntegerUnion; only_distinct_parts::Bool = false) where T <: IntegerUnion
  @req m >= 0 "m >= 0 required"
  @req n >= 0 "n >= 0 required"

  # Special cases
  if m == n
    return (p for p in [ partition(T[ 1 for i in 1:m], check = false) ])
  elseif m < n || n == 0
    return (p for p in Partition{T}[])
  elseif n == 1
    return (p for p in [ partition(T[m], check = false) ])
  end

  return (p for p in partitions(m, n, 1, m; only_distinct_parts = only_distinct_parts))
end

function partitions(m::T, n::IntegerUnion, v::Vector{T}, mu::Vector{S}) where {T <: IntegerUnion, S <: IntegerUnion}
  # Algorithm "partb" in [RJ76](@cite), de-gotoed from old ALGOL 60 code by E. Thiel.
  # The algorithm as published in the paper has several issues and we hope to have fixed
  # them all, see below for details. Some initial fixing was done by T. Schmit.

  @req m >= 0 "m >= 0 required"
  @req n >= 0 "n >= 0 required"
  @req length(mu) == length(v) "mu and v should have the same length"

  # Algorithm partb assumes that v is strictly increasing.
  # Added (and noticed) on Mar 22, 2023.
  @req all([v[i] < v[i + 1] for i in 1:length(v) - 1]) "v must be strictly increasing"

  # Parta allows v[1] = 0 but this is nonsense for entries of a partition
  @req all(>(0), v) "Entries of v must be positive"

  # For safety
  @req all(>(0), mu) "Entries of mu must be positive"

  # Special cases
  if n == 0
    # TODO: I don't understand this distinction here
    # (it also makes the function tabe instable)
    if m == 0
      return (p for p in [ Partition{T}[] ])
    else
      return (p for p in Partition{T}[])
    end
  end

  if isempty(mu)
    return (p for p in Partition{T}[])
  end

  #This will be the list of all partitions found.
  P = Partition{T}[]

  # Now, we get to the partb algorithm. This is a hell of an algorithm and the
  # published code has several issues.
  # The original ALGOL 60 code is electronically available at
  # https://gist.github.com/ulthiel/99de02994fc31fe614586ed0c930f744.
  # First, there were some issues with termination and indices for the arrays 
  # x, y, ii getting out of bounds. We had to introduce some additional checks 
  # and breaks to take care of this. Some initial fixing was done by T. Schmit.
  # An example showing the problem is 17, 3, [1, 4], [1, 4].

  # Initialize variables
  r = length(v)
  j = 1
  k = mu[1]
  ll = v[1]
  x = zeros(T, n)
  y = zeros(T, n)
  ii = zeros(T, n)

  # The algorithm has three goto labels b1, b2, b3.
  # b3 terminates the algorithm.
  # We introduce bools to indicate the jumps.
  gotob2 = false
  gotob1 = true

  # The first step in the algorithm is to initialize the arrays x and y
  # for the backtrack search. 
  # This fills the array x from right with the values from v from the left 
  # (the smallest) up to their specified multiplicity.
  # If this is larger than m, there cannot be a partition.
  for i = n:-1:1
    x[i] = ll
    y[i] = ll
    k = k - 1
    m = m - ll
    if k == 0
      if j == r
        # In the original algorithm there's a goto b3 here which means
        # the program will terminate, and thus return an empty list.
        # But this is wrong, an example is partitions(1, 1, [1], [1]).
        # What we need to do instead here is to break the loop.
        break
      end
      j = j + 1
      k = mu[j]
      ll = v[j]
    end
    if i == 1
      break
    else
      i = i - 1
    end
  end #for i

  lr = v[r]
  ll = v[1]

  # This is a necessary condition for existence of a partition
  if m < 0 || m > n * (lr - ll)
    return (p for p in P) #goto b3
  end

  # The following is a condition for when only a single partition
  # exists. We added the || i == 1 condition because without it the example
  # partitions(17, 7, [1, 4], [1, 4]) returns [0, 0, 4, 4, 4, 4, 1], which is 
  # nonsense. So, we have to make sure that all entries of x were modified in
  # the initial backtracking, which means i was counted down to 1.
  # Noticed on Mar 23, 2023.
  if m == 0 && x[1] != 0
    push!(P, partition(copy(x), check = false))
    return (p for p in P)
  end

  # Now, the actual algorithm starts
  i = 1
  m = m + y[1]

  # label b1
  while gotob1 == true
    if !gotob2
      for j = mu[r]:-1:1
        if m <= lr
          gotob2 = true
          break
        end
        x[i] = lr
        ii[i] = r - 1
        if i == n # Added, otherwise get out of bounds
          break
        end
        i = i + 1
        m = m - lr + y[i]
      end #for j

      if !gotob2
        r = r - 1
      end

      gotob2 = true
    end #if

    # label b2
    if gotob2
      while r > 0 && v[r] > m # Added additional r > 0, otherwise get out of bounds
        r = r - 1
      end

      if r == 0
        break
      end

      lr = v[r]
      if m == lr
        x[i] = lr
        if i <= n # Added, otherwise get out of bounds
          push!(P, partition(copy(x), check = false)) #need copy here!
        else
          break
        end

        r = r - 1
        if r == 0 # Added, otherwise get out of bounds
          break
        end
        lr = v[r]
      end #if

      k = y[i]
      # Here comes the most intricate mistake.
      # On "Knuth's" problem
      # 100, 7, [1, 2, 5, 10, 20, 50], [2, 2, 2, 2, 2, 2]
      # the published algorithm does not find the valid partition
      # 50 + 20 + 20 + 5 + 2 + 2 + 1.
      # But when we replace lr > k by lr >= k, everything works.
      # Finding this was a wild guess!
      # Found on Mar 27, 2023.
      if lr >= k && m - lr <= (n - i)*(lr - ll)
        gotob2 = false
        continue
      else
        x[i] = k
      end #if
      for i_0 = i - 1:-1:1 #this is to replace the for i = i - 1 in ALGOL code
        i = i_0
        r = ii[i]
        lr = v[r]
        m = m + x[i] - k
        k = y[i]
        if lr >= k && m - lr <= (n - i)*(lr - ll) # >= here as well (probably...)
          gotob2 = false
          break
        else
          x[i] = k
        end #if
      end #for
      if gotob2
        gotob1 = false
      end
    end #if gotob2
  end #while

  return (p for p in P)
end

function partitions(m::T, v::Vector{T}, mu::Vector{S}) where {T <: IntegerUnion, S <: IntegerUnion}

  @req m >= 0 "m >= 0 required"
  @req length(mu) == length(v) "mu and v should have the same length"

  # Algorithm partb assumes that v is strictly increasing.
  # Added (and noticed) on Mar 22, 2023.
  @req all([v[i] < v[i + 1] for i in 1:length(v) - 1]) "v must be strictly increasing"

  # Parta allows v[1] = 0 but this is nonsense for entries of a partition
  @req all(>(0), v) "Entries of v must be positive"

  # For safety
  @req all(>(0), mu) "Entries of mu must be positive"
  res = Partition{T}[]

  if isempty(v)
    return (p for p in res)
  end

  if m == 0
    # TODO: I don't understand this return (and it is type instable)
    return (p for p in [ Partition{T}[] ])
  end

  # We will loop over the number of parts.
  # We first determine the minimal and maximal number of parts.
  r = length(v)

  nmax = 0
  cursum = 0
  for i in 1:r, j in 1:mu[i]
    nmax += 1
    cursum += v[i]
    if cursum >= m
      break
    end
  end

  nmin = 0
  cursum = 0
  for i in r:-1:1, j in 1:mu[i]
    nmin += 1
    cursum += v[i]
    if cursum >= m
      break
    end
  end

  for n = nmin:nmax
    append!(res, partitions(m, n, v, mu))
  end

  return (p for p in res)
end

@doc raw"""
    partitions(m::T, v::Vector{T}) where T <: IntegerUnion
    partitions(m::T, v::Vector{T}, mu::Vector{<:IntegerUnion}) where T <: IntegerUnion
    partitions(m::T, n::IntegerUnion, v::Vector{T}, mu::Vector{<:IntegerUnion}) where T <: IntegerUnion

Return an iterator over all partitions of a non-negative integer `m` where each
part is an element in the vector `v` of positive integers.
It is assumed that the entries in `v` are strictly increasing.

If the optional vector `mu` is supplied, then each `v[i]` occurs a maximum of
`mu[i] > 0` times per partition.

If the optional integer `n >= 0` is supplied, the partitions will be into `n`
parts. In this case, the partitions are produced in lexicographically *decreasing*
order.

The implemented algorithm is "partb" in [RJ76](@cite).

# Example
The number of partitions of 100 where the parts are from {1, 2, 5, 10, 20, 50}:
```jldoctest
julia> length(partitions(100, [1, 2, 5, 10, 20, 50]))
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
function partitions(m::T, v::Vector{T}) where T <: IntegerUnion
  @req m >= 0 "m >= 0 required"
  @req all([v[i] < v[i + 1] for i in 1:length(v) - 1]) "v must be strictly increasing"
  @req all(>(0), v) "Entries of v must be positive"

  res = Partition{T}[]

  if isempty(v)
    return (p for p in res)
  end

  if m == 0
    # TODO: I don't understand this return (and it is type instable)
    return (p for p in [ Partition{T}[] ])
  end

  # We will loop over the number of parts.
  # We first determine the minimal and maximal number of parts.
  r = length(v)
  nmin = div(m, v[r])
  nmax = div(m, v[1])

  # Set the maximum multiplicity (equal to nmax above)
  mu = [ nmax for i in 1:r ]

  for n = nmin:nmax
    append!(res, partitions(m, n, v, mu))
  end

  return (p for p in res)
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
