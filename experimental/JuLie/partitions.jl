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

export Partition

export ascending_partitions
export conjugate
export dominates
export getindex_safe
export num_partitions
export partitions 


################################################################################
# Partition type
################################################################################
@doc raw"""
    Partition{T<:IntegerUnion} <: AbstractVector{T}

A **partition** of a non-negative integer ``n`` is a decreasing sequence ``λ₁ ≥ λ₂ ≥ … ≥
λᵣ`` of positive integers ``λᵢ`` such that ``n = λ₁ + … + λᵣ``. The ``λᵢ`` are called the
**parts** of the partition and ``r`` is called the **length**. 

A partition can be encoded as an array with elements ``λᵢ``. We provide the parametric type
`Partition{T}` which is a subtype of `AbstractVector{T}` where `T` can be any subtype of
`IntegerUnion`. All functions that can be used for vectors (1-dimensional arrays) can thus
be used for partitions as well. There is no performance impact by using an own type for
partitions rather than simply using arrays. The parametric type allows to increase
performance by using smaller integer types. For efficiency, the `Partition` constructor does
not check whether the given array is indeed a decreasing sequence.

# Examples
A partition can be created by either calling `Partition` on an array of integers or by
calling `Partition` with arguments being the sequence of parts.
```jldoctest
julia> P = Partition([6,4,4,2]) #The partition 6+4+4+2 of 16.
[6, 4, 4, 2]

julia> P = Partition(6,4,4,2) #Same as above but less to type
[6, 4, 4, 2]

julia> length(P)
4

julia> P[1]
6
```
Usually, ``|λ| ≔ n`` is called the **size** of ``λ``. In Julia, the function `size` for
arrays already exists and returns the *dimension* of an array. Instead, you can use the
Julia function `sum` to get the sum of the parts.
```jldoctest
julia> P = Partition(6,4,4,2)
[6, 4, 4, 2]

julia> sum(P) 
16
```
You can create partitions with smaller integer types as follows.
```jldoctest
julia> P = Partition{Int8}(6,4,4,2) #Or Partition(Int8[6,4,4,2])
Int8[6, 4, 4, 2]
```
There is a unique partition of 0, namely the **empty partition** (of length 0). It can be
created as follows.
```jldoctest
julia> P = Partition() #Or Partition([])
Int64[]
julia> sum(P)
0
julia> length(P)
0
julia> P = Partition{Int8}() #Or Partition(Int8[])
Int8[]
```

# References
1. [Ful97](@cite)
2. [Knu11](@cite), Section 7.2.1.4 (starting on page 390).
"""
struct Partition{T<:IntegerUnion} <: AbstractVector{T}
  p::Vector{T}
end

# The following are functions to make the Partition struct array-like.
function Base.show(io::IO, ::MIME"text/plain", P::Partition)
  print(io, P.p)
end

function Base.size(P::Partition)
  return size(P.p)
end

function Base.length(P::Partition)
  return length(P.p)
end

function Base.getindex(P::Partition, i::IntegerUnion)
  return getindex(P.p,Int(i))
end

function Base.setindex!(P::Partition, x::IntegerUnion, i::IntegerUnion)
  return setindex!(P.p,x,Int(i))
end

function Partition(parts::IntegerUnion...)
  return Partition(collect(Int, parts))
end

function Partition{T}(parts::IntegerUnion...) where T<:IntegerUnion
  return Partition(collect(T, parts))
end

# The empty array is of "Any" type, and this is not what we want. 
# We want it to be an array of integers of the default type Int64. 
function Partition(p::Vector{Any})
  return Partition(Vector{Int64}(p))
end

function Base.copy(P::Partition{T}) where T<:IntegerUnion
  return Partition{T}(copy(P.p))
end

@doc raw"""
    getindex_safe(P::Partition, i::IntegerUnion)

In algorithms involving partitions it is sometimes convenient to be able to access parts
beyond the length of the partition and then one wants to get the value zero instead of an
error. This function is a shortcut for
```
return (i>length(P.p) ? 0 : getindex(P.p,i))
```
If you are sure that `P[i]` exists, use `getindex` because this will be faster.

# Examples
```jldoctest
julia> P=Partition([3,2,1])
[3, 2, 1]

julia> getindex_safe(P, 3)
1

julia> getindex_safe(P, 4)
0
```
"""
function getindex_safe(P::Partition, i::IntegerUnion)
  return (i>length(P.p) ? 0 : getindex(P.p,Int(i)))
end


################################################################################
# Generating and counting unrestricted partitions
################################################################################

@doc raw"""
    num_partitions(n::IntegerUnion)

The number of integer partitions of the non-negative integer `n`. 

# Examples
```jldoctest
julia> num_partitions(1000)
24061467864032622473692149727991
```

# Algorithm

We use the function
[`arith_number_of_partitions`](http://flintlib.org/doc/arith.html?highlight=partitions#c.arith_number_of_partitions)
from [FLINT](@cite) which is very fast. It is based on the Hardy-Ramanujan-Rademacher
formula, see [Joh12](@cite) for details. 

# Further references
1. [Knu11](@cite), Section 7.2.1.4 (starting on page 395).
2. [OEIS](@cite), [A000041](https://oeis.org/A000041)
"""
function num_partitions(n::IntegerUnion)
  @req n >= 0 "n >= 0 required"
  n = ZZ(n)
  z = ZZ()
  ccall((:arith_number_of_partitions, Nemo.libflint), Cvoid, (Ref{ZZRingElem}, Culong), z, UInt(n))
  return z
end


@doc raw"""
    partitions(n::IntegerUnion)

A list of all partitions of a non-negative integer `n`, produced in lexicographically
*descending* order. This ordering is like in Sage, but opposite to GAP. You can apply the
function `reverse` to reverse the order. As usual, you may increase performance by using
smaller integer types.

# Algorithm
The algorithm used is "Algorithm ZS1" by [ZS98](@cite). This algorithm is also discussed in
[Knu11](@cite), Algorithm P (page 392).

# Examples
```jldoctest
julia> partitions(4) # Use partitions(Int8(4)) to use 8-bit integers
5-element Vector{Partition{Int64}}:
 [4]
 [3, 1]
 [2, 2]
 [2, 1, 1]
 [1, 1, 1, 1]
```
"""
function partitions(n::IntegerUnion)

  #Argument checking
  @req n >= 0 "n >= 0 required"

  # Use type of n
  T = typeof(n)

  # Some trivial cases
  if n == 0
    return Partition{T}[ Partition{T}([]) ]
  elseif n == 1
    return Partition{T}[ Partition{T}([1]) ]
  end

  # Now, the algorithm starts
  P = Partition{T}[]    #this will be the array of all partitions
  k = 1
  q = 1
  d = fill( T(1), n )
  d[1] = n
  push!(P, Partition{T}(d[1:1]))
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
    push!(P, Partition{T}(d[1:k]))
  end
  return P

end


@doc raw"""
    ascending_partitions(n::IntegerUnion;alg="ks")

Instead of encoding a partition of an integer ``n ≥ 0`` as a *descending*
sequence (which is our convention), one can also encode it as an *ascending*
sequence. In the papers Kelleher & O'Sullivan (2014) and Merca (2012) it is
said that generating the list of all ascending partitions is more efficient
than generating descending ones. To test this, I have implemented the
algorithms given in the papers:
1. "ks" (*default*) is the algorithm "AccelAsc" (Algorithm 4.1) in [KO14](@cite).
2. "m" is Algorithm 6 in [Mer12](@cite). This is actually similar to "ks".

The ascending partitions are stored here as arrays and are not of type
`Partition` since the latter are descending by our convention. I am using "ks"
as default since it looks slicker and I believe there is a tiny mistake in the
publication of "m" (which I fixed).

# Comparison

I don't see a significant speed difference to the descending encoding:
```julia-repl
julia> @btime partitions(Int8(90));
3.376 s (56634200 allocations: 6.24 GiB)

julia> @btime ascending_partitions(Int8(90),alg="ks");
3.395 s (56634200 allocations: 6.24 GiB)

julia> @btime ascending_partitions(Int8(90),alg="m");
3.451 s (56634200 allocations: 6.24 GiB)
```
"""
function ascending_partitions(n::IntegerUnion; alg="ks")

  #Argument checking
  @req n >= 0 "n >= 0 required"

  # Use type of n
  T = typeof(n)

  # Some trivial cases
  if n==0
    return Vector{T}[ [] ]
  elseif n==1
    return Vector{T}[ [1] ]
  end

  # Now, the algorithm starts
  if alg=="ks"
    P = Vector{T}[]    #this will be the array of all partitions
    a = zeros(T, n)
    k = 2
    y = n-1
    while k != 1
      k -= 1
      x = a[k] + 1
      while 2*x <= y
        a[k] = x
        y -= x
        k += 1
      end
      l = k + 1
      while x <= y
        a[k] = x
        a[l] = y
        push!(P, a[1:l])
        x += 1
        y -= 1
      end
      y += x - 1
      a[k] = y + 1
      push!(P, a[1:k])
    end
    return P

  elseif alg=="m"
    P = Vector{T}[]    #this will be the array of all partitions
    a = zeros(T, n)
    k = 1
    x = 1
    y = n-1
    while true
      while 3*x <= y
        a[k] = x
        y = y-x
        k += 1
      end
      t = k + 1
      u = k + 2
      while 2*x <= y
        a[k] = x
        a[t] = x
        a[u] = y - x
        push!(P, a[1:u])
        p = x + 1
        q = y - p
        while p <= q
          a[t] = p
          a[u] = q
          push!(P, a[1:u])
          p += 1
          q -= 1
        end
        a[t] = y
        push!(P, a[1:t])
        x += 1
        y -= 1
      end
      while x<=y
        a[k] = x
        a[t] = y
        push!(P, a[1:t])
        x += 1
        y -= 1
      end
      y += x-1
      a[k] = y+1
      push!(P, a[1:k])
      k -= 1

      # I think there's a mistake in the publication
      # because here k could be zero and then we access
      # a[k].
      # That's why I do a while true and check k > 0 here.
      if k == 0
        break
      else
        x = a[k] + 1
      end
    end
    return P
  else
    error("alg must be either ks or m")
  end

end

################################################################################
# Generating and counting restricted partitions
################################################################################

@doc raw"""
    num_partitions(n::IntegerUnion, k::IntegerUnion)

The number of integer partitions of the non-negative integer `n` into `k >= 0` parts. 

# Algorithm
We use the recurrence relation ``p_k(n) = p_{k-1}(n-1) + p_k(n-k)``, where ``p_k(n)``
denotes the number of partitions of ``n`` into ``k`` parts; see [Knu11](@cite), Section
7.2.1.4, Equation (39) on page 399.

# References
1. [OEIS](@cite), [A008284](https://oeis.org/A008284)
"""
function num_partitions(n::IntegerUnion, k::IntegerUnion)
  
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
    return num_partitions(n-k) #n-k>=0 holds since the case n<k was already handled

  # See https://oeis.org/A008284
  elseif n <= 2+3*k
    p = num_partitions(n-k) #n-k>=0 holds since the case n<k was already handled
    for i=0:Int(n)-2*Int(k)-1
      p = p - num_partitions(ZZ(i))
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
      for m = l+1:n-l+1
        p[m] = p[m] + p[m-l]
      end
    end
    return p[n-k+1]
  end

end

@doc raw"""
    partitions(m::T, n::IntegerUnion, l1::IntegerUnion, l2::IntegerUnion; only_distinct_parts::Bool = false) where T <: IntegerUnion

A list of all partitions of a non-negative integer `m` into `n >= 0` parts with lower bound
`l1 >= 0` and upper bound `l2` for the parts. There are two choices for the parameter
`only_distinct_parts`:
* `false`: no further restriction (*default*);
* `true`: only distinct parts. The partitions are produced in *decreasing* order.

# Examples
We compute all partitions of 7 into 3 parts where all parts are between 1 and 4:
```jldoctest
julia> partitions(7, 3, 1, 4)
3-element Vector{Partition{Int64}}:
 [4, 2, 1]
 [3, 3, 1]
 [3, 2, 2]
```

Same as above but requiring all parts to be distinct:
```jldoctest
julia> partitions(7, 3, 1, 4 ; only_distinct_parts=true)
1-element Vector{Partition{Int64}}:
 [4, 2, 1]
```
# Algorithm
The algorithm used is "parta" in [RJ76](@cite), de-gotoed from old ALGOL 60 code by E.
Thiel.
"""
function partitions(m::T, n::IntegerUnion, l1::IntegerUnion, l2::IntegerUnion; only_distinct_parts::Bool = false) where T <: IntegerUnion

  # Note that we are considering partitions of m here. I would switch m and n
  # but the algorithm was given like that and I would otherwise confuse myself
  # implementing it.

  #Argument checking
  @req m >= 0 "m >= 0 required"
  @req n >= 0 "n >= 0 required"
  @req l1 >= 0 "l1 >= 0 required"

  # Use type of n
  n = convert(T, n)

  # Some trivial cases
  if m == 0 && n == 0
    return Partition{T}[ Partition{T}([]) ]
  end

  if n == 0 || n > m
    return Partition{T}[]
  end

  if l2 < l1
    return Partition{T}[]
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
  j = only_distinct_parts*n*(n-1)
  m = m-n*l1-div(j,2)
  l2 = l2 - l1
  if 0 <= m <= n*l2-j

    for i = 1:n
      y[i] = x[i] = l1+only_distinct_parts*(n-i)
    end

    i = 1
    l2 = l2-only_distinct_parts*(n-1)

    while true
      while m > l2
        m -= l2
        x[i] = y[i] + l2
        i += 1
      end

      x[i] = y[i] + m
      push!(P, Partition{T}(x[1:n]))

      if i < n && m > 1
        m = 1
        x[i] = x[i] - 1
        i += 1
        x[i] = y[i] + 1
        push!(P, Partition{T}(x[1:n]))
      end

      lcycle = false
      for j = i-1:-1:1
        l2 = x[j] - y[j] - 1
        m = m + 1
        if m <= (n-j)*l2
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

  return P
end


@doc raw"""
    partitions(m::T, n::IntegerUnion) where T<:IntegerUnion

All partitions of a non-negative integer `m` into `n` parts (no further restrictions).

# Examples
```jldoctest
# All partitions of 7 into 3 parts.
julia> partitions(7, 3)
4-element Vector{Partition{Int64}}:
 [5, 1, 1]
 [4, 2, 1]
 [3, 3, 1]
 [3, 2, 2]
```
# Algorithm
This function is a shortcut for `partitions(m,n,1,m;only_distinct_parts=false)`.
"""
function partitions(m::T, n::IntegerUnion) where T<:IntegerUnion
  
  @req m >= 0 "m >= 0 required"
  @req n >= 0 "n >= 0 required"

  # Special cases
  if m == n
    return [ Partition(T[ 1 for i in 1:m]) ]
  elseif m < n || n == 0
    return Partition{T}[]
  elseif n == 1
    return [ Partition(T[m]) ]
  end

  return partitions(m,n,1,m; only_distinct_parts = false)
end


@doc raw"""
    partitions(m::T, n::IntegerUnion, v::Vector{T}, mu::Vector{S}) where {T<:IntegerUnion, S<:IntegerUnion}

All partitions of a non-negative integer `m` into `n >= 0` parts, where each part is an
element in the vector `v` of positive integers and each `v[i]` occurs a maximum of `mu[i] >
0` times. We assume (without loss of generality) that the entries in `v` are strictly
increasing. The partitions are produced in lexicographically *decreasing* order. 

# Examples
We compute the partitions of 100 into seven parts, where the parts are required to be
elements from {1, 2, 5, 10, 20, 50} and each part is allowed to occur at most twice.
```jldoctest
julia> partitions(100, 7, [1,2,5,10,20,50], [2,2,2,2,2,2])
1-element Vector{Partition{Int64}}:
 [50, 20, 20, 5, 2, 2, 1]
```

# Algorithm
The algorithm used is "partb" in [RJ76](@cite), de-gotoed from old ALGOL 60 code by E.
Thiel. The algorithm as published in the paper has several issues and we hope to have fixed
them all, see the code for details. Some initial fixing was done by T. Schmit.
"""
function partitions(m::T, n::IntegerUnion, v::Vector{T}, mu::Vector{S}) where {T<:IntegerUnion,S<:IntegerUnion}
  @req m >= 0 "m >= 0 required"
  @req n >= 0 "n >= 0 required"
  @req length(mu) == length(v) "mu and v should have the same length"

  # Algorithm partb assumes that v is strictly increasing. 
  # Added (and noticed) on Mar 22, 2023.
  @req all([v[i] < v[i+1] for i in 1:length(v)-1]) "v must be strictly increasing"

  # Parta allows v[1] = 0 but this is nonsense for entries of a partition
  @req all(>(0),v) "Entries of v must be positive"

  # For safety
  @req all(>(0),mu) "Entries of mu must be positive"

  # Special cases
  if n == 0
    if m == 0
      return [ Partition{T}[] ]
    else
      return Partition{T}[]
    end
  end

  if isempty(mu)
    return Partition{T}[]
  end
  
  #This will be the list of all partitions found.
  P = Partition{T}[] 

  # Now, we get to the partb algorithm. This is a hell of an algorithm and the
  # published code has several issues.
  # The original ALGOL 60 code is electronically available at
  # https://gist.github.com/ulthiel/99de02994fc31fe614586ed0c930f744.
  # First, there were some issues with termination and indices for the arrays 
  # x,y,ii getting out of bounds. We had to introduce some additional checks 
  # and breaks to take care of this. Some initial fixing was done by T. Schmit.
  # An example showing the problem is 17, 3, [1,4], [1,4].

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
  for i=n:-1:1
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
    return P #goto b3
  end

  # The following is a condition for when only a single partition
  # exists. We added the || i==1 condition because without it the example
  # partitions(17, 7, [1,4], [1,4]) returns [0, 0, 4, 4, 4, 4, 1], which is 
  # nonsense. So, we have to make sure that all entries of x were modified in
  # the initial backtracking, which means i was counted down to 1.
  # Noticed on Mar 23, 2023.
  if m == 0 && x[1] != 0
    push!(P,Partition{T}(copy(x)))
    return P
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
      while r > 0 && v[r] > m # Added additional r>0, otherwise get out of bounds
        r = r - 1
      end

      if r == 0
        break
      end

      lr = v[r]
      if m == lr
        x[i] = lr
        if i <= n # Added, otherwise get out of bounds
          push!(P, Partition{T}(copy(x))) #need copy here!
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
      # 100, 7, [1,2,5,10,20,50], [2,2,2,2,2,2] 
      # the published algorithm does not find the valid partition 
      # 50+20+20+5+2+2+1.
      # But when we replace lr > k by lr >= k, everything works.
      # Finding this was a wild guess!
      # Found on Mar 27, 2023.
      if lr >= k && m - lr <= (n-i)*(lr - ll)
        gotob2 = false
        continue
      else
        x[i] = k
      end #if
      for i_0 = i - 1:-1:1 #this is to replace the for i=i-1 in ALGOL code
        i = i_0
        r = ii[i]
        lr = v[r]
        m = m + x[i] - k
        k = y[i]
        if lr >= k && m - lr <= (n-i)*(lr-ll) # >= here as well (probably...)
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
  
  return P
end

@doc raw"""
    partitions(m::T, v::Vector{T}, mu::Vector{S}) where {T<:IntegerUnion,S<:IntegerUnion}
  
All partitions of a non-negative integer `m` where each part is an element in the vector `v`
of positive integers and each `v[i]` occurs a maximum of `mu[i] > 0` times. We assume
(without loss of generality) that the entries in `v` are strictly increasing.

# Example
We compute all partitions of 100 where the parts are from {1, 2, 5, 10, 20, 50} and each
part is allowed to occurr at most twice:
```jldoctest
julia> partitions(100, [1,2,5,10,20,50], [2,2,2,2,2,2])
6-element Vector{Partition{Int64}}:
 [50, 50]
 [50, 20, 20, 10]
 [50, 20, 20, 5, 5]
 [50, 20, 10, 10, 5, 5]
 [50, 20, 20, 5, 2, 2, 1]
 [50, 20, 10, 10, 5, 2, 2, 1]
```

# Algorithm
We use the function `partitions(m,n,v,mu)`, looping over the number of possible parts of
partitions.
"""
function partitions(m::T, v::Vector{T}, mu::Vector{S}) where {T<:IntegerUnion,S<:IntegerUnion}

  @req m >= 0 "m >= 0 required"
  @req length(mu) == length(v) "mu and v should have the same length"

  # Algorithm partb assumes that v is strictly increasing. 
  # Added (and noticed) on Mar 22, 2023.
  @req all([v[i] < v[i+1] for i in 1:length(v)-1]) "v must be strictly increasing"

  # Parta allows v[1] = 0 but this is nonsense for entries of a partition
  @req all(>(0),v) "Entries of v must be positive"

  # For safety
  @req all(>(0),mu) "Entries of mu must be positive"
  res = Partition{T}[]

  if isempty(v)
    return res
  end

  if m == 0
    return [ Partition{T}[] ]
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

  for n=nmin:nmax
    append!(res, partitions(m, n, v, mu))
  end

  return res

end

@doc raw"""
    function partitions(m::T, v::Vector{T}) where T<:IntegerUnion
  
All partitions of a non-negative integer `m` where each part is an element in the vector
`v`. We assume (without loss of generality) that the entries in `v` are strictly increasing. 

# Example
We compute the number of partitions of 100 where the parts are from {1, 2, 5, 10, 20, 50}:
```jldoctest
julia> length(partitions(100, [1,2,5,10,20,50]))
4562
```
# Algorithm
We use the function `partitions(m,n,v,mu)`, looping over the number of possible parts of
partitions.
"""
function partitions(m::T, v::Vector{T}) where T<:IntegerUnion

  @req m >= 0 "m >= 0 required"
  @req all([v[i] < v[i+1] for i in 1:length(v)-1]) "v must be strictly increasing"
  @req all(>(0),v) "Entries of v must be positive"

  res = Partition{T}[]

  if isempty(v)
    return res
  end

  if m == 0
    return [ Partition{T}[] ]
  end

  # We will loop over the number of parts.
  # We first determine the minimal and maximal number of parts.
  r = length(v)
  nmin = div(m, v[r])
  nmax = div(m, v[1])

  # Set the maximum multiplicity (equal to nmax above)
  mu = [ nmax for i in 1:r ]
  
  for n=nmin:nmax
    append!(res, partitions(m, n, v, mu))
  end

  return res

end

################################################################################
# Relations
################################################################################

@doc raw"""
    dominates(lambda::Partition, mu::Partition)

The **dominance order** on partitions is the partial order ``⊵`` defined by ``λ ⊵ μ`` if and
only if ``λ₁ + … + λᵢ ≥ μ₁ + … + μᵢ`` for all i. If ``λ ⊵ μ`` one says that ``λ``
**dominates** ``μ``. This function returns true if and only if `lambda` dominates `mu`.

Note that whereas the lexicographic ordering is a total ordering, the dominance ordering is
not.

# Examples
```jldoctest
julia> dominates( Partition(3,1), Partition(2,2) )
true

julia> dominates( Partition(4,1), Partition(3,3) )
false
```

# Remarks
[Knu11](@cite) says **majorizes** instead of **dominates** and uses the symbol ``⪰`` instead
of ``⊵``.

# References
1. [Ful97](@cite), page 26
2. [Knu11](@cite), Section 7.2.1.4, Exercise 54 (page 412)
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
# Operations
################################################################################
@doc raw"""
    conjugate(lambda::Partition{T}) where T<:IntegerUnion

The **conjugate** of a partition is obtained by considering its Young diagram
(see [Tableaux](@ref)) and then flipping it along its main diagonal.

# Examples
```jldoctest
julia> conjugate(Partition(8,8,8,7,2,1,1))
[7, 5, 4, 4, 4, 4, 4, 3]
```
# References
1. [Ful97](@cite), page 2
2. [Knu11](@cite), Section 7.2.1.4, page 394.
"""
function conjugate(lambda::Partition{T}) where T<:IntegerUnion
  if isempty(lambda)
    return copy(lambda)
  end

  mu = zeros(T, lambda[1])

  for i = 1:length(lambda)
    for j = 1:lambda[i]
      mu[j] += 1
    end
  end

  return Partition(mu)
end
