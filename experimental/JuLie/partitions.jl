################################################################################
# Partitions.
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
# by Ulrich Thiel and OSCAR-ified by Claudia He Yun and Matthias Zach.
################################################################################

export Partition

export ascending_partitions
export conjugate
export dominates
export getindex_safe
export num_partitions
export partitions

@doc Markdown.doc"""
    Partition{T} <: AbstractVector{T}

A **partition** of an integer ``n ≥ 0`` is a decreasing sequence
``λ=(λ₁,…,λᵣ)`` of positive integers ``λᵢ`` whose sum is equal to ``n``. The
``λᵢ`` are called the **parts** of the partition. We encode a partition as an
array with elements ``λᵢ``. You may increase performance by using smaller
integer types, see the examples below. For efficiency, the `Partition`
constructor does not check whether the given array is in fact a partition,
i.e. a decreasing sequence.

# Examples
```jldoctest
julia> P=Partition([3,2,1]) #The partition 3+2+1 of 6
[3, 2, 1]

julia> sum(P) #The sum of the parts.
6

julia> P[1] #First component
3

julia> P=Partition(Int8[3,2,1]) #Same partition but using 8 bit integers
Int8[3, 2, 1]
```

# Remarks

* Usually, ``|λ| ≔ n`` is called the **size** of ``λ``. In Julia, the function
`size` for arrays already exists and returns the *dimension* of an array.
Instead, you can use the Julia function `sum` to get the sum of the parts.

* There is no performance impact by using an own type for partitions rather
than simply using arrays—and this is of course much cleaner. The
implementation of a subtype of AbstractArray is explained in the [Julia
documentation](https://docs.julialang.org/en/v1/manual/interfaces/#man-
interface-array).

# References
1. Wikipedia, [Partition (number theory)](https://en.wikipedia.org/wiki/Partition_(number_theory))
"""
struct Partition{T} <: AbstractVector{T}
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

# The empty array is of "Any" type, and this is stupid. We want it here
# to get it into the default type Int64. This constructor is also called by
# MultiPartition, and this casts the whole array into "Any" whenever there's
# the empty partition inside.
function Partition(p::Vector{Any})
  return Partition(Vector{Int64}(p))
end

function Base.copy(P::Partition{T}) where T<:IntegerUnion
  return Partition{T}(copy(P.p))
end

@doc Markdown.doc"""
    getindex_safe(P::Partition, i::IntegerUnion)

In algorithms involving partitions it is sometimes convenient to be able to
access parts beyond the length of the partition, and then you want to get zero
instead of an error. This function is a shortcut for
```
return (i>length(P.p) ? 0 : getindex(P.p,i))
```
If you are sure that `P[i]` exists, use `getindex` because this will be faster.
"""
function getindex_safe(P::Partition, i::IntegerUnion)
  return (i>length(P.p) ? 0 : getindex(P.p,Int(i)))
end


@doc Markdown.doc"""
    num_partitions(n::IntegerUnion)

The number of integer partitions of the integer ``n ≥ 0``. Uses the function
from FLINT, which is very fast.

# References
1. The On-Line Encyclopedia of Integer Sequences, [A000041](https://oeis.org/A000041)
2. FLINT, [Number of partitions](http://flintlib.org/doc/arith.html?highlight=partitions#number-of-partitions)
"""
function num_partitions(n::IntegerUnion)
  @req n >= 0 "n >= 0 required"
  n = ZZ(n)
  z = ZZ()
  ccall((:arith_number_of_partitions, Nemo.libflint), Cvoid, (Ref{fmpz}, Culong), z, UInt(n))
  return z
end


@doc Markdown.doc"""
    num_partitions(n::IntegerUnion, k::IntegerUnion)

The number of integer partitions of the integer ``n ≥ 0`` into ``k ≥ 0``
parts. The implementation uses a recurrence relation.

# References
1. The On-Line Encyclopedia of Integer Sequences, [A008284](https://oeis.org/A008284)
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

    # Otherwise, use recurrence
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


@doc Markdown.doc"""
    partitions(n::IntegerUnion)

A list of all partitions of an integer ``n ≥ 0``, produced in
lexicographically *descending* order. This ordering is like in Sage, but
opposite to GAP. You can apply the function `reverse` to reverse the
order. As usual, you may increase performance by using smaller integer types.
The algorithm used is "Algorithm ZS1" by Zoghbi & Stojmenovic (1998); see
[ZS98](@cite).

# Examples
```jldoctest
julia> partitions(Int8(4))
5-element Vector{Partition{Int8}}:
 Int8[4]
 Int8[3, 1]
 Int8[2, 2]
 Int8[2, 1, 1]
 Int8[1, 1, 1, 1]
```
"""
function partitions(n::IntegerUnion)

  #Argument checking
  @req n >= 0 "n >= 0 required"

  # Use type of n
  T = typeof(n)

  # Some trivial cases
  if n==0
    return Partition{T}[ Partition{T}([]) ]
  elseif n==1
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


@doc Markdown.doc"""
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



@doc Markdown.doc"""
    partitions(m::IntegerUnion, n::IntegerUnion, l1::IntegerUnion, l2::IntegerUnion; only_distinct_parts::Bool = false)

A list of all partitions of an integer ``m ≥ 0`` into ``n ≥ 0`` parts with
lower bound ``l1 ≥ 0`` and upper bound ``l2 ≥ l1`` for the parts. There are
two choices for the parameter `only_distinct_parts`:
* `false`: no further restriction (*default*);
* `true`: only distinct parts.
The partitions are produced in *decreasing* order.

The algorithm used is "parta" in [RJ76](@cite), de-gotoed from old ALGOL code by E. Thiel!
"""
function partitions(m::T, n::IntegerUnion, l1::IntegerUnion, l2::IntegerUnion; only_distinct_parts::Bool = false) where T <: IntegerUnion

  # Note that we are considering partitions of m here. I would switch m and n
  # but the algorithm was given like that and I would otherwise confuse myself
  # implementing it.

  #Argument checking
  @req m >= 0 "m >= 0 required"
  @req n >= 0 "n >= 0 required"

  # Use type of n
  n = convert(T, n)

  # Some trivial cases
  if m == 0 && n == 0
    return Partition{T}[ Partition{T}([]) ]
  end

  if n == 0 || n > m
    return Partition{T}[]
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

      if i<n && m>1
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



@doc Markdown.doc"""
    partitions(m::IntegerUnion, n::IntegerUnion)

All partitions of an integer ``m ≥ 0`` into ``n ≥ 1`` parts (no further restrictions).
"""
function partitions(m::IntegerUnion, n::IntegerUnion)
  return partitions(m,n,1,m; only_distinct_parts = false)
end



@doc Markdown.doc"""
    partitions(mu::Vector{IntegerUnion}, m::IntegerUnion, v::Vector{IntegerUnion}, n::IntegerUnion)

All partitions of an integer ``m >= 0`` into ``n >= 1`` parts, where each part
is an element in ``v`` and each ``v[i]`` occurs a maximum of ``mu[i]`` times.
The partitions are produced in    *decreasing* order. The algorithm used is a
de-gotoed version (by E. Thiel!) of algorithm "partb" in [RJ76](@cite).

# Remark
The original algorithm lead to BoundsErrors, since r could get smaller than 1.
Furthermore x and y are handled as arrays with an infinite length. After
finding all valid partitions, the algorithm will continue searching for
partitions of length n+1. We thus had to add a few additional checks and
interruptions. Done by T. Schmit.
"""
function partitions(mu::Vector{S}, m::T, v::Vector{S}, n::IntegerUnion) where {T<:IntegerUnion, S<:IntegerUnion}
  @req length(mu)==length(v) "mu and v should have the same length"
  @req m >= 0 "m >= 0 required"
  @req n >= 1 "n >= 1 required"

  if isempty(mu)
    return Partition{T}[]
  end

  r = length(v)
  j = 1
  k = mu[1]
  ll = v[1]
  x = zeros(Int8, n+1)
  y = zeros(Int8, n+1)
  ii = zeros(Int8, n)
  i_1 = 0
  P = Partition{T}[]

  num = 0
  gotob2 = false
  gotob1 = true

  for i = n:-1:1
    x[i] = ll
    y[i] = ll
    k = k - 1
    m = m - ll
    if k == 0
      if j == r
        return P
      end
      j = j + 1
      k = mu[j]
      ll = v[j]
    end
  end #for i

  lr = v[r]
  ll = v[1]

  if m < 0 || m > n * (lr - ll)
    return P
  end

  if m == 0
    push!(P,Partition{T}(x[1:n]))
    return P
  end

  i = 1
  m = m + y[1]

  while gotob1 == true
    if !gotob2
      for j = mu[r]:-1:1
        if m<=lr
          gotob2 = true
          break
        end
        x[i] = lr
        ii[i] = r - 1
        i = i + 1
        m = m - lr + y[i]
      end #for j

      if !gotob2
        r = r - 1
      end

      gotob2 = true
    end #if

    if gotob2
      while r>0 && v[r] > m
        r = r - 1
      end

      if r == 0
        break
      end

      lr = v[r]
      if m == lr
        x[i] = lr
        if i<=n
          push!(P, Partition{T}(x[1:n]))
        else
          break
        end

        r = r - 1
        if r==0
          break
        end
        lr = v[r]
      end #if

      k = y[i]
      if lr > k && m - lr <= (n-i)*(lr - ll)
        gotob2 = false
        continue
      else
        x[i] = k
      end #if
      i_1 = i - 1
      for i_0 = i_1:-1:1
        i = i_0
        r = ii[i]
        lr = v[r]
        m = m + x[i] - k
        k = y[i]
        if lr > k && m - lr <= (n-i)*(lr-ll)
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



@doc Markdown.doc"""
    dominates(lambda::Partition, mu::Partition)

The **dominance order** on partitions is the partial order ``⊵`` defined by
``λ ⊵ μ`` if and only if ``λ₁ + … + λᵢ ≥ μ₁ + … + μᵢ`` for all i. This
function returns true if ``λ ⊵ μ``.

# References
1. Wikipedia, [Dominance order](https://en.wikipedia.org/wiki/Dominance_order)
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


@doc Markdown.doc"""
    conjugate(lambda::Partition{T}) where T<:IntegerUnion

The **conjugate** of a partition is obtained by considering its Young diagram
(see [Tableaux](@ref)) and then flipping it along its main diagonal.

# References
1. Wikipedia, [Partition (number theory)](https://en.wikipedia.org/wiki/Partition_(number_theory)#Conjugate_and_self-conjugate_partitions)
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
