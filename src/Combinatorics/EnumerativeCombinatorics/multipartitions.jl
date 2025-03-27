################################################################################
# Multipartitions.
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
# by Ulrich Thiel and OSCAR-ified by Morgan Rodgers.
################################################################################




"""
    multipartition(mp::Vector{Partition{T}}) where T <: IntegerUnion
    multipartition(mp::Vector{Vector{T}}) where T <: IntegerUnion

Return the multipartition given by the vector `mp` of integer sequences `mp`
(which are interpreted as integer partitions) as an object of type
`Multipartition{T}`.

The element type `T` may be optionally specified, see also the examples below.

# Examples
```jldoctest
julia> P = multipartition([[2,1], [], [3,2,1]])
Partition{Int64}[[2, 1], [], [3, 2, 1]]
julia> sum(P)
9
julia> P[2]
Empty partition
julia> P = multipartition(Vector{Int8}[[2,1], [], [3,2,1]])
Partition{Int8}[[2, 1], [], [3, 2, 1]]
```

# References
1. Wikipedia, [Multipartition](https://en.wikipedia.org/wiki/Multipartition)
"""
function multipartition(mp::Vector{Partition{T}}) where T <: IntegerUnion
  return Multipartition(mp)
end

function multipartition(mp::Vector{Vector{T}}) where T <: IntegerUnion
  return Multipartition([partition(p) for p in mp])
end

# This is only called when the empty array is part of mp (because then it's
# "Any" type and not of Integer type).
function multipartition(mp::Vector{Vector{Any}})
  return Multipartition([partition(p) for p in mp])
end

function Base.show(io::IO, ::MIME"text/plain", MP::Multipartition)
  print(io, MP.mp)
end

function Base.size(MP::Multipartition)
  return size(MP.mp)
end

function Base.length(MP::Multipartition)
  return length(MP.mp)
end

function Base.getindex(MP::Multipartition, i::Int)
  return getindex(MP.mp,i)
end


"""
    sum(P::Multipartition{T}) where T<:IntegerUnion

If `P` is a multipartition of the integer n, this function returns n.
"""
function Base.sum(MP::Multipartition{T}) where T<:IntegerUnion
  return sum(sum, MP.mp)
end


"""
    multipartitions(n::T, r::IntegerUnion)  where T<:IntegerUnion

A list of all `r`-component multipartitions of `n`, as elements of type Multipartitoon{T}.
The algorithm is recursive and  based on [`partitions(::IntegerUnion)`](@ref).

# Example
```jldoctest
julia> multipartitions(2,2)
5-element Vector{Multipartition{Int64}}:
 Partition{Int64}[[], [2]]
 Partition{Int64}[[], [1, 1]]
 Partition{Int64}[[1], [1]]
 Partition{Int64}[[2], []]
 Partition{Int64}[[1, 1], []]
```
"""
function multipartitions(n::T, r::IntegerUnion) where T<:IntegerUnion
  #Argument checking
  @req n >= 0 "n >= 0 required"
  @req r >= 1 "r >= 1 required"

  #This will be the list of multipartitions
  MP = Multipartition{T}[]

  #We will go through all compositions of n into r parts, and for each such #composition, we collect all the partitions for each of the components of
  #the composition.
  #We create the compositions here in place for efficiency.

  #recursively produces all Integer Vectors p of length r such that the sum of all the Elements equals n. Then calls recMultipartitions!
  function recP!(p::Vector{T}, i::T, n::T) #where T<:IntegerUnion
    if i==length(p) || n==0
      p[Int(i)] = n
      recMultipartitions!(fill(Partition(T[]),r), p, T(1))
    else
      for j=0:n
        p[i] = T(j)
        recP!(copy(p), T(i+1), T(n-j))
      end
    end
  end

  #recursively produces all multipartitions such that the i-th partition sums up to p[i]
  function recMultipartitions!(mp::Vector{Partition{T}}, p::Vector{T}, i::T) #where T<:IntegerUnion
    if i == length(p)
      for q in partitions(p[i])
        mp[i] = q
        push!(MP, Multipartition{T}(copy(mp)))
      end
    else
      for q in partitions(p[i])
        mp[i] = q
        recMultipartitions!(copy(mp), p, T(i+1))
      end
    end
  end

  recP!(zeros(T,r), T(1), n)
  return MP
end


@doc raw"""
    number_of_multipartitions(n::Int, k::Int)

The number of multipartitions of ``n`` into ``k`` parts is equal to
```math
\sum_{a=1}^k {k \choose a} \sum_{λ} p(λ₁) p(λ₂) ⋯ p(λ_a) \;,
```
where the second sum is over all [compositions](@ref Compositions) ``λ`` of ``n`` into ``a`` parts. I found this formula in the Proof of Lemma 2.4 in Craven (2006).

# References
1. Craven, D. (2006). The Number of t-Cores of Size n. [http://web.mat.bham.ac.uk/D.A.Craven/docs/papers/tcores0608.pdf](http://web.mat.bham.ac.uk/D.A.Craven/docs/papers/tcores0608.pdf)
"""
function number_of_multipartitions(n::Int, k::Int)

  # Special cases
  n == 0 && return ZZ(k)

  z = ZZ(0)
  for a in 1:k
    w = ZZ(0)
    for lambda in compositions(n,a)
      w += prod([number_of_partitions(lambda[i]) for i in 1:length(lambda)])
    end
    z += binomial(k,a)*w
  end

  return z

end
