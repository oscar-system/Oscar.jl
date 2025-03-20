################################################################################
# Multipartitions.
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
################################################################################

export Multipartition, multipartitions, num_multipartitions

"""
    Multipartition{T} <: AbstractArray{Partition{T},1}

Multipartitions are generalizations of partitions. An r-component **multipartition** of an integer n is an r-tuple of partitions λ¹, λ², …, λʳ where each λⁱ is a partition of some integer nᵢ ≥ 0 and the nᵢ sum to n.

Multipartitions are implemented as a subtype of 1-dimensional arrays of partitions. You can use smaller integer types to increase performance.

# Examples
```julia-repl
julia> P=Multipartition( [[2,1], [], [3,2,1]] )
Partition{Int64}[[2, 1], [], [3, 2, 1]]
julia> sum(P)
9
julia> P[2]
Int64[]
julia> P=Multipartition( Array{Int8,1}[[2,1], [], [3,2,1]] )
Partition{Int8}[[2, 1], [], [3, 2, 1]]
```

# References
1. Wikipedia, [Multipartition](https://en.wikipedia.org/wiki/Multipartition)
"""
struct Multipartition{T} <: AbstractArray{Partition{T},1}
    mp::Array{Partition{T},1}
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

# More convenient constructors
function Multipartition(mp::Array{Array{T,1},1}) where T<:Integer
    return Multipartition([Partition(p) for p in mp])
end

# This is only called when the empty array is part of mp (because then it's
# "Any" type and not of Integer type).
function Multipartition(mp::Array{Array{Any,1},1})
    return Multipartition([Partition(p) for p in mp])
end


"""
    sum(P::Multipartition{T}) where T<:Integer

If P is a multipartition of the integer n, this function returns n.
"""
function Base.sum(P::Multipartition{T}) where T<:Integer
    s = zero(T)
    for i=1:length(P)
        s += sum(P[i])
    end
    return s
end


"""
    multipartitions(n::T, r::Integer) where T<:Integer

A list of all r-component multipartitions of n. The algorithm is recursive and  based on [`partitions(::Integer)`](@ref).

# Example
```julia-repl
julia> multipartitions(2,2)
5-element Vector{Multipartition{Int64}}:
 Partition{Int64}[[], [2]]
 Partition{Int64}[[], [1, 1]]
 Partition{Int64}[[1], [1]]
 Partition{Int64}[[2], []]
 Partition{Int64}[[1, 1], []]
```
"""
function multipartitions(n::T, r::Integer) where T<:Integer
    #Argument checking
    n >= 0 || throw(ArgumentError("n >= 0 required"))
    r >= 1 || throw(ArgumentError("r >= 1 required"))

    #This will be the list of multipartitions
    MP = Multipartition{T}[]

    #We will go through all compositions of n into r parts, and for each such #composition, we collect all the partitions for each of the components of
    #the composition.
    #We create the compositions here in place for efficiency.

    #recursively produces all Integer Arrays p of length r such that the sum of all the Elements equals n. Then calls recMultipartitions!
    function recP!(p::Array{T,1}, i::T, n::T) #where T<:Integer
        if i==length(p) || n==0
            p[i] = n
            recMultipartitions!(fill(Partition(T[]),r), p, T(1))
        else
            for j=0:n
                p[i] = T(j)
                recP!(copy(p), T(i+1), T(n-j))
            end
        end
    end

    #recursively produces all multipartitions such that the i-th partition sums up to p[i]
    function recMultipartitions!(mp::Array{Partition{T},1}, p::Array{T,1}, i::T) #where T<:Integer
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
    num_multipartitions(n::Int, k::Int)

The number of multipartitions of ``n`` into ``k`` parts is equal to
```math
\sum_{a=1}^k {k \choose a} \sum_{λ} p(λ₁) p(λ₂) ⋯ p(λ_a) \;,
```
where the second sum is over all [compositions](@ref Compositions) ``λ`` of ``n`` into ``a`` parts. I found this formula in the Proof of Lemma 2.4 in Craven (2006).

# References
1. Craven, D. (2006). The Number of t-Cores of Size n. [http://web.mat.bham.ac.uk/D.A.Craven/docs/papers/tcores0608.pdf](http://web.mat.bham.ac.uk/D.A.Craven/docs/papers/tcores0608.pdf)
"""
function num_multipartitions(n::Int, k::Int)

    z = ZZ(0)

    # Special cases
    if n==0
        return ZZ(k)
    end

    for a=1:k
        w = ZZ(0)
        for λ in compositions(n,a)
            w += prod([num_partitions(λ[i]) for i=1:length(λ)])
        end
        z += binomial(k,a)*w
    end

    return z

end
