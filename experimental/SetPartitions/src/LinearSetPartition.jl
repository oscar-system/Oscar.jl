"""
    LinearSetPartition

LinearSetPartition represents a linear combination of set-partitions. See Chapter 5 in [Gro20](@cite).
"""
struct LinearSetPartition{S <: AbstractPartition, T <: RingElement}
    coefficients :: Dict{S, T}

    function LinearSetPartition{S, T}(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}
        key_types = Set(typeof(k) for k in keys(coeffs))
        @req length(key_types) <= 1 "keys in the coefficients dictionary must be of the same subtype of AbstractPartition"
        
        return new(simplify_operation_zero(coeffs))
    end
end 

"""
    linear_partition(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}

Initialize LinearSetPartition object, while simplifying the term.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d") 
(Univariate polynomial ring in d over QQ, d)

julia> linear_partition(Dict(set_partition([1, 2], [1, 1]) => S(4), set_partition([1, 1], [1, 1]) => 4*d))
LinearSetPartition{SetPartition, QQPolyRingElem}(Dict{SetPartition, QQPolyRingElem}(
    SetPartition([1, 2], [1, 1]) => 4, SetPartition([1, 1], [1, 1]) => 4*d))
"""
function linear_partition(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}
    return LinearSetPartition{S, T}(coeffs)
end

"""
    coefficients(p::LinearSetPartition)

Get the constructor field `coefficients`, the partition term, in form of a `Dict` from
`SetPartition` to coefficient.
""" 
function coefficients(p::LinearSetPartition)
    return p.coefficients
end

function hash(p::LinearSetPartition, h::UInt)
    return hash(coefficients(p), h) 
end

function ==(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
        return coefficients(p) == coefficients(q)
end

function deepcopy_internal(p::LinearSetPartition, stackdict::IdDict)
    if haskey(stackdict, p)
        return stackdict[p]
    end
    q = linear_partition(deepcopy_internal(coefficients(p), stackdict))
    stackdict[p] = q
    return q
end

"""
    linear_partition(term::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }

Return a `LinearSetPartition` generated from the Vector `term` of 2-tuples,
where the first element in the tuple is a `RingElement` and the second
a `SetPartition`. Furthermore simplify the term before initializing
the `LinearSetPartition` object with the corresponding dict.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in d over QQ, d)

julia> linear_partition([(set_partition([1, 1], [1, 1]), S(4)), (set_partition([1, 1], [1, 1]), 4*d)])
LinearSetPartition{SetPartition, QQPolyRingElem}(Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 1], [1, 1]) => 4*d + 4))
```
"""
function linear_partition(term::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }
    return linear_partition(Dict{S, T}(x[1] => x[2] for x in simplify_operation(term)))
end

function +(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement } 

    result = deepcopy(coefficients(p))

    for i in pairs(coefficients(q))
        result[i[1]] = get(result, i[1], 0) + i[2]
    end
    
    return linear_partition(result)
end

function *(a::RingElement, p::LinearSetPartition{S, T}) where
        { S <: AbstractPartition, T <: RingElement }
    result = Dict{S, T}()

    for (i, n) in pairs(coefficients(p))
        result[i] = a * n
    end
    return linear_partition(result)
end

"""
    compose(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}, d::T) where 
        { S <: AbstractPartition, T <: RingElement }

Multiply each coefficient from `p` with each coefficient from `q`
as well as perform a composition between each SetPartition part 
in `p` and `q` and return the result.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in d over QQ, d)

julia> a = linear_partition([(set_partition([1, 2], [1, 1]), S(4)), (set_partition([1, 1], [1, 1]), 4*d)])
LinearSetPartition{SetPartition, QQPolyRingElem}(Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 2], [1, 1]) => 4, SetPartition([1, 1], [1, 1]) => 4*d))

julia> compose(a, a, d)
LinearSetPartition{SetPartition, QQPolyRingElem}(Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 2], [1, 1]) => 16*d + 16, SetPartition([1, 1], [1, 1]) => 16*d^2 + 16*d))
```
"""
function compose(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}, d::T) where 
        { S <: AbstractPartition, T <: RingElement }
    result = Dict{S, T}()
    
    for i in pairs(coefficients(p))
        for ii in pairs(coefficients(q))
            (composition, loop) = compose_count_loops(i[1], ii[1])
            new_coefficient = i[2] * ii[2] * (d^loop)
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end
    return linear_partition(result)
end 

"""
    tensor_product(p::LinearSetPartition{P, R}, q::LinearSetPartition{P, R}) where 
        { P <: AbstractPartition, R <: RingElement }

Return the tensor product of `p` and `q`.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in d over QQ, d)

julia> a = linear_partition([(set_partition([1, 2], [1, 1]), S(4)), (set_partition([1, 1], [1, 1]), 4*d)])
LinearSetPartition{SetPartition, QQPolyRingElem}(Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 2], [1, 1]) => 4, SetPartition([1, 1], [1, 1]) => 4*d))

julia> tensor_product(a, a)
LinearSetPartition{SetPartition, QQPolyRingElem}(Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 1, 2, 2], [1, 1, 2, 2]) => 16*d^2, SetPartition([1, 2, 3, 3], [1, 1, 3, 3]) => 16*d, SetPartition([1, 2, 3, 4], [1, 1, 3, 3]) => 16, SetPartition([1, 1, 2, 3], [1, 1, 2, 2]) => 16*d))
```
"""
function tensor_product(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

    result = Dict{S, T}()
    
    for i in pairs(coefficients(p))
        for ii in pairs(coefficients(q))
            composition = tensor_product(i[1], ii[1])
            new_coefficient = i[2] * ii[2]
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end

    return linear_partition(result)
end 

function ⊗(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    return tensor_product(p, q)
end

function -(p::LinearSetPartition)
    return (-1 * p)
end

function -(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    return p + (-q)
end
