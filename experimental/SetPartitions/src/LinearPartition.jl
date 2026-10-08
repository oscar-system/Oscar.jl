"""
    LinearPartition{S <: AbstractPartition, T <: RingElem}

`LinearPartition` represents a linear combination of partitions of type `S`
with coefficients of type `T`. 

See for example Chapter 5 in [Gro20](@cite) for further information on linear 
combinations of partitions.
"""
struct LinearPartition{S <: AbstractPartition, T <: RingElem}
    base_ring :: Ring  # parent_type(T)
    coefficients :: Dict{S, T}

    function LinearPartition{S, T}(ring::Ring, coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElem}
        @req isconcretetype(S) "Linear combinations are only defined for concrete subtypes of AbstractPartition"
        @req ring isa parent_type(T) "Ring type is not the parent of the coefficient type"
        for c in values(coeffs)
            @req (parent(c) == ring) "Coefficient does not belong to the base ring"
        end
        return new(ring, simplify_operation_zero(coeffs))
    end
end 

"""
    linear_partition(ring::Ring, coeffs::Dict{S, <: Any}) where {S <: AbstractPartition}

Construct a linear partition over `ring` with coefficients `coeffs`.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, :d) 
(Univariate polynomial ring in d over QQ, d)

julia> linear_partition(S, Dict(set_partition([1, 2], [1, 1]) => 4, set_partition([1, 1], [1, 1]) => 4*d))
LinearPartition{SetPartition, QQPolyRingElem}(Univariate polynomial ring in d over QQ, Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 2], [1, 1]) => 4, SetPartition([1, 1], [1, 1]) => 4*d))
```
"""
function linear_partition(ring::Ring, coeffs::Dict{S, <: Any}) where {S <: AbstractPartition}
    ring_coeffs = Dict(p => ring(c) for (p, c) in coeffs)
    return LinearPartition{S, elem_type(ring)}(ring, ring_coeffs)
end

"""
    linear_partition(ring::Ring, terms::Vector{Tuple{S, <: Any}}) where {S <: AbstractPartition}

Return a `LinearPartition` generated from the vector `term` of 2-tuples,
where the first element in the tuple is an `AbstractPartition` and the second
a `RingElem`. The `ring` argument defines the base ring into which the second 
elements of the tuples are converted. Furthermore simplify the term before initializing
the `LinearPartition` object with the corresponding dict.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, :d)
(Univariate polynomial ring in d over QQ, d)

julia> linear_partition(S, [(set_partition([1, 1], [1, 1]), 4), (set_partition([1, 1], [1, 1]), 4*d)])
LinearPartition{SetPartition, QQPolyRingElem}(Univariate polynomial ring in d over QQ, Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 1], [1, 1]) => 4*d + 4))
```
"""
function linear_partition(ring::Ring, terms::Vector{Tuple{S, <: Any}}) where {S <: AbstractPartition}
    simpl_terms = simplify_operation([(p, ring(c)) for (p, c) in terms])
    return linear_partition(ring, Dict{S, elem_type(ring)}(simpl_terms))
end

"""
    base_ring(p::LinearPartition{S, T}) where {S <: AbstractPartition, T <: RingElem}

Return the underlying coefficient ring of `p`.
""" 
function base_ring(p::LinearPartition{S, T}) where {S <: AbstractPartition, T <: RingElem}
    return p.base_ring::parent_type(T) 
end

function base_ring_type(::Type{LinearPartition{S, T}}) where {S <: AbstractPartition, T <: RingElem} 
    return parent_type(T)
end

"""
    coefficients(p::LinearPartition)

Return the coefficients of `p` as dictionary from partitions to elements 
of the underlying ring.
""" 
function coefficients(p::LinearPartition)
    return p.coefficients
end

function hash(p::LinearPartition, h::UInt)
    return hash(base_ring(p), hash(coefficients(p), h)) 
end

function ==(p::LinearPartition, q::LinearPartition)
    return base_ring(p) == base_ring(q) && coefficients(p) == coefficients(q)
end

function deepcopy_internal(p::LinearPartition, stackdict::IdDict)
    if haskey(stackdict, p)
        return stackdict[p]
    end
    q = linear_partition(
            base_ring(p),
            deepcopy_internal(coefficients(p), stackdict))
    stackdict[p] = q
    return q
end

function +(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElem } 
    @req base_ring(p) == base_ring(q) "Linear partitions are defined over different base rings"
    result = deepcopy(coefficients(p))
    for i in pairs(coefficients(q))
        result[i[1]] = get(result, i[1], 0) + i[2]
    end
    return linear_partition(base_ring(p), result)
end

function *(a::RingElement, p::LinearPartition{S, T}) where
        { S <: AbstractPartition, T <: RingElem }
    a = base_ring(p)(a)
    result = Dict{S, T}()
    for (i, n) in pairs(coefficients(p))
        result[i] = a * n
    end
    return linear_partition(base_ring(p), result)
end

"""
    compose(p::LinearPartition{S, T}, q::LinearPartition{S, T}, d::T) where 
        { S <: AbstractPartition, T <: RingElem }

Return the composition between `p` and `q`.

The composition is obtained by multiplying each coefficient
of `p` with each coefficient of `q` and composing the corresponding
partitions. The `RingElem` parameter `d` multiplies each 
coefficient based on the number of loops identified during the  
composition.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, :d)
(Univariate polynomial ring in d over QQ, d)

julia> a = linear_partition(S, [(set_partition([1, 2], [1, 1]), 4), (set_partition([1, 1], [1, 1]), 4*d)])
LinearPartition{SetPartition, QQPolyRingElem}(Univariate polynomial ring in d over QQ, Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 2], [1, 1]) => 4, SetPartition([1, 1], [1, 1]) => 4*d))

julia> compose(a, a, d)
LinearPartition{SetPartition, QQPolyRingElem}(Univariate polynomial ring in d over QQ, Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 2], [1, 1]) => 16*d + 16, SetPartition([1, 1], [1, 1]) => 16*d^2 + 16*d))
```
"""
function compose(p::LinearPartition{S, T}, q::LinearPartition{S, T}, d::T) where 
        { S <: AbstractPartition, T <: RingElem }
    @req base_ring(p) == base_ring(q) "Linear partitions are defined over different base rings"
    result = Dict{S, T}()
    for i in pairs(coefficients(p))
        for ii in pairs(coefficients(q))
            (composition, loop) = compose_count_loops(i[1], ii[1])
            new_coefficient = i[2] * ii[2] * (d^loop)
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end
    return linear_partition(base_ring(p), result)
end 

"""
    tensor_product(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElem }

Return the tensor product of `p` and `q`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> S, d = polynomial_ring(QQ, :d)
(Univariate polynomial ring in d over QQ, d)

julia> a = linear_partition(S, [(set_partition([1, 2], [1, 1]), 4), (set_partition([1, 1], [1, 1]), 4*d)])
LinearPartition{SetPartition, QQPolyRingElem}(Univariate polynomial ring in d over QQ, Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 2], [1, 1]) => 4, SetPartition([1, 1], [1, 1]) => 4*d))

julia> tensor_product(a, a)
LinearPartition{SetPartition, QQPolyRingElem}(Univariate polynomial ring in d over QQ, Dict{SetPartition, QQPolyRingElem}(SetPartition([1, 1, 2, 2], [1, 1, 2, 2]) => 16*d^2, SetPartition([1, 2, 3, 3], [1, 1, 3, 3]) => 16*d, SetPartition([1, 2, 3, 4], [1, 1, 3, 3]) => 16, SetPartition([1, 1, 2, 3], [1, 1, 2, 2]) => 16*d))
```
"""
function tensor_product(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElem }
    @req base_ring(p) == base_ring(q) "Linear partitions are defined over different base rings"
    result = Dict{S, T}()
    for i in pairs(coefficients(p))
        for ii in pairs(coefficients(q))
            composition = tensor_product(i[1], ii[1])
            new_coefficient = i[2] * ii[2]
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end
    return linear_partition(base_ring(p), result)
end 

function ⊗(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElem }
    return tensor_product(p, q)
end

function -(p::LinearPartition)
    return (-1 * p)
end

function -(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElem }
    return p + (-q)
end
