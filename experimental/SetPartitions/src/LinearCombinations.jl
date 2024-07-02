"""
Fragen an Nicolas:
- name LinearSetPartition okay? (davor war er nur LinearPartition)
- Kommentar zeile 140
"""

"""
    LinearSetPartition{S <: AbstractPartition, T <: RingElement}

Initialize LinearSetPartition object, while simplifying
the term. 

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> LinearSetPartition(Dict(set_partition([1, 2], [1, 1]) => S(4), set_partition([1, 1], [1, 1]) => 4*d))
LinearSetPartition(Dict(set_partition([1, 2], [1, 1]) => S(4), set_partition([1, 1], [1, 1]) => 4*d))
```
"""
struct LinearSetPartition{S <: AbstractPartition, T <: RingElement}
    coefficients :: Dict{S, T}

    function LinearSetPartition{S, T}(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}
        new(simplify_operation_zero(coeffs))
    end
end 

function LinearSetPartition(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}
    LinearSetPartition{S, T}(coeffs)
end

function linear_partition(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}
    LinearSetPartition{S, T}(coeffs)
end

"""
    get_term(p::LinearSetPartition{S, T}) where { S <: AbstractPartition, T <: RingElement }

Get the constructor field `coefficients`, the partition term, in form of a `Dict` from
`SetPartition` to coefficient.
"""
function get_term(p::LinearSetPartition{S, T}) where { S <: AbstractPartition, T <: RingElement }
    p.coefficients
end

function hash(p::LinearSetPartition, h::UInt)
    return hash(p.coefficients, h)
end

function ==(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
        get_term(p) == get_term(q)
end

function deepcopy_internal(p::LinearSetPartition, stackdict::IdDict)
    if haskey(stackdict, p)
        return stackdict[p]
    end
    q = LinearSetPartition(deepcopy_internal(get_term(p), stackdict))
    stackdict[p] = q
    return q
end

"""
    LinearSetPartition(term::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }

Return LinearSetPartition object generated from the Vector `term` of 2-tuples,
where the first element in the tuple is a RingElement and the second
a SetPartition. Furthermore simplify the term before initializing
the LinearSetPartition object with the corresponding dict.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> LinearSetPartition([(set_partition([1, 1], [1, 1]), S(4)), (set_partition([1, 1], [1, 1]), 4*d)])
LinearSetPartition(Dict(set_partition([1, 1], [1, 1]) => 4 + 4*d))
```
"""
function LinearSetPartition(term::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }
    LinearSetPartition(Dict{S, T}(x[1] => x[2] for x in simplify_operation(term)))
end

function linear_partition(term::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }
    LinearSetPartition(term)
end



"""
    add(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

Return the addition between `p` and `q`.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = linear_partition([(set_partition([1, 2], [1, 1]), S(4)), 
(set_partition([1, 1], [1, 1]), 4*d)])
julia> add(a, a)
LinearSetPartition(Dict(set_partition([1, 2], [1, 1]) => S(8), 
SetPartition([1, 1], [1, 1]) => 8*d))
```
"""
function add(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement } 

    result = deepcopy(p)

    for i in pairs(get_term(q))
        get_term(result)[i[1]] = get(get_term(result), i[1], 0) + i[2]
    end
    
    LinearSetPartition(get_term(result))
end

function +(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    add(p, q)
end

"""
    scale(a::RingElement, p::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

Multiply each coefficient in `p` with `a` and return
the result.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> scale(0.5, linear_partition(set_partition([1, 2], [1, 1]) => S(8), 
set_partition([1, 1], [1, 1]) => 8*d))
LinearSetPartition(Dict(SetPartition([1, 2], [1, 1]) => S(4), 
SetPartition([1, 1], [1, 1]) => 4*d))
```
"""
function scale(a::RingElement, p::LinearSetPartition{S, T}) where # Doctest geht leider nicht weil QQPolyRing kein float erlaubt
        { S <: AbstractPartition, T <: RingElement }
    result = Dict{S, T}()

    for (i, n) in pairs(get_term(p))
        result[i] = a * n
    end
    LinearSetPartition(result)
end 

function *(a::RingElement, p::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    scale(a, p)
end

"""
    linear_composition(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}, d::T) where 
        { S <: AbstractPartition, T <: RingElement }

Multiply each coefficient from `p` with each coefficient from `q`
as well as perform a composition between each SetPartition part 
in `p` and `q` and return the result.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = linear_partition([(set_partition([1, 2], [1, 1]), S(4)), 
(SetPartition([1, 1], [1, 1]), 4*d)])
julia> linear_composition(a, a, d)
LinearSetPartition(Dict(SetPartition([1, 2], [1, 1]) => 16*d + 16, 
SetPartition([1, 1], [1, 1]) => 16*d^2 + 16*d))
```
"""
function linear_composition(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}, d::T) where 
        { S <: AbstractPartition, T <: RingElement }
    result = Dict{S, T}()
    
    for i in pairs(get_term(p))
        for ii in pairs(get_term(q))
            (composition, loop) = compose_count_loops(i[1], ii[1])
            new_coefficient = i[2] * ii[2] * (d^loop)
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end
    LinearSetPartition(result)
end 

"""
    linear_tensor_product(p::LinearSetPartition{P, R}, q::LinearSetPartition{P, R}) where 
        { P <: AbstractPartition, R <: RingElement }

Perform a tensor product similar

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = linear_partition([(set_partition([1, 2], [1, 1]), S(4)), 
(set_partition([1, 1], [1, 1]), 4*d)])
julia> linear_tensor_product(a, a)
LinearSetPartition(Dict(SetPartition([1, 1, 2, 2], [1, 1, 2, 2]) => 16*d^2, 
SetPartition([1, 2, 3, 3], [1, 1, 3, 3]) => 16*d, 
SetPartition([1, 2, 3, 4], [1, 1, 3, 3]) => S(16), 
SetPartition([1, 1, 2, 3], [1, 1, 2, 2]) => 16*d)))
```
"""
function linear_tensor_product(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

    result = Dict{S, T}()
    
    for i in pairs(get_term(p))
        for ii in pairs(get_term(q))
            composition = tensor_product(i[1], ii[1])
            new_coefficient = i[2] * ii[2]
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end

    LinearSetPartition(result)
    
end 

function ⊗(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    linear_tensor_product(p, q)
end

"""
    subtract(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

Perform a subtraction between `p` and `q`and return the result.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = linear_partition([(set_partition([1, 2], [1, 1]), S(4)), 
(SetPartition([1, 1], [1, 1]), 4*d)])
julia> subtract(a, a)
LinearSetPartition(Dict(SetPartition([1, 1], [1, 1]) => S(0)))
```
"""
function subtract(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    add(p, scale(-1, q))
end 

function -(p::LinearSetPartition{S, T}, q::LinearSetPartition{S, T}) where 
    { S <: AbstractPartition, T <: RingElement }
    subtract(p, q)
end 
