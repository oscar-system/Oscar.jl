"""
Sebvz todos:
- deepcopy
- dict zugriffe verbessern wenn möglich
- tests
- rename composition
"""


"""
    LinearPartition{S <: AbstractPartition, T <: RingElement}

Initialize LinearPartition object, while simplifying
the term. 

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> LinearPartition(Dict(SetPartition([1, 2], [1, 1]) => S(4), SetPartition([1, 1], [1, 1]) => 4*d))
LinearPartition(Dict(SetPartition([1, 2], [1, 1]) => S(4), SetPartition([1, 1], [1, 1]) => 4*d))
```
"""
struct LinearPartition{S <: AbstractPartition, T <: RingElement}
    coefficients :: Dict{S, T}

    function LinearPartition{S, T}(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}
        new(simplify_operation_zero(coeffs))
    end
end 

function LinearPartition(coeffs::Dict{S, T}) where {S <: AbstractPartition, T <: RingElement}
    LinearPartition{S, T}(coeffs)
end


# TODO hash

function ==(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    p == q
end

"""
    LinearPartition(terms::Vector{Tuple{S, T}}) where 
        { S <: AbstractPartition, T <: RingElement }

Return LinearPartition object generated from the Vector `terms` of 2-tuples,
where the first element in the tuple is a RingElement and the second
a SetPartition. Furthermore simplify the term before initializing
the LinearPartition object with the corresponding dict.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> LinearPartition([(SetPartition([1, 1], [1, 1]), S(4)), (SetPartition([1, 1], [1, 1]), 4*d)])
LinearPartition(Dict(SetPartition([1, 1], [1, 1]) => 4 + 4*d))
```
"""
function LinearPartition(terms::Vector{Tuple{S, T}}) where { S <: AbstractPartition, T <: RingElement }
    LinearPartition(Dict{S, T}(x[1] => x[2] for x in simplify_operation(terms)))
end

#function LinearPartition{T}(p::S) where { S <: AbstractPartition, T <: RingElement }
#    # LinearPartition([(p, one(T))])
#end

"""
    add(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

Return the addition between `p` and `q`.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = LinearPartition([(SetPartition([1, 2], [1, 1]), 4), 
(SetPartition([1, 1], [1, 1]), 4*d)])
julia> add(a, a)
LinearPartition(SetPartition([1, 2], [1, 1]) => 8, 
SetPartition([1, 1], [1, 1]) => 8*d)
```
"""
function add(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement } 

    result = copy(p.coefficients)

    for i in pairs(q.coefficients)
        result[i[1]] = get(result, i[1], 0) + i[2]
    end
    
    LinearPartition(result)
end

function +(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    add(p, q)
end

"""
    scale(a::RingElement, p::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

Multiply each coefficient in `p` with `a` and return
the result.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> scale(0.5, LinearPartition(SetPartition([1, 2], [1, 1]) => 8, 
SetPartition([1, 1], [1, 1]) => 8*d))
LinearPartition(SetPartition([1, 2], [1, 1]) => 4, 
SetPartition([1, 1], [1, 1]) => 4*d)
```
"""
function scale(a::RingElement, p::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    result = Dict{S, T}()

    for (i, n) in pairs(p.coefficients)
        result[i] = a * n
    end
    LinearPartition(result)
end 

function *(a::RingElement, p::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    scale(a, p)
end

"""
    composition(p::LinearPartition{P, R}, q::LinearPartition{P, R}, d::R) where 
        { P <: AbstractPartition, R <: RingElement }

Multiply each coefficient from `p` with each coefficient from `q`
as well as perform a composition between each SetPartition part 
in `p` and `q` and return the result.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = LinearPartition([(SetPartition([1, 2], [1, 1]), 4), 
(SetPartition([1, 1], [1, 1]), 4*d)])
julia> linear_composition(a, a)
LinearPartition(SetPartition([1, 2], [1, 1]) => 16*d + 16, 
SetPartition([1, 1], [1, 1]) => 16*d^2 + 16*d)
```
"""
function composition(p::LinearPartition{S, T}, q::LinearPartition{S, T}, d::RingElement) where 
        { S <: AbstractPartition, T <: RingElement }
    result = Dict{S, T}()
    
    for i in pairs(p.coefficients)
        for ii in pairs(q.coefficients)
            (composition, loop) = composition_loops(i[1], ii[1])
            new_coefficient = i[2] * ii[2] * (d^loop)
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end
    LinearPartition(result)
end 

function composition(p::LinearPartition{S, T}, q::LinearPartition{S, T}, d::RingElement) where 
        { S <: AbstractPartition, T <: RingElement }
    composition(p, q, d)
end

"""
    tensor_product(p::LinearPartition{P, R}, q::LinearPartition{P, R}) where 
        { P <: AbstractPartition, R <: RingElement }

Perform a tensor product similar

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = LinearPartition([(SetPartition([1, 2], [1, 1]), 4), 
(SetPartition([1, 1], [1, 1]), 4*d)])
julia> linear_tensor_product(a, a)
LinearPartition(SetPartition([1, 1, 2, 2], [1, 1, 2, 2]) => 8*d, 
SetPartition([1, 2, 3, 3], [1, 1, 3, 3]) => 4*d + 4, 
SetPartition([1, 2, 3, 4], [1, 1, 3, 3]) => 8, 
SetPartition([1, 1, 2, 3], [1, 1, 2, 2]) => 4*d + 4))
```
"""
function tensor_product(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

    result = Dict{S, T}()
    
    for i in pairs(p.coefficients)
        for ii in pairs(q.coefficients)
            composition = tensor_product(i[1], ii[1])
            new_coefficient = i[2] * ii[2]
            result[composition] = get(result, composition, 0) + new_coefficient
        end
    end

    LinearPartition(result)
    
end 

function ⊗(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    tensor_product(p, q)
end

# TODO think about involution

"""
    subtract(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }

Perform a subtraction between `p` and `q`and return the result.

# Examples
```jldoctest
julia> S, d = polynomial_ring(QQ, "d")
(Univariate polynomial ring in x over QQ, d)
julia> a = LinearPartition([(SetPartition([1, 2], [1, 1]), 4), 
(SetPartition([1, 1], [1, 1]), 4*d)])
julia> linear_tensor_product(a, a)
LinearPartition(SetPartition([1, 1], [1, 1]) => 0)
```
"""
function subtract(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
        { S <: AbstractPartition, T <: RingElement }
    add(p, scale(-1, q))
end 

function -(p::LinearPartition{S, T}, q::LinearPartition{S, T}) where 
    { S <: AbstractPartition, T <: RingElement }
    subtract(p, q)
end 
