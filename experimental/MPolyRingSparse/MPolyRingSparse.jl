"""
This is an approach to create a multivariate polynomial ring implementation
where the size of a monomial is independent of the number of variables in the ring.
Thus this implementation is to be used with low degree polynomials over rings
with lots of variables.
This strongly follows Generic.MPolyRingSparse from AbstractAlgebra.
"""

using AbstractAlgebra

import Base: deepcopy_internal, hash, isone, iszero, length, parent, +, -, *

import AbstractAlgebra: CacheDictType, get_cached!

import AbstractAlgebra: Ring, RingElement

import AbstractAlgebra: base_ring, check_parent, coeff, degree, elem_type, expressify, gen, gens, isconstant, isgen, isunit, nvars, parent_type, symbols, total_degree, vars

import AbstractAlgebra.Generic: ordering

export MPolyRingSparse, MPolySparse, PolynomialRingSparse

###############################################################################
#
#   MPolyRingSparse / MPolySparse
#
###############################################################################

# S is a Symbol which can take the values:
# :lex
# :revlex
# :deglex
# :degrevlex

mutable struct MPolyRingSparse{T <: RingElement} <: AbstractAlgebra.MPolyRing{T}
    base_ring::Ring
    S::Vector{Symbol}
    ord::Symbol
    num_vars::Int

    function MPolyRingSparse{T}(R::Ring, s::Vector{Symbol}, ord::Symbol,
                            cached::Bool = true) where T <: RingElement
        return get_cached!(MPolySparseID, (R, s, ord), cached) do
            new{T}(R, s, ord, length(s))
        end::MPolyRingSparse{T}
    end
end

const MPolySparseID = CacheDictType{Tuple{Ring, Vector{Symbol}, Symbol}, Ring}()

mutable struct MPolySparse{T <: RingElement} <: AbstractAlgebra.MPolyElem{T}
    coeffs::Vector{T}
    exps::Vector{Vector{Tuple{Int,Int}}}
    length::Int
    parent::MPolyRingSparse{T}

    function MPolySparse{T}(R::MPolyRingSparse) where T <: RingElement
        return new{T}(Vector{T}(undef, 0), Vector{Vector{Tuple{Int,Int}}}(undef, 0), 0, R)
    end

    MPolySparse{T}(R::MPolyRingSparse, a::Vector{T}, b::Vector{Vector{Tuple{Int,Int}}}) where T <: RingElement = new{T}(a, b, length(a), R)

    function MPolySparse{T}(R::MPolyRingSparse, a::T) where T <: RingElement
        return iszero(a) ? new{T}(Vector{T}(undef, 0), Vector{Vector{Tuple{Int,Int}}}(undef, 0), 0, R) :
                                            new{T}([a], [Vector{Tuple{Int,Int}}(undef, 0)], 1, R)
    end
end


###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent(a::MPolySparse{T}) where T <: RingElement = a.parent

parent_type(::Type{MPolySparse{T}}) where T <: RingElement = MPolyRingSparse{T}

elem_type(::Type{MPolyRingSparse{T}}) where T <: RingElement = MPolySparse{T}

base_ring(R::MPolyRingSparse{T}) where T <: RingElement = R.base_ring::parent_type(T)

symbols(a::MPolyRingSparse) = a.S

nvars(a::MPolyRingSparse) = a.num_vars

function gen(a::MPolyRingSparse{T}, i::Int, ::Type{Val{:lex}}) where {T <: RingElement}
    n = a.num_vars
    return a([one(base_ring(a))], [[(i,1)]])
end

function gen(a::MPolyRingSparse{T}, i::Int, ::Type{Val{:deglex}}) where {T <: RingElement}
    n = a.num_vars
    return a([one(base_ring(a))], [[(i,1)]])
end

function gen(a::MPolyRingSparse{T}, i::Int, ::Type{Val{:degrevlex}}) where {T <: RingElement}
    n = a.num_vars
    return a([one(base_ring(a))], [[(n-i+1,1)]])
end

function gens(a::MPolyRingSparse{T}) where {T <: RingElement}
    n = a.num_vars
    return [gen(a, i, Val{a.ord}) for i in 1:n]
end

function gen(a::MPolyRingSparse{T}, i::Int) where {T <: RingElement}
    return gen(a, i, Val{a.ord})
end

function vars(p::MPolySparse{T}) where {T <: RingElement}
    exps = p.exps
    gen_list = gens(p.parent)
    inds_in_p = sort!(unique([v[1] for v in exps[i] for i in 1:length(exps)]))
    
    return map(ind -> gen_list[ind], inds_in_p)
end

function ordering(a::MPolyRingSparse{T}) where {T <: RingElement}
    return a.ord
end

function check_parent(a::MPolySparse{T}, b::MPolySparse{T}, throw::Bool = true) where T <: RingElement
    b = parent(a) != parent(b)
    b & throw && error("Incompatible polynomial rings in polynomial operation")
    return !b
end


###############################################################################
#
#   Manipulating terms and monomials
#
###############################################################################

# TODO


###############################################################################
#
#   Monomial operations
#
###############################################################################

# TODO


###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(x::MPolySparse{T}, h::UInt) where {T <: RingElement}
    b = 0xf89bb54654fe9eae%UInt
    for i in 1:length(x)
        b = xor(b, xor(hash(x.coeffs[i], h), h))
        b = xor(b, xor(hash(x.exps[i], h), h))
        b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
    end
    return b
end

isunit(x::MPolySparse) = x.length == 1 && x.exps[1] == Vector{Tuple{Int,Int}}(undef, 0) && isunit(x.coeffs[1])

function isgen(x::MPolySparse{T}) where {T <: RingElement}
    if length(x) != 1
        return false
    end
    if !isone(x.coeff[1])
        return false
    end
    if length(x.exps[1]) != 1
        return false
    end
    return x.exps[1][1][2] == 1
end

# TODO
# function ishomogeneous(x::MPoly{T}) where {T <: RingElement}
#     last_deg = 0
#     is_first = true

#     for e in exponent_vectors(x)
#         d = sum(e)
#         if !is_first
#             if d != last_deg
#                 return false
#             else
#                 last_deg = d
#             end
#         else
#             is_first = false
#             last_deg = d
#         end
#     end
#     return true
# end

function coeff(x::MPolySparse, i::Int)
    return x.coeffs[i]
end

# TODO
# function monomial(x::MPoly, i::Int)
#     R = base_ring(x)
#     N = size(x.exps, 1)
#     exps = Matrix{UInt}(undef, N, 1)
#     monomial_set!(exps, 1, x.exps, i, N)
#     return parent(x)([one(R)], exps)
# end

# function monomial!(m::MPoly{T}, x::MPoly{T}, i::Int) where T <: RingElement
#     N = size(x.exps, 1)
#     fit!(m, 1)
#     monomial_set!(m.exps, 1, x.exps, i, N)
#     m.coeffs[1] = one(base_ring(x))
#     m.length = 1
#     return m
#  end

# function term(x::MPoly, i::Int)
#     R = base_ring(x)
#     N = size(x.exps, 1)
#     exps = Matrix{UInt}(undef, N, 1)
#     monomial_set!(exps, 1, x.exps, i, N)
#     return parent(x)([deepcopy(x.coeffs[i])], exps)
# end

function degree(f::MPolySparse{T}, i::Int) where T <: RingElement
    biggest = -1
    for j = 1:length(f)
        k = searchsortedfirst(f.exps[j], i; by=first)
        if 1 <= k <= length(f.exps[j]) && k == f.exps[j][k][1]
            d = f.exps[j][k][2]
        else
            d = 0
        end
        if d > biggest
            biggest = d
        end
    end
    return biggest
 end


function total_degree(f::MPolySparse{T}) where {T <: RingElement}
    ord = ordering(parent(f))
    if ord == :lex
        max_deg = -1
        for i = 1:length(f)
            sum_deg = sum(x -> x[2], f.exps[i])
            if sum_deg > max_deg
                max_deg = sum_deg
            end
        end
        return max_deg
    elseif ord == :deglex || ord == :degrevlex
        return length(f) == 0 ? -1 : sum(x -> x[2], f.exps[end])
    else
        error("total_degree is not implemented for this ordering.")
    end
end

length(x::MPolySparse) = x.length

isone(x::MPolySparse) = x.length == 1 && x.coeffs[1] == 1 && x.exps[1] == Vector{Tuple{Int,Int}}(undef, 0)

iszero(x::MPolySparse) = x.length == 0

isconstant(x::MPolySparse) = x.length == 0 || (x.length == 1 && x.exps[1] == Vector{Tuple{Int,Int}}(undef, 0))

function Base.deepcopy_internal(a::MPolySparse{T}, dict::IdDict) where {T <: RingElement}
    Re = deepcopy_internal(a.exps, dict)
    Rc = Array{T}(undef, a.length)
    for i = 1:a.length
        Rc[i] = deepcopy(a.coeffs[i])
    end
    return parent(a)(Rc, Re)
end


###############################################################################
#
#   Iterators
#
###############################################################################

# TODO


###############################################################################
#
#   Arithmetic functions
#
###############################################################################

# TODO
function -(a::MPolySparse{T}) where {T <: RingElement}
    return a
end

function +(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    return a
end

function -(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    return a
end

function *(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    return a
end

###############################################################################
#
#   Square root
#
###############################################################################

# TODO


###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

# TODO


###############################################################################
#
#   Comparison functions
#
###############################################################################

# TODO


###############################################################################
#
#   Ad hoc comparison functions
#
###############################################################################

# TODO


###############################################################################
#
#   Powering
#
###############################################################################

# TODO


###############################################################################
#
#   Inflation/deflation
#
###############################################################################

# TODO


###############################################################################
#
#   Exact division
#
###############################################################################

# TODO


###############################################################################
#
#   Euclidean division
#
###############################################################################

# TODO


###############################################################################
#
#   Evaluation
#
###############################################################################

# TODO


###############################################################################
#
#   GCD
#
###############################################################################

# TODO


###############################################################################
#
#   Build context
#
###############################################################################

# TODO


###############################################################################
#
#   GCD
#
###############################################################################

# TODO


###############################################################################
#
#   Unsafe functions
#
###############################################################################

# TODO


###############################################################################
#
#   Promotion rules
#
###############################################################################

# TODO


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::MPolyRingSparse{T})(b::RingElement) where {T <: RingElement}
    return a(base_ring(a)(b))
end

function (a::MPolyRingSparse{T})() where {T <: RingElement}
    z = MPolySparse{T}(a)
    return z
end

function (a::MPolyRingSparse{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
    z = MPolySparse{T}(a, base_ring(a)(b))
    return z
end

function (a::MPolyRingSparse{T})(b::T) where {T <: Union{Integer, Rational, AbstractFloat}}
    z = MPolySparse{T}(a, b)
    return z
end

function (a::MPolyRingSparse{T})(b::T) where {T <: RingElement}
    parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
    z = MPolySparse{T}(a, b)
    return z
end

function (a::MPolyRingSparse{T})(b::MPolySparse{T}) where {T <: RingElement}
    parent(b) != a && error("Unable to coerce polynomial")
    return b
end

function (a::MPolyRingSparse{T})(b::Vector{T}, e::Vector{Vector{Tuple{Int,Int}}}) where {T <: RingElement}
    if length(b) > 0 && isassigned(b, 1)
        parent(b[1]) != base_ring(a) && error("Unable to coerce to polynomial")
    end
    z = MPolySparse{T}(a, b, e)
    return z
end

# This is the main user interface for efficiently creating a polynomial. It accepts
# an array of coefficients and an array of exponent vectors. Sorting, coalescing of
# like terms and removal of zero terms is performed.
function (a::MPolyRingSparse{T})(b::Vector{T}, m::Vector{Vector{Int}}) where {T <: RingElement}
    if length(b) > 0 && isassigned(b, 1)
        parent(b[1]) != base_ring(a) && error("Unable to coerce to polynomial")
    end

    length(m) != length(b) && error("Exponent vector has length $(length(m)) (expected $(length(b)))")

    ord = ordering(a)
    # TODO

    # Pe = Matrix{UInt}(undef, N, length(m))

    # if ord == :lex
    #     for i = 1:length(m)
    #         for j = 1:N
    #             Pe[j, i] = UInt(m[i][N - j + 1])
    #         end
    #     end
    # elseif ord == :deglex
    #     for i = 1:length(m)
    #         for j = 1:N - 1
    #             Pe[j, i] = UInt(m[i][N - j])
    #         end
    #         Pe[N, i] = UInt(sum(m[i]))
    #     end
    # else # degrevlex
    #     for i = 1:length(m)
    #         for j = 1:N - 1
    #             Pe[j, i] = UInt(m[i][j])
    #         end
    #         Pe[N, i] = UInt(sum(m[i]))
    #     end
    # end

    # z = MPolySparse{T}(a, b, Pe)
    # z = sort_terms!(z)
    # z = combine_like_terms!(z)
    # return z
end


###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::MPolySparse, x = symbols(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    n = nvars(parent(a))
    for i in 1:length(a)
        prod = Expr(:call, :*)
        if !isone(a.coeffs[i])
            push!(prod.args, expressify(a.coeffs[i], context = context))
        end
        for (var, exp) in a.exps[i]
            if exp == 1
                push!(prod.args, x[var])
            else
                push!(prod.args, Expr(:call, :^, x[var], exp))
            end
        end
        push!(sum.args, prod)
    end
    return sum
end

###############################################################################
#
#   PolynomialRingSparse constructor
#
###############################################################################

function PolynomialRingSparse(R::Ring, s::Vector{Symbol}; cached::Bool = true, ordering::Symbol = :lex)
    T = elem_type(R)
    parent_obj = MPolyRingSparse{T}(R, s, ordering, cached)

    return parent_obj, gens(parent_obj)
end

function PolynomialRingSparse(R::Ring, s::Vector{String}; cached::Bool = true, ordering::Symbol = :lex)
    return PolynomialRingSparse(R, [Symbol(v) for v in s]; cached=cached, ordering=ordering)
end

function PolynomialRingSparse(R::Ring, s::Vector{Char}; cached::Bool = true, ordering::Symbol = :lex)
    return PolynomialRingSparse(R, [Symbol(v) for v in s]; cached=cached, ordering=ordering)
end

function PolynomialRingSparse(R::Ring, n::Int, s::Symbol=:x; cached::Bool = false, ordering::Symbol = :lex)
    return PolynomialRingSparse(R, [Symbol(s, i) for i=1:n], cached = cached, ordering = ordering)
end

function PolynomialRingSparse(R::Ring, n::Int, s::String; cached::Bool = false, ordering::Symbol = :lex)
    return PolynomialRingSparse(R, n, Symbol(s); cached=cached, ordering=ordering)
end

function PolynomialRingSparse(R::Ring, n::Int, s::Char; cached::Bool = false, ordering::Symbol = :lex)
    return PolynomialRingSparse(R, n, Symbol(s); cached=cached, ordering=ordering)
end
