"""
This is an approach to create a multivariate polynomial ring implementation
where the size of a monomial is independent of the number of variables in the ring.
Thus this implementation is to be used with low degree polynomials over rings
with lots of variables.
This strongly follows Generic.MPolyRingSparse from AbstractAlgebra.
"""

import Base: deepcopy_internal, divrem, hash, isless, isone, iszero, 
             length, parent, sqrt, +, -, *, ^, ==

import AbstractAlgebra: CacheDictType, get_cached!, internal_power

import AbstractAlgebra: Ring, RingElement

import AbstractAlgebra: base_ring, change_base_ring, change_coefficient_ring,
                        check_parent, coeff, combine_like_terms!,
                        degree, divexact, divides, elem_type, exponent, exponent_vector,
                        expressify, fit!, gen, gens,
                        isconstant, isgen, ishomogeneous, issquare, isunit,
                        map_coefficients, monomial, monomial!,
                        nvars, parent_type, setcoeff!, set_exponent_vector!,
                        sort_terms!, symbols, term, total_degree, vars, zero!

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
    return a([one(base_ring(a))], [[(i,1)]])
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
    inds_in_p = sort!(unique([v[1] for i in 1:length(p) for v in exps[i]]))
    
    return map(ind -> gen_list[ind], inds_in_p)
end

function var_index(x::MPolySparse{T}) where {T <: RingElement}
    !isgen(x) && error("Not a variable in var_index")
    return x.exps[1][1][1]
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

function transform_exps(e::Vector{Int})
    exp = Vector{Tuple{Int,Int}}(undef, 0)
    for (i,j) in enumerate(e)
        if !iszero(j)
            push!(exp, (i,j))
        end
    end
    return exp
end

function transform_exps2(e::Vector{Tuple{Int,Int}}, N::Int)
    exp = zeros(Int, N)
    for (i,j) in e
        exp[i] = j
    end
    return exp
end

function exponent_vector(a::MPolySparse{T}, i::Int) where T <: RingElement
    return transform_exps2(a.exps[i], nvars(parent(a)))
end

function exponent(a::MPolySparse{T}, i::Int, j::Int) where T <: RingElement
    k = searchsortedfirst(a.exps[i], j; by=first)
    if 1 <= k <= length(a.exps[i])
        return a.exps[i][k][2]
    else
        return 0
    end
end

function coeff(a::MPolySparse, i::Int)
    return a.coeffs[i]
end

function coeff(a::MPolySparse{T}, e::Vector{Tuple{Int,Int}}) where T <: RingElement
    monomial_lt = (x, y) -> monomial_isless(x, y, parent(a))
    k = searchsortedfirst(a.exps, e; lt = monomial_lt)
    if 1 <= k <= length(a)
        return a.coeffs[k]
    else
        return parent(a)()
    end
end

function coeff(a::MPolySparse{T}, e::Vector{Int}) where T <: RingElement
    return coeff(a, transform_exps(e))
end

###############################################################################
#
#   Monomial operations
#
###############################################################################

function monomial_isless(a::Vector{Tuple{Int,Int}}, b::Vector{Tuple{Int,Int}}, R::MPolyRing{T}) where {T <: RingElement}
    return monomial_cmp(a, b, R) < 0
 end

function monomial_cmp(a::Vector{Tuple{Int,Int}}, b::Vector{Tuple{Int,Int}}, R::MPolyRing{T}) where {T <: RingElement}
    if R.ord == :lex
        for i in 1:min(length(a), length(b))
            if a[i][1] != b[i][1]
                return b[i][1] - a[i][1]
            elseif a[i][2] != b[i][2]
                return a[i][2] - b[i][2]
            end
        end
        return length(a) - length(b)
    elseif R.ord == :deglex
        dega = sum(x -> x[2], a; init=0)
        degb = sum(x -> x[2], b; init=0)
        if dega != degb
            return dega - degb
        end
        for i in 1:min(length(a), length(b))
            if a[i][1] != b[i][1]
                return b[i][1] - a[i][1]
            elseif a[i][2] != b[i][2]
                return a[i][2] - b[i][2]
            end
        end
        return length(a) - length(b)
    elseif R.ord == :degrevlex
        dega = sum(x -> x[2], a; init=0)
        degb = sum(x -> x[2], b; init=0)
        if dega != degb
            return dega - degb
        end
        for i in 0:min(length(a), length(b))-1
            if a[length(a)-i][1] != b[length(b)-i][1]
                return b[length(b)-i][1] - a[length(a)-i][1]
            elseif a[length(a)-i][2] != b[length(b)-i][2]
                return b[length(b)-i][2] - a[length(a)-i][2]
            end
        end
        return length(a) - length(b)
    end
end

function merge_exps(a::Vector{Tuple{Int,Int}}, b::Vector{Tuple{Int,Int}})
    r = Tuple{Int,Int}[]
    i = 1
    j = 1
    while i <= length(a) && j <= length(b)
        if a[i][1] < b[j][1]
            push!(r, a[i])
            i += 1
        elseif a[i][1] == b[j][1]
            push!(r, (a[i][1], a[i][2] + b[j][2]))
            i += 1
            j += 1
        else
            push!(r, b[j])
            j += 1
        end
    end
    while i <= length(a)
        push!(r, a[i])
        i += 1
    end
    while j <= length(b)
        push!(r, b[j])
        j += 1
    end
    return r
end


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
    if !isone(x.coeffs[1])
        return false
    end
    if length(x.exps[1]) != 1
        return false
    end
    return x.exps[1][1][2] == 1
end

function ishomogeneous(x::MPolySparse{T}) where {T <: RingElement}
    last_deg = 0
    is_first = true
 
    for e in x.exps
       d = sum(x -> x[2], e; init=0)
       if !is_first
          if d != last_deg
             return false
          else
             last_deg = d
          end
       else
          is_first = false
          last_deg = d
       end
    end
    return true
 end

function monomial(x::MPolySparse{T}, i::Int) where {T <: RingElement}
    return parent(x)([one(base_ring(x))], [x.exps[i]])
end

function monomial!(m::MPolySparse{T}, x::MPolySparse{T}, i::Int) where {T <: RingElement}
    fit!(m, 1)
    m.exps[1] = x.exps[i]
    m.coeffs[1] = one(base_ring(x))
    m.length = 1
    return m
 end


function term(x::MPolySparse{T}, i::Int) where {T <: RingElement}
    return parent(x)([x.coeffs[i]], [x.exps[i]])
end


function degree(f::MPolySparse{T}, i::Int) where {T <: RingElement}
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
            sum_deg = sum(x -> x[2], f.exps[i]; init=0)
            if sum_deg > max_deg
                max_deg = sum_deg
            end
        end
        return max_deg
    elseif ord == :deglex || ord == :degrevlex
        return length(f) == 0 ? -1 : sum(x -> x[2], f.exps[1]; init = 0)
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

function -(a::MPolySparse{T}) where {T <: RingElement}
    r = zero(a)
    fit!(r, length(a))
    for i in 1:length(a)
        r.coeffs[i] = -a.coeffs[i]
        r.exps[i] = deepcopy(a.exps[i])
    end
    r.length = a.length
    return r
end

function +(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    r = zero(a)
    fit!(r, length(a) + length(b))
    i = 1
    j = 1
    k = 1
    while i <= length(a) && j <= length(b)
        cmpexp = monomial_cmp(a.exps[i], b.exps[j], parent(a))
        if cmpexp > 0
            r.coeffs[k] = a.coeffs[i]
            r.exps[k] = deepcopy(a.exps[i])
            i += 1
        elseif cmpexp == 0
            c = a.coeffs[i] + b.coeffs[j]
            if !iszero(c)
                r.coeffs[k] = c
                r.exps[k] = deepcopy(a.exps[i])
            else
                k -= 1
            end
            i += 1
            j += 1
        else
            r.coeffs[k] = b.coeffs[j]
            r.exps[k] = deepcopy(b.exps[j])
            j += 1
        end
        k += 1
    end
    while i <= length(a)
        r.coeffs[k] = a.coeffs[i]
        r.exps[k] = deepcopy(a.exps[i])
        i += 1
        k += 1
    end
    while j <= length(b)
        r.coeffs[k] = b.coeffs[j]
        r.exps[k] = deepcopy(b.exps[j])
        j += 1
        k += 1
    end
    r.length = k - 1
    return r
end

function -(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    r = zero(a)
    fit!(r, length(a) + length(b))
    i = 1
    j = 1
    k = 1
    while i <= length(a) && j <= length(b)
        cmpexp = monomial_cmp(a.exps[i], b.exps[j], parent(a))
        if cmpexp > 0
            r.coeffs[k] = a.coeffs[i]
            r.exps[k] = deepcopy(a.exps[i])
            i += 1
        elseif cmpexp == 0
            c = a.coeffs[i] - b.coeffs[j]
            if !iszero(c)
                r.coeffs[k] = c
                r.exps[k] = deepcopy(a.exps[i])
            else
                k -= 1
            end
            i += 1
            j += 1
        else
            r.coeffs[k] = -b.coeffs[j]
            r.exps[k] = deepcopy(b.exps[j])
            j += 1
        end
        k += 1
    end
    while i <= length(a)
        r.coeffs[k] = a.coeffs[i]
        r.exps[k] = deepcopy(a.exps[i])
        i += 1
        k += 1
    end
    while j <= length(b)
        r.coeffs[k] = -b.coeffs[j]
        r.exps[k] = deepcopy(b.exps[j])
        j += 1
        k += 1
    end
    r.length = k - 1
    return r
end

function *(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    r = zero(a)
    fit!(r, length(a) * length(b))
    for i in 1:length(a), j in 1:length(b)
        r.coeffs[(i-1)*length(b) + j] = a.coeffs[i]*b.coeffs[j]
        r.exps[(i-1)*length(b) + j] = merge_exps(a.exps[i], b.exps[j])
    end
    r.length = length(a) * length(b)
    return combine_like_terms!(sort_terms!(r))
end


###############################################################################
#
#   Transformation to/from dense polynomials
#
###############################################################################

struct SparseToDenseData{T <: RingElement}
    vmap::Vector{Int}
    sparse_R::MPolyRingSparse{T}
    dense_R::MPolyRing{T}
end

function sparse_to_dense(as::MPolySparse{T}...) where {T <: RingElement}
    for a in as
        check_parent(a, as[1], true)
    end
    
    sparse_R = parent(as[1])

    vmap = sort!(unique!(map(var_index, reduce(vcat, map(vars, as))))) # vars_indices occurring in as
    N = length(vmap) > 0 ? length(vmap) : 1

    dense_R, _ = PolynomialRing(base_ring(sparse_R), N; ordering=sparse_R.ord)
    dense_as = map(function (a)
        ctx = MPolyBuildCtx(dense_R)
        for i in 1:length(a)
            expv = zeros(Int, N)
            for (j,k) in a.exps[i]
                expv[findfirst(==(j), vmap)] = k
            end
            push_term!(ctx, coeff(a, i), expv)
        end
        return finish(ctx)
    end, as)
    return dense_as, SparseToDenseData{T}(vmap, sparse_R, dense_R)
end

function dense_to_sparse(data::SparseToDenseData{T}, dense_as::MPolyElem{T}...) where {T <: RingElement}
    sparse_R = data.sparse_R
    vmap = data.vmap
    as = map(function (da)
        a = zero(sparse_R)
        fit!(a, length(da))
        for i in 1:length(da)
            expv = Vector{Tuple{Int,Int}}(undef, 0)
            for (j, k) in enumerate(exponent_vector(da, i))
                if !iszero(k)
                    push!(expv, (vmap[j],k))
                end
            end
            sort!(expv; by=first)
            set_exponent_vector!(a, i, expv)
            setcoeff!(a, i, coeff(da, i))
        end
        return combine_like_terms!(sort_terms!(a))
    end, dense_as)
    return as
end


###############################################################################
#
#   Square root
#
###############################################################################

function Base.sqrt(a::MPolySparse{T}; check::Bool=true) where {T <: RingElement}
    (da,), data = sparse_to_dense(a)
    dq = sqrt(da)
    (q,) = dense_to_sparse(data, dq)
    return q
end

function issquare(a::MPolySparse{T}) where {T <: RingElement}
    (da,), _ = sparse_to_dense(a)
    return issquare(da)
end


###############################################################################
#
#   Comparison functions
#
###############################################################################

function ==(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    fl = check_parent(a, b, false)
    !fl && return false
    if length(a) != length(b)
        return false
    end
    for i = 1:length(a)
        if a.exps[i] != b.exps[i]
            return false
        end
        if a.coeffs[i] != b.coeffs[i]
            return false
        end
    end
    return true
end

function Base.isless(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    check_parent(a, b)
    (!ismonomial(a) || !ismonomial(b)) && error("Not monomials in comparison")
    return monomial_isless(a.exps[1], b.exps[1], parent(a))
end


###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::MPolySparse{T}, b::Int) where {T <: RingElement}
    b < 0 && throw(DomainError(b, "exponent must be >= 0"))
    # special case powers of x for constructing polynomials efficiently
    if iszero(a)
        if iszero(b)
            return one(a)
        else
            return zero(a)
        end
    elseif b == 0
        return one(a)
    elseif b == 1
        return deepcopy(a)
    elseif b == 2
        return a*a
    else
        return internal_power(a, b)
    end
end

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

function divides(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    (da, db), data = sparse_to_dense(a, b)
    (flag, dq) = divides(da, db)
    (q,) = dense_to_sparse(data, dq)
    return flag, q
end

function divexact(a::MPolySparse{T}, b::MPolySparse{T}; check::Bool=true) where {T <: RingElement}
    d, q = divides(a, b)
    check && d == false && error("Not an exact division in divexact")
    return q
end

function divexact(a::MPolySparse, n::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
    (da,), data = sparse_to_dense(a)
    (dq) = divexact(da, n; check=check)
    (q,) = dense_to_sparse(data, dq)
    return q
end

function divexact(a::MPolySparse{T}, n::T; check::Bool=true) where {T <: RingElem}
    (da,), data = sparse_to_dense(a)
    (dq) = divexact(da, n; check=check)
    (q,) = dense_to_sparse(data, dq)
    return q
end


###############################################################################
#
#   Euclidean division
#
###############################################################################

function Base.divrem(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    (da, db), data = sparse_to_dense(a, b)
    (dq, dr) = divrem(da, db)
    (q, r) = dense_to_sparse(data, dq, dr)
    return q, r
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function (a::MPolySparse{T})(vals::Union{NCRingElem, RingElement}...) where T <: RingElement
    length(vals) != nvars(parent(a)) && error("Number of variables does not match number of values")
    R = base_ring(a)
    # The best we can do here is to cache previously used powers of the values
    # being substituted, as we cannot assume anything about the relative
    # performance of powering vs multiplication. The function should not try
    # to optimise computing new powers in any way.
    # Note that this function accepts values in a non-commutative ring, so operations
    # must be done in a certain order.
    powers = [Dict{Int, Any}() for i in 1:length(vals)]
    # First work out types of products
    r = R()
    c = zero(R)
    U = Vector{Any}(undef, length(vals))
    for j = 1:length(vals)
        W = typeof(vals[j])
        if ((W <: Integer && W != BigInt) ||
            (W <: Rational && W != Rational{BigInt}))
            c = c*zero(W)
            U[j] = parent(c)
        else
            U[j] = parent(vals[j])
            c = c*zero(parent(vals[j]))
        end
    end
    cvzip = zip(coefficients(a), exponent_vectors(a))
    for (c, v) in cvzip
        t = c
        for j = 1:length(vals)
            exp = v[j]
            if !haskey(powers[j], exp)
                powers[j][exp] = (U[j](vals[j]))^exp
            end
            t = t*powers[j][exp]
        end
        r += t
    end
    return r
end


###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    (da, db), data = sparse_to_dense(a, b)
    (dg) = gcd(da, db)
    (g,) = dense_to_sparse(data, dg)
    return g
end

function lcm(a::MPolySparse{T}, b::MPolySparse{T}) where {T <: RingElement}
    (da, db), data = sparse_to_dense(a, b)
    (dl) = lcm(da, db)
    (l,) = dense_to_sparse(data, dl)
    return l
end


###############################################################################
#
#  Change base ring
#
###############################################################################

function map_coefficients(f, p::MPolySparse; cached = true, parent::MPolyRingSparse = _change_mpoly_ring(parent(f(zero(base_ring(p)))), parent(p), cached))
    return _map(f, p, parent)
end

function _map(g, p::MPolySparse, Rx)
    cvzip = zip(coefficients(p), exponent_vectors(p))
    M = MPolyBuildCtx(Rx)
    for (c, v) in cvzip
        push_term!(M, g(c), v)
    end

    return finish(M)
end

function _change_mpoly_ring(R, Rx, cached)
    P, _ = PolynomialRingSparse(R, map(string, symbols(Rx)), ordering = ordering(Rx), cached = cached)
    return P
 end

function change_base_ring(R::Ring, p::MPolySparse{T}; cached = true, parent::MPolyRingSparse = _change_mpoly_ring(R, parent(p), cached)) where {T <: RingElement}
    base_ring(parent) != R && error("Base rings do not match.")
    return _map(R, p, parent)
end

function change_coefficient_ring(R::Ring, p::MPolySparse{T}; cached = true, parent::MPolyRingSparse = _change_mpoly_ring(R, parent(p), cached)) where {T <: RingElement}
   return change_base_ring(R, p, cached = cached, parent = parent)
end


###############################################################################
#
#   Unsafe functions
#
###############################################################################

function fit!(a::MPolySparse{T}, n::Int) where {T <: RingElement}
    if length(a) < n
        resize!(a.coeffs, n)
        resize!(a.exps, n)
    end
    return nothing
end

function zero!(a::MPolySparse{T}) where {T <: RingElement}
    a.length = 0
    return a
end

function set_exponent_vector!(a::MPolySparse{T}, i::Int, e::Vector{Tuple{Int,Int}}) where T <: RingElement
    n = nvars(parent(a))
    all(x -> 0 < x[1] <= n, e) || error("variable index out of range")
    fit!(a, i)
    a.exps[i] = filter(x -> !iszero(x[2]), e)
    if i > length(a)
        a.length = i
    end
    return a
end

function set_exponent_vector!(a::MPolySparse{T}, i::Int, e::Vector{Int}) where T <: RingElement
    return set_exponent_vector!(a, i, transform_exps(e))
end

for T in [RingElem, Integer, Rational, AbstractFloat]
    @eval begin
        function setcoeff!(a::MPolySparse{S}, i::Int, c::S) where {S <: $T}
            fit!(a, i)
            a.coeffs[i] = c
            if i > length(a)
                a.length = i
            end
            return a
        end
    end
end

function sort_terms!(a::MPolySparse{T}) where {T <: RingElement}
    n = length(a)
    monomial_lt = (x, y) -> monomial_isless(x, y, parent(a))
    if n > 1
       p = sortperm(view(a.exps, 1:n), lt = monomial_lt, rev = true)
       a.coeffs = [a.coeffs[p[i]] for i in 1:n]
       a.exps = [a.exps[p[i]] for i in 1:n]
    end
    return a
end

function combine_like_terms!(a::MPolySparse{T}) where {T <: RingElement}
    o = 0
    i = 1
    while i <= length(a)
        if o > 0 && monomial_cmp(a.exps[o], a.exps[i], parent(a)) == 0
            a.coeffs[o] += a.coeffs[i]
        else
            o += (o < 1 || !iszero(a.coeffs[o]))
            a.exps[o] = a.exps[i]
            a.coeffs[o] = a.coeffs[i]
        end
        i += 1
    end
    o += (o < 1 || !iszero(a.coeffs[o]))
    a.length = o - 1
    return a
end


###############################################################################
#
#   Promotion rules
#
###############################################################################

AbstractAlgebra.promote_rule(::Type{MPolySparse{T}}, ::Type{MPolySparse{T}}) where T <: RingElement = MPolySparse{T}

function AbstractAlgebra.promote_rule(::Type{MPolySparse{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
    AbstractAlgebra.promote_rule(T, U) == T ? MPolySparse{T} : Union{}
end


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
    z = MPolySparse{T}(a, b, map(x -> filter(y -> !iszero(y[2]), x), e))
    return z
end

# This is the main user interface for efficiently creating a polynomial. It accepts
# an array of coefficients and an array of exponent vectors. Sorting, coalescing of
# like terms and removal of zero terms is performed.
function (R::MPolyRingSparse{T})(b::Vector{T}, e::Vector{Vector{Int}}) where {T <: RingElement}
    if length(b) > 0 && isassigned(b, 1)
        parent(b[1]) != base_ring(R) && error("Unable to coerce to polynomial")
    end

    length(e) != length(b) && error("Exponent vector has length $(length(e)) (expected $(length(b)))")

    z = MPolySparse{T}(R, b, map(transform_exps, e))
    z = sort_terms!(z)
    z = combine_like_terms!(z)
    return z
end


###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::MPolySparse, x = symbols(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    for i in 1:length(a)
        prod = Expr(:call, :*)
        if !isone(a.coeffs[i])
            push!(prod.args, expressify(a.coeffs[i], context = context))
        end
        for (var, exp) in a.exps[i]
            if exp == 1
                push!(prod.args, x[var])
            elseif exp > 1
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
