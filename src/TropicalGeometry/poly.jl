################################################################################
#
#  Tropical polynomials
#
################################################################################


################################################################################
#
#  Tropical alternatives to generic functions
#
################################################################################

###
# Alternatives to generic functions that use R(1) and R(0) instead of one(R) and zero(R)
###
one(R::Generic.PolyRing{<:TropicalSemiringElem}) = R(one(base_ring(R)))
zero(R::Generic.PolyRing{<:TropicalSemiringElem}) = R(zero(base_ring(R)))
one(R::MPolyRing{<:TropicalSemiringElem}) = R(one(base_ring(R)))
zero(R::MPolyRing{<:TropicalSemiringElem}) = R(zero(base_ring(R)))

function polynomial_ring(R::TropicalSemiring, s::Symbol; cached::Bool = true)
   T = elem_type(R)
   parent_obj = Oscar.Generic.PolyRing{T}(R, s, cached)
   return parent_obj, parent_obj([zero(R), one(R)])
end


###
# Alternative to generic prints that display zero sums as 0
###
function AbstractAlgebra.expressify(@nospecialize(a::PolyRingElem{<:TropicalSemiringElem}), x = var(parent(a)); context = nothing)
    if iszero(a)
        return expressify(zero(base_ring(a)), context = context)
    end
    sum = Expr(:call, :+)
    for k in degree(a):-1:0
        c = coeff(a, k)
        if !iszero(c)
            xk = k < 1 ? expressify(one(base_ring(a)), context = context) : k == 1 ? x : Expr(:call, :^, x, k)
            if isone(c)
                push!(sum.args, Expr(:call, :*, xk))
            else
                push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
            end
        end
    end
    return sum
end
function AbstractAlgebra.expressify(a::MPolyRingElem{<:TropicalSemiringElem}, x = symbols(parent(a)); context = nothing)
    if iszero(a)
        return expressify(zero(base_ring(a)), context = context)
    end
    sum = Expr(:call, :+)
    n = nvars(parent(a))
    for (c, v) in zip(coefficients(a), exponents(a))
        prod = Expr(:call, :*)
        if !isone(c)
            push!(prod.args, expressify(c, context = context))
        end
        for i in 1:n
            if v[i] > 1
                push!(prod.args, Expr(:call, :^, x[i], v[i]))
            elseif v[i] == 1
                push!(prod.args, x[i])
            end
        end
        # Capture empty products
        if length(prod.args) == 1
            prod = expressify(one(base_ring(a)), context = context)
        end
        push!(sum.args, prod)
    end
    return sum
end


###
#  Disabling ideals in tropical polyomial rings
###
function MPolyIdeal(::Ring, ::Vector{Generic.MPoly{TropicalSemiringElem{minOrMax}}}) where {minOrMax <: Union{typeof(min),typeof(max)}}
    error("Ideals over tropical semirings not supported")
end


###
#  Roots of univariate tropical polynomials
###
"""
    roots(f::PolyRingElem{<:TropicalSemiringElem})

Return the tropical roots of a univariate tropical polynomial `f`, i.e., the variable values where
the min or max is attained at least twice.

# Examples

```jldoctest
julia> R,x = polynomial_ring(tropical_semiring(min),:x)
(Univariate polynomial ring in x over min tropical semiring, x)

julia> f = 1*x^2+x+0
(1)*x^2 + x + (0)

julia> roots(f)
2-element Vector{QQFieldElem}:
 0
 -1

julia> R,x = polynomial_ring(tropical_semiring(max),:x)
(Univariate polynomial ring in x over max tropical semiring, x)

julia> f = 1*x^2+x+0
(1)*x^2 + x + (0)

julia> roots(f)
1-element Vector{QQFieldElem}:
 -1//2

```
"""
function roots(f::PolyRingElem{TropicalSemiringElem{minOrMax}}) where {minOrMax <: Union{typeof(min),typeof(max)}}
    # Construct either the lower (min) or upper (max) convex hull of degrees and coefficients
    # Example: for min(1+2x,x,0) it is conv({(2,1), (1,0), (0,0)}) + RR_{>=0}*(0,1)
    #          for max(1+2x,x,0) it is conv({(2,1), (1,0), (0,0)}) + RR_{>=0}*(0,-1)
    # Note: univariate polynomials also enumerate over zero coefficients, necessitating !iszero(c) below
    valsAndExps = [ QQFieldElem[first(d),QQ(c)] for (d,c) in zip(exponents(f),coefficients(f)) if !iszero(c) ]

    if length(valsAndExps)<2
        return QQFieldElem[] # return empty set if f not at least binomial
    end

    s = minOrMax==typeof(min) ? 1 : -1 # +1 if min, -1 if max
    hullDirection = s*[0 1]
    newtonPolygon = convex_hull(valsAndExps,hullDirection)

    # The negated slopes of the lower (min) or upper (max) edges are when function value is attained at least twice
    # Example: for min(1+2x,x,0) it is -1 and 0, for max(1+2x,x,0) it is -1/2
    F = affine_inequality_matrix(facets(newtonPolygon))
    negatedSlopes = [ F[i,2]/F[i,3] for i in 1:nrows(F) if s*F[i,3]<0 ]
    return negatedSlopes
end



################################################################################
#
#  @tropical macro
#
################################################################################

"""
    @tropical(expr)

Translate the expression in the tropical world.

# Examples

```jldoctest
julia> T = tropical_semiring(min);

julia> Tx, x = polynomial_ring(T, :x => 1:3);

julia> @tropical min(1, x[1], x[2], 2*x[3])
x[3]^2 + x[1] + x[2] + (1)
```
"""
macro tropical(expr)
    e = _tropicalize(expr)
    return quote
        $(esc(e))
    end
end

_tropicalize(x::Symbol) = x

_tropicalize(x::Int) = x

function _tropicalize(x::Expr)
    if x.head == :call
        if x.args[1] == :min
            x.args[1] = :(+)
        elseif x.args[1] == :(*)
            length(x.args) <= 3 || error("Cannot convert")
            x.args[1] = :(Oscar._tropical_mul)
        elseif x.args[1] == :(+)
            x.args[1] = :*
                else
            error("Cannot convert")
        end
        for i in 2:length(x.args)
            x.args[i] = _tropicalize(x.args[i])
        end
    else
        return x
    end
    return x
end

function _tropical_mul(x, y)
    if x isa Union{Integer, Rational, QQFieldElem, ZZRingElem}
        if x isa Rational || x isa QQFieldElem
            _x = ZZ(x)
            return y^_x
        else
            return y^x
        end
    elseif y isa Union{Integer, Rational, QQFieldElem, ZZRingElem}
        if y isa Rational | y isa QQFieldElem
            _y = ZZ(y)
            return x^_y
        else
            return x^y
        end
    else
        error("Cannot convert ", x, " * ", y)
    end
end



################################################################################
#
#  Conversion to tropical polynomial
#
################################################################################

@doc raw"""
    tropical_polynomial(f::Union{<:MPolyRingElem,<:PolyRingElem},nu::TropicalSemiringMap)

Given a polynomial `f` and a tropical semiring map `nu`,
return the tropicalization of `f` as a polynomial over the tropical semiring.

# Examples
```jldoctest
julia> R, (x,y) = polynomial_ring(QQ,[:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> nu = tropical_semiring_map(QQ,7)
Map into Min tropical semiring encoding the 7-adic valuation on Rational field

julia> f = 7*x+y+49
7*x + y + 49

julia> tropical_polynomial(f,nu)
(1)*x + y + (2)
```
"""
function tropical_polynomial(f::Union{<:MPolyRingElem,<:PolyRingElem}, nu::Union{Nothing,TropicalSemiringMap}=nothing)

    # if unspecified, set nu to be the trivial valuation + min convention
    isnothing(nu) && (nu = tropical_semiring_map(coefficient_ring(f)))

    T = tropical_semiring(nu)
    Tx,x = polynomial_ring(T,[repr(x) for x in gens(parent(f))])
    tropf = inf(T)

    for (c,alpha) in zip(coefficients(f), exponents(f))
        tropf = tropf + nu(c)*monomial(Tx,alpha)
    end

    return tropf
end



################################################################################
#
#  Basic functions
#
################################################################################

@doc raw"""
    newton_subdivision(f::Generic.MPoly{TropicalSemiringElem{minOrMax}}) where minOrMax<:Union{typeof(min),typeof(max)}

Return the dual subdivision on `exponents(f)` with weights `coefficients(f)` (min-convention) or `-coefficients(f)` (max-convention).  It is dual to `tropical_hypersurface(f)`.

# Examples
```jldoctest
julia> _, (x,y) = polynomial_ring(tropical_semiring(),[:x, :y]);

julia> f = 1+x+y+x^2;

julia> Deltaf = newton_subdivision(f)
Subdivision of points in ambient dimension 2

julia> points(Deltaf)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [2, 0]
 [1, 0]
 [0, 1]
 [0, 0]

julia> maximal_cells(Deltaf)
2-element SubObjectIterator{Vector{Int64}}:
 [2, 3, 4]
 [1, 2, 3]

```
"""
function newton_subdivision(f::Generic.MPoly{TropicalSemiringElem{minOrMax}}) where minOrMax<:Union{typeof(min),typeof(max)}
    # preserve_ordering=true, since weights of regular subdivisions have to be in min-convention,
    # e.g., [0 0; 1 0; 0 1; 2 0] decomposed into [0 0; 1 0; 0 1] and [1 0; 0 1; 2 0] has min_weight [+1,0,0,0]
    # which is dual to the tropical hypersurface of min(+1, x, y, 2*x) or max(-1, x, y, 2*x)
    weights = QQ.(coefficients(f); preserve_ordering=true)
    points = matrix(QQ,collect(exponents(f)))
    return subdivision_of_points(points,weights)
end


@doc raw"""
    newton_subdivision(f::MPolyRingElem, nu::Union{Nothing,TropicalSemiringMap}=nothing)

Return the dual subdivision on `exponents(f)` with weights `nu.(coefficients(f))` (min-convention) or `-nu.(coefficients(f))` (max-convention).  It is dual to `tropical_hypersurface(f,nu)`.

# Examples
```jldoctest
julia> _, (x,y) = QQ[:x, :y];

julia> nu = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> f = 2+x+y+x^2;

julia> Deltaf = newton_subdivision(f,nu)
Subdivision of points in ambient dimension 2

julia> points(Deltaf)
4-element SubObjectIterator{PointVector{QQFieldElem}}:
 [2, 0]
 [1, 0]
 [0, 1]
 [0, 0]

julia> maximal_cells(Deltaf)
2-element SubObjectIterator{Vector{Int64}}:
 [2, 3, 4]
 [1, 2, 3]

```
"""
function newton_subdivision(f::MPolyRingElem, nu::Union{Nothing,TropicalSemiringMap}=nothing)
    tropf = tropical_polynomial(f,nu)
    return newton_subdivision(tropf)
end

################################################################################
#
#  Outdated code
#
################################################################################

# # Disabled due to potential clash with
# # tropical_polynomial(f::MPolyRingElem, nu::Union{Nothing,TropicalSemiringMap}=nothing)
# @doc raw"""
#     tropical_polynomial(f::MPolyRingElem,M::Union{typeof(min),typeof(max)}=min)

# Given a polynomial `f` over a field with an intrinsic valuation (i.e., a field
# on which a function `valuation` is defined such as `PadicField(7,2)`),
# return the tropicalization of `f` as a polynomial over the min tropical semiring
# (default) or the max tropical semiring.

# # Examples
# ```jldoctest
# julia> K = PadicField(7, 2)
# Field of 7-adic numbers

# julia> Kxy, (x,y) = K[:x, :y]
# (Multivariate polynomial ring in 2 variables over QQ_7, AbstractAlgebra.Generic.MPoly{PadicFieldElem}[x, y])

# julia> f = 7*x+y+49
# (7^1 + O(7^3))*x + y + 7^2 + O(7^4)

# julia> tropical_polynomial(f,min)
# (1)*x + y + (2)

# julia> tropical_polynomial(f,max)
# (-1)*x + y + (-2)
# ```
# """
# function tropical_polynomial(f::MPolyRingElem, M::Union{typeof(min),typeof(max)}=min)
#   T = tropical_semiring(M)
#   if M==min
#     s=1
#   else
#     s=-1
#   end

#   Tx,x = polynomial_ring(T,[repr(x) for x in gens(parent(f))])
#   tropf = inf(T)

#   if base_ring(parent(f)) isa NonArchLocalField
#     for (c,alpha) in zip(coefficients(f), exponents(f))
#       tropf = tropf + T(s*valuation(c))*monomial(Tx,alpha)
#     end
#   else
#     for alpha in exponents(f)
#       tropf = tropf + T(0)*monomial(Tx,alpha)
#     end
#   end

#   return tropf
# end
