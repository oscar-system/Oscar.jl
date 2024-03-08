###############################################################################
# FreeModElem constructors
###############################################################################

@doc raw"""
    FreeModElem(c::SRow{T}, parent::FreeMod{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
FreeModElem(c::SRow{T}, parent::FreeMod{T}) where T = FreeModElem{T}(c, parent)

@doc raw"""
    FreeModElem(c::Vector{T}, parent::FreeMod{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
function FreeModElem(c::Vector{T}, parent::FreeMod{T}) where T
  @assert length(c) == rank(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:rank(parent)), c)
  return FreeModElem{T}(sparse_coords,parent)
end

#@doc raw"""
#    (F::FreeMod{T})(c::SRow{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#"""
function (F::FreeMod{T})(c::SRow{T}) where T
  return FreeModElem(c, F)
end

#@doc raw"""
#    (F::FreeMod{T})(c::Vector{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#
## Examples
#```jldoctest
#julia> R, (x,y) = polynomial_ring(QQ, ["x", "y"])
#(Multivariate Polynomial Ring in x, y over Rational Field, QQMPolyRingElem[x, y])
#
#julia> F = FreeMod(R,3)
#Free module of rank 3 over Multivariate Polynomial Ring in x, y over Rational Field
#
#julia> V = [x, zero(R), y]
#3-element Vector{QQMPolyRingElem}:
# x
# 0
# y
#
#julia> f = F(V)
#x*e[1] + y*e[3]
#```
#"""
function (F::FreeMod{T})(c::Vector{T}) where T 
 return FreeModElem(c, F)
end

function (F::AbstractFreeMod{T})(v::AbstractFreeModElem{T}) where T
  @assert parent(v) === F
  return v
end

function in(v::AbstractFreeModElem, F::AbstractFreeMod)
  return parent(v) === F
end

function in(v::AbstractFreeModElem, M::SubquoModule)
  return represents_element(v, M)
end

@doc raw"""
    coordinates(v::AbstractFreeModElem)

Return the entries (with respect to the standard basis) of `v` as a sparse row.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> F = FreeMod(R,3)
Free module of rank 3 over Multivariate polynomial ring in 2 variables over QQ

julia> f = x*gen(F,1)+y*gen(F,3)
x*e[1] + y*e[3]

julia> coordinates(f)
Sparse row with positions [1, 3] and values QQMPolyRingElem[x, y]
```
"""
function coordinates(v::AbstractFreeModElem)
  return v.coords
end

#########################################################

@doc raw"""
    repres(v::AbstractFreeModElem)

Return just `v`. This function exists for compatibility (with subquotient elements) reasons.
"""
function repres(v::AbstractFreeModElem)
  return v
end

function getindex(v::AbstractFreeModElem, i::Int)
  if isempty(coordinates(v))
    return zero(base_ring(parent(v)))
  end
  return coordinates(v)[i]
end

elem_type(::Type{FreeMod{T}}) where {T} = FreeModElem{T}
parent_type(::Type{FreeModElem{T}}) where {T} = FreeMod{T}

function generator_symbols(F::FreeMod)
  return F.S
end

function expressify(e::AbstractFreeModElem; context = nothing)
  sum = Expr(:call, :+)
  for (pos, val) in e.coords
     # assuming generator_symbols(parent(e)) is an array of strings/symbols
     push!(sum.args, Expr(:call, :*, expressify(val, context = context), generator_symbols(parent(e))[pos]))
  end
  return sum
end
@enable_all_show_via_expressify FreeModElem

@doc raw"""
    Vector(e::FreeModElem)

Return the coordinates of `e` as a Vector.
"""
function Vector(e::FreeModElem)
   return [e[i] for i in 1:rank(parent(e))]
end

@doc raw"""
    basis(F::AbstractFreeMod)

Return the standard basis of `F`.
"""
function basis(F::AbstractFreeMod)
  bas = elem_type(F)[]
  for i=1:dim(F)
    s = Hecke.sparse_row(base_ring(F), [(i, base_ring(F)(1))])
    push!(bas, F(s))
  end
  return bas
end


@doc raw"""
    gens(F::AbstractFreeMod)

Return the (canonical) generators of the free module `F`.
"""
gens(F::AbstractFreeMod) = basis(F)

@doc raw"""
    basis(F::AbstractFreeMod, i::Int)

    gen(F::AbstractFreeMod, i::Int)

Return the `i`th basis vector of `F`, that is, return the `i`th standard unit vector.
"""
function basis(F::AbstractFreeMod, i::Int)
  @assert 0 < i <= ngens(F)
  s = Hecke.sparse_row(base_ring(F), [(i, base_ring(F)(1))])
  return F(s)
end
gen(F::AbstractFreeMod, i::Int) = basis(F,i)

@doc raw"""
    base_ring(F::AbstractFreeMod)

Return the underlying ring of `F`.
"""
base_ring(F::FreeMod) = (F.R)::base_ring_type(F)

#TODO: Parent - checks everywhere!!!

# the negative of a free module element
-(a::AbstractFreeModElem) = parent(a)(-coordinates(a))

# Addition of free module elements
function +(a::AbstractFreeModElem, b::AbstractFreeModElem)
   check_parent(a, b)
   return parent(a)(coordinates(a)+coordinates(b))
end

# Subtraction of free module elements
function -(a::AbstractFreeModElem, b::AbstractFreeModElem)
    check_parent(a,b)
    return parent(a)(coordinates(a)-coordinates(b))
end

# Equality of free module elements
function (==)(a::AbstractFreeModElem, b::AbstractFreeModElem) 
  if parent(a) !== parent(b)
    return false
  end
  return a.coords == b.coords
end

function hash(a::AbstractFreeModElem, h::UInt)
  b = 0xaa2ba4a32dd0b431 % UInt
  h = hash(typeof(a), h)
  h = hash(parent(a), h)
  return xor(h, b)
end

# A special method for the case where we can safely assume 
# that the coordinates of elements allow hashing.
function hash(a::AbstractFreeModElem{<:MPolyRingElem{<:FieldElem}}, h::UInt)
  b = 0xaa2ba4a32dd0b431 % UInt
  h = hash(typeof(a), h)
  h = hash(parent(a), h)
  h = hash(coordinates(a), h)
  return xor(h, b)
end

# Once we have a suitable implementation of an in-place version 
# of simplify! for sparse vectors of quotient ring elements, we 
# should probably enable the hash method below.
# function hash(a::AbstractFreeModElem{<:MPolyQuoRingElem}, h::UInt)
#   simplify!(a)
#   b = 0xaa2ba4a32dd0b431 % UInt
#   h = hash(typeof(a), h)
#   h = hash(parent(a), h)
#   h = hash(coordinates(a), h)
#   return xor(h, b)
# end

simplify(a::AbstractFreeModElem{<:MPolyQuoRingElem}) = parent(a)(map_entries(simplify, coordinates(a)))

function Base.deepcopy_internal(a::AbstractFreeModElem, dict::IdDict)
  return parent(a)(deepcopy_internal(coordinates(a), dict))
end

# scalar multiplication with polynomials, integers
function *(a::MPolyDecRingElem, b::AbstractFreeModElem)
  @req parent(a) === base_ring(parent(b)) "elements not compatible"
  return parent(b)(a*coordinates(b))
end

function *(a::MPolyRingElem, b::AbstractFreeModElem) 
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  return parent(b)(a*coordinates(b))
end

function *(a::RingElem, b::AbstractFreeModElem) 
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  return parent(b)(a*coordinates(b))
end

*(a::Int, b::AbstractFreeModElem) = parent(b)(a*coordinates(b))
*(a::Integer, b::AbstractFreeModElem) = parent(b)(base_ring(parent(b))(a)*coordinates(b))
*(a::QQFieldElem, b::AbstractFreeModElem) = parent(b)(base_ring(parent(b))(a)*coordinates(b))

@doc raw"""
    zero(F::AbstractFreeMod)

Return the zero element of  `F`.
"""
zero(F::AbstractFreeMod) = F(sparse_row(base_ring(F), Tuple{Int, elem_type(base_ring(F))}[]))

@doc raw"""
    parent(a::AbstractFreeModElem)

Return the free module where `a` lives in.
"""
parent(a::AbstractFreeModElem) = a.parent

@doc raw"""
    is_zero(f::AbstractFreeModElem)

Return `true` if `f` is zero, `false` otherwise.
"""
is_zero(f::AbstractFreeModElem) = iszero(coordinates(f))

simplify!(a::FreeModElem) = a
simplify(a::FreeModElem) = a
