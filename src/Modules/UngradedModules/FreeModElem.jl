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

function FreeModElem(i::Int, x::T, parent::FreeMod{T}) where T
  @assert 1 <= i <= rank(parent)
  sparse_coords = sparse_row(base_ring(parent), [i], [x])
  return FreeModElem{T}(sparse_coords, parent)
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
#julia> R, (x,y) = polynomial_ring(QQ, [:x, :y])
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
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> F = FreeMod(R,3)
Free module of rank 3 over R

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
  for (pos, val) in coordinates(e)
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
  bas = Vector{elem_type(F)}(undef, dim(F))
  e = one(base_ring(F))
  for i=1:dim(F)
    s = sparse_row(base_ring(F), [i], [e])
    bas[i] = F(s)
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
  e = one(base_ring(F))
  s = sparse_row(base_ring(F), [i], [e])
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
  return coordinates(a) == coordinates(b)
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
*(a::Any, b::AbstractFreeModElem) = parent(b)(base_ring(parent(b))(a)*coordinates(b))

function *(a::T, b::AbstractFreeModElem{T}) where {T <: AdmissibleModuleFPRingElem}
  @assert is_left(parent(b)) "left multiplication is not defined for non-left module $(parent(b))"
  parent(a) === base_ring(parent(b)) && return parent(b)(a*coordinates(b))
  return parent(b)(base_ring(parent(b))(a)*coordinates(b))
end

function *(b::AbstractFreeModElem{T}, a::T) where {T <: AdmissibleModuleFPRingElem}
  @assert is_right(parent(b)) "right multiplication not defined for non-right module $(parent(b))"
  error("right multiplication is not supported at the moment")
end

function *(b::AbstractFreeModElem{T}, a::Any) where {T <: RingElem}
  error("scalar multiplication from the right is not yet supported")
end

# Methods to determine whether a module is a left-, right-, or bi-module. 
# We plan to have flags set for this. But for the moment the generic code only supports left-multiplication, 
# so we can decide this from the type alone. How we do it in the long run is not yet decided, but in either case
# we want to use these functions to decide as they are already there for ideals. 
is_left(M::ModuleFP) = is_left(typeof(M))
is_left(::Type{T}) where {RET<:RingElem, T<:ModuleFP{RET}} = true
is_left(::Type{T}) where {RET<:AdmissibleModuleFPRingElem, T<:ModuleFP{RET}} = true # Left multiplication is generically supported

is_right(M::ModuleFP) = is_right_module(typeof(M))
is_right(::Type{T}) where {RET<:RingElem, T<:ModuleFP{RET}} = true
is_right(::Type{T}) where {RET<:AdmissibleModuleFPRingElem, T<:ModuleFP{RET}} = false # Right multiplication is not supported by the generic code at the moment, but we plan to do so eventually. 

is_two_sided(M::ModuleFP) = is_right_module(typeof(M))
is_two_sided(::Type{T}) where {RET<:RingElem, T<:ModuleFP{RET}} = true
is_two_sided(::Type{T}) where {RET<:AdmissibleModuleFPRingElem, T<:ModuleFP{RET}} = false # see above


@doc raw"""
    zero(F::AbstractFreeMod)

Return the zero element of  `F`.
"""
zero(F::AbstractFreeMod) = F(sparse_row(base_ring(F)))

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
