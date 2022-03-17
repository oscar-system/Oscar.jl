# module Tropical

export tropical_semiring,
       @tropical,
       convention,
       det

# using Reexport

# @reexport using Oscar

# import Oscar: expressify, Ring, @enable_all_show_via_expressify, zero, one, iszero, isone, PolynomialRing

################################################################################
#
#  Types
#
################################################################################

# We use T to record whether we are in the min/max case
# T is either typeof(min) or typeof(max)
struct TropicalSemiring{T} <: Field
end

# We use the flag isinf to denote +/- infinity
# todo: should this be called TropicalNumber instead?
#   I have no preference, except that it should be consistent with the other libraries,
#   e.g. what are elements of p-adic number rings called?
mutable struct TropicalSemiringElem{T} <: FieldElem
  parent::TropicalSemiring{T}
  isinf::Bool
  data::fmpq

  function TropicalSemiringElem(R::TropicalSemiring{T}, isinf::Bool) where {T}
    @assert isinf
    z = new{T}()
    z.isinf = true
    z.parent = R
    return z
  end

  function TropicalSemiringElem(R::TropicalSemiring{T}, x::RingElem) where {T}
    return new{T}(R, false, x)
  end
end

# Type gymnastics

Oscar.elem_type(::Type{TropicalSemiring{T}}) where {T} = TropicalSemiringElem{T}

Oscar.parent_type(::Type{TropicalSemiringElem{T}}) where {T} = TropicalSemiring{T}

################################################################################
#
#  Constructors for the tropical semiring
#
################################################################################

@doc Markdown.doc"""

    tropical_semiring(M::Union{typeof(min),typeof(max)}=min)

The tropical semiring with min (default) or max.

!!! warning
    There is no substraction in the tropical semiring. Any substraction of two tropical numbers will yield an error.

# Examples
Basic arithmetics with tropical numbers:
```jldoctest
julia> T = tropical_semiring() # = tropical_semiring(min)
Tropical semiring (min)

julia> T = tropical_semiring(max)
Tropical semiring (max)

julia> 0*T(3) + 1*T(1)^2 + inf(T) # = max(0+3,1+2*1,-∞)
(3)

julia> T(0) == 0    # checks whether the tropical number is 0
true

julia> iszero(T(0)) # checks whether the tropical number is neutral element of addition
false
```

Tropical polynomials:
```jldoctest
julia> T = tropical_semiring()
Tropical semiring (min)

julia> Tx,(x1,x2) = PolynomialRing(T,3)
(Multivariate Polynomial Ring in x1, x2, x3 over Tropical semiring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}[x1, x2, x3])

julia> f = x1 + -1*x2 + 0
x1 + (-1)*x2 + (0)

julia> evaluate(f,[T(-1//2),T(1//2)]) # warning: omitting T(0) gives an error
(-1//2)
```

Tropical matrices:
```jldoctest
julia> T = tropical_semiring()
Tropical semiring (min)

julia> A = [T(0) inf(T); inf(T) T(0)] # = tropical identity matrix
2×2 Matrix{Oscar.TropicalSemiringElem{typeof(min)}}:
 (0)  ∞
 ∞    (0)

julia> 2*A
2×2 Matrix{Oscar.TropicalSemiringElem{typeof(min)}}:
 (2)  ∞
 ∞    (2)

julia> A*A
2×2 Matrix{Oscar.TropicalSemiringElem{typeof(min)}}:
 (0)  ∞
 ∞    (0)

julia> det(A)
(0)
```
"""
tropical_semiring() = TropicalSemiring{typeof(min)}()
tropical_semiring(::typeof(max)) = TropicalSemiring{typeof(max)}()
tropical_semiring(::typeof(min)) = TropicalSemiring{typeof(min)}()

################################################################################
#
#  Constructors for elements
#
################################################################################

function (R::TropicalSemiring)(u::TropicalSemiringElem)
  @assert parent(u) === R
  return u
end

function (R::TropicalSemiring)(u::RingElem)
  v = QQ(u)
  @assert parent(v) === QQ
  return TropicalSemiringElem(R, v)
end

function (R::TropicalSemiring)(u::Union{Integer, Rational})
  return TropicalSemiringElem(R, QQ(u))
end

inf(R::TropicalSemiring) = TropicalSemiringElem(R, true)

Oscar.zero(T::TropicalSemiring) = inf(T)

Oscar.one(R::TropicalSemiring) = R(zero(QQ))

Oscar.zero(x::TropicalSemiringElem) = zero(parent(x))

(R::TropicalSemiring)() = zero(R)

################################################################################
#
#  Basic access
#
################################################################################

# The underlying rational number. This is undefined for inf.
data(x::TropicalSemiringElem) = x.data

# Test if something is inf.
isinf(x::TropicalSemiringElem) = x.isinf

Oscar.parent(x::TropicalSemiringElem{T}) where T = TropicalSemiring{T}()

@doc Markdown.doc"""
    convention(T::TropicalSemiring)

Returns `min` if `T` is the min tropical semiring,
returns `max` if `T` is the max tropical semiring.

# Examples
```jldoctest
julia> T = tropical_semiring(min)
Tropical semiring (min)

julia> convention(T)
min (generic function with 12 methods)

julia> T = tropical_semiring(max)
Tropical semiring (max)

julia> convention(T)
max (generic function with 14 methods)
```
"""
convention(x::TropicalSemiring{typeof(min)}) = min
convention(x::TropicalSemiring{typeof(max)}) = max

################################################################################
#
#  Printing
#
################################################################################

# Hook into the fancy printing

# We use (x) for finite values and ±∞ for infinity.
function AbstractAlgebra.expressify(x::TropicalSemiringElem{T}; context = nothing) where {T}
  if isinf(x)
    return T === typeof(min) ? "∞" : "-∞"
  end
  return Expr(:call, "", expressify(data(x), context = context))
end

AbstractAlgebra.expressify(R::TropicalSemiring{typeof(min)}; context = nothing) = "Tropical semiring (min)"

AbstractAlgebra.expressify(R::TropicalSemiring{typeof(max)}; context = nothing) = "Tropical semiring (max)"

@enable_all_show_via_expressify TropicalSemiringElem

@enable_all_show_via_expressify TropicalSemiring


################################################################################
#
#  Predicates
#
################################################################################

Oscar.iszero(x::TropicalSemiringElem) = isinf(x)

Oscar.isone(x::TropicalSemiringElem) = !isinf(x) && iszero(data(x))

################################################################################
#
#  Equality
#
################################################################################

function Base.:(==)(x::TropicalSemiringElem, y::TropicalSemiringElem)
  (isinf(x) && isinf(y)) && return true
  ((isinf(x) && !isinf(y)) || (!isinf(x) && isinf(y))) && return false
  return data(x) == data(y)
end

################################################################################
#
#  Comparison
#    * min is ordered as usual: -inf < ... < -1 < 0 < 1 < ...
#    * max is ordered reversely:       ... > -1 > 0 > 1 > ... > inf
#    (see Section 2.7 in Joswig: "Essentials of Tropical Combinatorics"
#
################################################################################

function isless(x::TropicalSemiringElem{typeof(min)}, y::TropicalSemiringElem{typeof(min)})
  if isinf(x)
    return false # x=-inf, no y is smaller
  end
  if isinf(y)
    return true # y=-inf, smaller than all x except x=-inf, which was handled above
  end
  return data(x) < data(y)
end

function isless(x::TropicalSemiringElem{typeof(max)}, y::TropicalSemiringElem{typeof(max)})
  if isinf(x)
    return false # x=inf, no y is smaller
  end
  if isinf(y)
    return true # y=inf, smaller than all x except x=inf, which was handled above
  end
  return data(x) > data(y)
end

################################################################################
#
#  Copying
#
################################################################################

Base.copy(a::TropicalSemiringElem) = a

function Base.deepcopy_internal(x::TropicalSemiringElem, dict::IdDict)
  if !isinf(x)
    return TropicalSemiringElem(x.parent, Base.deepcopy_internal(data(x), dict))
  else
    return inf(parent(x))
  end
end

################################################################################
#
#  Addition
#
################################################################################

function Base.:(+)(x::TropicalSemiringElem{T}, y::TropicalSemiringElem{T}) where {T}
  isinf(x) && return deepcopy(y)
  isinf(y) && return deepcopy(x)
  return parent(x)(convention(parent(x))(data(x), data(y)))
end

################################################################################
#
#  Multiplication
#
################################################################################

function Base.:(*)(x::TropicalSemiringElem{T}, y::TropicalSemiringElem{T}) where {T}
  if isinf(x)
    return x
  else
    if isinf(y)
      return y
    else
      return parent(x)(data(x) + data(y))
    end
  end
end

################################################################################
#
#  Division / Inversion
#
################################################################################

function divexact(a::TropicalSemiringElem{T}, b::TropicalSemiringElem{T}) where T <: FieldElement
    if iszero(b)
        error("dividing by (tropical) zero")
    end
    if iszero(a)
        return a
    end
    return parent(a)(data(a)-data(b))
end

function inv(a::TropicalSemiringElem{T}) where T <: FieldElem
    if iszero(a)
        error("inverting (tropical) zero")
    end
    return parent(a)(-data(a))
end

function Base.:(//)(a::TropicalSemiringElem{T}, b::TropicalSemiringElem{T}) where {T}
    if iszero(b)
        error("dividing by (tropical) zero")
    end
    if iszero(a)
        return a
    end
    return parent(a)(data(a)-data(b))
end

################################################################################
#
#  Powering
#
################################################################################

function Base.:(^)(a::TropicalSemiringElem, n::Int)
  return Base.power_by_squaring(a, n)
end

################################################################################
#
#  Disable subtraction
#
################################################################################

function Base.:(-)(x::TropicalSemiringElem, y::TropicalSemiringElem...)
  error("Tropical subtraction not defined (use tropical division for classical substraction)")
end

################################################################################
#
#  Unsafe operations
#
################################################################################

Oscar.mul!(x::TropicalSemiringElem, y::TropicalSemiringElem, z::TropicalSemiringElem) = y * z

Oscar.addeq!(y::TropicalSemiringElem, z::TropicalSemiringElem) = y + z

################################################################################
#
#  Conversions
#
################################################################################

# converts an element x of the tropical semiring to Int
#   throws an error if x is infinity / -infinity
#   preserves ordering if preserve_ordering==true (default: false),
#   meaning that given x in the max-plus semiring returns -x
#   meaning that the actual valuation of x is returned
function (::Type{Int})(x::TropicalSemiringElem{typeof(min)}; preserve_ordering::Bool=false)
  @assert !iszero(x)
  return Int(ZZ(data(x)))
end
function (::Type{Int})(x::TropicalSemiringElem{typeof(max)}; preserve_ordering::Bool=false)
  @assert !iszero(x)
  if preserve_ordering
    return -Int(ZZ(data(x)))
  end
  return Int(ZZ(data(x)))
end

################################################################################
#
#  Determinant
#
################################################################################

function det(x::AbstractAlgebra.Generic.MatSpaceElem{Oscar.TropicalSemiringElem{T}}) where {T}
  R = base_ring(x)
  S = AbstractAlgebra.SymmetricGroup(nrows(x))
  res = zero(R)
  for s in S
    o = one(R)
    for i in 1:nrows(x)
      o = o * x[i, s[i]]
    end
    res = res + o
  end
  return res
end

function det(x::Matrix{Oscar.TropicalSemiringElem{T}}) where {T}
  if length(x)==0
    if T == typeof(min)
      return tropical_semiring(min)(1)
    else
      return tropical_semiring(max)(1)
    end
  end

  return det(matrix(parent(x[1,1]),x))
end

################################################################################
#
#  Polynomials
#
################################################################################

# The generic functions use R(1) and R(0), which is bad.
one(R::AbstractAlgebra.Generic.PolyRing{<:TropicalSemiringElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.PolyRing{TropicalSemiringElem{S}}) where {S} = R(zero(base_ring(R)))

function Oscar.PolynomialRing(R::TropicalSemiring, s::Symbol; cached::Bool = true)
   T = elem_type(R)
   parent_obj = Oscar.Generic.PolyRing{T}(R, s, cached)

   return parent_obj, parent_obj([zero(R), one(R)])
end

# Oscar will print zero sums as 0, which we do not want.
# So we have to adjust the printing code for polynomials
function AbstractAlgebra.expressify(@nospecialize(a::PolyElem{<:TropicalSemiringElem}),
                                    x = var(parent(a)); context = nothing)
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

# As above, now for multivariate polynomials
function AbstractAlgebra.expressify(a::MPolyElem{<:TropicalSemiringElem}, x = symbols(parent(a)); context = nothing)
  if iszero(a)
    return expressify(zero(base_ring(a)), context = context)
  end
  sum = Expr(:call, :+)
  n = nvars(parent(a))
  for (c, v) in zip(coefficients(a), exponent_vectors(a))
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

one(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalSemiringElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalSemiringElem}) = R(zero(base_ring(R)))

################################################################################
#
#  @tropical macro
#
################################################################################

"""
    @tropical(expr)

Translates the expression in the tropical world.

# Examples

```jlexample
julia> T = tropical_semiring(min);

julia> Tx, x = Tropical.PolynomialRing(T, "x" => 1:3);

julia> @tropical min(1, x[1], x[2], 2*x[3])
x[1] + x[2] + x[3]^2 + (1)
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
      x.args[1] = :(Tropical._tropical_mul)
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
  if x isa Union{Integer, Rational, fmpq, fmpz}
    if x isa Rational || x isa fmpq
      _x = ZZ(x)
      return y^_x
    else
      return y^x
    end
  elseif y isa Union{Integer, Rational, fmpq, fmpz}
    if y isa Rational | y isa fmpq
      _y = ZZ(y)
      return x^_y
    else
      return x^y
    end
  else
    error("Cannot convert ", x, " * ", y)
  end
end

# end # module
