# module Tropical

export tropical_numbers,
       @tropical,
       determinant,
       permanent

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
mutable struct TropicalNumbers{T} <: Ring
end

# We use the flag isinf to denote +/- infinity
# todo: should this be called TropicalNumber instead?
#   I have no preference, except that it should be consistent with the other libraries,
#   e.g. what are elements of p-adic number rings called?
mutable struct TropicalNumbersElem{T} <: RingElem
  parent::TropicalNumbers{T}
  isinf::Bool
  data::fmpq

  function TropicalNumbersElem(R::TropicalNumbers{T}, isinf::Bool) where {T}
    @assert isinf
    z = new{T}()
    z.isinf = true
    z.parent = R
    return z
  end

  function TropicalNumbersElem(R::TropicalNumbers{T}, x::RingElem) where {T}
    return new{T}(R, false, x)
  end
end

# Type gymnastics

Oscar.elem_type(::Type{TropicalNumbers{T}}) where {T} = TropicalNumbersElem{T}

Oscar.parent_type(::Type{TropicalNumbersElem{T}}) where {T} = TropicalNumbers{T}

################################################################################
#
#  Constructors for the tropical semiring
#
################################################################################

# Invoke via tropical_numbers(max)
tropical_numbers(::typeof(max)) = TropicalNumbers{typeof(max)}()

# Invoke via tropical_numbers(min)
tropical_numbers(::typeof(min)) = TropicalNumbers{typeof(min)}()

################################################################################
#
#  Constructors for elements
#
################################################################################

function (R::TropicalNumbers)(u::TropicalNumbersElem)
  @assert parent(u) === R
  return u
end

function (R::TropicalNumbers)(u::RingElem)
  v = QQ(u)
  @assert parent(v) === QQ
  return TropicalNumbersElem(R, v)
end

function (R::TropicalNumbers)(u::Union{Integer, Rational})
  return TropicalNumbersElem(R, QQ(u))
end

inf(R::TropicalNumbers) = TropicalNumbersElem(R, true)

Oscar.zero(T::TropicalNumbers) = inf(T)

Oscar.one(R::TropicalNumbers) = R(zero(QQ))

Oscar.zero(x::TropicalNumbersElem) = zero(parent(x))

(R::TropicalNumbers)() = zero(R)

################################################################################
#
#  Basic access
#
################################################################################

# The underlying rational number. This is undefined for inf.
data(x::TropicalNumbersElem) = x.data

# Test if something is inf.
isinf(x::TropicalNumbersElem) = x.isinf

Oscar.parent(x::TropicalNumbersElem) = x.parent

# get the underlyling min/max function
fun(x::TropicalNumbers{typeof(min)}) = min

fun(x::TropicalNumbers{typeof(max)}) = max

################################################################################
#
#  Printing
#
################################################################################

# Hook into the fancy printing

# We use (x) for finite values and ±∞ for infinity.
function AbstractAlgebra.expressify(x::TropicalNumbersElem{T}; context = nothing) where {T}
  if isinf(x)
    return T === typeof(min) ? "∞" : "-∞"
  end
  return Expr(:call, "", expressify(data(x), context = context))
end

AbstractAlgebra.expressify(R::TropicalNumbers{typeof(min)}; context = nothing) = "Tropical ring (min)"

AbstractAlgebra.expressify(R::TropicalNumbers{typeof(max)}; context = nothing) = "Tropical ring (max)"

@enable_all_show_via_expressify TropicalNumbersElem

@enable_all_show_via_expressify TropicalNumbers


################################################################################
#
#  Predicates
#
################################################################################

Oscar.iszero(x::TropicalNumbersElem) = isinf(x)

Oscar.isone(x::TropicalNumbersElem) = !isinf(x) && iszero(data(x))

################################################################################
#
#  Equality
#
################################################################################

function Base.:(==)(x::TropicalNumbersElem, y::TropicalNumbersElem)
  (isinf(x) && isinf(y)) && return true
  ((isinf(x) && !isinf(y)) || (!isinf(x) && isinf(y))) && return false
  return data(x) == data(y)
end

################################################################################
#
#  Copying
#
################################################################################

Base.copy(a::TropicalNumbersElem) = a

function Base.deepcopy_internal(x::TropicalNumbersElem, dict::IdDict)
  if !isinf(x)
    return TropicalNumbersElem(x.parent, Base.deepcopy_internal(data(x), dict))
  else
    return inf(parent(x))
  end
end

################################################################################
#
#  Addition
#
################################################################################

function Base.:(+)(x::TropicalNumbersElem{T}, y::TropicalNumbersElem{T}) where {T}
  if isinf(x)
    return deepcopy(y)
  else
    if isinf(y)
      return deepcopy(x)
    else
      return parent(x)(fun(parent(x))(data(x), data(y)))
    end
  end
end

################################################################################
#
#  Multiplication
#
################################################################################

function Base.:(*)(x::TropicalNumbersElem{T}, y::TropicalNumbersElem{T}) where {T}
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
#  Powering
#
################################################################################

function Base.:(^)(a::TropicalNumbersElem, n::Int)
  return Base.power_by_squaring(a, n)
end

################################################################################
#
#  Disable subtraction
#
################################################################################

function Base.:(-)(x::TropicalNumbersElem, y::TropicalNumbersElem...)
  error("Computer says no!")
end

################################################################################
#
#  Unsafe operations
#
################################################################################

Oscar.mul!(x::TropicalNumbersElem, y::TropicalNumbersElem, z::TropicalNumbersElem) = y * z

Oscar.addeq!(y::TropicalNumbersElem, z::TropicalNumbersElem) = y + z

################################################################################
#
#  Permanent
#
################################################################################

# todo: maybe this should be called tropical determinant
#   lest it might crash with the non-tropical notion of determinant
function determinant(x)
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

function permanent(x)
  return determinant(x)
end

################################################################################
#
#  Polynomials
#
################################################################################

# The generic functions use R(1) and R(0), which is bad.
one(R::AbstractAlgebra.Generic.PolyRing{<:TropicalNumbersElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.PolyRing{TropicalNumbersElem{S}}) where {S} = R(zero(base_ring(R)))

function Oscar.PolynomialRing(R::TropicalNumbers, s::Symbol; cached::Bool = true)
   T = elem_type(R)
   parent_obj = Oscar.Generic.PolyRing{T}(R, s, cached)

   return parent_obj, parent_obj([zero(R), one(R)])
end

# Oscar will print zero sums as 0, which we do not want.
# So we have to adjust the printing code for polynomials
function AbstractAlgebra.expressify(@nospecialize(a::PolyElem{<:TropicalNumbersElem}),
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
function AbstractAlgebra.expressify(a::MPolyElem{<:TropicalNumbersElem}, x = symbols(parent(a)); context = nothing)
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

one(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalNumbersElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalNumbersElem}) = R(zero(base_ring(R)))

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
julia> T = tropical_numbers(min);

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
