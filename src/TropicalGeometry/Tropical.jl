# module Tropical

export tropical_ring,
       @tropical,
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
mutable struct TropicalRing{T} <: Ring
end

# We use the flag isinf to denote +/- infinity
mutable struct TropicalRingElem{T} <: RingElem
  parent::TropicalRing{T}
  isinf::Bool
  data::fmpq

  function TropicalRingElem(R::TropicalRing{T}, isinf::Bool) where {T}
    @assert isinf
    z = new{T}()
    z.isinf = true
    z.parent = R
    return z
  end

  function TropicalRingElem(R::TropicalRing{T}, x::RingElem) where {T}
    return new{T}(R, false, x)
  end
end

# Type gymnastics

Oscar.elem_type(::Type{TropicalRing{T}}) where {T} = TropicalRingElem{T}

Oscar.parent_type(::Type{TropicalRingElem{T}}) where {T} = TropicalRing{T}

################################################################################
#
#  Constructors for the tropical semiring
#
################################################################################

# Invoke via tropical_ring(max)
tropical_ring(::typeof(max)) = TropicalRing{typeof(max)}()

# Invoke via tropical_ring(min)
tropical_ring(::typeof(min)) = TropicalRing{typeof(min)}()

################################################################################
#
#  Constructors for elements
#
################################################################################

function (R::TropicalRing)(u::TropicalRingElem)
  @assert parent(u) === R
  return u
end

function (R::TropicalRing)(u::RingElem)
  v = QQ(u)
  @assert parent(v) === QQ
  return TropicalRingElem(R, v)
end

function (R::TropicalRing)(u::Union{Integer, Rational})
  return TropicalRingElem(R, QQ(u))
end

inf(R::TropicalRing) = TropicalRingElem(R, true)

Oscar.zero(T::TropicalRing) = inf(T)

Oscar.one(R::TropicalRing) = R(zero(QQ))

Oscar.zero(x::TropicalRingElem) = zero(parent(x))

(R::TropicalRing)() = zero(R)

################################################################################
#
#  Basic access
#
################################################################################

# The underlying rational number. This is undefined for inf.
data(x::TropicalRingElem) = x.data

# Test if something is inf.
isinf(x::TropicalRingElem) = x.isinf

Oscar.parent(x::TropicalRingElem) = x.parent

# get the underlyling min/max function
fun(x::TropicalRing{typeof(min)}) = min

fun(x::TropicalRing{typeof(max)}) = max

################################################################################
#
#  Printing
#
################################################################################

# Hook into the fancy printing

# We use (x) for finite values and ±∞ for infinity.
function AbstractAlgebra.expressify(x::TropicalRingElem{T}; context = nothing) where {T}
  if isinf(x) 
    return T === typeof(min) ? "∞" : "-∞"
  end
  return Expr(:call, "", expressify(data(x), context = context))
end

AbstractAlgebra.expressify(R::TropicalRing{typeof(min)}; context = nothing) = "Tropical ring (min)"

AbstractAlgebra.expressify(R::TropicalRing{typeof(max)}; context = nothing) = "Tropical ring (max)"

@enable_all_show_via_expressify TropicalRingElem

@enable_all_show_via_expressify TropicalRing


################################################################################
#
#  Predicates
#
################################################################################

Oscar.iszero(x::TropicalRingElem) = isinf(x)

Oscar.isone(x::TropicalRingElem) = !isinf(x) && iszero(data(x))

################################################################################
#
#  Equality
#
################################################################################

function Base.:(==)(x::TropicalRingElem, y::TropicalRingElem)
  (isinf(x) && isinf(y)) && return true
  ((isinf(x) && !isinf(y)) || (!isinf(x) && isinf(y))) && return false
  return data(x) == data(y)
end

################################################################################
#
#  Copying
#
################################################################################

Base.copy(a::TropicalRingElem) = a

function Base.deepcopy_internal(x::TropicalRingElem, dict::IdDict)
  if !isinf(x)
    return TropicalRingElem(x.parent, Base.deepcopy_internal(data(x), dict))
  else
    return inf(parent(x))
  end
end

################################################################################
#
#  Addition
#
################################################################################

function Base.:(+)(x::TropicalRingElem{T}, y::TropicalRingElem{T}) where {T}
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

function Base.:(*)(x::TropicalRingElem{T}, y::TropicalRingElem{T}) where {T}
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

function Base.:(^)(a::TropicalRingElem, n::Int)
  return Base.power_by_squaring(a, n)
end

################################################################################
#
#  Disable subtraction
#
################################################################################

function Base.:(-)(x::TropicalRingElem, y::TropicalRingElem...)
  error("Computer says no!")
end

################################################################################
#
#  Unsafe operations
#
################################################################################

Oscar.mul!(x::TropicalRingElem, y::TropicalRingElem, z::TropicalRingElem) = y * z

Oscar.addeq!(y::TropicalRingElem, z::TropicalRingElem) = y + z

################################################################################
#
#  Permanent
#
################################################################################

function permanent(x)
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

################################################################################
#
#  Polynomials
#
################################################################################

# The generic functions use R(1) and R(0), which is bad.
one(R::AbstractAlgebra.Generic.PolyRing{<:TropicalRingElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.PolyRing{TropicalRingElem{S}}) where {S} = R(zero(base_ring(R)))

function Oscar.PolynomialRing(R::TropicalRing, s::Symbol; cached::Bool = true)
   T = elem_type(R)
   parent_obj = Oscar.Generic.PolyRing{T}(R, s, cached)

   return parent_obj, parent_obj([zero(R), one(R)])
end

# Oscar will print zero sums as 0, which we do not want.
# So we have to adjust the printing code for polynomials
function AbstractAlgebra.expressify(@nospecialize(a::PolyElem{<:TropicalRingElem}),
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
function AbstractAlgebra.expressify(a::MPolyElem{<:TropicalRingElem}, x = symbols(parent(a)); context = nothing)
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

one(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalRingElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalRingElem}) = R(zero(base_ring(R)))

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
julia> T = tropical_ring(min);

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
