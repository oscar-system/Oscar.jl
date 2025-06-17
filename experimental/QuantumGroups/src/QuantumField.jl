@doc raw"""
    q_integer(n::Int, q::QuantumFieldElem) -> QuantumFieldElem

Return the q-integer $[n]_q$.

# Examples
```jldoctest
julia> q_integer(3, q)
(q^6 + q^4 + q^2 + 1)/q^3

julia> q_integer(2, q^2)
(q^4 + 1)/q^2
"""
function q_integer(n::Int, q::QuantumFieldElem)
  return div(q^n - q^-n, q - q^-1)
end

@doc raw"""
    q_factorial(n::Int, q::QuantumFieldElem) -> QuantumFieldElem

Return the q-factorial $[n]_q!$.

# Examples
```jldoctest
julia> q_factorial(3, q)
(q^6 + q^4 + q^2 + 1)/q^3

julia> q_integer(2, q^2)
(q^4 + 1)/q^2
"""
function q_factorial(n::Int, q::QuantumFieldElem)
  return get!(
    parent(q).q_factorial::Dict{Tuple{Int,QuantumFieldElem},QuantumFieldElem}, (n, q)
  ) do
    prod(div!(q^k - q^-k, q - q^-1) for k in 2:n; init=one(q))
  end
end

@doc raw"""
    q_binomial(n::Int, k::Int, q::QuantumFieldElem) -> QuantumFieldElem

Return the q-binomial coefficient.

# Examples
```jldoctest
julia> q_binomial(4, 2, q)
(q^8 + q^6 + 2*q^4 + q^2 + 1)//q^4

julia> q_binomial(5, 3, q^2)
(q^24 + q^20 + 2*q^16 + 2*q^12 + 2*q^8 + q^4 + 1)//q^12
"""
function q_binomial(n::Int, k::Int, q::QuantumFieldElem)
  z = one(q)
  for i in 0:(k - 1)
    z = mul!(z, q^(n - i) - q^(i - n))
  end
  for i in 1:k
    z = div!(z, q^i - q^-i)
  end
  return z
end

function image!(z::QuantumFieldElem, ::Type{_BarAutomorphism}, x::QuantumFieldElem)
  s = degree(x.d.d.num) - degree(x.d.d.den)
  if s < 0
    z.d.d.num = reverse!(z.d.d.num, x.d.d.num, length(x.d.d.num) - s)
    z.d.d.den = reverse!(z.d.d.den, x.d.d.den, length(x.d.d.den))
  else
    z.d.d.num = reverse!(z.d.d.num, x.d.d.num, length(x.d.d.num))
    z.d.d.den = reverse!(z.d.d.den, x.d.d.den, length(x.d.d.den) + s)
  end
  return z
end

function image!(z::QuantumFieldElem, ::_BarAutomorphism, x::QuantumFieldElem)
  return image!(z, _BarAutomorphism, x)
end

###############################################################################
#
#   Constructors
#
###############################################################################

function quantum_field()
  return quantum_field(:q)
end

@doc raw"""
    quantum_field(s::Symbol) -> QuantumField, QuantumFieldElem

# Examples
```jldoctest
julia> QF, v = quantum_field(:v)
(Rational function field over QQ, v)
```
"""
function quantum_field(s::Symbol)
  QQq, _ = rational_function_field(QQ, s; cached=false)
  QF = QuantumField(QQq, Dict{Tuple{Int,QuantumFieldElem},QuantumFieldElem}())
  return QF, gen(QF)
end

function (QF::QuantumField)()
  return zero(QF)
end

function (QF::QuantumField)(n::Integer)
  return QuantumFieldElem(QF, QF.d(n))
end

###############################################################################
#
#   Accessors
#
###############################################################################

function parent(x::QuantumFieldElem)
  return x.parent
end

###############################################################################
#
#   Copying / Hashing
#
###############################################################################

function Base.hash(x::QuantumFieldElem, h::UInt)
  b = 0xf37b7f4ddd4dbe54 % UInt
  h = hash(parent(x), h)
  h = hash(x.d, h)

  return xor(h, b)
end

function deepcopy_internal(x::QuantumFieldElem, dict::IdDict)
  return get!(dict, x) do
    QuantumFieldElem(
      x.parent,
      deepcopy_internal(x.d, dict),
    )
  end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function gen(QF::QuantumField)
  return QuantumFieldElem(QF, gen(QF.d))
end

function one(QF::QuantumField)
  return QuantumFieldElem(QF, one(QF.d))
end

function one!(x::QuantumFieldElem)
  x.d = one!(x.d)
  return x
end

function zero(QF::QuantumField)
  return QuantumFieldElem(QF, zero(QF.d))
end

function zero!(x::QuantumFieldElem)
  x.d = zero!(x.d)
  return x
end

function inv(x::QuantumFieldElem)
  return QuantumFieldElem(x.parent, inv(x.d))
end

function inv!(z::QuantumFieldElem, x::QuantumFieldElem)
  z.d = inv!(z.d, x.d)
  return z
end

function neg!(z::QuantumFieldElem, x::QuantumFieldElem)
  z.d = neg!(z.d, x.d)
  return z
end

function set!(x::QuantumFieldElem, y::QuantumFieldElem)
  x.d.d.num = set!(x.d.d.num, y.d.d.num)
  x.d.d.den = set!(x.d.d.den, y.d.d.den)
  return x
end

function characteristic(QF::QuantumField)
  return characteristic(QF.d)
end

function is_unit(x::QuantumFieldElem)
  return is_unit(x.d)
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function Base.:+(x::QuantumFieldElem, y::QuantumFieldElem)
  check_parent(x, y)
  return add!(zero(x), x, y)
end

function Base.:(/)(x::QuantumFieldElem, y::QuantumFieldElem)
  return divexact(x, y)
end

function Base.:(//)(x::QuantumFieldElem, y::QuantumFieldElem)
  return divexact(x, y)
end

function Base.:*(x::QuantumFieldElem, y::QuantumFieldElem)
  check_parent(x, y)
  return mul!(zero(x), x, y)
end

function Base.:-(x::QuantumFieldElem, y::QuantumFieldElem)
  check_parent(x, y)
  return sub!(zero(x), x, y)
end

function Base.:-(x::QuantumFieldElem)
  return neg!(zero(x), x)
end

function Base.:^(x::QuantumFieldElem, n::Int)
  return QuantumFieldElem(parent(x), x.d^n)
end

function divexact(x::QuantumFieldElem, y::QuantumFieldElem)
  check_parent(x, y)
  return div!(zero(x), x, y)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(x::QuantumFieldElem, y::QuantumFieldElem)
  return x.d == y.d
end

function isone(x::QuantumFieldElem)
  return isone(x.d)
end

function iszero(x::QuantumFieldElem)
  return iszero(x.d)
end

###############################################################################
#
#   Printing
#
###############################################################################

function Base.show(io::IO, QF::QuantumField)
  print(io, "Rational function field over QQ")
end

function expressify(x::QuantumFieldElem; context=nothing)
  return expressify(x.d; context=context)
end

@enable_all_show_via_expressify QuantumFieldElem

###############################################################################
#
#   Type system
#
###############################################################################

function elem_type(::Type{QuantumField})
  return QuantumFieldElem
end

function parent_type(::Type{QuantumFieldElem})
  return QuantumField
end

function is_domain_type(::Type{QuantumFieldElem})
  return true
end

###############################################################################
#
#   Unsafe methods
#
###############################################################################

function add!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = add!(z.d, x.d, y.d)
  return z
end

function addmul!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = addmul!(z.d, x.d, y.d)
  return z
end

function div!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = div!(z.d, x.d, y.d)
  return z
end

function mul!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = mul!(z.d, x.d, y.d)
  return z
end

function pow!(z::QuantumFieldElem, x::QuantumFieldElem, n::Int)
  z.d = pow!(z.d, x.d, n)
  return z
end

function sub!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = sub!(z.d, x.d, y.d)
  return z
end

function submul!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = submul!(z.d, x.d, y.d)
  return z
end

###############################################################################
#
#   Conformance 
#
###############################################################################

function ConformanceTests.generate_element(QF::QuantumField)
  return QuantumFieldElem(ConformanceTests.generate_element(QF.d))
end
