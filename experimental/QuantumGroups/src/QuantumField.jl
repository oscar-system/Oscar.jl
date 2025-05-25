function q_factorial(n::Int, q::QuantumFieldElem)
  return get!(
    parent(q).q_factorial::Dict{Tuple{Int,QuantumFieldElem},QuantumFieldElem}, (n, q)
  ) do
    prod(div!(q^k - q^-k, q - q^-1) for k in 2:n; init=one(q))
  end
end

function q_binomial(n::Int, k::Int, q::QuantumFieldElem)
  z = one(q)
  for i in 0:(k - 1)
    z *= q^(n - i) - q^(i - n)
  end
  for i in 1:k
    z /= q^i - q^-i
  end
  return z
end

function _bar!(z::LaurentPolyWrap, x::LaurentPolyWrap)
  reverse!(z.poly, x.poly, length(x.poly))
  z.mindeg = -degree(x.poly) - x.mindeg
  return z
end

function bar!(z::QuantumFieldElem, x::QuantumFieldElem)
  check_parent(z, x)
  z.d.num = _bar!(z.d.num, x.d.num)
  z.d.den = _bar!(z.d.den, x.d.den)
  return z
end

###############################################################################
#
#   Constructors
#
###############################################################################

function quantum_field()
  A, _ = laurent_polynomial_ring(ZZ, :q; cached=false)
  QF = QuantumField(fraction_field(A), Dict{Tuple{Int,QuantumFieldElem},QuantumFieldElem}())
  return QF, gen(QF)
end

function (QF::QuantumField)()
  return zero(QF)
end

function (QF::QuantumField)(n::Int)
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
#   Arithmetic
#
###############################################################################

function Base.hash(x::QuantumFieldElem, h::UInt)
  return hash(x.d, h)
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
#   Arithmetic
#
###############################################################################

function gen(QF::QuantumField)
  return QuantumFieldElem(QF, gen(QF.d))
end

function one(QF::QuantumField)
  return QuantumFieldElem(QF, one(QF.d))
end

function zero(QF::QuantumField)
  return QuantumFieldElem(QF, zero(QF.d))
end

function one!(x::QuantumFieldElem)
  x.d = one!(x.d)
  return x
end

function zero!(x::QuantumFieldElem)
  x.d = zero!(x.d)
  return x
end

function inv(x::QuantumFieldElem)
  return QuantumFieldElem(x.parent, inv(x.d))
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function Base.:+(x::QuantumFieldElem, y::QuantumFieldElem)
  return add!(zero(x), x, y)
end

function Base.:(//)(x::QuantumFieldElem, y::QuantumFieldElem)
  return div!(zero(x), x, y)
end

function Base.:*(x::QuantumFieldElem, y::QuantumFieldElem)
  return mul!(zero(x), x, y)
end

function Base.:-(x::QuantumFieldElem, y::QuantumFieldElem)
  return sub!(zero(x), x, y)
end

function Base.:-(x::QuantumFieldElem)
  return neg!(zero(x), x)
end

function Base.:^(x::QuantumFieldElem, n::Int)
  return QuantumFieldElem(parent(x), x.d^n)
end

function div(x::QuantumFieldElem, y::QuantumFieldElem)
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

function sub!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = sub!(z.d, x.d, y.d)
  return z
end

function submul!(z::QuantumFieldElem, x::QuantumFieldElem, y::QuantumFieldElem)
  z.d = submul!(z.d, x.d, y.d)
  return z
end
