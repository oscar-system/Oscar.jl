
################################################################################
#
#  unary operators
#
################################################################################

Base.:-(x::CliffordAlgebraElem) = parent(x)(map(y -> -1 * y, coefficients(x)))
Base.:+(x::CliffordOrderElem) = x
Base.:-(x::CliffordOrderElem) = parent(x)(map(y -> -1 * y, coefficients(x)))
Base.:+(x::ZZCliffordOrderElem) = x
Base.:-(x::ZZCliffordOrderElem) = parent(x)(map(y -> -1 * y, coefficients(x)))

function inv(x::CliffordOrderElem)
  CO = parent(x)
  CA = ambient_algebra(CO)
  xinv = inv(CA(x))
  if !(xinv in CO)
    error("Element is not invertible")
  end
  return CO(xinv)
end

function inv(x::ZZCliffordOrderElem)
  CO = parent(x)
  CA = ambient_algebra(CO)
  xinv = inv(CA(x))
  if !(xinv in CO)
    error("Element is not invertible")
  end
  return CO(xinv)
end

################################################################################
#
#  binary operators
#
################################################################################

##### Algebra #####
function Base.:+(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:FieldElem}
  check_parent(x, y)
  return parent(x)(coefficients(x) .+ coefficients(y))
end

Base.:-(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:FieldElem} = x + -y

function Base.:*(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:FieldElem}
  check_parent(x, y)
  xcoeffs, ycoeffs = copy(coefficients(x)), copy(coefficients(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram_matrix(parent(x)), 1))
end

@doc raw"""
    divexact(x::CliffordAlgebraElem, a::T) where {T<:RingElement} -> CliffordAlgebraElem

Return the element `y` in the Clifford algebra containing $x$ such that $ay = x$,
if it exists. Otherwise an error is raised.
"""
divexact(x::CliffordAlgebraElem, a::T) where {T<:RingElement} =
  parent(x)(divexact.(coefficients(x), a))

##### Order #####

function Base.:+(x::CliffordOrderElem{T, CliffordAlgebra{U, V}},
                  y::CliffordOrderElem{T, CliffordAlgebra{U, V}}) where {T, U<:NumFieldElem, V}
  check_parent(x, y)
  return parent(x)(coefficients(x) .+ coefficients(y))
end

Base.:-(x::CliffordOrderElem{T, CliffordAlgebra{U, V}},
          y::CliffordOrderElem{T, CliffordAlgebra{U, V}}) where {T, U<:NumFieldElem, V} = x + -y

function Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U, V}},
                  y::CliffordOrderElem{T, CliffordAlgebra{U, V}}) where {T, U<:NumFieldElem, V}
  check_parent(x, y)
  xcoeffs, ycoeffs = copy(coefficients(x)), copy(coefficients(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram_matrix(parent(x)), 1))
end

Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U,V}}, a::W) where {T, U<:NumFieldElem, V, W<:FieldElem} = parent(x)(a * algebra(parent(x))(x))
Base.:*(a::W, x::CliffordOrderElem{T, CliffordAlgebra{U,V}}) where {T, U<:NumFieldElem, V, W<:FieldElem} = parent(x)(a * algebra(parent(x))(x))
Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U,V}}, a::Rational{Int}) where {T, U<:NumFieldElem, V} = parent(x)(a .* coefficients(x))
Base.:*(a::Rational{Int}, x::CliffordOrderElem{T, CliffordAlgebra{U,V}}) where {T, U<:NumFieldElem, V} = parent(x)(a .* coefficients(x))

@doc raw"""
    divexact(x::CliffordOrderElem, a::T) where {T<:RingElement} -> CliffordOrderElem

Return the element `y` in the Clifford order containing $x$ such that $ay = x$,
if it exists. Otherwise an error is raised.
"""
function divexact(x::CliffordOrderElem, elt::T) where {T<:RingElement}
  ambalg = algebra(parent(x))
  res = divexact(ambalg(x), elt)
  @req res in parent(x) "Not an exact division"
  return parent(x)(res)
end

### ZZ ###
function Base.:+(x::ZZCliffordOrderElem, y::ZZCliffordOrderElem)
  check_parent(x, y)
  return parent(x)(coefficients(x) .+ coefficients(y))
end

Base.:-(x::ZZCliffordOrderElem, y::ZZCliffordOrderElem) = x + -y

function Base.:*(x::ZZCliffordOrderElem, y::ZZCliffordOrderElem)
  check_parent(x, y)
  xcoeffs, ycoeffs = copy(coefficients(x)), copy(coefficients(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram_matrix(parent(x)), 1))
end

Base.:*(x::ZZCliffordOrderElem, a::QQFieldElem) = parent(x)(a .* coefficients(x))
Base.:*(a::QQFieldElem, x::ZZCliffordOrderElem) = parent(x)(a .* coefficients(x))
Base.:*(x::ZZCliffordOrderElem, a::Rational{Int}) = parent(x)(a .* coefficients(x))
Base.:*(a::Rational{Int}, x::ZZCliffordOrderElem) = parent(x)(a .* coefficients(x))

@doc raw"""
    divexact(x::ZZCliffordOrderElem, a::RingElement) -> ZZCliffordOrderElem

Return the element `y` in the Clifford order containing $x$ such that $ay = x$,
if it exists. Otherwise an error is raised.
"""
function divexact(x::ZZCliffordOrderElem, a::T) where {T<:RingElement}
  ambalg = algebra(parent(x))
  res = divexact(ambalg(x), a)
  @req res in parent(x) "Not an exact division"
  return parent(x)(res)
end

################################################################################
#
#  Equality and hash
#
################################################################################

##### Algebra ######
Base.:(==)(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T} = parent(x) === parent(y) && coefficients(x) == coefficients(y)

function Base.hash(x::CliffordAlgebraElem, h::UInt)
  b = 0x1c4629b4de23b24c % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

##### Order #####
Base.:(==)(x::CliffordOrderElem{T}, y::CliffordOrderElem{T}) where {T} = parent(x) === parent(y) && coefficients(x) == coefficients(y)

function Base.hash(x::CliffordOrderElem, h::UInt)
  b = 0x8f3a7c2b1d4e5f6a % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

### ZZ ###
Base.:(==)(x::ZZCliffordOrderElem, y::ZZCliffordOrderElem) = parent(x) === parent(y) && coefficients(x) == coefficients(y)

function Base.hash(x::ZZCliffordOrderElem, h::UInt)
  b = 0x924e7a492844d7a0 % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

