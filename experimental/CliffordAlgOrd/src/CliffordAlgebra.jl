export CliffordAlgebra,
  CliffordAlgebraElem,
  clifford_algebra,
  even_coefficients,
  odd_coefficients,
  even_part,
  odd_part,
  basis_of_center,
  basis_of_centroid,
  quadratic_discriminant,
  disq

################################################################################
#
#  Datatypes
#
################################################################################

### Algebra ###
# Data structure for Clifford algebras. The type variable 'T' represents the element type
# of the base ring, i.e. it is usually QQFieldElem or AbsSimpleNumFieldElem.
# The type variable 'S' represents the type of the Gram matrix of the underlying quadratic space,
# i.e. it is usually QQMatrix or AbstactAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}.
mutable struct CliffordAlgebra{T,S} <: Hecke.AbstractAssociativeAlgebra{T}
  base_ring::Ring
  space::Hecke.QuadSpace{K,S} where {K} # K = parent_type(T), e.g. T is of type AbsSimpleNumFieldElem and K is of type AbsSimpleNumField
  gram::S
  dim::Int
  basis_of_centroid::Any # Vector{elem_type(C)}, with C an instance of CliffordAlgebra
  disq::T
  basis_of_center::Any # Vector{elem_type(C)}, with C an instance of CliffordAlgebra

  #Return the Clifford algebra of the quadratic space 'qs' 
  function CliffordAlgebra{T,S}(qs::Hecke.QuadSpace{K,S}) where {T,S,K}
    gram = gram_matrix(qs)
    return new{T,S}(base_ring(qs), qs, gram, 2^ncols(gram))
  end
end

### Elements ###
# Data structure for the elements of a Clifford algebra. The type variables serve the
# same purpose as they do for Clifford algebras.
mutable struct CliffordAlgebraElem{T,S} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::CliffordAlgebra{T, S}
  coeffs::Vector{T}
  even_coeffs::Vector{T}
  odd_coeffs::Vector{T}

  #Return the 0-element of the Clifford algebra C
  function CliffordAlgebraElem{T,S}(C::CliffordAlgebra{T,S}) where {T,S}
    newelt = new{T,S}(C, fill(C.base_ring(), C.dim))
    _set_even_odd_coefficients!(newelt)
    return newelt
  end

  CliffordAlgebraElem(C::CliffordAlgebra) =
    CliffordAlgebraElem{elem_type(C.base_ring),typeof(C.gram)}(C)

  #Return the element in the Clifford algebra C with coefficient vector coeff wrt. the canonical basis
  function CliffordAlgebraElem{T,S}(
    C::CliffordAlgebra{T,S}, coeff::Vector{R}
  ) where {T,S,R<:FieldElem}
    @req length(coeff) == C.dim "invalid length of coefficient vector"
    newelt = new{T,S}(C, coeff)
    _set_even_odd_coefficients!(newelt)
    return newelt
  end

  CliffordAlgebraElem(C::CliffordAlgebra{T,S}, coeff::Vector{T}) where {T,S} =
    CliffordAlgebraElem{elem_type(C.base_ring),typeof(C.gram)}(C, coeff)

  function CliffordAlgebraElem(C::CliffordAlgebra{T,S}, coeff::Vector{R}) where {T,S,R}
    K = C.base_ring
    @req __can_convert_coefficients(coeff, K) "entries of coefficient vector are not contained in $(K)"
    return CliffordAlgebraElem{elem_type(C.base_ring),typeof(C.gram)}(
      C, K.(coeff)
    )
  end
end

elem_type(::Type{CliffordAlgebra{T, S}}) where {T, S} = CliffordAlgebraElem{T, S}

parent_type(::Type{CliffordAlgebraElem{T, S}}) where {T, S} = CliffordAlgebra{T, S}

base_ring_type(::Type{CliffordAlgebra{T, S}}) where {T, S} = parent_type(T)

is_domain_type(::Type{CliffordAlgebraElem{T, S}}) where {T, S} = false

is_exact_type(::Type{CliffordAlgebraElem{T, S}}) where {T, S} = true

################################################################################
#
#  Construction
#
################################################################################

### Algebra ###

@doc raw"""
    clifford_algebra(qs::QuadSpace) -> CliffordAlgebra

Return the Clifford algebra of the quadratic space 'qs'.
"""
clifford_algebra(qs::Hecke.QuadSpace) =
  CliffordAlgebra{elem_type(base_ring(qs)),typeof(gram_matrix(qs))}(qs)

### Elements ###
(C::CliffordAlgebra)() = CliffordAlgebraElem(C)

function (C::CliffordAlgebra)(a::R) where {R<:RingElem}
  res = fill(zero(a), dim(C))
  res[1] = a
  return CliffordAlgebraElem(C, res)
end

function (C::CliffordAlgebra)(a::R) where {R<:Number}
  res = fill(zero(a), dim(C))
  res[1] = a
  return CliffordAlgebraElem(C, res)
end

function (C::CliffordAlgebra)(a::CliffordAlgebraElem)
  @req parent(a) === C "The element does not lie in the Clifford algebra"
  return a
end

(C::CliffordAlgebra)(coeff::Vector{R}) where {R} = CliffordAlgebraElem(C, coeff)

################################################################################
#
#  Basic field access
#
################################################################################

### Algebra ###
@doc raw"""
    base_ring(C::CliffordAlgebra) -> Ring

Return the base ring of the Clifford algebra $C$.
"""
base_ring(C::CliffordAlgebra) = C.base_ring::base_ring_type(typeof(C))

@doc raw"""
    space(C::CliffordAlgebra) -> QuadSpace

Return the underlying quadratic space of the Clifford algebra $C$.
"""
space(C::CliffordAlgebra{T, S}) where {T, S} = C.space::Hecke.QuadSpace{parent_type(T), S}

@doc raw"""
    gram_matrix(C::CliffordAlgebra) -> MatElem

Return the Gram matrix of the free quadratic module with Clifford algebra $C$.
"""
gram_matrix(C::CliffordAlgebra) = C.gram

@doc raw"""
    dim(C::CliffordAlgebra) -> Int

Return the dimension of the Clifford algebra $C$ over its base ring.
"""
dim(C::CliffordAlgebra) = C.dim

### Elements ###
@doc raw"""
    parent(x::CliffordAlgebraElem) -> CliffordAlgebra

Return the Clifford algebra containing $x$.
"""
parent(x::CliffordAlgebraElem) = x.parent

@doc raw"""
    coefficients(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical basis of its parent Clifford algebra.  
"""
coefficients(x::CliffordAlgebraElem) = x.coeffs

@doc raw"""
    even_coefficients(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical basis of its parent Clifford algebra,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coefficients(x::CliffordAlgebraElem)
  if isdefined(x, :even_coeffs)
    return x.even_coeffs
  end
  _set_even_odd_coefficients!(x)
  return x.even_coeffs
end

@doc raw"""
    odd_coefficients(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical basis of its parent Clifford algebra,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coefficients(x::CliffordAlgebraElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs
  end
  _set_even_odd_coefficients!(x)
  return x.odd_coeffs
end

################################################################################
#
#  String I/O
#
################################################################################

### Algebra ###
function Base.show(io::IO, ::MIME"text/plain", C::CliffordAlgebra)
  io = pretty(io)
  print(io, "Clifford algebra of quadratic space with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent(), "\ndefined over $(base_ring(C))")
end

Base.show(io::IO, C::CliffordAlgebra) = print(io, "Clifford algebra over $(base_ring(C))")

### Elements ###
function Base.show(io::IO, x::CliffordAlgebraElem)
  print(io, "[")
  foreach(y -> print(io, "$y "), coefficients(x)[1:(end - 1)])
  print(io, "$(coefficients(x)[end])]")
end

################################################################################
#
#  basic functionality
#
################################################################################

@doc raw"""
    zero(C::CliffordAlgebra) -> CliffordAlgebraElem

Return the additive identity of the Clifford algebra $C$.
"""
zero(C::CliffordAlgebra) = C()

@doc raw"""
    one(C::CliffordAlgebra) -> CliffordAlgebraElem

Return the multiplicative identity of the Clifford algebra $C$.
"""
function one(C::CliffordAlgebra)
  res = C()
  res[1] = base_ring(C)(1)
  return res
end

@doc raw"""
    basis(C::CliffordAlgebra, i::Int) -> CliffordAlgebraElem

Return the $i$-th canonical basis vector of the Clifford algebra $C$.
"""
function basis(C::CliffordAlgebra, i::Int)
  res = C()
  res[i] = base_ring(C)(1)
  return res
end

@doc raw"""
    gen(C::CliffordAlgebra, i::Int) -> CliffordAlgebraElem

Return the $i$-th canonical algebra generator of the Clifford
algebra $C$. This is just the image of the $i$-th basis vector
of the underlying quadratic space under the canonical embedding
into the Clifford algebra.
"""
function gen(C::CliffordAlgebra, i::Int)
  res = C()
  if i <= 0
    throw(BoundsError(res, i))
  end
  res[2^(i - 1) + 1] = base_ring(C)(1)
  return res
end

@doc raw"""
    basis(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the canonical basis of the Clifford algebra $C$.
"""
basis(C::CliffordAlgebra) = map(i -> basis(C, i), 1:dim(C))

@doc raw"""
    gens(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the vector of canonical algebra generators of the Clifford algebra $C$, i.e.,
if `gram_matrix(C)` is the Gram matrix of the underlying quadratic space with respect
to the basis (e_1,...,e_n) then the vector of images under the canonical embedding into
the Clifford algebra is returned.
"""
gens(C::CliffordAlgebra) = map(i -> gen(C, i), 1:_dim_qf(C))

@doc raw"""
    is_commutative(C::CliffordAlgebra) -> Bool

Return `true` if $C$ is commutative and `false` otherwise.
"""
is_commutative(C::CliffordAlgebra) = dim(C) == 1 || dim(C) == 2

################################################################################
#
#  Element Access
#
################################################################################

@doc raw"""
    coeff(x::CliffordAlgebraElem, i::Int) -> FieldElem

Return the `i`-th coefficient of the element `x`.
"""
coeff(x::CliffordAlgebraElem, i::Int64) = coefficients(x)[i]

getindex(x::CliffordAlgebraElem, i::Int64) = coefficients(x)[i]

function setindex!(x::CliffordAlgebraElem, newentry::RingElement, i::Int64) 
  coefficients(x)[i] = newentry
end

################################################################################
#
#  Other functionality
#
################################################################################

@doc raw"""
    basis_of_centroid(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return a basis of the centroid of $C$. Unless `dim(space(C)) = 0`, it consists of
two elements. The first one is the multiplicative identity of $C$. The square of
the second basis element, if present, equals `quadratic_discriminant(C)`.
"""
function basis_of_centroid(C::CliffordAlgebra)
  if isdefined(C, :basis_of_centroid)
    return C.basis_of_centroid::Vector{elem_type(C)}
  end
  n = _dim_qf(C)
  if n == 0
    C.basis_of_centroid = [one(C)]
    return C.basis_of_centroid::Vector{elem_type(C)}
  end
  T = orthogonal_basis(space(C))
  orth_elt = prod(map(i -> sum(map(j -> gen(C, j) * T[i, j], 1:n)), 1:n))
  orth_elt *= denominator(orth_elt)
  C.disq = (orth_elt^2)[1]
  C.basis_of_centroid = [one(C), orth_elt]::Vector{elem_type(C)}
end

@doc raw"""
    quadratic_discriminant(C::CliffordAlgebra) -> FieldElem

Return the quadratic discriminant of $C$ as an element of `base_ring(C)`.
"""
function quadratic_discriminant(C::CliffordAlgebra)
  if isdefined(C, :disq)
    return C.disq
  end
  if dim(space(C)) == 0
    C.disq = one(base_ring(C))
  else
    C.disq = (basis_of_centroid(C)[2]^2)[1]
  end
end

@doc raw"""
    disq(C::CliffordAlgebra) -> FieldElem

Alias for `quadratic_discriminant`.
"""
disq(C::CliffordAlgebra) = quadratic_discriminant(C)

@doc raw"""
    basis_of_center(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return a basis of the center of $C$. It equals `basis_of_centroid(C)`, if and only
if `dim(space(C))` is odd. Otherwise it contains only the
multiplicative identity of $C$.
"""
function basis_of_center(C::CliffordAlgebra)
  if isdefined(C, :basis_of_center)
    return C.basis_of_center::Vector{elem_type(C)}
  end
  if is_odd(dim(space(C)))
    C.basis_of_center = basis_of_centroid(C)::Vector{elem_type(C)}
  else
    C.basis_of_center = [one(C)]::Vector{elem_type(C)}
  end
end

@doc raw"""
    representation_matrix(x::CliffordAlgebraElem, action::Symbol = :left) -> MatElem

Return the representation matrix of the element $x$ with respect to `basis(parent(x))`. The multiplication is from
the left if `action == :left` and from the right if `action == :right`.
"""
function representation_matrix(x::CliffordAlgebraElem, action::Symbol = :left)
  @req (action == :left) || (action == :right) "The action is either $(:left) or $(:right)."
  C = parent(x)
  n = dim(C)
  res = zero_matrix(base_ring(C), n, n)
  if action == :left
    for i in 1:n
      res[i, :] = coefficients(x * basis(C, i))
    end
  elseif action == :right
    for i in 1:n
      res[i, :] = coefficients(basis(C, i) * x)
    end
  end
  return res
end

################################################################################
#
#  unary operators
#
################################################################################

Base.:-(x::CliffordAlgebraElem) = parent(x)(map(y -> -1 * y, coefficients(x)))

################################################################################
#
#  binary operators
#
################################################################################

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

################################################################################
#
#  Equality and hash
#
################################################################################

Base.:(==)(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T} = parent(x) === parent(y) && coefficients(x) == coefficients(y)

function Base.hash(x::CliffordAlgebraElem, h::UInt)
  b = 0x1c4629b4de23b24c % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

################################################################################
#
#  Graded parts
#
################################################################################

@doc raw"""
    even_part(x::CliffordAlgebraElem) -> CliffordAlgebraElem

Return the projection of $x$ onto the even Clifford algebra
"""
even_part(x::CliffordAlgebraElem) = parent(x)(even_coefficients(x))

@doc raw"""
    odd_part(x::CliffordAlgebraElem) -> CliffordAlgebraElem

Return the projection of $x$ onto the odd Clifford algebra.
"""
odd_part(x::CliffordAlgebraElem) = parent(x)(odd_coefficients(x))

################################################################################
#
#  Auxillary functions
#
################################################################################

function _set_even_odd_coefficients!(x::CliffordAlgebraElem)
  R = base_ring(parent(x))
  d = dim(parent(x))
  x.even_coeffs = [ iseven(count_ones(y - 1)) ? x.coeffs[y] : R() for y in 1:d]
  x.odd_coeffs = x.coeffs - x.even_coeffs
  return x
end

function __can_convert_coefficients(coeff::Vector{R}, K::Field) where {R}
  if length(coeff) == 0
    return true
  end
  try
    K.(coeff)  # Try converting the coefficient vector to K
    true
  catch
    false
  end
end

_dim_qf(C::CliffordAlgebra) = ncols(C.gram)

function _mul_aux(x::Vector{T}, y::Vector{T}, gram::MatElem{T}, i::Int) where {T<:RingElement}
  if length(y) == 1
    return x .* y[1]
  elseif is_zero(y)
    return fill(y[1], 2^ncols(gram))
  end
  return _mul_aux(_mul_with_gen(x, i, gram), y[2:2:end], gram, i + 1) +
         _mul_aux(x, y[1:2:end], gram, i + 1)
end

#Implements the right multiplication with 'i'-th generator e_i of the Clifford algebra containing x.
#The result is returned as a coefficient vector for further computations.
_mul_with_gen(x::Vector{T}, i::Int, gram::MatElem{T}) where {T<:RingElement} = sum(
  map(
    char ->
      x[char] .*
      _mul_baseelt_with_gen(char, i, gram),
      1:2^ncols(gram),
  ),
)

#Shift all entries of the vector 'X' to the right by 's' entries without changing the length of X.
function _shift_entries!(X::Vector, s::Int)
  X[(s + 1) : end] = X[1 : (end - s)]
  X[1 : s] = fill(parent(X[1])(), s)
  return X
end

#Right multiplication of the basis element represented by 'char' with the
#'i'-th generator of the Clifford algebra/order containing said basis vector
function _mul_baseelt_with_gen(char::Int, i::Int, gram::MatElem)
  R = base_ring(gram)
  res = fill(R(), 2^ncols(gram))
  if char == 1
    res[char + 2^(i - 1)] = R(1)
    return res
  end
  j = 8 * sizeof(Int) - leading_zeros(char - 1)
  if j < i
    res[char + 2^(i - 1)] = R(1)
    return res
  end
  if j == i
    res[char - 2^(i - 1)] = R(divexact(gram[i, i], 2))
    return res
  end
  res[char - 2^(j - 1)] = gram[i, j]
  res -= _shift_entries!(_mul_baseelt_with_gen(char - 2^(j - 1), i, gram), 2^(j - 1))
end
