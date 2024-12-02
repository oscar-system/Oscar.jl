export CliffordAlgebra,
  CliffordAlgebraElem,
  clifford_algebra,
  even_coeff,
  odd_coeff,
  even_part,
  odd_part,
  centroid,
  quadratic_discriminant,
  disq

################################################################################
#
#  Datatypes
#
################################################################################

### Algebra ###
mutable struct CliffordAlgebra{T,S} <: Hecke.AbstractAssociativeAlgebra{T}
  base_ring::Ring
  space::Hecke.QuadSpace{K,S} where {K}
  gram::S
  dim::Int
  centroid::Any
  disq::T
  center::Any

  #Return the Clifford algebra of the quadratic space 'qs' 
  function CliffordAlgebra{T,S}(qs::Hecke.QuadSpace{K,S}) where {T,S,K}
    gram = gram_matrix(qs)
    return new{T,S}(base_ring(qs), qs, gram, 2^ncols(gram))
  end
end

### Elements ###
mutable struct CliffordAlgebraElem{T,S,P} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::P
  coeffs::Vector{T}
  even_coeffs::Vector{T}
  odd_coeffs::Vector{T}

  #Return the 0-element of the Clifford algebra C
  function CliffordAlgebraElem{T,S,P}(C::CliffordAlgebra{T,S}) where {T,S,P}
    newelt = new{T,S,CliffordAlgebra{T,S}}(C, fill(C.base_ring(), C.dim))
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  CliffordAlgebraElem(C::CliffordAlgebra) =
    CliffordAlgebraElem{elem_type(C.base_ring),typeof(C.gram),typeof(C)}(C)

  #Return the element in the Clifford algebra C with coefficient vector coeff wrt. the canonical basis
  function CliffordAlgebraElem{T,S,P}(
    C::CliffordAlgebra{T,S}, coeff::Vector{R}
  ) where {T,S,P,R<:FieldElem}
    @req length(coeff) == C.dim "invalid length of coefficient vector"
    newelt = new{T,S,CliffordAlgebra{T,S}}(C, coeff)
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  CliffordAlgebraElem(C::CliffordAlgebra{T,S}, coeff::Vector{T}) where {T,S} =
    CliffordAlgebraElem{elem_type(C.base_ring),typeof(C.gram),typeof(C)}(C, coeff)

  function CliffordAlgebraElem(C::CliffordAlgebra{T,S}, coeff::Vector{R}) where {T,S,R}
    K = C.base_ring
    @req __can_convert_coeff(coeff, K) "entries of coefficient vector are not contained in $(K)"
    return CliffordAlgebraElem{elem_type(C.base_ring),typeof(C.gram),typeof(C)}(
      C, K.(coeff)
    )
  end
end

elem_type(::Type{CliffordAlgebra{T,S}}) where {T,S} =
  CliffordAlgebraElem{T,S,CliffordAlgebra{T,S}}

parent_type(::Type{CliffordAlgebraElem{T,S,P}}) where {T,S,P} = P

base_ring_type(C::CliffordAlgebra{T,S}) where {T,S} = typeof(base_ring(C))

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
base_ring(C::CliffordAlgebra) = C.base_ring

@doc raw"""
    space(C::CliffordAlgebra) -> QuadSpace

Return the underlying quadratic space of the Clifford algebra $C$.
"""
space(C::CliffordAlgebra{T,S}) where {T,S} = C.space::Hecke.QuadSpace{parent_type(T),S}

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
  coeff(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical basis of its parent Clifford algebra.  
"""
coeff(x::CliffordAlgebraElem) = x.coeffs

@doc raw"""
    even_coeff(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical basis of its parent Clifford algebra,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coeff(x::CliffordAlgebraElem)
  if isdefined(x, :even_coeffs)
    return x.even_coeffs
  end
  _set_even_odd_coeff!(x)
  return x.even_coeffs
end

@doc raw"""
    odd_coeff(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical basis of its parent Clifford algebra,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coeff(x::CliffordAlgebraElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs
  end
  _set_even_odd_coeff!(x)
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
  foreach(y -> print(io, "$y "), coeff(x)[1:(end - 1)])
  print(io, "$(coeff(x)[end])]")
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
  v = fill(base_ring(C)(), dim(C))
  v[1] = base_ring(C)(1)
  return C(v)
end

@doc raw"""
    basis(C::CliffordAlgebra, i::Int) -> CliffordAlgebraElem

Return the $i$-th canonical basis vector of the Clifford algebra $C$.
"""
function basis(C::CliffordAlgebra, i::Int)
  res = CliffordAlgebraElem(C)
  res.coeffs[i] = base_ring(C)(1)
  _set_even_odd_coeff!(res)
  return res
end

@doc raw"""
    gen(C::CliffordAlgebra, i::Int) -> CliffordAlgebraElem

Return the $i$-th canonical multiplicative generator of the Clifford
algebra $C$. This is just the $i$-th basis vector of the underlying
quadratic module.
"""
function gen(C::CliffordAlgebra, i::Int)
  res = CliffordAlgebraElem(C)
  if i <= 0
    res.coeffs[i] #Throws a BoundsError instead of a DomainError for consistency
  end
  res.coeffs[2^(i - 1) + 1] = base_ring(C)(1)
  _set_even_odd_coeff!(res)
  return res
end

@doc raw"""
    basis(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the canonical basis of the Clifford algebra $C$.
"""
basis(C::CliffordAlgebra) = map(i -> basis(C, i), 1:dim(C))

@doc raw"""
    gens(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the vector of canonical multiplicative generators of the Clifford algebra $C$.
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

getindex(x::CliffordAlgebraElem, i::Int64) = coeff(x)[i]

################################################################################
#
#  Other functionality
#
################################################################################

@doc raw"""
    centroid(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the centroid of $C$. Unless `dim(space(C)) = 0`, it is always two-dimensional,
so it is returned as a vector containing the basis elements. The first one is the
multiplicative identity of $C$. The square of the second basis element, if present,
equals `quadratic_discriminant(C)`.
"""
function centroid(C::CliffordAlgebra)
  if isdefined(C, :centroid)
    return C.centroid
  end
  n = _dim_qf(C)
  if n == 0
    C.centroid = [one(C)]
    return C.centroid
  end
  T = orthogonal_basis(space(C))
  orth_elt = prod(map(i -> sum(map(j -> gen(C, j) * T[i, j], 1:n)), 1:n))
  orth_elt *= denominator(orth_elt)
  C.disq = coeff(orth_elt^2)[1]
  C.centroid = [one(C), orth_elt]
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
    C.disq = coeff(centroid(C)[2]^2)[1]
  end
end

@doc raw"""
    disq(C::CliffordAlgebra) -> FieldElem

Alias for `quadratic_discriminant`.
"""
disq(C::CliffordAlgebra) = quadratic_discriminant(C)

@doc raw"""
    center(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the center of $C$. It equals `centroid(C)`, if and only if `dim(space(C))`
is odd. Otherwise it is one-dimensional and returned as a tuple of length one
containing the multiplicative identity of $C$.
"""
function center(C::CliffordAlgebra)
  if isdefined(C, :center)
    return C.center
  end
  if is_odd(dim(space(C)))
    C.center = centroid(C)
  else
    C.center = [one(C)]
  end
end

@doc raw"""
    representation_matrix(x::CliffordAlgebraElem) -> MatElem

Return the representation matrix of the element $x$ with respect to `basis(parent(x))`.
"""
function representation_matrix(x::CliffordAlgebraElem, action::Symbol=:left)
  @req (action == :left) || (action == :right) "The action is either $(:left) or $(:right)."
  C = parent(x)
  n = dim(C)
  res = zero_matrix(base_ring(C), n, n)
  if action == :left
    for i in 1:n
      res[:, i] = coeff(x * basis(C, i))
    end
  elseif action == :right
    for i in 1:n
      res[:, i] = coeff(basis(C, i) * x)
    end
  end
  return res
end

################################################################################
#
#  unary operators
#
################################################################################


Base.:-(x::CliffordAlgebraElem) = parent(x)(map(y -> -1 * y, coeff(x)))

################################################################################
#
#  binary operators
#
################################################################################

function Base.:+(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:FieldElem}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  return parent(x)(coeff(x) .+ coeff(y))
end

Base.:-(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:FieldElem} = x + -y

function Base.:*(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:FieldElem}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  xcoeffs, ycoeffs = copy(coeff(x)), copy(coeff(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram_matrix(parent(x)), 1))
end

@doc raw"""
    divexact(x::CliffordAlgebraElem, a::T) where {T<:Union{RingElem, Number}} -> CliffordAlgebraElem

Return the element `y` in the Clifford algebra containing $x$ such that $ay = x$,
if it exists. Otherwise an error is raised.
"""
divexact(x::CliffordAlgebraElem, a::T) where {T<:RingElem} =
  parent(x)(divexact.(coeff(x), a))

divexact(x::CliffordAlgebraElem, a::T) where {T<:Number} = parent(x)(divexact.(coeff(x), a))

################################################################################
#
#  Equality and hash
#
################################################################################

function Base.:(==)(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  return coeff(x) == coeff(y)
end

function Base.hash(x::CliffordAlgebraElem, h::UInt)
  b = 0x1c4629b4de23b24c % UInt
  h = hash(parent(x), h)
  h = hash(coeff(x), h)
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
even_part(x::CliffordAlgebraElem) = parent(x)(even_coeff(x))

@doc raw"""
    odd_part(x::CliffordAlgebraElem) -> CliffordAlgebraElem

Return the projection of $x$ onto the odd Clifford algebra.
"""
odd_part(x::CliffordAlgebraElem) = parent(x)(odd_coeff(x))

################################################################################
#
#  Auxillary functions
#
################################################################################

function _set_even_odd_coeff!(x::CliffordAlgebraElem)
  x.even_coeffs = map(
    y -> if sum(digits(y - 1; base=2, pad=_dim_qf(x.parent))) % 2 == 0
      x.coeffs[y]
    else
      x.parent.base_ring()
    end, 1:(x.parent.dim)
  )
  x.odd_coeffs = x.coeffs - x.even_coeffs
  return x
end

function __can_convert_coeff(coeff::Vector{R}, K::Field) where {R}
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

function _mul_aux(x::Vector{T}, y::Vector{T}, gram::MatElem{T}, i::Int) where {T}
  if length(y) == 1 #|| is_zero(y)
    return x .* y[1]
  end
  return _mul_aux(_mul_with_gen(x, i, gram), y[2:2:end], gram, i + 1) +
         _mul_aux(x, y[1:2:end], gram, i + 1)
end

#Implements the right multiplication with 'i'-th generator e_i of the Clifford algebra containing x.
#The result is returned as a coefficient vector for further computations.
_mul_with_gen(x::Vector{T}, i::Int, gram::MatElem{T}) where {T} = sum(
  map(
    j ->
      x[j] .*
      _mul_baseelt_with_gen(BitVector(digits(j - 1; base=2, pad=ncols(gram))), i, gram),
    1:length(x),
  ),
)

#In the following, bitvectors of length n are used, where 2^n is the Dimension of the Clifford algebra
#containing x. A bitvector will represent the characteristic function of a subset I of {1,...,n}, thus
#representing the element e_I of the canonical basis of the Clifford algebra. Here, e_I is the ordered
#product of the elements of I.

#Converts a bitvector to an integer. The first entry corresponds to 2^0, the second one to 2^1 and so on.
function _bitvec_to_int(x::BitVector)
  pow_of_two = 1
  res = 0
  for i in view(x, 1:length(x))
    res += pow_of_two * i
    pow_of_two <<= 1
  end
  return res
end

#Return the index of e_I in the canonical basis of the Clifford algebra.
_get_basis_index(x::BitVector) = _bitvec_to_int(x) + 1

#Return the index of the last non-zero entry
function _get_last(x::BitVector)
  res = findlast(isone, x)
  if res == nothing
    return 0
  end
  return res
end

function _shift_entries!(X::Vector, s::Int)
  X[(s + 1):end] = X[1:(end - s)]
  X[1:s] = fill(parent(X[1])(), s)
  return X
end

#Right multiplication of the basis element represented by the bitvector 'x' with the
#'i'-th generator of the Clifford algebra containing said basis vector.
function _mul_baseelt_with_gen(x::BitVector, i::Int, gram::MatElem)
  R = base_ring(gram)
  j = _get_last(x)
  res = fill(R(), 2^length(x))
  if j < i
    res[_get_basis_index(x) + 2^(i - 1)] = R(1)
    return res
  end
  if j == i
    res[_get_basis_index(x) - 2^(i - 1)] = R(divexact(gram[i, i], 2))
    return res
  end
  xj_to_zero = copy(x)
  xj_to_zero[j] = 0
  res[_get_basis_index(xj_to_zero)] = gram[i, j]
  res -= _shift_entries!(_mul_baseelt_with_gen(xj_to_zero, i, gram), 2^(j - 1))
end
