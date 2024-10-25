export CliffordAlgebra,
  CliffordAlgebraElem,
  clifford_algebra,
  dim_qf,
  even_coeff,
  odd_coeff,
  even_part,
  odd_part,
  coeff,
  centroid,
  center,
  change_basis,
  change_basis!,
  max_orth_elt

################################################################################
#
#  Datatypes
#
################################################################################

### Algebra ###
mutable struct CliffordAlgebra{T} <: Hecke.AbstractAssociativeAlgebra{T}
  base_ring::Ring
  gram::MatElem{T}
  dim_qf::Int
  dim::Int
  centroid::Any
  max_orth_elt::Any

  #Returns the Clifford algebra with the quadratic form defined by the Gram matrix 'gram' 
  function CliffordAlgebra{T}(gram::MatElem{T}) where {T}
    @req T <: RingElem "The base ring is not commutative"
    @req is_symmetric(gram) "The Gram matrix must be symmetric"
    @req all(map(x -> is_divisible_by(x, base_ring(gram)(2)), diagonal(gram))) "The Gram matrix must define a quadratic form with values in $(base_ring(gram)), that is, all its diagonal entries must be divisible by 2."
    return new{T}(base_ring(gram), gram, ncols(gram), 2^ncols(gram))
  end
end

### Elements ###
mutable struct CliffordAlgebraElem{T,S} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::S
  coeffs::Vector{T}
  even_coeffs::Vector{T}
  odd_coeffs::Vector{T}

  #Returns the 0-element of the Clifford algebra C
  function CliffordAlgebraElem{T,S}(C::CliffordAlgebra{T}) where {T,S}
    z = new{T,CliffordAlgebra{T}}(C, fill(C.base_ring(), C.dim))
    set_even_odd_coeff!(z)
    return z
  end

  CliffordAlgebraElem(C::CliffordAlgebra) = CliffordAlgebraElem{elem_type(C.base_ring),typeof(C)}(C)

  #Returns the element in the Clifford algebra C with coefficient vector coeff wrt. the canonical basis
  function CliffordAlgebraElem{T,S}(C::CliffordAlgebra{T}, coeff::Vector{T}) where {T,S}
    @req length(coeff) == C.dim "invalid length of coefficient vector"
    z = new{T,CliffordAlgebra{T}}(C, coeff)
    set_even_odd_coeff!(z)
    return z
  end

  CliffordAlgebraElem(C::CliffordAlgebra{T}, coeff::Vector{T}) where {T} =
    CliffordAlgebraElem{elem_type(base_ring(C)),typeof(C)}(C, coeff)
end

elem_type(::Type{CliffordAlgebra{T}}) where {T} = CliffordAlgebraElem{T,CliffordAlgebra{T}}

parent_type(::Type{CliffordAlgebraElem{T,S}}) where {T,S} = S

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
    gram(C::CliffordAlgebra) -> MatElem

Return the Gram matrix of the free quadratic module with Clifford algebra $C$.
"""
gram(C::CliffordAlgebra) = C.gram

@doc raw"""
    dim(C::CliffordAlgebra) -> Int

Return the dimension of the Clifford algebra $C$ over its base ring.
"""
dim(C::CliffordAlgebra) = C.dim

@doc raw"""
    dim_qf(C::CliffordAlgebra) -> Int

Return the dimension of the free quadratic module with Clifford algebra $C$.
"""
dim_qf(C::CliffordAlgebra) = C.dim_qf

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

function set_even_odd_coeff!(x::CliffordAlgebraElem)
  x.even_coeffs = map(
    y -> if sum(digits(y - 1; base=2, pad=x.parent.dim_qf)) % 2 == 0
      x.coeffs[y]
    else
      x.parent.base_ring()
    end, 1:(x.parent.dim)
  )
  x.odd_coeffs = x.coeffs - x.even_coeffs
  return x
end

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
  set_even_odd_coeff!(x)
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
  set_even_odd_coeff!(x)
  return x.odd_coeffs
end

################################################################################
#
#  Construction
#
################################################################################

### Algebra ###

@doc raw"""
    clifford_algebra(G::MatElem) -> CliffordAlgebra

Return the Clifford algebra of the quadratic $R$-module $(R^n,q)$ where
$R$ is the base ring of the Gram matrix $G$, $n$ is the number of rows
or columns of G and $q$ is the quadratic form, defined by $G$.
"""
clifford_algebra(gram::MatElem) = CliffordAlgebra{elem_type(base_ring(gram))}(gram)

### Elements ###
(C::CliffordAlgebra)() = CliffordAlgebraElem(C)

#The implemention for AbsAlgElem should already work.
#(C::CliffordAlgebra{T})(a::T) where {T<:RingElem} = CliffordAlgebraElem(C, a .* coeff(one(C)))

(C::CliffordAlgebra{T})(a::Int) where {T<:RingElem} = C(base_ring(C)(a))

(C::CliffordAlgebra{T})(v::Vector{T}) where {T<:RingElem} = CliffordAlgebraElem(C, v)

(C::CliffordAlgebra)(v::Vector) = CliffordAlgebraElem(C, base_ring(C).(v))

################################################################################
#
#  String I/O
#
################################################################################

### Algebra ###
function Base.show(io::IO, ::MIME"text/plain", C::CliffordAlgebra)
  io = pretty(io)
  print(io, "Clifford algebra of quadratic form with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram(C))
  print(io, Dedent(), "\ndefined over $(base_ring(C))")
end

Base.show(io::IO, C::CliffordAlgebra) = print(io, "Clifford algebra over $(base_ring(C))")

### Elements ###
function Base.show(io::IO, x::CliffordAlgebraElem)
  print(io, "(")
  foreach(y -> print(io, "$y "), coeff(x)[1:(end - 1)])
  print(io, "$(coeff(x)[end]))")
end

################################################################################
#
#  basic functionality
#
################################################################################

zero(C::CliffordAlgebra) = C()

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
  set_even_odd_coeff!(res)
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
    res.coeffs[i]
  end
  res.coeffs[2^(i - 1) + 1] = base_ring(C)(1)
  set_even_odd_coeff!(res)
  return res
end

basis(C::CliffordAlgebra) = map(i -> basis(C, i), 1:dim(C))

gens(C::CliffordAlgebra) = map(i -> gen(C, i), 1:dim_qf(C))

is_commutative(C::CliffordAlgebra) = dim(C) == 1

################################################################################
#
#  other functionality
#
################################################################################

@doc raw"""
    centroid(C::CliffordAlgebra) -> ResidueRing, Map

Return the centroid of $C$, together with an embedding into $C$. The centroid is
the centraliser of the even Clifford algebra and always a two-dimensional étale
algebra over the base ring.
"""
function centroid(C::CliffordAlgebra{T}) where {T<:FieldElem}
  
  return 0
end

function centroid(C::CliffordAlgebra{T}) where {T<:RingElem}

  return 0
end

function centroid(C::CliffordAlgebra{ZZRingElem})
  T = orthogonal_basis(quadratic_space(QQ, gram(C)))
  n = ncols(T)
  for i in 2:n
    T[i,:] *= denominator(T[i,:])
  end
  orth = prod(map(i -> sum(map(j -> gen(C, j) * ZZ(T[i, j]), 1:n)), 1:n))
  orth = divexact(orth, gcd(coeff(orth)))

end

function max_orth_elt(C::CliffordAlgebra{ZZRingElem})
  T = orthogonal_basis(quadratic_space(QQ, gram(C)))
  n = ncols(T)
  for i in 2:n
    T[i,:] *= denominator(T[i,:])
  end
  orth = prod(map(i -> sum(map(j -> gen(C, j) * ZZ(T[i, j]), 1:n)), 1:n))
  return divexact(orth, gcd(coeff(orth)))
end

@doc raw"""
    center(C::CliffordAlgebra) -> ResidueRing, Map

Return the center of $C$, together with an embedding into $C$. The center is
equal to the base ring of $C$, if $dim_qf(C)$ is even and equal to the centroid
otherwise.
"""
function center(C::CliffordAlgebra{T}) where {T<:FieldElem}
  if dim(C) % 2 == 1
    if isdefined(C, :centroid)
      return C.centroid
    end
    return centroid(C)
  end
  R, x = polynomial_ring(base_ring(C), 'x')
  qR, canmap = quo(R, x - 1)
  invcanmap = canmap.section
  return qR, map_from_func(p -> (invcanmap(p))(1), qR, C)
end

function center(C::CliffordAlgebra{T}) where {T<:NumFieldOrderElem}
  if dim(C) % 2 == 1
    if isdefined(C, :centroid)
      return C.centroid
    end
    return centroid(C)
  end
  R, x = polynomial_ring(base_ring(C), 'x')
  qR, canmap = quo(R, x - 1)
  return qR, map_from_func(p -> (canmap.section(p))(1), qR, C)
end

function center(C::CliffordAlgebra{ZZRingElem})
  if dim(C) % 2 == 1
    if isdefined(C, :centroid)
      return C.centroid
    end
    return centroid(C)
  end
  R, x = polynomial_ring(ZZRing(), 'x')
  qR, canmap = quo(R, x - 1)
  return qR, map_from_func(p -> C((canmap.section(p))(1)), qR, C)
end

@doc raw"""
    change_basis(C::CliffordAlgebra{T}, X::MatElem{T}) where {T :< RingElem} -> CliffordAlgebra

Return the Clifford algebra of the underlying quadratic module that
is obtained after performing the base change given by $X$, viewing
its rows as the coefficient vectors of the new basis elements.
"""
function change_basis(C::CliffordAlgebra{T}, X::MatElem{T}) where {T<:RingElem}
  @req is_invertible(X) "The transformation matrix must be invertible."
  return clifford_algebra(X * gram(C) * transpose(X))
end

@doc raw"""
    change_basis!(C::CliffordAlgebra{T}, X::MatElem{T}) where {T <: RingElem} -> CliffordAlgebra

Return the Clifford algebra of the underlying quadratic module that
is obtained after performing the base change given by $X$ in-place,
viewing its rows as the coefficient vectors of the new basis elements. 
"""
function change_basis!(C::CliffordAlgebra{T}, X::MatElem{T}) where {T<:RingElem}
  @req is_invertible(X) "The transformation matrix must be invertible."
  C.gram = X * gram(C) * transpose(X)
  return C
end

################################################################################
#
#  unary operators
#
################################################################################

Base.:+(x::CliffordAlgebraElem) = x

Base.:-(x::CliffordAlgebraElem) = parent(x)(map(y -> -1 * y, coeff(x)))

################################################################################
#
#  binary operators
#
################################################################################

function Base.:+(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:RingElem}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  return parent(x)(coeff(x) .+ coeff(y))
end

Base.:-(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:RingElem} = x + -y

function Base.:*(x::CliffordAlgebraElem{T}, y::CliffordAlgebraElem{T}) where {T<:RingElem}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  xcoeffs, ycoeffs = copy(coeff(x)), copy(coeff(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram(parent(x)), 1))
end

#The implementation from AbsAlgElem should work already in this case!
#Base.:*(x::CliffordAlgebraElem{T}, elt::T) where {T <: RingElem} = parent(x)(elt .* coeff(x))

Base.:*(x::CliffordAlgebraElem{T}, elt::Union{Int,Rational{Int}}) where {T<:RingElem} =
  base_ring(parent(x))(elt) * x

Base.:*(elt::Union{Int,Rational{Int}}, x::CliffordAlgebraElem{T}) where {T<:RingElem} = x * elt

@doc raw"""
    divexact(x::CliffordAlgebraElem{T}, a::T) where {T<:RingElem} -> CliffordAlgebraElem{T}

Return the element 'y' in the given Clifford algebra such that $ay = x$,
if it exists. Otherwise an error is raised.
"""
divexact(x::CliffordAlgebraElem{T}, elt::T) where {T<:RingElem} =
  parent(x)(divexact.(coeff(x), elt))

divexact(x::CliffordAlgebraElem{T}, elt::Int) where {T<:RingElem} =
  divexact(x, base_ring(parent(x))(elt))

################################################################################
#
#  equality and hash
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

Return the of the projection of $x$ onto the even Clifford algebra
"""
even_part(x::CliffordAlgebraElem) = parent(x)(even_coeff(x))

@doc raw"""
    odd_part(x::CliffordAlgebraElem) -> CliffordAlgebraElem

Return the projection of $x$ onto the odd Clifford algebra.
"""
odd_part(x::CliffordAlgebraElem) = parent(x)(odd_coeff(x))

################################################################################
#
#  Auxillary functions for multiplication
#
################################################################################

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

#Returns the index of e_I in the canonical basis of the Clifford algebra.
_get_basis_index(x::BitVector) = _bitvec_to_int(x) + 1

#=
#Returns the index of the first non-zero entry
function _get_first(x::BitVector)
  res = findfirst(isone, x)
  if res == nothing
    return length(x) + 1
  end
  return res
end
=#

#Returns the index of the last non-zero entry
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
