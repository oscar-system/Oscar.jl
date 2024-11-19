export
  CliffordOrder,
  ZZCliffordOrder,
  CliffordOrderElem,
  ZZCliffordOrderElem,
  clifford_order,
  ambient_algebra,
  coefficient_ideals,
  even_coeff,
  odd_coeff,
  even_part,
  odd_part,
  center,
  centroid,
  quadratic_discriminant,
  disq

#############################################################
#
#  Datatypes
#
#############################################################

mutable struct CliffordOrder{T, C} <: Hecke.AbstractAssociativeAlgebra{T}

  base_ring::Ring
  ambient_algebra::C
  rank::Int
  lattice::QuadLat
  gram::MatElem
  coefficient_ideals::Vector{NumFieldOrderFractionalIdeal}
  centroid::Any
  disq::Any
  center::Any

  function CliffordOrder{T, C}(ls::QuadLat{S, M}) where {T, C, S<:NumField, M<:MatElem}
    @req is_integral(fractional_ideal(base_ring(ls), base_field(ls)(1//2)) * norm(ls)) "The given lattice is not even!"
    qs = rational_span(ls)
    coeff_ids = _set_coefficient_ideals!(ls)
    return new{T, C}(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs), coeff_ids)
  end

end


mutable struct ZZCliffordOrder <: Hecke.AbstractAssociativeAlgebra{ZZRingElem}

  base_ring::ZZRing
  ambient_algebra::CliffordAlgebra{QQFieldElem, QQMatrix}
  rank::Int
  lattice::ZZLat
  gram::QQMatrix
  centroid::Any
  disq::ZZRingElem
  center::Any

  function ZZCliffordOrder(ls::ZZLat)
    @req is_integral(1//2 * norm(ls)) "The given lattice is not even!"
    qs = rational_span(ls) 
    return new(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs))
  end

end

### Elements ###
mutable struct CliffordOrderElem{T, C, P} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::P
  coeffs::Vector{S} where S<:NumFieldElem
  even_coeffs::Vector{S} where S<:NumFieldElem
  odd_coeffs::Vector{S} where S<:NumFieldElem

  #Returns the 0-element of the Clifford order C
  function CliffordOrderElem{T, C, P}(CO::CliffordOrder{T, C}) where {T, C, P}
    ambalg = CO.ambient_algebra
    newelt = new{T, C, P}(CO, fill(ambalg.base_ring(), CO.rank))
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  CliffordOrderElem(CO::CliffordOrder) = CliffordOrderElem{elem_type(CO.base_ring), typeof(CO.ambient_algebra), typeof(CO)}(CO)

  #Returns the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function CliffordOrderElem{T, C, P}(CO::CliffordOrder{T, C}, coeff::Vector{S}) where {T, C, P, S<:NumFieldElem}
    @req length(coeff) == CO.rank "invalid length of coefficient vector"
    
    for i in 1:CO.rank
      @req coeff[i] in CO.coefficient_ideals[i] "The element does not lie in the Clifford order."
    end

    newelt = new{T, C, CliffordOrder{T, C}}(CO, coeff)
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  function CliffordOrderElem(CO::CliffordOrder{T, C}, coeff::Vector{S}) where {T, C, S}
    K = CO.ambient_algebra.base_ring
    @req _can_convert_coeff(coeff, K) "entries of coefficient vector are not contained in $(K)"
    return CliffordOrderElem{elem_type(CO.base_ring), typeof(CO.ambient_algebra), typeof(CO)}(CO, K.(coeff))
  end

end

### Elements ###
mutable struct ZZCliffordOrderElem <: Hecke.AbstractAssociativeAlgebraElem{ZZRingElem}
  parent::ZZCliffordOrder
  coeffs::Vector{QQFieldElem}
  even_coeffs::Vector{QQFieldElem}
  odd_coeffs::Vector{QQFieldElem}

  #Returns the 0-element of the Clifford order C
  function ZZCliffordOrderElem(CO::ZZCliffordOrder)
    newelt = new(CO, fill(QQ(), CO.rank))
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  #Returns the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function ZZCliffordOrderElem(CO::ZZCliffordOrder, coeff::Vector{QQFieldElem})
    @req length(coeff) == CO.rank "invalid length of coefficient vector"
    for i in 1:CO.rank
      @req is_integer(coeff[i]) "The element does not lie in the Clifford order."
    end
    newelt = new(CO, coeff)
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  function ZZCliffordOrderElem(CO::ZZCliffordOrder, coeff::Vector{S}) where {S}
    @req _can_convert_coeff(coeff, QQ) "entries of coefficient vector are not contained in $(QQField)"
    return ZZCliffordOrderElem(CO, QQ.(coeff))
  end

end

#############################################################
#
#  Type computation
#
#############################################################

elem_type(::Type{CliffordOrder{T, C}}) where {T, C} = CliffordOrderElem{T, C, CliffordOrder{T, C}}

elem_type(::Type{ZZCliffordOrder}) = ZZCliffordOrderElem

parent_type(::Type{CliffordOrderElem{T, C, P}}) where {T, C, P} = P

parent_type(::Type{ZZCliffordOrderElem}) = ZZCliffordOrder

base_ring_type(CO::CliffordOrder{T, C}) where {T, C} = typeof(base_ring(CO))

base_ring_type(CO::ZZCliffordOrder) = ZZRingElem

#############################################################
#
#  Construction
#
#############################################################

### Order ###
@doc raw"""
    clifford_order(ls::QuadLat) -> CliffordOrder

Return the Clifford order of the even quadratic lattice $ls$.
"""

clifford_order(ls::QuadLat) = CliffordOrder{elem_type(base_ring(ls)), typeof(clifford_algebra(rational_span(ls)))}(ls)
@doc raw"""
    clifford_order(ls::ZZLat) -> ZZCliffordOrder

Return the Clifford order of the even quadratic lattice $ls$.
"""
clifford_order(ls::ZZLat) = ZZCliffordOrder(ls)

### Elements ###
(C::CliffordOrder)() = CliffordOrderElem(C)

function (C::CliffordOrder)(a::S) where {S<:RingElem}
  res = fill(zero(a), rank(C))
  res[1] = a
  return CliffordOrderElem(C, res)
end

function (C::CliffordOrder)(a::S) where {S<:Number}
  res = fill(zero(a), rank(C))
  res[1] = a
  return CliffordOrderElem(C, res)
end

(C::CliffordOrder)(coeff::Vector{S}) where {S} = CliffordOrderElem(C, coeff)

(C::ZZCliffordOrder)() = ZZCliffordOrderElem(C)

function (C::ZZCliffordOrder)(a::S) where {S<:RingElem}
  res = fill(zero(a), rank(C))
  res[1] = a
  return ZZCliffordOrderElem(C, res)
end

function (C::ZZCliffordOrder)(a::S) where {S<:Number}
  res = fill(zero(a), rank(C))
  res[1] = a
  return ZZCliffordOrderElem(C, res)
end

(C::ZZCliffordOrder)(coeff::Vector{S}) where {S} = ZZCliffordOrderElem(C, coeff)

################################################################################
#
#  Element conversion to ambient algebra and vice versa
#
################################################################################

function (C::CliffordAlgebra)(elt::CliffordOrderElem)
  @req C === ambient_algebra(parent(elt)) "The Clifford algebra provided and the
    ambient algebra of the Clifford order containing the element provided are not
    identical."
  return C(coeff(elt))
end

function (C::CliffordOrder)(elt::CliffordAlgebraElem)
  @req ambient_algebra(C) === parent(elt) "The ambient algebra of the Clifford order
    provided and the parent object of the element provided are not identical."
    return C(coeff(elt))
end

function (C::CliffordAlgebra)(elt::ZZCliffordOrderElem)
  @req C === ambient_algebra(parent(elt)) "The Clifford algebra provided and the
    ambient algebra of the Clifford order containing the element provided are not
    identical."
  return C(coeff(elt))
end

function (C::ZZCliffordOrder)(elt::CliffordAlgebraElem)
  @req ambient_algebra(C) === parent(elt) "The ambient algebra of the Clifford order
    provided and the parent object of the element provided are not identical."
    return C(coeff(elt))
end

################################################################################
#
#  Containment
#
################################################################################

@doc raw"""
    in(x::CliffordAlgebraElem, C::CliffordOrder)

Return true if the element $x$ is contained in the Clifford order $C$.
"""
function Base.in(x::CliffordAlgebraElem, C::CliffordOrder)
  if ambient_algebra(C) === parent(x)
    return false
  end
  coe, coeids = coeff(x), coefficient_ideals(C)
  for i in 1:rank(C)
    if !(coe[i] in coeids[i])
      return false
    end
  end
  true 
end

@doc raw"""
    in(x::CliffordAlgebraElem, C::ZZCliffordOrder)

Return true if the element $x$ is contained in the Clifford order $C$.
"""
function Base.in(x::CliffordAlgebraElem, C::ZZCliffordOrder)
  if !(ambient_algebra(C) === parent(x))
    return false
  end
  coe = coeff(x)
  for i in 1:rank(C)
    if !(is_integral(coe[i]))
      return false
    end
  end
  true 
end

################################################################################
#
#  Basic field access
#
################################################################################

### Order  ###
@doc raw"""
    base_ring(C::CliffordOrder) -> Ring

Return the base ring of the Clifford order $C$.
"""
base_ring(C::CliffordOrder) = C.base_ring

@doc raw"""
    ambient_algebra(C::CliffordOrder) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of $C$.
"""
ambient_algebra(C::CliffordOrder) = C.ambient_algebra

@doc raw"""
    rank(C::CliffordOrder) -> Int

Return the rank of the Clifford order $C$ over its base ring.
"""
rank(C::CliffordOrder) = C.rank

@doc raw"""
    lattice(C::CliffordOrder) -> QuadLat

Return the underlying even quadratic lattice of the Clifford order $C$.
"""
lattice(C::CliffordOrder) = C.lattice

@doc raw"""
    gram_matrix(C::CliffordOrder) -> MatElem

Return the Gram matrix with respect to the fixed pseudo-basis of the Clifford order $C$.
"""
gram_matrix(C::CliffordOrder) = C.gram

@doc raw"""
    coefficient_ideals(C::CliffordOrder) -> Vector{NumFieldOrderFractionalIdeal}

Return the vector of coefficient ideals of the canonical pseudo-basis of $C$.
"""
coefficient_ideals(C::CliffordOrder) = C.coefficient_ideals


@doc raw"""
    base_ring(C::ZZCliffordOrder) -> ZZRing

Return the base ring of the Clifford order $C$.
"""
base_ring(C::ZZCliffordOrder) = C.base_ring

@doc raw"""
    ambient_algebra(C::ZZCliffordOrder) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of $C$.
"""
ambient_algebra(C::ZZCliffordOrder) = C.ambient_algebra

@doc raw"""
    rank(C::ZZCliffordOrder) -> Int

Return the rank of the Clifford order $C$ over its base ring.
"""
rank(C::ZZCliffordOrder) = C.rank

@doc raw"""
    lattice(C::ZZCliffordOrder) -> ZZLat

Return the underlying even quadratic lattice of the Clifford order $C$.
"""
lattice(C::ZZCliffordOrder) = C.lattice

@doc raw"""
    gram_matrix(C::ZZCliffordOrder) -> QQMatrix

Return the Gram matrix with respect to the fixed basis of the Clifford order $C$.
"""
gram_matrix(C::ZZCliffordOrder) = C.gram



### Elements ###
@doc raw"""
    parent(x::CliffordOrderElem) -> CliffordOrder

Return the Clifford order containing $x$.
"""
parent(x::CliffordOrderElem) = x.parent

@doc raw"""
  coeff(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ with respect to the
canonical pseudo-basis of its parent Clifford order.
"""
coeff(x::CliffordOrderElem) = x.coeffs

@doc raw"""
    even_coeff(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coeff(x::CliffordOrderElem)
  if isdefined(x, :even_coeffs)
    return x.even_coeffs
  end
  _set_even_odd_coeff!(x)
  return x.even_coeffs
end

@doc raw"""
    odd_coeff(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coeff(x::CliffordOrderElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs
  end
  _set_even_odd_coeff!(x)
  return x.odd_coeffs
end

@doc raw"""
    parent(x::ZZCliffordOrderElem) -> ZZCliffordOrder

Return the Clifford order containing $x$.
"""
parent(x::ZZCliffordOrderElem) = x.parent

@doc raw"""
  coeff(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of $x$ with respect to the
canonical basis of its parent Clifford order.
"""
coeff(x::ZZCliffordOrderElem) = x.coeffs

@doc raw"""
    even_coeff(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of $x$ with respect to the
canonical basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coeff(x::ZZCliffordOrderElem)
  if isdefined(x, :even_coeffs)
    return x.even_coeffs
  end
  _set_even_odd_coeff!(x)
  return x.even_coeffs
end

@doc raw"""
    odd_coeff(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of $x$ with respect to the
canonical basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coeff(x::ZZCliffordOrderElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs
  end
  _set_even_odd_coeff!(x)
  return x.odd_coeffs
end

#############################################################
#
#  String I/O
#
#############################################################

### Order ###
function Base.show(io::IO, ::MIME"text/plain", C::CliffordOrder)
  io = pretty(io)
  print(io, "Clifford order of even lattice over $(base_ring(C)) with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent(), "\nand coefficient ideals\n")
  print(io, Indent())
  show(io, "text/plain", _coefficient_ideals_of_lattice(lattice(C)))
  print(io, Dedent())
end

Base.show(io::IO, C::CliffordOrder) = print(io, "Clifford order over $(base_ring(C))")

function Base.show(io::IO, ::MIME"text/plain", C::ZZCliffordOrder)
  io = pretty(io)
  print(io, "Clifford order of even integer lattice with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent())
end

Base.show(io::IO, C::Union{ZZCliffordOrder, CliffordOrder}) = print(io, "Clifford order over $(base_ring(C))")

### Elements ###
function Base.show(io::IO, x::Union{ZZCliffordOrderElem, CliffordOrderElem})
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
    zero(C::CliffordOrder) -> CliffordOrderElem

Return the additive identity of the Clifford order $C$.
"""
zero(C::CliffordOrder) = C()

@doc raw"""
    one(C::CliffordOrder) -> CliffordOrderElem

Return the multiplicative identity of the Clifford order $C$.
"""
function one(C::CliffordOrder)
  v = fill(base_ring(ambient_algebra(C))(), rank(C))
  v[1] = base_ring(ambient_algebra(C))(1)
  return C(v)
end

@doc raw"""
    pseudo_basis(C::CliffordOrder, i::Int) -> PMat{NumFieldElem, NumFieldOrderFractionalIdeal}

Return the $i$-th canonical pseudo-basis element of $C$.
"""
function pseudo_basis(C::CliffordOrder, i::Int)
  ba = base_ring(ambient_algebra(C))
  res = zero_matrix(ba, 1, rank(C))
  res[1, i] = ba(1)
  return pseudo_matrix(res, [coefficient_ideals(C)[i]])
end

@doc raw"""
    gen(C::CliffordOrder, i::Int) -> PMat{NumFieldElem, NumFieldOrderFractionalIdeal}

Return the $i$-th pseudo-element of the canonical multiplicative pseudo-generating set
of the Clifford order $C$.
"""
function gen(C::CliffordOrder, i::Int)
  ba = base_ring(ambient_algebra(C))
  res = zero_matrix(ba, 1, rank(C))
  if i<= 0
    res[1, i] #Throw a BoundsError for consistency
  end
  res[1, 2^(i - 1) + 1] = ba(1)
  return pseudo_matrix(res, [coefficient_ideals(lattice(C))[i]])
end

@doc raw"""
    pseudo_basis(C::CliffordOrder) -> PMat{NumFieldElem, NumFieldOrderFractionalIdeal}

Return the canonical pseudo-basis of the Clifford order $C$.
"""
pseudo_basis(C::CliffordOrder) = pseudo_matrix(identity_matrix(base_ring(ambient_algebra(C)), rank(C)), coefficient_ideals(C))

@doc raw"""
    gens(C::CliffordOrder) -> PMat{NumFieldElem, NumFieldOrderFractionalIdeal}

Return the canonical multiplicative pseudo-generating set of the Clifford order $C$.
"""
function gens(C::CliffordOrder) 
  ba = base_ring(ambient_algebra(C))
  n = rank(lattice(C))
  res = zero_matrix(ba, n, rank(C))
  for i in 1:n
    res[i, 2^(i - 1) + 1] = ba(1)
  end
  return pseudo_matrix(res, coefficient_ideals(lattice(C)))
end

@doc raw"""
    is_commutative(C::CliffordOrder) -> Bool

Return `true` if $C$ is commutative and `false` otherwise.
"""
is_commutative(C::CliffordOrder) = rank(C) == 1 || rank(C) == 2

#######

@doc raw"""
    zero(C::ZZCliffordOrder) -> ZZCliffordOrderElem

Return the additive identity of the Clifford order $C$.
"""
zero(C::ZZCliffordOrder) = C()

@doc raw"""
    one(C::ZZCliffordOrder) -> ZZCliffordOrderElem

Return the multiplicative identity of the Clifford order $C$.
"""
function one(C::ZZCliffordOrder)
  coeff_res = fill(QQ(), rank(C))
  coeff_res[1] = QQ(1)
  return C(coeff_res)
end

@doc raw"""
    basis(C::ZZCliffordOrder, i::Int) -> ZZCliffordOrderElem

Return the $i$-th canonical basis element of $C$.
"""
function basis(C::ZZCliffordOrder, i::Int)
  coeff_res = fill(QQ(), rank(C))
  coeff_res[i] = QQ(1)
  return C(coeff_res)
end

@doc raw"""
    gen(C::CliffordOrder, i::Int) -> ZZCliffordOrderElem

Return the $i$-th element of the canonical multiplicative generating set
of the Clifford order $C$.
"""
function gen(C::ZZCliffordOrder, i::Int)
  coeff_res = fill(QQ(), rank(C))
  if i <= 0
    coeff_res[i] #Throw a BoundsError for consistency
  end
  coeff_res[2^(i - 1) + 1] = QQ(1)
  return C(coeff_res)
end

@doc raw"""
    basis(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return the canonical basis of the Clifford order $C$.
"""
basis(C::ZZCliffordOrder) = map(i -> basis(C, i), 1:rank(C))

@doc raw"""
    gens(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return the canonical multiplicative generating set of the Clifford order $C$.
"""
gens(C::ZZCliffordOrder) = map(i -> gen(C, i), 1:rank(lattice(C)))

@doc raw"""
    is_commutative(C::ZZCliffordOrder) -> Bool

Return `true` if $C$ is commutative and `false` otherwise.
"""
is_commutative(C::ZZCliffordOrder) = rank(C) == 1 || rank(C) == 2

################################################################################
#
#  Other functionality
#
################################################################################

@doc raw"""
    centroid(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return the centroid of $C$. Unless `rank(lattice(C)) = 0`, it is always free of rank
two, so it is returned as a vector containing the basis elements. The first one is the
multiplicative identity of $C$. The square of the second basis element, if present,
equals `quadratic_discriminant(C)`.
"""
function centroid(C::ZZCliffordOrder)
  if isdefined(C, :centroid)
    return C.centroid
  end
  n = rank(lattice(C))
  if n == 0
    C.centroid = [one(C)]
    return C.centroid
  end
  
  T = orthogonal_basis(space(ambient_algebra(C)))
  for i in 1:n
    T[i,:] *= 1//gcd(T[i,:])
  end

  orth_elt = prod(map(i -> sum(map(j -> gen(C, j) * T[i, j], 1:n)), 1:n))
  orth_elt *= 1//gcd(coeff(orth_elt))
  C.disq = ZZ(coeff(orth_elt^2)[1])
  C.centroid = [one(C), orth_elt]
end


@doc raw"""
    quadratic_discriminant(C::ZZCliffordOrder) -> ZZRingElem

Returns the quadratic discriminant of $C$ as an integer.
"""
function quadratic_discriminant(C::ZZCliffordOrder) 
  if isdefined(C, :disq)
    return C.disq
  end
  if rank(lattice(C)) == 0
    C.disq = one(base_ring(C))
  else
    C.disq = ZZ(coeff(centroid(C)[2]^2)[1])
  end
end

@doc raw"""
    disq(C::ZZCliffordOrder) -> ZZRingElem

Alias for `quadratic_discriminant`.
"""
disq(C::ZZCliffordOrder) = quadratic_discriminant(C)

@doc raw"""
    center(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Returns the center of $C$. It equals `centroid(C)`, if and only if `rank(lattice(C))`
is odd. Otherwise it free of rank one, generated by the multiplicative identity of $C$.
"""
function center(C::ZZCliffordOrder)
  if isdefined(C, :center)
    return C.center
  end
  if is_odd(rank(lattice(C)))
    C.center = centroid(C)
  else 
    C.center = [one(C)]
  end
end

###

@doc raw"""
    centroid(C::CliffordOrder) -> PMat{NumFieldElem, NumFieldOrderFractionalIdeal}

Return the centroid of $C$. Unless `rank(lattice(C)) = 0`, it is a lattice of rank two,
so it is returned as a pseudo-matrix with two rows. The first row is the coefficient
vector of the multiplicative identity of $C$ with coefficients in `base_ring(C)`.
"""
function centroid(C::CliffordOrder)
  if isdefined(C, :centroid)
    return C.centroid
  end
  _set_centroid_and_disq!(C)
  return C.centroid
end

@doc raw"""
    quadratic_discriminant(C::CliffordOrder) -> Tuple{NumFieldElem, NumFieldOrderFractionalIdeal}

Returns the quadratic discriminant of $C$.
"""
function quadratic_discriminant(C::CliffordOrder) 
  if isdefined(C, :disq)
    return C.disq
  end
  _set_centroid_and_disq!(C)
  return C.disq
end

@doc raw"""
    disq(C::CliffordOrder) -> Tuple{NumFieldElem, NumFieldOrderFractionalIdeal}

Alias for `quadratic_discriminant`.
"""
disq(C::CliffordOrder) = quadratic_discriminant(C)

@doc raw"""
    center(C::CliffordOrder) -> PMat{NumFieldElem, NumFieldOrderFractionalIdeal}

Returns the center of $C$. It equals `centroid(C)`, if and only if `rank(lattice(C))`
is odd. Otherwise it free of rank one, generated by the multiplicative identity of $C$.
"""
function center(C::CliffordOrder)
  if isdefined(C, :center)
    return C.center
  end
  if is_odd(rank(lattice(C)))
    C.center = centroid(C)
  else 
    C.center = pseudo_basis(C, 1)
  end
end

################################################################################
#
#  unary operators
#
################################################################################

Base.:+(x::CliffordOrderElem) = x

Base.:-(x::CliffordOrderElem) = parent(x)(map(y -> -1 * y, coeff(x)))

Base.:+(x::ZZCliffordOrderElem) = x

Base.:-(x::ZZCliffordOrderElem) = parent(x)(map(y -> -1 * y, coeff(x)))

################################################################################
#
#  binary operators
#
################################################################################

function Base.:+(x::CliffordOrderElem{T, CliffordAlgebra{U, V}},
                  y::CliffordOrderElem{T, CliffordAlgebra{U, V}}) where {T, U<:NumFieldElem, V}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  return parent(x)(coeff(x) .+ coeff(y))
end

Base.:-(x::CliffordOrderElem{T, CliffordAlgebra{U, V}},
          y::CliffordOrderElem{T, CliffordAlgebra{U, V}}) where {T, U<:NumFieldElem, V} = x + -y

function Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U, V}},
                  y::CliffordOrderElem{T, CliffordAlgebra{U, V}}) where {T, U<:NumFieldElem, V}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  xcoeffs, ycoeffs = copy(coeff(x)), copy(coeff(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram_matrix(parent(x)), 1))
end

Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U,V}}, a::W) where {T, U<:NumFieldElem, V, W<:FieldElem} = parent(x)(a * ambient_algebra(parent(x))(x))

Base.:*(a::W, x::CliffordOrderElem{T, CliffordAlgebra{U,V}}) where {T, U<:NumFieldElem, V, W<:FieldElem} = parent(x)(a * ambient_algebra(parent(x))(x))

Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U,V}}, a::Rational{Int}) where {T, U<:NumFieldElem, V} = parent(x)(a .* coeff(x))

Base.:*(a::Rational{Int}, x::CliffordOrderElem{T, CliffordAlgebra{U,V}}) where {T, U<:NumFieldElem, V} = parent(x)(a .* coeff(x))


function Base.:+(x::ZZCliffordOrderElem, y::ZZCliffordOrderElem)
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  return parent(x)(coeff(x) .+ coeff(y))
end

Base.:-(x::ZZCliffordOrderElem, y::ZZCliffordOrderElem) = x + -y

function Base.:*(x::ZZCliffordOrderElem, y::ZZCliffordOrderElem)
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford algebra"
  xcoeffs, ycoeffs = copy(coeff(x)), copy(coeff(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram_matrix(parent(x)), 1))
end

Base.:*(x::ZZCliffordOrderElem, a::QQFieldElem) = parent(x)(a .* coeff(x))

Base.:*(a::QQFieldElem, x::ZZCliffordOrderElem) = parent(x)(a .* coeff(x))

Base.:*(x::ZZCliffordOrderElem, a::Rational{Int}) = parent(x)(a .* coeff(x))

Base.:*(a::Rational{Int}, x::ZZCliffordOrderElem) = parent(x)(a .* coeff(x))

#=
@doc raw"""
    divexact(x::CliffordAlgebraElem, a::RingElem) -> CliffordAlgebraElem

Return the element 'y' in the given Clifford algebra such that $ay = x$,
if it exists. Otherwise an error is raised.
"""
divexact(x::CliffordAlgebraElem{T}, elt::T) where {T<:FieldElem} =
  parent(x)(divexact.(coeff(x), elt))

divexact(x::CliffordAlgebraElem{T}, elt::QQFieldElem) where {T<:NumFieldElem} =
  parent(x)(divexact.(coeff(x), elt))

divexact(x::CliffordAlgebraElem{T}, elt::Rational{Int}) where {T<:FieldElem} =
  parent(x)(divexact.(coeff(x), elt))

divexact(x::CliffordAlgebraElem{T}, elt::Int) where {T<:FieldElem} =
  divexact(x, base_ring(parent(x))(elt))

divexact(x::CliffordAlgebraElem{T}, elt::ZZRingElem) where {T<:FieldElem} =
  divexact(x, base_ring(parent(x))(elt))
=# 
################################################################################
#
#  Equality and hash
#
################################################################################

function Base.:(==)(x::CliffordOrderElem{T}, y::CliffordOrderElem{T}) where {T}
  @req parent(x) === parent(y) "The inputs must lie in the same Clifford order."
  return coeff(x) == coeff(y)
end

function Base.hash(x::CliffordOrderElem, h::UInt)
  b = 0x8f3a7c2b1d4e5f6a % UInt
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
    even_part(x::CliffordOrderElem) -> CliffordOrderElem

Return the of the projection of $x$ onto the even Clifford order.
"""
even_part(x::CliffordOrderElem) = parent(x)(even_coeff(x))

@doc raw"""
    odd_part(x::CliffordOrderElem) -> CliffordOrderElem

Return the projection of $x$ onto the odd Clifford order.
"""
odd_part(x::CliffordOrderElem) = parent(x)(odd_coeff(x))

@doc raw"""
    even_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the of the projection of $x$ onto the even Clifford order.
"""
even_part(x::ZZCliffordOrderElem) = parent(x)(even_coeff(x))

@doc raw"""
    odd_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the projection of $x$ onto the odd Clifford order.
"""
odd_part(x::ZZCliffordOrderElem) = parent(x)(odd_coeff(x))

#############################################################
#
#  Auxillary functions
#
#############################################################

# Computes the ideals present in the canonical pseudo-basis of the Clifford order
# of the lattice 'ls' in its fixed pseudo-basis.
function _set_coefficient_ideals!(ls::QuadLat)
  coeff_ids = coefficient_ideals(ls)
  res = [fractional_ideal(base_ring(ls), base_ring(ls)(1))]
  for i in 1:rank(ls)
    tmpres = map(x -> x*coeff_ids[i], res)
    res = vcat(res, tmpres)
  end
  return res
end

function _set_even_odd_coeff!(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
  x.even_coeffs = map(
                      y -> if sum(digits(y - 1; base=2, pad=rank(x.parent.lattice))) % 2 == 0
      x.coeffs[y]
    else
      x.parent.ambient_algebra.base_ring()
    end, 1:(x.parent.rank)
  )
  x.odd_coeffs = x.coeffs - x.even_coeffs
  return x
end

function _can_convert_coeff(coeff::Vector{S}, K::Field) where {S}
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

#Returns the coefficient_ideals of the given lattice 'ls'
_coefficient_ideals_of_lattice(ls::QuadLat) = coefficient_ideals(ls)::Vector{<:NumFieldOrderFractionalIdeal}

function _set_centroid_and_disq!(C::CliffordOrder)
  n = rank(lattice(C))
  if n == 0
    pb1 = pseudo_basis(C,1)
    C.centroid = pb1
    C.disq = (coeff(pb1[1]), coefficient_ideals(pb1)[1])
  end
  
  br = base_ring(C)
  z_elt = coeff(centroid(ambient_algebra(C))[2])
  lambda_empt = z_elt[1]
  
  #compute ideal intersection
  ideal_array = map(i -> !is_zero(z_elt[i]) ? (z_elt[i])^(-1) * coefficient_ideals(C)[i] : fractional_ideal(br, zero(br)), 2:length(z_elt)) 
  filter!(!is_zero, ideal_array)
  len = length(ideal_array)
  c_ideal = ideal_array[1]
  for i in 2:len
    c_ideal = lcm(c_ideal, ideal_array[i])
  end
  
  if is_zero(lambda_empt)
    C.disq = tuple(disq(ambient_algebra(C)) * c_ideal^2, disq(ambient_algebra(C)))
    C.centroid = pseudo_matrix(matrix(base_ring(ambient_algebra(C)), 2, 2^n, vcat(coeff(one(C)), z_elt)), [fractional_ideal(br, one(br)), c_ideal])
  else
    z_elt = (2 * lambda_empt)^(-1) .* z_elt
    c_ideal *= (2 * lambda_empt)
    b_ideal = lcm(fractional_ideal(br, br(2)), c_ideal)
    C.disq = tuple(disq(ambient_algebra(C)) * b_ideal^2, disq(ambient_algebra(C)))
    z_elt[1] += base_ring(ambient_algebra(C))(1//2)
    b_ideal = lcm(fractional_ideal(br, one(br)), c_ideal)
    C.centroid = pseudo_matrix(matrix(base_ring(ambient_algebra(C)), 2, 2^n, vcat(coeff(one(C)), z_elt)), [fractional_ideal(br, one(br)), b_ideal])
  end
end
