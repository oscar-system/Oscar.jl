export
  CliffordOrder,
  ZZCliffordOrder,
  CliffordOrderElem,
  ZZCliffordOrderElem,
  clifford_order,
  algebra,
  coefficient_ideals,
  even_coefficients,
  odd_coefficients,
  even_part,
  odd_part,
  basis_of_center,
  pseudo_basis_of_center,
  basis_of_centroid,
  pseudo_basis_of_centroid,
  pseudo_gen,
  pseudo_gens,
  quadratic_discriminant,
  disq

#############################################################
#
#  Datatypes
#
#############################################################

# Data structure for Clifford orders over rings distinct from the integers. The type variable 'T' represents the element type
# of the base ring, i.e. it is usually AbsSimpleNumFieldOrderElem. The type variable 'C' represents the type of its ambient algebra,
# i.e. it is usually CliffordAlgebra{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}.
mutable struct CliffordOrder{T, C} <: Hecke.AbstractAssociativeAlgebra{T}

  base_ring::Ring
  algebra::C
  rank::Int
  lattice::QuadLat
  gram::MatElem

  # In the 4 lines below let CO be an instance of CliffordOrder and R = base_ring(CO), the base ring of CO
  coefficient_ideals::Any # Vector{Hecke.fractional_ideal_type(base_ring_type(C))}
  pseudo_basis_of_centroid::Any # Vector{Tuple{elem_type(CO), Hecke.fractional_ideal_type(base_ring_type(CO))}}
  disq::Any # Tuple{Hecke.fractional_ideal_type(base_ring_type(CO)), elem_type(base_ring(algebra(CO)))}
  pseudo_basis_of_center::Any # Vector{Tuple{elem_type(CO), Hecke.fractional_ideal_type(base_ring_type(CO))}}

  function CliffordOrder{T, C}(ls::QuadLat{S, M}) where {T, C, S<:NumField, M<:MatElem}
    if !is_zero(rank(ls))
      @req is_integral(fractional_ideal(base_ring(ls), base_field(ls)(1//2)) * norm(ls)) "The given lattice is not even!"
    end
    qs = rational_span(ls)
    coeff_ids = _set_coefficient_ideals!(ls)
    return new{T, C}(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs), coeff_ids)
  end
end

mutable struct ZZCliffordOrder <: Hecke.AbstractAssociativeAlgebra{ZZRingElem}

  base_ring::ZZRing
  algebra::CliffordAlgebra{QQFieldElem, QQMatrix}
  rank::Int
  lattice::ZZLat
  gram::QQMatrix
  basis_of_centroid::Any #Always of type Vector{ZZCliffordOrderElem} 
  disq::ZZRingElem
  basis_of_center::Any #Always of type Vector{ZZCliffordOrderElem}

  function ZZCliffordOrder(ls::ZZLat)
    if !is_zero(rank(ls))
      @req is_even(ls) "The given lattice is not even!"
    end
    qs = rational_span(ls) 
    return new(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs))
  end

end

##### Elements #####
# Data structure for elements of Clifford orders over rings distinct from the integers. The type
# variables serve the same purpose as they do for Clifford orders.
mutable struct CliffordOrderElem{T, C} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::CliffordOrder{T, C}
  coeffs::Any 
  even_coeffs::Any
  odd_coeffs::Any

  #Return the 0-element of the Clifford order C
  function CliffordOrderElem{T, C}(CO::CliffordOrder{T, C}) where {T, C}
    ambalg = CO.algebra
    newelt = new{T, C}(CO, fill(ambalg.base_ring(), CO.rank))
    _set_even_odd_coefficients!(newelt)
    return newelt
  end

  CliffordOrderElem(CO::CliffordOrder) = CliffordOrderElem{elem_type(CO.base_ring), typeof(CO.algebra)}(CO)

  #Return the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function CliffordOrderElem{T, C}(CO::CliffordOrder{T, C}, coeff::Vector{S}) where {T, C, S<:NumFieldElem}
    @req length(coeff) == rank(CO) "invalid length of coefficient vector"
    
    for i in 1:rank(CO)
      @req coeff[i] in coefficient_ideals(CO)[i] "The element does not lie in the Clifford order."
    end

    newelt = new{T, C}(CO, coeff)
    _set_even_odd_coefficients!(newelt)
    return newelt
  end

  function CliffordOrderElem(CO::CliffordOrder{T, C}, coeff::Vector{S}) where {T, C, S}
    K = base_ring(algebra(CO))
    @req _can_convert_coefficients(coeff, K) "entries of coefficient vector are not contained in $(K)"
    return CliffordOrderElem{elem_type(base_ring(CO)), typeof(algebra(CO))}(CO, K.(coeff))
  end

end

##### Elements #####
mutable struct ZZCliffordOrderElem <: Hecke.AbstractAssociativeAlgebraElem{ZZRingElem}
  parent::ZZCliffordOrder
  coeffs::Vector{QQFieldElem}
  even_coeffs::Vector{QQFieldElem}
  odd_coeffs::Vector{QQFieldElem}

  #Return the 0-element of the Clifford order CO
  function ZZCliffordOrderElem(CO::ZZCliffordOrder)
    newelt = new(CO, fill(QQ(), CO.rank))
    _set_even_odd_coefficients!(newelt)
    return newelt
  end

  #Return the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function ZZCliffordOrderElem(CO::ZZCliffordOrder, coeff::Vector{QQFieldElem})
    @req length(coeff) == CO.rank "invalid length of coefficient vector"
    for i in 1:CO.rank
      @req is_integer(coeff[i]) "The element does not lie in the Clifford order."
    end
    newelt = new(CO, coeff)
    _set_even_odd_coefficients!(newelt)
    return newelt
  end

  function ZZCliffordOrderElem(CO::ZZCliffordOrder, coeff::Vector{S}) where {S}
    @req _can_convert_coefficients(coeff, QQ) "entries of coefficient vector are not contained in $(QQField)"
    return ZZCliffordOrderElem(CO, QQ.(coeff))
  end

end

#############################################################
#
#  Type computation
#
#############################################################

elem_type(::Type{CliffordOrder{T, C}}) where {T, C} = CliffordOrderElem{T, C}

elem_type(::Type{ZZCliffordOrder}) = ZZCliffordOrderElem

parent_type(::Type{CliffordOrderElem{T, C}}) where {T, C} = CliffordOrder{T, C}

parent_type(::Type{ZZCliffordOrderElem}) = ZZCliffordOrder

base_ring_type(::Type{CliffordOrder{T, C}}) where {T, C} = parent_type(T)

base_ring_type(::Type{ZZCliffordOrder}) = ZZRing

is_domain_type(::Type{CliffordOrderElem{T, C}}) where {T, C} = false

is_domain_type(::Type{ZZCliffordOrderElem}) = false

is_exact_type(::Type{CliffordOrderElem{T, C}}) where {T, C} = true

is_exact_type(::Type{ZZCliffordOrderElem}) = true

#############################################################
#
#  Construction
#
#############################################################

##### Order #####
@doc raw"""
    clifford_order(ls::QuadLat) -> CliffordOrder

Return the Clifford order of the even lattice $ls$.

# Examples

```jldoctest
julia> K,a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(lattice(quadratic_space(K, 2*identity_matrix(K, 2))))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5
with basis [1, sqrt(-5)] with Gram matrix
  [2   0]
  [0   2]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   1//1 * <1, 1>
  Norm: 1
  principal generator 1
  two normal wrt: 1
   1//1 * <1, 1>
  Norm: 1
  principal generator 1
  two normal wrt: 1
```
"""
clifford_order(ls::QuadLat) = CliffordOrder{elem_type(base_ring(ls)), typeof(clifford_algebra(rational_span(ls)))}(ls)
@doc raw"""
    clifford_order(ls::ZZLat) -> ZZCliffordOrder

Return the Clifford order of the even integer lattice $ls$.

# Examples

```jldoctest
julia> C = clifford_order(root_lattice(:E,8))
Clifford order of even integer lattice with Gram matrix
  [ 2   -1    0    0    0    0    0    0]
  [-1    2   -1    0    0    0    0    0]
  [ 0   -1    2   -1    0    0    0   -1]
  [ 0    0   -1    2   -1    0    0    0]
  [ 0    0    0   -1    2   -1    0    0]
  [ 0    0    0    0   -1    2   -1    0]
  [ 0    0    0    0    0   -1    2    0]
  [ 0    0   -1    0    0    0    0    2]
```
"""
clifford_order(ls::ZZLat) = ZZCliffordOrder(ls)

##### Elements #####
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

function (C::CliffordOrder)(a::CliffordOrderElem)
  @req parent(a) === C "The element does not lie in the Clifford order"
  return a
end

### ZZ ###

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

function (C::ZZCliffordOrder)(a::ZZCliffordOrderElem)
  @req parent(a) === C "The element does not lie in the Clifford order"
  return a
end

################################################################################
#
#  Element conversion to ambient algebra and vice versa
#
################################################################################

function (C::CliffordAlgebra)(elt::CliffordOrderElem)
  @req C === algebra(parent(elt)) "The Clifford algebra provided and the
    ambient algebra of the Clifford order containing the element provided are not
    identical."
  return C(coefficients(elt))
end

function (C::CliffordOrder)(elt::CliffordAlgebraElem)
  @req algebra(C) === parent(elt) "The ambient algebra of the Clifford order
    provided and the parent object of the element provided are not identical."
    return C(coefficients(elt))
end

function (C::CliffordAlgebra)(elt::ZZCliffordOrderElem)
  @req C === algebra(parent(elt)) "The Clifford algebra provided and the
    ambient algebra of the Clifford order containing the element provided are not
    identical."
  return C(coefficients(elt))
end

function (C::ZZCliffordOrder)(elt::CliffordAlgebraElem)
  @req algebra(C) === parent(elt) "The ambient algebra of the Clifford order
    provided and the parent object of the element provided are not identical."
    return C(coefficients(elt))
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
  if !(algebra(C) === parent(x))
    return false
  end
  coe, coeids = coefficients(x), coefficient_ideals(C)
  for i in 1:rank(C)
    if !(coe[i] in coeids[i])
      return false
    end
  end
  return true 
end

### ZZ ###

@doc raw"""
    in(x::CliffordAlgebraElem, C::ZZCliffordOrder)

Return true if the element $x$ is contained in the Clifford order $C$.
"""
function Base.in(x::CliffordAlgebraElem, C::ZZCliffordOrder)
  if !(algebra(C) === parent(x))
    return false
  end
  coe = coefficients(x)
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

##### Order  #####
@doc raw"""
    base_ring(C::CliffordOrder) -> Ring

Return the base ring of the Clifford order $C$.
"""
base_ring(C::CliffordOrder) = C.base_ring::base_ring_type(typeof(C))

@doc raw"""
    algebra(C::CliffordOrder) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of $C$.
"""
algebra(C::CliffordOrder) = C.algebra

@doc raw"""
    rank(C::CliffordOrder) -> Int

Return the rank of the Clifford order $C$ over its base ring.
"""
rank(C::CliffordOrder) = C.rank

@doc raw"""
    dim(C::CliffordOrder) -> Int

Alias for `rank`
"""
dim(C::CliffordOrder) = rank(C)

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
coefficient_ideals(C::CliffordOrder) = C.coefficient_ideals::Vector{Hecke.fractional_ideal_type(base_ring_type(C))}

### ZZ ###

@doc raw"""
    base_ring(C::ZZCliffordOrder) -> ZZRing

Return the base ring of the Clifford order $C$.
"""
base_ring(C::ZZCliffordOrder) = C.base_ring::ZZRing

@doc raw"""
    algebra(C::ZZCliffordOrder) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of $C$.
"""
algebra(C::ZZCliffordOrder) = C.algebra

@doc raw"""
    rank(C::ZZCliffordOrder) -> Int

Return the rank of the Clifford order $C$ over its base ring.
"""
rank(C::ZZCliffordOrder) = C.rank

@doc raw"""
    dim(C::ZZCliffordOrder) -> Int

Alias for `rank`
"""
dim(C::ZZCliffordOrder) = rank(C)

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


##### Elements #####
@doc raw"""
    parent(x::CliffordOrderElem) -> CliffordOrder

Return the Clifford order containing $x$.
"""
parent(x::CliffordOrderElem) = x.parent

@doc raw"""
  coefficients(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ with respect to the
canonical pseudo-basis of its parent Clifford order.
"""
coefficients(x::CliffordOrderElem) = x.coeffs::Vector{elem_type(base_ring(algebra(parent(x))))}

@doc raw"""
    even_coefficients(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coefficients(x::CliffordOrderElem)
  if isdefined(x, :even_coefficientss)
    return x.even_coeffs::Vector{elem_type(base_ring(algebra(parent(x))))}
  end
  _set_even_odd_coefficients!(x)
  return x.even_coeffs::Vector{elem_type(base_ring(algebra(parent(x))))}
end

@doc raw"""
    odd_coefficients(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coefficients(x::CliffordOrderElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs::Vector{elem_type(base_ring(algebra(parent(x))))}
  end
  _set_even_odd_coefficients!(x)
  return x.odd_coeffs::Vector{elem_type(base_ring(algebra(parent(x))))}
end

### ZZ ###

@doc raw"""
    parent(x::ZZCliffordOrderElem) -> ZZCliffordOrder

Return the Clifford order containing $x$.
"""
parent(x::ZZCliffordOrderElem) = x.parent

@doc raw"""
  coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of $x$ with respect to the
canonical basis of its parent Clifford order.
"""
coefficients(x::ZZCliffordOrderElem) = x.coeffs

@doc raw"""
    even_coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of $x$ with respect to the
canonical basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coefficients(x::ZZCliffordOrderElem)
  if isdefined(x, :even_coeffs)
    return x.even_coeffs
  end
  _set_even_odd_coefficients!(x)
  return x.even_coeffs
end

@doc raw"""
    odd_coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of $x$ with respect to the
canonical basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coefficients(x::ZZCliffordOrderElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs
  end
  _set_even_odd_coefficients!(x)
  return x.odd_coeffs
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
  res = C()
  res[1] = base_ring(algebra(C))(1)
  return res
end

@doc raw"""
    pseudo_basis(C::CliffordOrder, i::Int) -> Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}

Return the $i$-th canonical pseudo-basis element of $C$.

# Examples

```jldoctest
julia> K,a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(lattice(quadratic_space(K, 2*identity_matrix(K, 2))));

julia> pseudo_basis(C,1)
([1, 0, 0, 0], 1//1 * <1, 1>
Norm: 1
principal generator 1
two normal wrt: 1)

julia> pseudo_basis(C,3)
([0, 0, 1, 0], 1//1 * <1, 1>
Norm: 1
principal generator 1
two normal wrt: 1)
```
"""
function pseudo_basis(C::CliffordOrder, i::Int)
  res = C()
  res[i] = base_ring(algebra(C))(1)
  return (res, coefficient_ideals(C)[i])
end

@doc raw"""
    pseudo_gen(C::CliffordOrder, i::Int) -> Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}

Return the $i$-th pseudo-element of the canonical algebra pseudo-generating set
of the Clifford order $C$.

# Examples

```jldoctest
julia> K,a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(lattice(quadratic_space(K, 2*identity_matrix(K, 2))));

julia> pseudo_gen(C,1)
([0, 1, 0, 0], 1//1 * <1, 1>
Norm: 1
principal generator 1
two normal wrt: 1)

julia> pseudo_gen(C,2)
([0, 0, 1, 0], 1//1 * <1, 1>
Norm: 1
principal generator 1
two normal wrt: 1)
```
"""
function pseudo_gen(C::CliffordOrder, i::Int)
  res = C()
  if i<= 0
    res[i] #Throw a BoundsError instead of a DomainError for consistency
  end
  res[2^(i - 1) + 1] = base_ring(algebra(C))(1)
  return (res, coefficient_ideals(lattice(C))[i])
end

@doc raw"""
    pseudo_basis(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return the canonical pseudo-basis of the Clifford order $C$.
"""
pseudo_basis(C::CliffordOrder) = map(i -> pseudo_basis(C, i), 1:rank(C))

@doc raw"""
    pseudo_gens(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return the canonical algebra pseudo-generating set of the Clifford order $C$.
"""
pseudo_gens(C::CliffordOrder) = map(i -> pseudo_gen(C, i), 1:ncols(gram_matrix(C))) 

@doc raw"""
    is_commutative(C::CliffordOrder) -> Bool

Return `true` if $C$ is commutative and `false` otherwise.
"""
is_commutative(C::CliffordOrder) = rank(C) == 1 || rank(C) == 2

### ZZ ###

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
  res = C()
  res[1] = QQ(1)
  return res
end

@doc raw"""
    basis(C::ZZCliffordOrder, i::Int) -> ZZCliffordOrderElem

Return the $i$-th canonical basis element of $C$.
"""
function basis(C::ZZCliffordOrder, i::Int)
  res = C()
  res[i] = QQ(1)
  return res
end

@doc raw"""
    gen(C::CliffordOrder, i::Int) -> ZZCliffordOrderElem

Return the $i$-th element of the canonical multiplicative generating set
of the Clifford order $C$.
"""
function gen(C::ZZCliffordOrder, i::Int)
  res = C()
  if i <= 0
    res[i] #Throw a BoundsError for consistency
  end
  res[2^(i - 1) + 1] = QQ(1)
  return res
end

@doc raw"""
    basis(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return the canonical basis of the Clifford order $C$.
"""
basis(C::ZZCliffordOrder) = map(i -> basis(C, i), 1:rank(C))

@doc raw"""
    gens(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return the vector of canonical algebra generators of the Clifford order $C$, i.e.,
if `gram_matrix(C)` is the Gram matrix of the underlying quadratic lattice with respect
to the basis (e_1,...,e_n) then the vector of images under the canonical embedding into
the Clifford order is returned.
"
"""
gens(C::ZZCliffordOrder) = map(i -> gen(C, i), 1:rank(lattice(C)))

@doc raw"""
    is_commutative(C::ZZCliffordOrder) -> Bool

Return `true` if $C$ is commutative and `false` otherwise.
"""
is_commutative(C::ZZCliffordOrder) = rank(C) == 1 || rank(C) == 2

################################################################################
#
#  Element access
#
################################################################################

@doc raw"""
    coeff(x::CliffordOrderElem, i::Int) -> RingElem

Return the `i`-th coefficient of the element `x`.
"""
coeff(x::CliffordOrderElem, i::Int) = coefficients(x)[i]

getindex(x::CliffordOrderElem, i::Int) = coefficients(x)[i]

function setindex!(x::CliffordOrderElem, newentry::RingElement, i::Int64) 
  coefficients(x)[i] = newentry
end

@doc raw"""
    coeff(x::ZZCliffordOrderElem, i::Int) -> ZZRingElem

Return the `i`-th coefficient of the element `x`.
"""
coeff(x::ZZCliffordOrderElem, i::Int64) = coefficients(x)[i]

getindex(x::ZZCliffordOrderElem, i::Int64) = coefficients(x)[i]

function setindex!(x::ZZCliffordOrderElem, newentry::RingElement, i::Int64) 
  coefficients(x)[i] = newentry
end

################################################################################
#
#  Other functionality
#
################################################################################

@doc raw"""
    pseudo_basis_of_centroid(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return a pseudo-basis of the centroid of $C$. Unless `rank(lattice(C)) = 0`, it contains
two pseudo-elements
"""
function pseudo_basis_of_centroid(C::CliffordOrder)
  if isdefined(C, :pseudo_basis_of_centroid)
    return C.pseudo_basis_of_centroid::Vector{Tuple{elem_type(C), Hecke.fractional_ideal_type(base_ring_type(C))}}
  end
  _set_centroid_and_disq!(C)
  return C.pseudo_basis_of_centroid::Vector{Tuple{elem_type(C), Hecke.fractional_ideal_type(base_ring_type(C))}}
end

@doc raw"""
    quadratic_discriminant(C::CliffordOrder) -> Tuple{NumFieldOrderFractionalIdeal, NumFieldElem}

Return the quadratic discriminant of $C$.
"""
function quadratic_discriminant(C::CliffordOrder)
  if isdefined(C, :disq)
    return C.disq::Tuple{Hecke.fractional_ideal_type(base_ring_type(C)), elem_type(base_ring(algebra(C)))} 
  end
  _set_centroid_and_disq!(C)
  return C.disq::Tuple{Hecke.fractional_ideal_type(base_ring_type(C)), elem_type(base_ring(algebra(C)))} 
end

@doc raw"""
    disq(C::CliffordOrder) -> Tuple{NumFieldOrderFractionalIdeal, NumFieldElem}

Alias for `quadratic_discriminant`.
"""
disq(C::CliffordOrder) = quadratic_discriminant(C)

@doc raw"""
    pseudo_basis_of_center(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return a pseudo-basis of the center of $C$. It equals `pseudo_basis_of_centroid(C)`, if and only if
`rank(lattice(C))` is odd. Otherwise it is trivial. 
"""
function pseudo_basis_of_center(C::CliffordOrder)
  if isdefined(C, :pseudo_basis_of_center)
    return C.pseudo_basis_of_center::Vector{Tuple{elem_type(C), Hecke.fractional_ideal_type(base_ring_type(C))}}
  end
  if is_odd(rank(lattice(C)))
    C.pseudo_basis_of_center = pseudo_basis_of_centroid(C)::Vector{Tuple{elem_type(C), Hecke.fractional_ideal_type(base_ring_type(C))}}
  else 
    C.pseudo_basis_of_center = [pseudo_basis(C, 1)]::Vector{Tuple{elem_type(C), Hecke.fractional_ideal_type(base_ring_type(C))}}
  end
end

### ZZ ###

@doc raw"""
    basis_of_centroid(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return a basis of the centroid of $C$. Unless `rank(lattice(C)) = 0`, it is always free of rank
two, so it is returned as a vector containing the basis elements. The first one is the
multiplicative identity of $C$. The square of the second basis element, if present,
equals `quadratic_discriminant(C)`.
"""
function basis_of_centroid(C::ZZCliffordOrder)
  if isdefined(C, :basis_of_centroid)
    return C.basis_of_centroid::Vector{ZZCliffordOrderElem}
  end
  n = rank(lattice(C))
  if n == 0
    C.basis_of_centroid = [one(C)]
    C.disq = ZZ(1)
    return C.basis_of_centroid::Vector{ZZCliffordOrderElem}
  end
  
  T = orthogonal_basis(space(algebra(C)))
  for i in 1:n
    T[i,:] *= 1//gcd(T[i,:])
  end

  orth_elt = prod(map(i -> sum(map(j -> gen(C, j) * T[i, j], 1:n)), 1:n))
  orth_elt *= 1//gcd(coefficients(orth_elt))
  C.disq = ZZ(coefficients(orth_elt^2)[1])
  z_elt = 1//2*algebra(C)(1 + orth_elt)
  
  if z_elt in C 
    C.basis_of_centroid = [one(C), C(z_elt)]
  else
    C.basis_of_centroid = [one(C), orth_elt]
  end
  return C.basis_of_centroid::Vector{ZZCliffordOrderElem}
end

@doc raw"""
    quadratic_discriminant(C::ZZCliffordOrder) -> ZZRingElem

Return the quadratic discriminant of $C$ as an integer.
"""
function quadratic_discriminant(C::ZZCliffordOrder) 
  if isdefined(C, :disq)
    return C.disq
  end
  basis_of_centroid(C)
  return C.disq
end

@doc raw"""
    disq(C::ZZCliffordOrder) -> ZZRingElem

Alias for `quadratic_discriminant`.
"""
disq(C::ZZCliffordOrder) = quadratic_discriminant(C)

@doc raw"""
    basis_of_center(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return a basis of the center of $C$. It equals `basis_of_centroid(C)`, if and only if
`rank(lattice(C))` is odd. Otherwise it contains only the multiplicative identity of $C$.
"""
function basis_of_center(C::ZZCliffordOrder)
  if isdefined(C, :basis_of_center)
    return C.basis_of_center::Vector{ZZCliffordOrderElem}
  end
  if is_odd(rank(lattice(C)))
    C.basis_of_center = basis_of_centroid(C)
  else 
    C.basis_of_center = [one(C)]
  end
  return C.basis_of_center::Vector{ZZCliffordOrderElem}
end

################################################################################
#
#  Inherited functionality for elements
#
################################################################################

# Note: This small section only contains methods on elements that rely on identical methods that exist for the ambient algebra.
# In fact, they convert the given element into an element of the ambient algebra then apply the existing method and interprete
# the result accordingly.


@doc raw"""
    representation_matrix(x::CliffordOrderElem, action::Symbol = :left) -> MatElem

Return the representation matrix of $x$ as an element of `algebra(parent(x))` with respect to `basis(algebra(parent(x)))`. The multiplication is from
the left if `action == :left` and from the right if `action == :right`.
"""
representation_matrix(x::CliffordOrderElem, action::Symbol = :left) = representation_matrix(algebra(parent(x))(x), action)

minimal_polynomial(x::CliffordOrderElem) = minimal_polynomial(algebra(parent(x))(x))

characteristic_polynomial(x::CliffordOrderElem) = characteristic_polynomial(algebra(parent(x))(x))

# Compute a/b if action is :right and b\a if action is :left (and if this is possible)
function divexact(a::CliffordOrderElem, b::CliffordOrderElem, action::Symbol = :left; check::Bool=true)
  check_parent(a, b)
  CO = parent(a)
  CA = algebra(CO)
  res = divexact(CA(a), CA(b), action; check = check)
  if !(res in CO)
    error("Divison not possible")
  end
  return CO(res)
end

@doc raw"""
    divexact_right(a::CliffordOrderElem, b::CliffordOrderElem) -> CliffordOrderElem

Return an element $c$ such that $a = c \cdot b$.
"""
divexact_right(a::CliffordOrderElem, b::CliffordOrderElem; check::Bool=true) = divexact(a, b, :right; check=check)

@doc raw"""
    divexact_left(a::CliffordOrderElem, b::CliffordOrderElem) -> CliffordOrderElem

Return an element $c$ such that $a = b \cdot c$.
"""
divexact_left(a::CliffordOrderElem, b::CliffordOrderElem; check::Bool=true) = divexact(a, b, :left; check=check)

function inv(x::CliffordOrderElem)
  CO = parent(x)
  CA = algebra(CO)
  xinv = inv(CA(x))
  if !(xinv in CO)
    error("Element is not invertible")
  end
  return CO(xinv)
end

### ZZ ###

@doc raw"""
    representation_matrix(x::ZZCliffordOrderElem, action::Symbol = :left) -> MatElem

Return the representation matrix of $x$ as an element of `algebra(parent(x))` with respect to `basis(algebra(parent(x)))`. The multiplication is from
the left if `action == :left` and from the right if `action == :right`.
"""
representation_matrix(x::ZZCliffordOrderElem, action::Symbol = :left) = representation_matrix(algebra(parent(x))(x), action)

minimal_polynomial(x::ZZCliffordOrderElem) = minimal_polynomial(algebra(parent(x))(x))

characteristic_polynomial(x::ZZCliffordOrderElem) = characteristic_polynomial(algebra(parent(x))(x))

# Compute a/b if action is :right and b\a if action is :left (and if this is possible)
function divexact(a::ZZCliffordOrderElem, b::ZZCliffordOrderElem, action::Symbol = :left; check::Bool=true)
  check_parent(a, b)
  CO = parent(a)
  CA = algebra(CO)
  res = divexact(CA(a), CA(b), action; check = check)
  if !(res in CO)
    error("Divison not possible")
  end
  return CO(res)
end

@doc raw"""
    divexact_right(a::ZZCliffordOrderElem, b::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return an element $c$ such that $a = c \cdot b$.
"""
divexact_right(a::ZZCliffordOrderElem, b::ZZCliffordOrderElem; check::Bool=true) = divexact(a, b, :right; check=check)

@doc raw"""
    divexact_left(a::ZZCliffordOrderElem, b::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return an element $c$ such that $a = b \cdot c$.
"""
divexact_left(a::ZZCliffordOrderElem, b::ZZCliffordOrderElem; check::Bool=true) = divexact(a, b, :left; check=check)

function inv(x::ZZCliffordOrderElem)
  CO = parent(x)
  CA = algebra(CO)
  xinv = inv(CA(x))
  if !(xinv in CO)
    error("Element is not invertible")
  end
  return CO(xinv)
end

################################################################################
#
#  unary operators
#
################################################################################

Base.:+(x::CliffordOrderElem) = x

Base.:-(x::CliffordOrderElem) = parent(x)(map(y -> -1 * y, coefficients(x)))

Base.:+(x::ZZCliffordOrderElem) = x

Base.:-(x::ZZCliffordOrderElem) = parent(x)(map(y -> -1 * y, coefficients(x)))

################################################################################
#
#  binary operators
#
################################################################################

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

################################################################################
#
#  Graded parts
#
################################################################################

@doc raw"""
    even_part(x::CliffordOrderElem) -> CliffordOrderElem

Return the projection of $x$ onto the even Clifford order.
"""
even_part(x::CliffordOrderElem) = parent(x)(even_coefficients(x))

@doc raw"""
    odd_part(x::CliffordOrderElem) -> CliffordOrderElem

Return the projection of $x$ onto the odd Clifford order.
"""
odd_part(x::CliffordOrderElem) = parent(x)(odd_coefficients(x))

@doc raw"""
    even_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the projection of $x$ onto the even Clifford order.
"""
even_part(x::ZZCliffordOrderElem) = parent(x)(even_coefficients(x))

@doc raw"""
    odd_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the projection of $x$ onto the odd Clifford order.
"""
odd_part(x::ZZCliffordOrderElem) = parent(x)(odd_coefficients(x))

#############################################################
#
#  Auxillary functions
#
#############################################################

# Compute the ideals present in the canonical pseudo-basis of the Clifford order
# of the lattice 'ls' in its fixed pseudo-basis.
function _set_coefficient_ideals!(ls::QuadLat)
  coeff_ids = coefficient_ideals(ls)
  res = [fractional_ideal(base_ring(ls), base_ring(ls)(1))]
  for i in 1:rank(ls)
    tmpres = map(x -> x*coeff_ids[i], res)
    res = vcat(res, tmpres)
  end
  return res::Vector{typeof(fractional_ideal(base_ring(ls), one(base_ring(ls))))}
end

function _set_even_odd_coefficients!(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
  R = base_ring(algebra(parent(x)))
  d = rank(parent(x))
  x.even_coeffs = [ iseven(count_ones(y - 1)) ? x.coeffs[y] : R() for y in 1:d]
  x.odd_coeffs = x.coeffs - x.even_coeffs
  return x
end


function _can_convert_coefficients(coeff::Vector{S}, K::Field) where {S}
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

#Return the coefficient_ideals of the given lattice 'ls'
_coefficient_ideals_of_lattice(ls::QuadLat) = coefficient_ideals(ls)::Vector{<:NumFieldOrderFractionalIdeal}

function _set_centroid_and_disq!(C::CliffordOrder)
  n = rank(lattice(C))
  if n == 0
    pb1 = pseudo_basis(C, 1)
    C.pseudo_basis_of_centroid = [pb1]
    C.disq = (pb1[2], disq(algebra(C)))
  else 
    br = base_ring(C)
    z_elt = coefficients(basis_of_centroid(algebra(C))[2])
    lambda_empt = z_elt[1]

    #compute ideal intersection
    ideal_array = map(i -> !is_zero(z_elt[i]) ? (z_elt[i])^(-1) * coefficient_ideals(C)[i] : fractional_ideal(br, zero(br)), 2:length(z_elt)) 
    filter!(!is_zero, ideal_array)
    len = length(ideal_array)
    c_ideal = ideal_array[1]
    for i in 2:len
      c_ideal = simplify(lcm(c_ideal, ideal_array[i]))
    end
    if is_zero(lambda_empt)
      C.disq = tuple(disq(algebra(C)) * c_ideal^2, disq(algebra(C)))
      C.pseudo_basis_of_centroid = pseudo_matrix(matrix(base_ring(algebra(C)), 2, 2^n, vcat(coefficients(one(C)), z_elt)), [fractional_ideal(br, one(br)), c_ideal])
    else
      z_elt = (2 * lambda_empt)^(-1) .* z_elt
      c_ideal *= (2 * lambda_empt)
      b_ideal = simplify(lcm(fractional_ideal(br, br(2)), c_ideal))
      C.disq = tuple(simplify((2*lambda_empt)^(-2) * disq(algebra(C)) * b_ideal^2), disq(algebra(C)))
      z_elt[1] += base_ring(algebra(C))(1//2)
      b_ideal = simplify(lcm(fractional_ideal(br, one(br)), c_ideal))
      C.pseudo_basis_of_centroid = [pseudo_basis(C, 1), (C(z_elt), b_ideal)]
    end
  end
end
