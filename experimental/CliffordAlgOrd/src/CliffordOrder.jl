
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

Return the Clifford order of the even lattice `ls`. If the lattice is not even an error is raised.

# Examples

```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(lattice(quadratic_space(K, 2*identity_matrix(K, 2))))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [2   0]
  [0   2]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <1>//1
   <1>//1

julia> ql = quadratic_lattice(K, pseudo_matrix(K[1 0 ;0 1], [ideal(OK, OK(1)), ideal(OK, [OK(2), a + 1])]); gram = K[0 1; 1 0]);

julia> clifford_order(ql)
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   1]
  [1   0]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <1>//1
   <2, sqrt(-5) + 1>//1
```
"""
clifford_order(ls::QuadLat) = CliffordOrder{elem_type(base_ring(ls)), typeof(clifford_algebra(rational_span(ls)))}(ls)

@doc raw"""
    clifford_order(ls::ZZLat) -> ZZCliffordOrder

Return the Clifford order of the even integer lattice `ls`.

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

# for disambiguation
function (C::CliffordOrder)(a::ZZRingElem)
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

# for disambiguation
function (C::ZZCliffordOrder)(a::ZZRingElem)
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
  @req C === ambient_algebra(parent(elt)) "The Clifford algebra provided and the
    ambient algebra of the Clifford order containing the element provided are not
    identical."
  return C(coefficients(elt))
end

function (C::CliffordOrder)(elt::CliffordAlgebraElem)
  @req ambient_algebra(C) === parent(elt) "The ambient algebra of the Clifford order
    provided and the parent object of the element provided are not identical."
    return C(coefficients(elt))
end

function (C::CliffordAlgebra)(elt::ZZCliffordOrderElem)
  @req C === ambient_algebra(parent(elt)) "The Clifford algebra provided and the
    ambient algebra of the Clifford order containing the element provided are not
    identical."
  return C(coefficients(elt))
end

function (C::ZZCliffordOrder)(elt::CliffordAlgebraElem)
  @req ambient_algebra(C) === parent(elt) "The ambient algebra of the Clifford order
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

Return `true` if the element `x` is contained in the Clifford order `C`.
"""
function Base.in(x::CliffordAlgebraElem, C::CliffordOrder)
  ambient_algebra(C) === parent(x) || return false
  coe, coeids = coefficients(x), coefficient_ideals(C)
  for i in 1:rank(C)
    coe[i] in coeids[i] || return false
  end
  return true 
end

### ZZ ###

@doc raw"""
    in(x::CliffordAlgebraElem, C::ZZCliffordOrder)

Return `true` if the element `x` is contained in the Clifford order `C`.
"""
function Base.in(x::CliffordAlgebraElem, C::ZZCliffordOrder)
  ambient_algebra(C) === parent(x) || return false
  coe = coefficients(x)
  for i in 1:rank(C)
    is_integral(coe[i]) || return false
  end
  return true 
end

################################################################################
#
#  Basic field access
#
################################################################################

##### Order  #####
@doc raw"""
    base_ring(C::CliffordOrder) -> Ring

Return the base ring of the Clifford order `C`.
"""
base_ring(C::CliffordOrder) = C.base_ring::base_ring_type(typeof(C))

@doc raw"""
    ambient_algebra(C::CliffordOrder) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of `C`.
"""
ambient_algebra(C::CliffordOrder) = C.ambient_algebra

@doc raw"""
    algebra(C::CliffordOrder) -> CliffordAlgebra

Alias for [`ambient_algebra(C)`](@ref amient_algebra(C::CliffordOrder)).
"""
algebra(C::CliffordOrder) = ambient_algebra(C)

@doc raw"""
    rank(C::CliffordOrder) -> Int

Return the rank of the Clifford order `C` over its base ring.
"""
rank(C::CliffordOrder) = C.rank

@doc raw"""
    dim(C::CliffordOrder) -> Int

Alias for `rank`
"""
dim(C::CliffordOrder) = rank(C)

@doc raw"""
    lattice(C::CliffordOrder) -> QuadLat

Return the underlying even quadratic lattice of the Clifford order `C`.
"""
lattice(C::CliffordOrder) = C.lattice

@doc raw"""
    gram_matrix(C::CliffordOrder) -> MatElem

Return the Gram matrix with respect to the fixed pseudo-basis of the Clifford order `C`.
"""
gram_matrix(C::CliffordOrder) = C.gram

@doc raw"""
    coefficient_ideals(C::CliffordOrder) -> Vector{NumFieldOrderFractionalIdeal}

Return the vector of coefficient ideals of the canonical pseudo-basis of `C`.
"""
coefficient_ideals(C::CliffordOrder) = C.coefficient_ideals::Vector{Hecke.fractional_ideal_type(base_ring_type(C))}

### ZZ ###

@doc raw"""
    base_ring(C::ZZCliffordOrder) -> ZZRing

Return the base ring of the Clifford order `C`.
"""
base_ring(C::ZZCliffordOrder) = C.base_ring::ZZRing

@doc raw"""
    ambient_algebra(C::ZZCliffordOrder) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of `C`.
"""
ambient_algebra(C::ZZCliffordOrder) = C.ambient_algebra

@doc raw"""
    algebra(C::ZZCliffordOrder) -> CliffordAlgebra

Alias for [`ambient_algebra(C)`](@ref amient_algebra(C::ZZCliffordOrder)).
"""
algebra(C::ZZCliffordOrder) = ambient_algebra(C)

@doc raw"""
    rank(C::ZZCliffordOrder) -> Int

Return the rank of the Clifford order `C` over its base ring.
"""
rank(C::ZZCliffordOrder) = C.rank

@doc raw"""
    dim(C::ZZCliffordOrder) -> Int

Alias for `rank`
"""
dim(C::ZZCliffordOrder) = rank(C)

@doc raw"""
    lattice(C::ZZCliffordOrder) -> ZZLat

Return the underlying even quadratic lattice of the Clifford order `C`.
"""
lattice(C::ZZCliffordOrder) = C.lattice

@doc raw"""
    gram_matrix(C::ZZCliffordOrder) -> QQMatrix

Return the Gram matrix with respect to the fixed basis of the Clifford order `C`.
"""
gram_matrix(C::ZZCliffordOrder) = C.gram


##### Elements #####
@doc raw"""
    parent(x::CliffordOrderElem) -> CliffordOrder

Return the Clifford order containing `x`.
"""
parent(x::CliffordOrderElem) = x.parent

@doc raw"""
  coefficients(x::CliffordOrderElem) -> Vector

Return the coefficient vector of `x` with respect to the
canonical pseudo-basis of its parent Clifford order.
"""
coefficients(x::CliffordOrderElem) = x.coeffs::Vector{elem_type(base_ring(ambient_algebra(parent(x))))}

@doc raw"""
    even_coefficients(x::CliffordOrderElem) -> Vector

Return the coefficient vector of `x` wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coefficients(x::CliffordOrderElem)
  if isdefined(x, :even_coefficientss)
    return x.even_coeffs::Vector{elem_type(base_ring(ambient_algebra(parent(x))))}
  end
  _set_even_odd_coefficients!(x)
  return x.even_coeffs::Vector{elem_type(base_ring(ambient_algebra(parent(x))))}
end

@doc raw"""
    odd_coefficients(x::CliffordOrderElem) -> Vector

Return the coefficient vector of `x` wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coefficients(x::CliffordOrderElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs::Vector{elem_type(base_ring(ambient_algebra(parent(x))))}
  end
  _set_even_odd_coefficients!(x)
  return x.odd_coeffs::Vector{elem_type(base_ring(ambient_algebra(parent(x))))}
end

### ZZ ###

@doc raw"""
    parent(x::ZZCliffordOrderElem) -> ZZCliffordOrder

Return the Clifford order containing `x`.
"""
parent(x::ZZCliffordOrderElem) = x.parent

@doc raw"""
  coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of `x` with respect to the
canonical basis of its parent Clifford order.
"""
coefficients(x::ZZCliffordOrderElem) = x.coeffs

@doc raw"""
    even_coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of `x` with respect to the
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

Return the coefficient vector of `x` with respect to the
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

Return the additive identity of the Clifford order `C`.
"""
zero(C::CliffordOrder) = C()

@doc raw"""
    one(C::CliffordOrder) -> CliffordOrderElem

Return the multiplicative identity of the Clifford order `C`.
"""
one(C::CliffordOrder) = C(1)

@doc raw"""
    pseudo_basis(C::CliffordOrder, i::Int) -> Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}

Return the `i`-th canonical pseudo-basis element of `C`.

# Examples

```jldoctest
julia> K,a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(lattice(quadratic_space(K, 2*identity_matrix(K, 2))));

julia> pseudo_basis(C,1)
([1, 0, 0, 0], <1>//1)

julia> pseudo_basis(C,3)
([0, 0, 1, 0], <1>//1)
```
"""
function pseudo_basis(C::CliffordOrder, i::Int)
  R = base_ring(ambient_algebra(C))
  coeffs = fill(zero(R), rank(C))
  coeffs[i] = one(R)
  return (C(coeffs), coefficient_ideals(C)[i])
end

@doc raw"""
    pseudo_gen(C::CliffordOrder, i::Int) -> Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}

Return the `i`-th pseudo-element of the canonical algebra pseudo-generating set
of the Clifford order `C`.

# Examples

```jldoctest
julia> K,a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(lattice(quadratic_space(K, 2*identity_matrix(K, 2))));

julia> pseudo_gen(C,1)
([0, 1, 0, 0], <1>//1)

julia> pseudo_gen(C,2)
([0, 0, 1, 0], <1>//1)
```
"""
function pseudo_gen(C::CliffordOrder, i::Int)
  R = base_ring(ambient_algebra(C))
  coeffs = fill(zero(R), rank(C))
  i > 0 || throw(BoundsError(coeffs, i))
  coeffs[2^(i - 1) + 1] = one(R)
  return (C(coeffs), coefficient_ideals(C)[i])
end

@doc raw"""
    pseudo_basis(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return the canonical pseudo-basis of the Clifford order `C`.
"""
pseudo_basis(C::CliffordOrder) = map(i -> pseudo_basis(C, i), 1:rank(C))

@doc raw"""
    pseudo_gens(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return the canonical algebra pseudo-generating set of the Clifford order `C`.
"""
pseudo_gens(C::CliffordOrder) = map(i -> pseudo_gen(C, i), 1:ncols(gram_matrix(C))) 

@doc raw"""
    is_commutative(C::CliffordOrder) -> Bool

Return `true` if `C` is commutative and `false` otherwise.
"""
is_commutative(C::CliffordOrder) = rank(C) == 1 || rank(C) == 2

### ZZ ###

@doc raw"""
    zero(C::ZZCliffordOrder) -> ZZCliffordOrderElem

Return the additive identity of the Clifford order `C`.
"""
zero(C::ZZCliffordOrder) = C()

@doc raw"""
    one(C::ZZCliffordOrder) -> ZZCliffordOrderElem

Return the multiplicative identity of the Clifford order `C`.
"""
one(C::ZZCliffordOrder) = C(1)

@doc raw"""
    basis(C::ZZCliffordOrder, i::Int) -> ZZCliffordOrderElem

Return the `i`-th canonical basis element of `C`.
"""
function basis(C::ZZCliffordOrder, i::Int)
  coeffs = fill(QQ(0), rank(C))
  coeffs[i] = QQ(1)
  return C(coeffs)
end

@doc raw"""
    gen(C::CliffordOrder, i::Int) -> ZZCliffordOrderElem

Return the `i`-th element of the canonical multiplicative generating set
of the Clifford order `C`.
"""
function gen(C::ZZCliffordOrder, i::Int)
  coeffs = fill(QQ(0), rank(C))
  i > 0 || throw(BoundsError(coeffs, i)) 
  coeffs[2^(i - 1) + 1] = QQ(1)
  return C(coeffs)
end

@doc raw"""
    basis(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return the canonical basis of the Clifford order `C`.
"""
basis(C::ZZCliffordOrder) = map(i -> basis(C, i), 1:rank(C))

@doc raw"""
    gens(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return the vector of canonical algebra generators of the Clifford order `C`, i.e.,
if `gram_matrix(C)` is the Gram matrix of the underlying quadratic lattice with respect
to the basis (e_1,...,e_n) then the vector of images under the canonical embedding into
the Clifford order is returned.
"
"""
gens(C::ZZCliffordOrder) = map(i -> gen(C, i), 1:rank(lattice(C)))

@doc raw"""
    is_commutative(C::ZZCliffordOrder) -> Bool

Return `true` if `C` is commutative and `false` otherwise.
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
  coefficients(x)[i] = number_field(base_ring(x))(newentry)
  _set_even_odd_coefficients!(x)
end

@doc raw"""
    coeff(x::ZZCliffordOrderElem, i::Int) -> ZZRingElem

Return the `i`-th coefficient of the element `x`.
"""
coeff(x::ZZCliffordOrderElem, i::Int64) = coefficients(x)[i]

getindex(x::ZZCliffordOrderElem, i::Int64) = coefficients(x)[i]

function setindex!(x::ZZCliffordOrderElem, newentry::RingElement, i::Int64) 
  coefficients(x)[i] = QQ(newentry)
  _set_even_odd_coefficients!(x)
end

################################################################################
#
#  Other functionality
#
################################################################################

@doc raw"""
    pseudo_basis_of_centroid(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return a pseudo-basis of the centroid of `C`. Unless `rank(lattice(C)) = 0`, it contains
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

Return the quadratic discriminant of `C`.
"""
function quadratic_discriminant(C::CliffordOrder)
  if isdefined(C, :disq)
    return C.disq::Tuple{Hecke.fractional_ideal_type(base_ring_type(C)), elem_type(base_ring(ambient_algebra(C)))} 
  end
  _set_centroid_and_disq!(C)
  return C.disq::Tuple{Hecke.fractional_ideal_type(base_ring_type(C)), elem_type(base_ring(ambient_algebra(C)))} 
end

@doc raw"""
    disq(C::CliffordOrder) -> Tuple{NumFieldOrderFractionalIdeal, NumFieldElem}

Alias for `quadratic_discriminant`.
"""
disq(C::CliffordOrder) = quadratic_discriminant(C)

@doc raw"""
    pseudo_basis_of_center(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return a pseudo-basis of the center of `C`. It equals `pseudo_basis_of_centroid(C)`, if and only if
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

Return a basis of the centroid of `C`. Unless `rank(lattice(C)) = 0`, it is always free of rank
two, so it is returned as a vector containing the basis elements. The first one is the
multiplicative identity of `C`. The square of the second basis element, if present,
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
  
  T = orthogonal_basis(space(ambient_algebra(C)))
  for i in 1:n
    T[i,:] *= 1//gcd(T[i,:])
  end

  orth_elt = prod(map(i -> sum(map(j -> gen(C, j) * T[i, j], 1:n)), 1:n))
  orth_elt *= 1//gcd(coefficients(orth_elt))
  C.disq = ZZ(coefficients(orth_elt^2)[1])
  z_elt = 1//2*ambient_algebra(C)(1 + orth_elt)
  
  if z_elt in C 
    C.basis_of_centroid = [one(C), C(z_elt)]
  else
    C.basis_of_centroid = [one(C), orth_elt]
  end
  return C.basis_of_centroid::Vector{ZZCliffordOrderElem}
end

@doc raw"""
    quadratic_discriminant(C::ZZCliffordOrder) -> ZZRingElem

Return the quadratic discriminant of `C` as an integer.
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

Return a basis of the center of `C`. It equals `basis_of_centroid(C)`, if and only if
`rank(lattice(C))` is odd. Otherwise it contains only the multiplicative identity of `C`.
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

Return the representation matrix of `x` as an element of `ambient_algebra(parent(x))` with respect to `basis(ambient_algebra(parent(x)))`. The multiplication is from
the left if `action == :left` and from the right if `action == :right`.
"""
representation_matrix(x::CliffordOrderElem, action::Symbol = :left) = representation_matrix(ambient_algebra(parent(x))(x), action)

minimal_polynomial(x::CliffordOrderElem) = minimal_polynomial(ambient_algebra(parent(x))(x))

characteristic_polynomial(x::CliffordOrderElem) = characteristic_polynomial(ambient_algebra(parent(x))(x))

### ZZ ###

@doc raw"""
    representation_matrix(x::ZZCliffordOrderElem, action::Symbol = :left) -> MatElem

Return the representation matrix of `x` as an element of `ambient_algebra(parent(x))` with respect to `basis(ambient_algebra(parent(x)))`. The multiplication is from
the left if `action == :left` and from the right if `action == :right`.
"""
representation_matrix(x::ZZCliffordOrderElem, action::Symbol = :left) = representation_matrix(ambient_algebra(parent(x))(x), action)

minimal_polynomial(x::ZZCliffordOrderElem) = minimal_polynomial(ambient_algebra(parent(x))(x))

characteristic_polynomial(x::ZZCliffordOrderElem) = characteristic_polynomial(ambient_algebra(parent(x))(x))

################################################################################
#
#  Graded parts
#
################################################################################

@doc raw"""
    even_part(x::CliffordOrderElem) -> CliffordOrderElem

Return the projection of `x` onto the even Clifford order.
"""
even_part(x::CliffordOrderElem) = parent(x)(even_coefficients(x))

@doc raw"""
    is_even(x::CliffordOrderElem) -> Bool

Return 'true' if 'x' is even, i.e. if [`even_part(x)`](@ref even_part(x::CliffordOrderElem))
and `x` coincide. Otherwise, return `false`.
"""
is_even(x::CliffordOrderElem) = (even_part(x) == x)

@doc raw"""
    odd_part(x::CliffordOrderElem) -> CliffordOrderElem

Return the projection of `x` onto the odd Clifford order.
"""
odd_part(x::CliffordOrderElem) = parent(x)(odd_coefficients(x))

@doc raw"""
    is_odd(x::CliffordOrderElem) -> Bool

Return 'true' if 'x' is odd, i.e. if [`odd_part(x)`](@ref odd_part(x::CliffordOrderElem))
and `x` coincide. Otherwise, return `false`.
"""
is_odd(x::CliffordOrderElem) = (odd_part(x) == x)

@doc raw"""
    even_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the projection of `x` onto the even Clifford order.
"""
even_part(x::ZZCliffordOrderElem) = parent(x)(even_coefficients(x))

@doc raw"""
    is_even(x::ZZCliffordOrderElem) -> Bool

Return 'true' if 'x' is even, i.e. if [`even_part(x)`](@ref even_part(x::ZZCliffordOrderElem))
and `x` coincide. Otherwise, return `false`.
"""
is_even(x::ZZCliffordOrderElem) = (even_part(x) == x)

@doc raw"""
    odd_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the projection of `x` onto the odd Clifford order.
"""
odd_part(x::ZZCliffordOrderElem) = parent(x)(odd_coefficients(x))

@doc raw"""
    is_odd(x::ZZCliffordOrderElem) -> Bool

Return 'true' if 'x' is odd, i.e. if [`odd_part(x)`](@ref odd_part(x::ZZCliffordOrderElem))
and `x` coincide. Otherwise, return `false`.
"""
is_odd(x::ZZCliffordOrderElem) = (odd_part(x) == x)

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
  R = base_ring(ambient_algebra(parent(x)))
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
    C.disq = (pb1[2], disq(ambient_algebra(C)))
  else 
    br = base_ring(C)
    z_elt = coefficients(basis_of_centroid(ambient_algebra(C))[2])
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
      C.disq = tuple(disq(ambient_algebra(C)) * c_ideal^2, disq(ambient_algebra(C)))
      C.pseudo_basis_of_centroid = pseudo_matrix(matrix(base_ring(ambient_algebra(C)), 2, 2^n, vcat(coefficients(one(C)), z_elt)), [fractional_ideal(br, one(br)), c_ideal])
    else
      z_elt = (2 * lambda_empt)^(-1) .* z_elt
      c_ideal *= (2 * lambda_empt)
      b_ideal = simplify(lcm(fractional_ideal(br, br(2)), c_ideal))
      C.disq = tuple(simplify((2*lambda_empt)^(-2) * disq(ambient_algebra(C)) * b_ideal^2), disq(ambient_algebra(C)))
      z_elt[1] += base_ring(ambient_algebra(C))(1//2)
      b_ideal = simplify(lcm(fractional_ideal(br, one(br)), c_ideal))
      C.pseudo_basis_of_centroid = [pseudo_basis(C, 1), (C(z_elt), b_ideal)]
    end
  end
end
