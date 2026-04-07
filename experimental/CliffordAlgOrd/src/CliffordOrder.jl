
#############################################################
#
#  Type computation
#
#############################################################

elem_type(::Type{CliffordOrder{T, C}}) where {T, C} = CliffordOrderElem{T, C, elem_type(base_ring_type(C))}
elem_type(::Type{ZZCliffordOrder}) = ZZCliffordOrderElem

parent_type(::Type{CliffordOrderElem{T, C, S}}) where {T, C, S} = CliffordOrder{T, C}
parent_type(::Type{ZZCliffordOrderElem}) = ZZCliffordOrder

base_ring_type(::Type{CliffordOrder{T, C}}) where {T, C} = parent_type(T)
base_ring_type(::Type{ZZCliffordOrder}) = ZZRing

is_domain_type(::Type{CliffordOrderElem{T, C, S}}) where {T, C, S} = false
is_domain_type(::Type{ZZCliffordOrderElem}) = false

is_exact_type(::Type{CliffordOrderElem{T, C, S}}) where {T, C, S} = true
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

Return the Clifford order of the even integer lattice `ls`. If the lattice is not even, an error is raised.

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

function (C::CliffordAlgebra)(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
  @req C === ambient_algebra(parent(x)) "The provided Clifford algebra and the
    ambient algebra of the Clifford order containing the provided element are not
    identical."
  return C(coefficients(x))
end

function (C::Union{CliffordOrder, ZZCliffordOrder})(x::CliffordAlgebraElem)
  @req ambient_algebra(C) === parent(x) "The ambient algebra of the provided Clifford
    order and the parent object of the element provided are not identical."
    return C(coefficients(x))
end

################################################################################
#
#  Containment
#
################################################################################

@doc raw"""
    in(x::CliffordAlgebraElem, C::CliffordOrder) -> Bool
    in(x::CliffordAlgebraElem, C::ZZCliffordOrder) -> Bool

Return `true` if the element `x` is contained in the Clifford order `C`.

# Examples
For Clifford orders over rings of integers, it is tested if all coefficients of a given element
lie in fractional ideal associated to the entry.
```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(K[1 0; 0 1], [ideal(OK, OK(2)), ideal(OK, OK(3))]); gram = K[0 1; 1 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   1]
  [1   0]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1

julia> x = pseudo_basis(C, 2)[1]
[0, 1, 0, 0]

julia> x in C
false

julia> 2*x in C
true
```
For Clifford order over the integers, it is tested of all coefficients of a given elements are integers.
```jldoctest
julia> C = clifford_order(root_lattice(:A, 2))
Clifford order of even integer lattice with Gram matrix
  [ 2   -1]
  [-1    2]

julia> CA = ambient_algebra(C); (x, y) = CA([6//3, 0, -4//2, 6]), CA([1, 4, 5//2, -2])
([2, 0, -2, 6], [1, 4, 5//2, -2])

julia> x in C
true

julia> y in C
false
```
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
    base_ring(C::ZZCliffordOrder) -> ZZRing

Return the base ring of the Clifford order `C`.
"""
base_ring(C::CliffordOrder) = C.base_ring::base_ring_type(typeof(C))
base_ring(C::ZZCliffordOrder) = C.base_ring::ZZRing

@doc raw"""
    ambient_algebra(C::Union{CliffordOrder, ZZCliffordOrder}) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of `C`.
"""
ambient_algebra(C::Union{CliffordOrder, ZZCliffordOrder}) = C.ambient_algebra

@doc raw"""
    algebra(C::CliffordOrder) -> CliffordAlgebra

Alias for [`ambient_algebra(C)`](@ref ambient_algebra(C::Union{CliffordOrder, ZZCliffordOrder})).
"""
algebra(C::Union{CliffordOrder, ZZCliffordOrder}) = ambient_algebra(C)

@doc raw"""
    rank(C::Union{CliffordOrder, ZZCliffordOrder}) -> Int

Return the rank of the Clifford order `C` over its base ring.
"""
rank(C::Union{CliffordOrder, ZZCliffordOrder}) = C.rank

@doc raw"""
    lattice(C::CliffordOrder) -> QuadLat
    lattice(C::ZZCliffordOrder) -> ZZLat

Return the underlying even quadratic lattice of the Clifford order `C`.
"""
lattice(C::Union{CliffordOrder, ZZCliffordOrder}) = C.lattice

@doc raw"""
    gram_matrix(C::CliffordOrder) -> MatElem
    gram_matrix(C::ZZCliffordOrder) -> QQMatrix

Return the Gram matrix with respect to the fixed pseudo-basis of the Clifford order `C`.
"""
gram_matrix(C::Union{CliffordOrder, ZZCliffordOrder}) = C.gram

@doc raw"""
    coefficient_ideals(C::CliffordOrder) -> Vector{NumFieldOrderFractionalIdeal}

Return the vector of coefficient ideals of the canonical pseudo-basis of `C`.
"""
coefficient_ideals(C::CliffordOrder) = C.coefficient_ideals::Vector{Hecke.fractional_ideal_type(base_ring_type(C))}

##### Elements #####
@doc raw"""
    parent(x::CliffordOrderElem) -> CliffordOrder
    parent(x::ZZCliffordOrderElem) -> ZZCliffordOrder

Return the Clifford order containing `x`.
"""
parent(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) = x.parent

@doc raw"""
    coefficients(x::CliffordOrderElem) -> Vector{NumFieldElem}
    coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of `x` with respect to the
canonical (pseudo-)basis of its parent Clifford order.
"""
coefficients(x::CliffordOrderElem) = x.coeffs::Vector{elem_type(base_ring(ambient_algebra(parent(x))))}
coefficients(x::ZZCliffordOrderElem) = x.coeffs

################################################################################
#
#  basic functionality
#
################################################################################

@doc raw"""
    zero(C::CliffordOrder) -> CliffordOrderElem
    zero(C::ZZCliffordOrder) -> ZZCliffordOrderElem

Return the additive identity of the Clifford order `C`.
"""
zero(C::Union{CliffordOrder, ZZCliffordOrder}) = C()

@doc raw"""
    one(C::CliffordOrder) -> CliffordOrderElem
    one(C::ZZCliffordOrder) -> ZZCliffordOrderElem

Return the multiplicative identity of the Clifford order `C`.
"""
one(C::Union{CliffordOrder, ZZCliffordOrder}) = C(1)

@doc raw"""
    pseudo_basis(C::CliffordOrder, i::Int) -> Tuple{CliffordAlgebraElem, NumFieldOrderFractionalIdeal}

Return the `i`-th canonical pseudo-basis element of `C`. The first coordinate is returned as an element
of `ambient_algebra(C)`.

# Examples

```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(K[1 0; 0 1], [ideal(OK, OK(2)), ideal(OK, OK(3))]); gram = K[0 1; 1 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   1]
  [1   0]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1

julia> pseudo_basis(C, 1)
([1, 0, 0, 0], <1>//1)

julia> pseudo_basis(C, 2)
([0, 1, 0, 0], <2>//1)

julia> pseudo_basis(C, 3)
([0, 0, 1, 0], <3>//1)

julia> pseudo_basis(C, 4)
([0, 0, 0, 1], <36, 3162>//1)
```
"""
pseudo_basis(C::CliffordOrder, i::Int) = (basis(ambient_algebra(C), i), coefficient_ideals(C)[i])

@doc raw"""
    pseudo_gen(C::CliffordOrder, i::Int) -> Tuple{CliffordAlgebraElem, NumFieldOrderFractionalIdeal}

Return the `i`-th pseudo-element of the canonical algebra pseudo-generating set
of the Clifford order `C`. The first coordinate is returned as an element of 
`ambient_algebra(C)`.

# Examples

```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(identity_matrix(K, 3), [ideal(OK, OK(2)), ideal(OK, OK(3)), ideal(OK, OK(5))]); gram = K[0 0 1; 0 2 0; 1 0 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   0   1]
  [0   2   0]
  [1   0   0]
and coefficient ideals of the lattice
  3-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1
   <5>//1

julia> pseudo_gen(C, 1)
([0, 1, 0, 0, 0, 0, 0, 0], <2>//1)

julia> pseudo_gen(C, 2)
([0, 0, 1, 0, 0, 0, 0, 0], <3>//1)

julia> pseudo_gen(C, 3)
([0, 0, 0, 0, 1, 0, 0, 0], <5>//1)
```
"""
pseudo_gen(C::CliffordOrder, i::Int) = (gen(ambient_algebra(C), i), coefficient_ideals(lattice(C))[i])

@doc raw"""
    pseudo_basis(C::CliffordOrder) -> Vector{Tuple{CliffordAlgebraElem, NumFieldOrderFractionalIdeal}}

Return the canonical pseudo-basis of the Clifford order `C`.

# Examples

```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(K[1 0; 0 1], [ideal(OK, OK(2)), ideal(OK, OK(3))]); gram = K[0 1; 1 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   1]
  [1   0]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1

julia> pseudo_basis(C)
4-element Vector{Tuple{CliffordAlgebraElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}, AbsSimpleNumFieldOrderFractionalIdeal}}:
 ([1, 0, 0, 0], <1>//1)
 ([0, 1, 0, 0], <2>//1)
 ([0, 0, 1, 0], <3>//1)
 ([0, 0, 0, 1], <36, 3162>//1)
```
"""
pseudo_basis(C::CliffordOrder) = map(i -> pseudo_basis(C, i), 1:rank(C))

@doc raw"""
    pseudo_gens(C::CliffordOrder) -> Vector{Tuple{CliffordAlgebraElem, NumFieldOrderFractionalIdeal}}

Return the canonical algebra pseudo-generating set of the Clifford order `C`.

# Examples

```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(identity_matrix(K, 3), [ideal(OK, OK(2)), ideal(OK, OK(3)), ideal(OK, OK(5))]); gram = K[0 0 1; 0 2 0; 1 0 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   0   1]
  [0   2   0]
  [1   0   0]
and coefficient ideals of the lattice
  3-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1
   <5>//1

julia> pseudo_gens(C)
3-element Vector{Tuple{CliffordAlgebraElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}, AbsSimpleNumFieldOrderFractionalIdeal}}:
 ([0, 1, 0, 0, 0, 0, 0, 0], <2>//1)
 ([0, 0, 1, 0, 0, 0, 0, 0], <3>//1)
 ([0, 0, 0, 0, 1, 0, 0, 0], <5>//1)
```
"""
pseudo_gens(C::CliffordOrder) = map(i -> pseudo_gen(C, i), 1:ncols(gram_matrix(C))) 

@doc raw"""
    is_commutative(C::Union{CliffordOrder, ZZCliffordOrder}) -> Bool

Return `true` if `C` is commutative and `false` otherwise.
"""
is_commutative(C::Union{CliffordOrder, ZZCliffordOrder}) = rank(C) == 1 || rank(C) == 2

### ZZ ###

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
    gen(C::ZZCliffordOrder, i::Int) -> ZZCliffordOrderElem

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

################################################################################
#
#  Element access
#
################################################################################

@doc raw"""
    coeff(x::CliffordOrderElem, i::Int) -> RingElem
    coeff(x::ZZCliffordOrderElem, i::Int) -> ZZRingElem

Return the `i`-th coefficient of the element `x`.
"""
coeff(x::Union{CliffordOrderElem, ZZCliffordOrderElem}, i::Int) = coefficients(x)[i]

getindex(x::Union{CliffordOrderElem, ZZCliffordOrderElem}, i::Int) = coefficients(x)[i]

function setindex!(x::CliffordOrderElem, newentry::RingElement, i::Int)
  ne = number_field(base_ring(x))(newentry)
  ci = coefficient_ideals(parent(x))[i]
  @req ne in ci "$newentry is not contained in the associated fractional ideal $ci"
  coefficients(x)[i] = ne
end

function setindex!(x::ZZCliffordOrderElem, newentry::RingElement, i::Int)
  @req is_integer(newentry) "$newentry is not an integer"
  coefficients(x)[i] = QQ(newentry)
end

################################################################################
#
#  Other functionality
#
################################################################################

@doc raw"""
    pseudo_basis_of_centroid(C::CliffordOrder) -> Vector{Tuple{CliffordAlgebraElem, NumFieldOrderFractionalIdeal}}

Return a pseudo-basis of the centroid of `C`. Unless `rank(lattice(C)) = 0`, it contains
two pseudo-elements, so it is returned as a vector containing these.

# Examples
```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(K[1 0; 0 1], [ideal(OK, OK(2)), ideal(OK, OK(3))]); gram = K[0 1; 1 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   1]
  [1   0]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1

julia> pseudo_basis_of_centroid(C)
2-element Vector{Tuple{CliffordAlgebraElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}, AbsSimpleNumFieldOrderFractionalIdeal}}:
 ([1, 0, 0, 0], <1>//1)
 ([1, 0, 0, -1], <6, 6>//1)
```
"""
function pseudo_basis_of_centroid(C::CliffordOrder)
  rank(lattice(C)) == 0 && return [pseudo_basis(C, 1)]
  z_amb, c_ideal_raw = _compute_raw_orth_data!(C)
  is_zero(coeff(z_amb, 1)) && return [pseudo_basis(C, 1), (z_amb, c_ideal_raw)]
  z_cent = deepcopy(coefficients(z_amb))
  z_cent[1] += 1//2
  br = base_ring(C)
  b_ideal = simplify(lcm(fractional_ideal(br, one(br)), c_ideal_raw))
  return [pseudo_basis(C, 1), (ambient_algebra(C)(z_cent), b_ideal)]
end

@doc raw"""
    pseudo_basis_of_max_orth_suborder_of_centroid(C::CliffordOrder) -> Vector{Tuple{CliffordAlgebraElem, NumFieldOrderFractionalIdeal}}

Return a pseudo-basis of the maximal orthogonal suborder of the centroid of `C`. Unless
`rank(lattice(C)) = 0`, it consists of two pseudo-elements, so it is returned as a vector
containing these. The first one is the multiplicative identity of `C` paired with the
ring of integers as fractional ideal. The square of the first coordinate of the second
pseudo-element, if present, equals the second coordinate of `quadratic_discriminant(C)`
if `rank(lattice(C))` is odd, and a quarter of it if the rank is even.

# Examples

```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(K[1 0; 0 1], [ideal(OK, OK(2)), ideal(OK, OK(3))]); gram = K[0 1; 1 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   1]
  [1   0]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1

julia> pseudo_basis_of_max_orth_suborder_of_centroid(C)
2-element Vector{Tuple{CliffordAlgebraElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}, AbsSimpleNumFieldOrderFractionalIdeal}}:
 ([1, 0, 0, 0], <1>//1)
 ([1//2, 0, 0, -1], <6, 6*sqrt(-5)>//1)

julia> x = pseudo_basis_of_max_orth_suborder_of_centroid(C)[2][1]; (4 * x^2)[1] == disq(C)[2]
true
```
"""
function pseudo_basis_of_max_orth_suborder_of_centroid(C::CliffordOrder)
  rank(lattice(C)) == 0 && return [pseudo_basis(C, 1)]
  return [pseudo_basis(C, 1), _max_orth_pseudo_elt(C)]
end

@doc raw"""
    quadratic_discriminant(C::CliffordOrder) -> Tuple{NumFieldOrderFractionalIdeal, NumFieldElem}

Return the quadratic discriminant of `C` as a tuple `(disq_ideal, alg_disq)` where:
* `disq_ideal` is the discriminant ideal of the maximal orthogonal suborder of the centroid of `C`. 
  Specifically, if this suborder is written as $R \cdot 1 \oplus \mathfrak{a} \cdot x$ 
  (where $R$ is the base ring of `C` and $x$ is an orthogonal generator), then `disq_ideal` 
  equals the fractional ideal $\mathfrak{a}^2 x^2$.
* `alg_disq` is the quadratic discriminant of the ambient Clifford algebra. It lies in the 
  same $K$-square class as $x^2$ (rather the scalar coordinate $x^2[1]$), where $K$ is the base field.
"""
function quadratic_discriminant(C::CliffordOrder)
  isdefined(C, :disq) && return C.disq::Tuple{Hecke.fractional_ideal_type(base_ring_type(C)), elem_type(base_ring(ambient_algebra(C)))}
  _compute_raw_orth_data!(C)
  return C.disq::Tuple{Hecke.fractional_ideal_type(base_ring_type(C)), elem_type(base_ring(ambient_algebra(C)))}
end

@doc raw"""
    disq(C::CliffordOrder) -> Tuple{NumFieldOrderFractionalIdeal, NumFieldElem}

Alias for [`quadratic_discriminant`](@ref quadratic_discriminant(::CliffordOrder)).
"""
disq(C::CliffordOrder) = quadratic_discriminant(C)

@doc raw"""
    pseudo_basis_of_center(C::CliffordOrder) -> Vector{Tuple{CliffordOrderElem, NumFieldOrderFractionalIdeal}}

Return a pseudo-basis of the center of `C`. It equals `pseudo_basis_of_centroid(C)`, if and only if
`rank(lattice(C))` is odd. Otherwise it is trivial. 

# Examples

In case that the underlying lattice has even rank, the center of the Clifford Order equals its base ring.
```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> C = clifford_order(quadratic_lattice(K, pseudo_matrix(K[1 0; 0 1], [ideal(OK, OK(2)), ideal(OK, OK(3))]); gram = K[0 1; 1 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   1]
  [1   0]
and coefficient ideals of the lattice
  2-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1

julia> pseudo_basis_of_center(C)
1-element Vector{Tuple{CliffordAlgebraElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}, AbsSimpleNumFieldOrderFractionalIdeal}}:
 ([1, 0, 0, 0], <1>//1)
```
If the underlying lattice has odd rank, then the center and the centroid of the Clifford algebra coincide.
```jldoctest
julia> K, a = quadratic_field(-5); OK = maximal_order(K);

julia> D = clifford_order(quadratic_lattice(K, pseudo_matrix(identity_matrix(K, 3), [ideal(OK, OK(2)), ideal(OK, OK(3)), ideal(OK, OK(5))]); gram = K[0 0 1; 0 2 0; 1 0 0]))
Clifford order of even lattice over maximal order of imaginary quadratic field defined by x^2 + 5 with Gram matrix
  [0   0   1]
  [0   2   0]
  [1   0   0]
and coefficient ideals of the lattice
  3-element Vector{AbsSimpleNumFieldOrderFractionalIdeal}:
   <2>//1
   <3>//1
   <5>//1

julia> pseudo_basis_of_center(D)
2-element Vector{Tuple{CliffordAlgebraElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}, AbsSimpleNumFieldOrderFractionalIdeal}}:
 ([1, 0, 0, 0, 0, 0, 0, 0], <1>//1)
 ([0, 0, -1, 0, 0, 0, 0, -2], <15, 15*sqrt(-5)>//1)

julia> pseudo_basis_of_center(D) == pseudo_basis_of_centroid(D)
true
```
"""
function pseudo_basis_of_center(C::CliffordOrder)
  is_odd(rank(lattice(C))) && return pseudo_basis_of_centroid(C)
  return [pseudo_basis(C, 1)]
end

### ZZ ###

@doc raw"""
    basis_of_centroid(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return a basis of the centroid of `C`. Unless `rank(lattice(C)) = 0`, it is always free of rank
two, so it is returned as a vector containing the basis elements. The first one is the
multiplicative identity of `C`. 

# Examples

```jldoctest
julia> C = clifford_order(root_lattice(:A, 2))
Clifford order of even integer lattice with Gram matrix
  [ 2   -1]
  [-1    2]

julia> basis_of_centroid(C)
2-element Vector{ZZCliffordOrderElem}:
 [1, 0, 0, 0]
 [1, 0, 0, 1]
```
"""
function basis_of_centroid(C::ZZCliffordOrder) 
  rank(lattice(C)) == 0 && return [one(C)]
  moe = _max_orth_elt(C)
  z_elt = 1//2*ambient_algebra(C)(1 + moe)
  z_elt in C && return [one(C), C(z_elt)]
  return [one(C), moe] 
end

@doc raw"""
    basis_of_max_orth_suborder_of_centroid(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return a basis of the maximal orthogonal suborder of the centroid of `C`. Unless
`rank(lattice(C)) = 0`, it is always free of rank two, so it is returned as a vector
containing the basis elements. The first one is the multiplicative identity of `C`.
The square of the second basis element, if present, equals `quadratic_discriminant(C)`.

# Examples
```jldoctest
julia> C = clifford_order(root_lattice(:A, 2))
Clifford order of even integer lattice with Gram matrix
  [ 2   -1]
  [-1    2]

julia> basis_of_max_orth_suborder_of_centroid(C)
2-element Vector{ZZCliffordOrderElem}:
 [1, 0, 0, 0]
 [1, 0, 0, 2]

julia> (C([1, 0, 0, 2])^2)[1] == disq(C)
true
```
"""
function basis_of_max_orth_suborder_of_centroid(C::ZZCliffordOrder)
  rank(lattice(C)) == 0 && return [one(C)]
  return [one(C), _max_orth_elt(C)]
end

@doc raw"""
    quadratic_discriminant(C::ZZCliffordOrder) -> ZZRingElem

Return the quadratic discriminant of `C` as an integer.

# Examples

```jldoctest
julia> C = clifford_order(root_lattice(:A, 2))
Clifford order of even integer lattice with Gram matrix
  [ 2   -1]
  [-1    2]

julia> quadratic_discriminant(C)
-3
```
"""
function quadratic_discriminant(C::ZZCliffordOrder) 
  isdefined(C, :disq) && return C.disq
  _max_orth_elt(C)
  return C.disq
end

@doc raw"""
    disq(C::ZZCliffordOrder) -> ZZRingElem

Alias for [`quadratic_discriminant`](@ref quadratic_discriminant(::ZZCliffordOrder)).
"""
disq(C::ZZCliffordOrder) = quadratic_discriminant(C)

@doc raw"""
    basis_of_center(C::ZZCliffordOrder) -> Vector{ZZCliffordOrderElem}

Return a basis of the center of `C`. It equals `basis_of_centroid(C)`, if and only if
`rank(lattice(C))` is odd. Otherwise it contains only the multiplicative identity of `C`.
"""
function basis_of_center(C::ZZCliffordOrder)
  is_odd(rank(lattice(C))) && return basis_of_centroid(C)
  return [one(C)]
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
    representation_matrix(x::ZZCliffordOrderElem, action::Symbol = :left) -> MatElem

Return the representation matrix of `x` as an element of `ambient_algebra(parent(x))` with respect to `basis(ambient_algebra(parent(x)))`. The multiplication is from
the left if `action == :left` and from the right if `action == :right`.
"""
representation_matrix(x::Union{CliffordOrderElem, ZZCliffordOrderElem}, action::Symbol = :left) = representation_matrix(ambient_algebra(parent(x))(x), action)

minimal_polynomial(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) = minimal_polynomial(ambient_algebra(parent(x))(x))

characteristic_polynomial(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) = characteristic_polynomial(ambient_algebra(parent(x))(x))

################################################################################
#
#  Graded parts
#
################################################################################

@doc raw"""
    even_coefficients(x::CliffordOrderElem) -> Vector{NumFieldElem}
    even_coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of `x` with respect to the
canonical (pseudo-)basis of its parent Clifford order,
but with all its coefficients that correspond to basis
elements with odd grading set to zero.
"""
even_coefficients(x::CliffordOrderElem) = [iseven(count_ones(y - 1)) ? coefficients(x)[y] : zero(base_ring(ambient_algebra(parent(x)))) for y in 1:rank(parent(x))]
even_coefficients(x::ZZCliffordOrderElem) = [iseven(count_ones(y - 1)) ? coefficients(x)[y] : QQ() for y in 1:rank(parent(x))]

@doc raw"""
    even_part(x::CliffordOrderElem) -> CliffordOrderElem
    even_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the projection of `x` onto the even Clifford order.
"""
even_part(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) = parent(x)(even_coefficients(x))

@doc raw"""
    is_even(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) -> Bool

Return 'true' if 'x' is even, i.e. if [`even_part(x)`](@ref even_part(x::Union{CliffordOrderElem, ZZCliffordOrderElem}))
and `x` coincide. Otherwise, return `false`.
"""
is_even(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) = (even_part(x) == x)

@doc raw"""
    odd_coefficients(x::CliffordOrderElem) -> Vector{NumFieldElem}
    odd_coefficients(x::ZZCliffordOrderElem) -> Vector{QQFieldElem}

Return the coefficient vector of `x` with respect to the
canonical (pseudo-)basis of its parent Clifford order,
but with all its coefficients that correspond to basis elements
with even grading set to zero.
"""
odd_coefficients(x::CliffordOrderElem) = [isodd(count_ones(y - 1)) ? coefficients(x)[y] : zero(base_ring(ambient_algebra(parent(x)))) for y in 1:rank(parent(x))]
odd_coefficients(x::ZZCliffordOrderElem) = [isodd(count_ones(y - 1)) ? coefficients(x)[y] : QQ() for y in 1:rank(parent(x))]

@doc raw"""
    odd_part(x::CliffordOrderElem) -> CliffordOrderElem
    odd_part(x::ZZCliffordOrderElem) -> ZZCliffordOrderElem

Return the projection of `x` onto the odd Clifford order.
"""
odd_part(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) = parent(x)(odd_coefficients(x))

@doc raw"""
    is_odd(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) -> Bool

Return 'true' if 'x' is odd, i.e. if [`odd_part(x)`](@ref odd_part(x::Union{CliffordOrderElem, ZZCliffordOrderElem}))
and `x` coincide. Otherwise, return `false`.
"""
is_odd(x::Union{CliffordOrderElem, ZZCliffordOrderElem}) = (odd_part(x) == x)

#############################################################
#
#  Promotion rules 
#
#############################################################

function AbstractAlgebra.promote_rule(::Type{COE}, ::Type{V}) where
        {T<:RingElement, S, COE<:CliffordOrderElem{T, S}, V<:RingElement}
    AbstractAlgebra.promote_rule(T, V) == T ? COE : Union{}
end

AbstractAlgebra.promote_rule(::Type{COE}, ::Type{COE}) where
        {T<:RingElement, S, COE<:CliffordOrderElem{T, S}} = COE

### ZZ ###
function AbstractAlgebra.promote_rule(::Type{ZZCliffordOrderElem}, ::Type{V}) where {V<:RingElement}
    AbstractAlgebra.promote_rule(ZZRingElem, V) == ZZRingElem ? ZZCliffordOrderElem : Union{}
end

AbstractAlgebra.promote_rule(::Type{ZZCliffordOrderElem}, ::Type{ZZCliffordOrderElem}) = ZZCliffordOrderElem

#############################################################
#
#  Random generation 
#
#############################################################

function rand(rng::Random.AbstractRNG, CO::CliffordOrder, v...)
  ids = coefficient_ideals(CO)
  coeffs = [rand(rng, ids[i], v...) for i in 1:CO.rank]
  return CO(coeffs)
end

rand(CO::CliffordOrder, v...) = rand(Random.default_rng(), CO, v...)

### ZZ ###
function rand(rng::Random.AbstractRNG, CO::ZZCliffordOrder, v...)
  coeffs = [rand(rng, base_ring(CO), v...) for _ in 1:CO.rank]
  return CO(coeffs)
end

rand(CO::ZZCliffordOrder, v...) = rand(Random.default_rng(), CO, v...)

#############################################################
#
#  Conformance test element generation 
#
#############################################################

function ConformanceTests.generate_element(CO::CliffordOrder)
  K = base_ring(ambient_algebra(CO))
  coeffs = map(coefficient_ideals(CO)) do id
    sum((K(rand(-3:3)) * b for b in basis(id)), init=zero(K))
  end
  return CO(coeffs)
end

### ZZ ###
ConformanceTests.generate_element(CO::ZZCliffordOrder) = CO([QQ(rand(-3:3)) for _ in 1:rank(CO)])

#############################################################
#
#  Auxillary functions
#
#############################################################

# Compute the ideals present in the canonical pseudo-basis of the Clifford order
# of the lattice 'ls' in its fixed pseudo-basis.
function _coefficient_ideals_of_CO(ls::QuadLat)
  coeff_ids = coefficient_ideals(ls)
  res = [fractional_ideal(base_ring(ls), base_ring(ls)(1))]
  for i in 1:rank(ls)
    tmpres = map(x -> x*coeff_ids[i], res)
    res = vcat(res, tmpres)
  end
  return res::Vector{typeof(fractional_ideal(base_ring(ls), one(base_ring(ls))))}
end

#Return the coefficient_ideals of the given lattice 'ls'
_coefficient_ideals_of_lattice(ls::QuadLat) = coefficient_ideals(ls)::Vector{<:NumFieldOrderFractionalIdeal}

function _max_orth_pseudo_elt(C::CliffordOrder)
  z_amb, c_ideal_raw = _compute_raw_orth_data!(C)
  rank(lattice(C)) == 0 && return (z_amb, c_ideal_raw)
  is_zero(coeff(z_amb, 1)) && return (z_amb, c_ideal_raw)
  br = base_ring(C)
  return (z_amb, simplify(lcm(fractional_ideal(br, br(2)), c_ideal_raw)))
end

function _compute_raw_orth_data!(C::CliffordOrder)
  isdefined(C, :_raw_orth_data) && return C._raw_orth_data::Tuple{elem_type(ambient_algebra(C)), Hecke.fractional_ideal_type(base_ring_type(C))}
    
  ambalg = ambient_algebra(C)
  br = base_ring(C)
  n = rank(lattice(C))
    
  if n == 0
    pb1 = pseudo_basis(C, 1)
    C.disq = tuple(pb1[2], disq(ambalg))
    C._raw_orth_data = pb1
    return C._raw_orth_data::Tuple{elem_type(ambalg), Hecke.fractional_ideal_type(base_ring_type(C))}
  end
    
  z_elt = coefficients(basis_of_centroid(ambalg)[2])
  lambda_empt = z_elt[1]
    
  coeff_ideals = coefficient_ideals(C) 
  valid_ideals = ((z_elt[i])^(-1) * coeff_ideals[i] for i in 2:length(z_elt) if !is_zero(z_elt[i]))
  c_ideal = reduce((x, y) -> simplify(lcm(x, y)), valid_ideals)
    
  if is_zero(lambda_empt)
    C.disq = tuple(disq(ambalg) * c_ideal^2, disq(ambalg))
    C._raw_orth_data = (ambalg(z_elt), c_ideal)
  else
    z_elt = (2 * lambda_empt)^(-1) .* z_elt
    c_ideal = c_ideal * (2 * lambda_empt)
        
    b_ideal = simplify(lcm(fractional_ideal(br, br(2)), c_ideal))
    C.disq = tuple(simplify((2*lambda_empt)^(-2) * disq(ambalg) * b_ideal^2), disq(ambalg))
    C._raw_orth_data = (ambalg(z_elt), c_ideal)
  end
    
  return C._raw_orth_data::Tuple{elem_type(ambalg), Hecke.fractional_ideal_type(base_ring_type(C))}
end

function _max_orth_elt(C::ZZCliffordOrder)
  isdefined(C, :max_orth_elt) && isdefined(C, :disq) && return C.max_orth_elt::ZZCliffordOrderElem
  n = rank(lattice(C))
  if n == 0
    C.max_orth_elt = C(1)
    C.disq = ZZ(1)
    return C.max_orth_elt::ZZCliffordOrderElem
  end
  
  T = orthogonal_basis(space(ambient_algebra(C)))
  for i in 1:n 
    T[i,:] *= 1//gcd(T[i,:])
  end

  orth_gens = [sum(gen(C, j) * T[i, j] for j in 1:n) for i in 1:n]
  orth_elt = prod(orth_gens)
  scale = 1//gcd(coefficients(orth_elt))
  C.max_orth_elt = orth_elt * scale

  if !isdefined(C, :disq)
    sign = n % 4 in (0, 1) ? 1 : -1
    C.disq = ZZ(sign * prod(coefficients(x^2)[1] for x in orth_gens) * scale^2)
  end

  return C.max_orth_elt::ZZCliffordOrderElem
end

