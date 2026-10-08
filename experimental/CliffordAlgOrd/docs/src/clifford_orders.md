```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Clifford orders over rings of integers](@id cliffordorders)

Throughout this page we make the additional assumption that the Dedekind domain ``R`` is the
ring of integers ``\mathcal{O}_K`` of the algebraic number field ``K``. For Clifford orders over rings of
integers, excluding ``\mathbb{Z}``, we introduce the following new types:
- `CliffordOrder{T, C} <: Hecke.AbstractAssociativeAlgebra{T}` for Clifford orders
- `CliffordOrderElem{T, C, S} <: Hecke.AbstractAssociativeAlgebraElem{T}` for elements of Clifford orders
Here, `T` is the element type of the base ring ``\mathcal{O}_K`` and `C` is the type of the Clifford
algebra of the ambient space of the underlying even ``\mathcal{O}_K``-lattice, respectively. Moreover,
``S`` is the element type of ``K``. It is used for storing coefficients of elements with respect
to some pseudo-basis of the Clifford order, as these need not lie in ``\mathcal{O}_K``.

In addition to that we introduce the following new types for Clifford orders over ``\mathbb{Z}``:
- `ZZCliffordOrder <: Hecke.AbstractAssociativeAlgebra{ZZRingElem}` for Clifford orders
- `ZZCliffordOrderElem <: Hecke.AbstractAssociativeAlgebraElem{ZZRingElem}` for elements of Clifford orders

Depending on whether the given even lattice is defined over ``\mathcal{O}_K`` or ``\mathbb{Z}``, we
provide a corresponding constructor to return its Clifford order.

```@docs
clifford_order(ls::QuadLat)
clifford_order(ls::ZZLat)
```

## [Element construction](@id elementconstruction_cliffordo)

Let ``L`` be an even ``R``-lattice with ambient space ``V = KL`` and suppose that the Clifford order
``C(L)`` was constructed with respect to the implicitly chosen pseudo-basis ``(e_i, \mathfrak{a}_i)_{i \in \underline{n}}``
of ``L``. An element ``a`` in the Clifford order ``C(L)`` is stored as it would be as an element of ``C(V)``;
as a vector of length ``2^n`` with entries in ``K``. This vector is the coefficient vector of ``a``
with respect to the ``K``-basis ``(e_I \mid I \subseteq \underline{n})`` of ``C(V)``. Upon creation and
whenever performing any arithmetic operations with ``a \in C(L)``, checks are performed to ensure
that the coefficient of ``e_I`` in this representation of ``a`` is actually contained in ``\mathfrak{a}_I``.

!!! note "Ordering of the basis elements of the Clifford order"
    We use the exact same ordering of the basis elements ``e_I``: They are enumerated in increasing order
    of the integer represented by the characteristic vectors of subsets ``I``, read as a binary number
    with the first entry being the least significant bit.

We provide multiple ways of directly constructing elements of a given Clifford order

- `(C::Union{CliffordOrder, ZZCliffordOrder})()` returns the zero element of the Clifford order `C`
- `(C::Union{CliffordOrder, ZZCliffordOrder})(a::T)` returns `a` as an element of `C`, if possible
- `(C::Union{CliffordOrder, ZZCliffordOrder})(coeffs::Vector)` returns the element in `C` with coefficient
    vector `coeffs`

## [(Pseudo-)basis and (pseudo-)generators](@id basisandgenerators_cliffordo)

Recall that over the integers, every even ``\mathbb{Z}``-lattice is free, so the theory is very similar
to the one over number fields. In particular, its Clifford order is free of rank ``2^n``, has the
``\mathbb{Z}``-basis ``(e_I \mid I \subseteq \underline{n})`` and is generated as a ``\mathbb{Z}``-algebra
by ``\lbrace e_i \mid i \in \underline{n}\rbrace``. We provide the following methods to access these canonical
basis elements and generators of the Clifford order.

```@docs
basis(C::ZZCliffordOrder, i::Int)
basis(C::ZZCliffordOrder)
gen(C::ZZCliffordOrder, i::Int)
gens(C::ZZCliffordOrder)
```

Over the ring of integers ``\mathcal{O}_K \neq \mathbb{Z}``, the Clifford order has the pseudo-basis
``(e_I, \mathfrak{a}_I)_{I \subseteq \underline{n}}`` and has the pseudo-generating system
``\lbrace e_i, \mathfrak{a}_i \rbrace_{i \in \underline{n}}`` as ``\mathcal{O}_K``-algebra. We provide
the following methods to access these canonical pseudo-basis elements and generators of the Clifford order.

```@docs
pseudo_basis(C::CliffordOrder, i::Int)
pseudo_basis(C::CliffordOrder)
pseudo_gen(C::CliffordOrder, i::Int)
pseudo_gens(C::CliffordOrder)
```

!!! note "Return type of pseudo-elements"
    In the methods above, the ``e_I`` or ``e_i`` are returned as elements of the ambient algebra. This is because they
    need not lie in the Clifford order themselves. In fact, this is the case if and only if ``1`` lies in the associated
    fractional ideal ``\mathfrak{a}_I`` or ``\mathfrak{a}_i``.
    In contrast, the methods for Clifford orders over the integers will always return the basis elements and generators
    as elements of the Clifford order, not the ambient algebra.


## [Basic methods for Clifford orders](@id basicmethodsforcliffordorders)

```@docs
zero(C::Union{CliffordOrder, ZZCliffordOrder})
one(C::Union{CliffordOrder, ZZCliffordOrder})
```
```@docs
base_ring(C::CliffordOrder)
lattice(C::Union{CliffordOrder, ZZCliffordOrder})
gram_matrix(C::Union{CliffordOrder, ZZCliffordOrder})
rank(C::Union{CliffordOrder, ZZCliffordOrder})
is_commutative(C::Union{CliffordOrder, ZZCliffordOrder})
```

## [Basic methods for elements of a Clifford order](@id basicmethodsforelementsofacliffordorder)

```@docs
parent(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
coeff(x::Union{CliffordOrderElem, ZZCliffordOrderElem}, i::Int)
coefficients(x::CliffordOrderElem)
```

In addition to the above `coeff`-method, you can also use standard indexing syntax to conveniently access
and modify the coefficients of an element of a Clifford order. When modifying coefficients this way, it
is verified that the new coefficient is valid. For `CliffordOrderElem`, it is checked that the value lies
within the associated fractional ideal of that entry. For `ZZCliffordOrderElem`, it is checked that the
value is an integer.

```julia
#Access the third entry of an element x of a Clifford order
x[3]

#Modify the third entry of an element x of a Clifford order
x[3] = 2
```

### [Containment and conversion](@id containmentandconversion_cliffordo)

Elements of the ambient Clifford algebra ``C(V)`` can be tested for membership in the Clifford order ``C(L)``.

```@docs
Base.in(x::CliffordAlgebraElem, C::CliffordOrder)
```

Additionally, we provide the following constructors to transition elements between the algebra and the order, raising an
error if an algebra element is not contained in the order:

- (C::CliffordAlgebra)(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
- (C::Union{CliffordOrder, ZZCliffordOrder})(x::CliffordAlgebraElem)


### [Functionality for graded parts](@id functionalityforgradedparts_cliffordo)

The Clifford order ``C(L)`` inherits the ``\mathbb{Z}/2\mathbb{Z}``-grading ``C(L) = C_0(L) \oplus C_1(L)`` from
the ``\mathbb{Z}/2\mathbb{Z}``-grading [``C(V) = C_0(V) \oplus C_1(V)``](@ref functionalityforgradedparts_clifford)
via ``C_i(L) = C_i(V) \cap C(L)``, for ``i \in \lbrace 0, 1 \rbrace``. A pseudo-generating set of ``C_i(L)`` is
given by the set of pseudo-basis elements ``(e_I, \mathfrak{a}_I)`` of ``C(L)`` with ``I \: \mathrm{mod} \: 2 = i``.
Thus, ``C_0(L)`` and ``C_1(L)`` are ``\mathcal{O}_K``-submodules of ``C(L)`` of rank ``2^{n-1}``. Moreover, ``C_0(L)``
is an ``\mathcal{O}_K``-suborder of ``C(L)``, called the *even Clifford order* of ``L`` and the *odd part*
``C_1(L)`` is a ``C_0(L)``-bimodule.

Currently, we only provide the following basic functionality on elements for the graded parts of the Clifford order:

```@docs
even_coefficients(x::CliffordOrderElem)
even_part(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
is_even(x::Union{CliffordOrderElem, ZZCliffordOrderElem})

odd_coefficients(x::CliffordOrderElem)
odd_part(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
is_odd(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
```

## [Center, centroid and quadratic discriminant](@id centerandcentroid_cliffordo)

Similar to the [field case](@ref centerandcentroid_clifford), the *centroid* of the even lattice ``L = (L,q)``
(or of ``C = C(L)``), denoted by ``\mathcal{Z}(L)``, is the centraliser of ``C_0(L)`` in ``C``. Thus, we have
``\mathcal{Z}(L) = \mathcal{Z}(V) \cap C(L)``. Unless ``n = \mathrm{dim}(V) = \mathrm{rank}(L) = 0``, in which
case ``C(L) = C_0(L) = \mathcal{Z}(L) \cong \mathcal{O}_K``, the centroid is always an ``\mathcal{O}_K``-order
in the separable quadratic ``K``-algebra ``\mathcal{Z}(V)``. This means that
```math
\mathcal{Z}(L) = \mathcal{O}_K\cdot 1_C \oplus \mathfrak{a}x
```
with a fractional ideal ``\mathfrak{a}`` and some root ``x \in C(V)`` of the polynomial ``X^2 - tX + n \in K[X]``,
where ``t \in \mathfrak{a}^{-1}`` and ``n \in \mathfrak{a}^{-2}``. Every such quadratic ``\mathcal{O}_K``-order
``\Lambda`` contains a unique maximal orthogonal suborder, denoted by ``\Lambda^o``, where orthogonal means that
one can choose ``t = 0`` in the above representation and maximal is meant with respect to set inclusion. 
If we write ``\Lambda^o = \mathcal{O}_K \cdot 1_\Lambda \oplus \mathfrak{b}z`` with ``z \in K\Lambda``
satisfying ``z^2 - d \in K^\times``, then ``\mathcal{Z}(L)^o`` is determined as an ``\mathcal{O}_K``-algebra
up to isomorphism by the fractional ideal ``\mathfrak{b}^2 n`` and the ``K``-square class of ``d``. The pair
``(\mathfrak{b}^2 n, d(K^\times)^2)`` is called the *quadratic discriminant* of ``\Lambda``.

The *quadratic discriminant* of the lattice ``L`` (or of ``C``), denoted by ``\mathrm{disq}(L)``, is defined as
the quadratic discriminant of its centroid ``\mathcal{Z}(L)``.

!!! note "Quadratic discriminant over principal ideal domains"
    Unlike the field case where the quadratic discriminant is simply a square class, for an even lattice it is
    a pair living in ``\mathcal{I}(\mathcal{O}_K) \times K^\times/(K^\times)^2``, where ``\mathcal{I}(\mathcal{O}_K)``
    is the group of fractional ideals of ``\mathcal{O}_K``.

    If ``\mathcal{O}_K`` is a principal ideal domain, the quadratic discriminant can be identified with the
    unique square class ``d(\mathcal{O}_K^\times)^2`` for some ``d \in \mathcal{O}_K`` such that
    ``\mathcal{Z}(L) = \mathcal{O}_K[X]/(X^2 - d)``. In particular, if ``\mathcal{O}_K = \Z``,
    then the quadratic discriminant can be regarded as an integer.

Just as in the field case, the behavior of the center and the centroid heavily depends on the parity of the rank ``n`` of the lattice:

- If the rank ``n`` is even, the centroid is entirely contained within the even Clifford order, i.e. ``\mathcal{Z}(L) \subseteq C_0(L)``.
  Moreover, the center is trivial, ``Z(C) = \mathcal{O}_K``.
- If the rank ``n`` is odd, the centroid is strictly orthogonal, meaning ``\mathcal{Z}(L) = \mathcal{Z}(L)^o``. 
  In this case, the generator of the centroid can be chosen from ``C_1(L)`` such that it has a trace of zero. Furthermore,
  the centroid coincides with the center of the Clifford order, ``Z(C) = \mathcal{Z}(L)``.

For Clifford orders over ``\mathcal{O}_K \neq \Z`` we provide the following methods:

```@docs
pseudo_basis_of_center(::CliffordOrder)
pseudo_basis_of_centroid(::CliffordOrder)
pseudo_basis_of_max_orth_suborder_of_centroid(::CliffordOrder)
quadratic_discriminant(::CliffordOrder)
disq(::CliffordOrder)
```

For Clifford orders over the integers we provide the following methods:

```@docs
basis_of_center(::ZZCliffordOrder)
basis_of_centroid(::ZZCliffordOrder)
basis_of_max_orth_suborder_of_centroid(::ZZCliffordOrder)
quadratic_discriminant(::ZZCliffordOrder)
disq(::ZZCliffordOrder)
```

