```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Clifford Algebras over fields](@id cliffordalgebras) 

For Clifford algebras over fields, we introduce the following new types:
- `CliffordAlgebra{T, S} <: Hecke.AbstractAssociativeAlgebra{T}` for Clifford algebras
- `CliffordAlgebraElem{T, S} <: Hecke.AbstractAssociativeAlgebraElem{T}` for elements of Clifford algebras
Here, `T` is the element type of the base field and `S` is the element type of the Gram matrix of the
underlying quadratic space, respectively. For example, something like `CliffordAlgebra{QQFieldElem, QQMatrix}`
and `CliffordAlgebraElem{QQFieldElem, QQMatrix}` would be expected here.

We provide a constructor that, given a [`quadratic space`](@ref quadratic_space(::NumField, ::Int)) over some algebraic number
field, including the rationals, returns its Clifford algebra.

```@docs
clifford_algebra(qs::Hecke.QuadSpace)
```

## [Element construction](@id elementconstruction_clifford)

The Clifford algebra ``C`` of a quadratic ``K``-space ``V`` of dimension ``n`` is a free ``K``-algebra
of dimension ``2^n``. Since in this project, ``C`` is constructed based on a Gram matrix ``G \in K^{n \times n}``
with respect to some implicitly chosen ``K``-basis ``(e_i \mid I \subseteq \underline{n})`` of ``V``, we
represent an element ``a \in C`` as its coefficient vector with respect to the ``K``-basis
``(e_I \mid I \subseteq \underline{n})`` of ``C``.

!!! note "Ordering of the basis elements of the Clifford algebra"
    The basis elements ``e_I`` are enumerated in increasing order of the integer represented by the
    characteristic vectors of subsets ``I``, read as a binary number with the first entry being the
    least significant bit. For example, if ``n = 3``, then this order is ``e_\emptyset, e_1, e_2, e_{1,2},
    e_{3}, e_{1,3}, e_{2,3}, e_{1,2,3}``, which directly corresponds to the increasing sequence of binary
    numbers ``000, 100, 010, 110, 001, 101, 011, 111``.

    As a consequence, in our data structure the element ``2 - 2e_2 + 3 e_{2,3} - 5e_{1,2,3}`` is
    represented by the vector `[2, 0, -2, 0, 0, 0, 3, -5]`.

We provide multiple ways of directly constructing elements of a given Clifford algebra:
- `(C::ActionPolyRing)()` returns the zero element of the Clifford algebra `C`
- `(C::ActionPolyRing)(a::T)` returns `a` as an element of `C`, if possible
- `(C::CliffordAlgebra)(coeffs::Vector)` returns the element in `C` with coefficient
    vector `coeffs`

## [Basis and generators](@id basisandgenerators_clifford)

In addition to constructing elements of a Clifford algebra directly, we provide the methods below
to allow for direct access to the canonical basis elements ``e_I`` and algebra generators ``e_i``. One
can then use the usual arithmetic operators `+,-,*` to obtain arbitrary linear combinations and products
of these elements.

```@docs
basis(C::CliffordAlgebra, i::Int)
basis(C::CliffordAlgebra)
gen(C::CliffordAlgebra, i::Int)
gens(C::CliffordAlgebra)
```

## [Basic methods for Clifford algebras](@id basicmethodsforcliffordalgebras)

```@docs
zero(C::CliffordAlgebra)
one(C::CliffordAlgebra)
```
```@docs
base_ring(C::CliffordAlgebra)
space(C::CliffordAlgebra)
gram_matrix(C::CliffordAlgebra)
dim(C::CliffordAlgebra)
is_commutative(::CliffordAlgebra)
```

## [Basic methods for elements of a Clifford algebra](@id basicmethodsforelementsofacliffordalgebra)

```@docs
parent(x::CliffordAlgebraElem)
coeff(x::CliffordAlgebraElem, i::Int)
coefficients(x::CliffordAlgebraElem)
```

In addition to the above `coeff`-method, you can also use standard indexing syntax to conveniently access
and modify the coefficients of an element of a Clifford algebra.
```julia
#Access the third entry of an element x of a Clifford algebra
x[3]

#Modify the third entry of an element x of a Clifford algebra
x[3] = 2
```
### [Functionality for graded parts](@id functionalityforgradedparts_clifford)

The Clifford algebra ``C = C(V)`` carries a natural ``\mathbb{Z}/2\mathbb{Z}``-grading ``C = C_0(V) \oplus C_1(V)``, where
``C_i = C_i(V)`` is generated as a ``K``-space by the set of basis elements ``e_I`` of ``C`` with ``I \: \mathrm{mod} \: 2 = i``,
for ``i \in \lbrace 0, 1\rbrace``. In particular, both ``C_0`` and ``C_1`` are ``2^{n-1}``-dimensional subspaces
of ``C``. Moreover, ``C_0`` is a subalgebra of ``C``, called the *even Clifford algebra* of ``V`` and the *odd part*
``C_1`` is a ``C_0``-bimodule.

Currently, we only provide the following basic functionality on elements for the graded parts of the Clifford algebra:

```@docs
even_coefficients(::CliffordAlgebraElem)
even_part(::CliffordAlgebraElem)
is_even(::CliffordAlgebraElem)

odd_coefficients(::CliffordAlgebraElem)
odd_part(::CliffordAlgebraElem)
is_odd(::CliffordAlgebraElem)
```

## [Center, centroid and quadratic discriminant](@id centerandcentroid_clifford)

The *centroid* of ``V = (V,q)`` (or of ``C = C(V)``), denoted by ``\mathcal{Z}(V)``, is defined as the centraliser
of ``C_0(V)`` in ``C``. Unless ``n = \mathrm{dim}(V) = 0``, in which case ``C(V) = C_0(V) = \mathcal{Z}(V) \cong K``,
the centroid is always a separable quadratic ``K``-algebra. This means that ``\mathcal{Z}(V) \cong K[x]/(x^2 - d)``,
for some non-zero ``d \in K``. The square class ``d(K^\times)^2`` is uniquely determined and called the *quadratic
discriminant* of ``V`` (or of ``C``), denoted by ``\mathrm{disq}(V)``. If ``n = 0``, we put
``\mathrm{disq}(V) \coloneqq 1(K^\times)^2``.

Note that there are also the usual notions of the center of ``C`` (not. ``Z(C)``) and the discriminant of ``V``
(not. ``\mathrm{disc}(V)``) as a quadratic ``K``-space; see [`discriminant`](@ref discriminant(::Hecke.QuadSpace)).
The four concepts of the discriminant and quadratic discriminant of ``V`` and the center and centroid of ``C``
are connected as follows:

- If the dimension ``n`` is even, then ``\mathrm{disc}(V) = \mathrm{disq}(V)`` and ``C`` is a central simple ``K``-algebra, so ``Z(C) = K``.
  Moreover, the centroid is entirely contained in the even Clifford algebra.
- If the dimension ``n`` is odd, then ``\mathrm{disc}(V) = 2\mathrm{disq}(V)`` and ``Z(C) = \mathcal{Z}(V)``.

```@docs
basis_of_center(::CliffordAlgebra)
basis_of_centroid(::CliffordAlgebra)
quadratic_discriminant(::CliffordAlgebra)
disq(::CliffordAlgebra)
```
