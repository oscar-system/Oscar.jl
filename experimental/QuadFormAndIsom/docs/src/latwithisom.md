# Lattices with isometry

We call *lattice with isometry* any pair $(L, f)$ consisting of an integer
lattice $L$ together with an isometry $f \in O(L)$. We refer to the section
about [Integer Lattices](@ref) of the documentation for new users.

In Oscar, such a pair is encoded in the type called `ZZLatWithIsom`:

```@docs
ZZLatWithIsom
```

It is seen as a quadruple $(Vf, L, f, n)$ where $Vf = (V, f_a)$ consists of
the ambient rational quadratic space $V$ of $L$, and an isometry $f_a$ of $V$
preserving $L$ and inducing $f$ on $L$. The integer $n$ is the order of $f$,
which is a divisor of the order of the isometry $f_a\in O(V)$.

Given a lattice with isometry $(L, f)$, we provide the following accessors to
the elements of the previously described quadruple:

```@docs
ambient_isometry(::ZZLatWithIsom)
ambient_space(::ZZLatWithIsom)
isometry(::ZZLatWithIsom)
lattice(::ZZLatWithIsom)
order_of_isometry(::ZZLatWithIsom)
```

Note that for some computations, it is more convenient to work either with the
isometry of the lattice itself, or with the fixed isometry of the ambient
quadratic space inducing it on the lattice.

## Constructors

We provide two ways to construct a pair $Lf = (L,f)$ consisting of an integer
lattice endowed with an isometry. One way to construct an object of type
`ZZLatWithIsom` is through the methods `integer_lattice_with_isometry`. These
two methods do not require as input an ambient quadratic space with isometry.

```@docs
integer_lattice_with_isometry(::ZZLat, ::QQMatrix)
integer_lattice_with_isometry(::ZZLat)
```

By default, the first constructor will always check whether the matrix
defines an isometry of the lattice, or its ambient space. We recommend not to
disable this parameter to avoid any further issues. Note that as in the case of
quadratic spaces with isometry, both isometries of integer lattices of *finite
order* and *infinite order* are supported.

Another way of constructing such lattices with isometry is by fixing an ambient
quadratic space with isometry, of type `QuadSpaceWithIsom`, and specifying a
basis for an integral lattice in that space. If this lattice is preserved by
the fixed isometry of the quadratic space considered, then we endow it with the
induced action.

```@docs
lattice(::QuadSpaceWithIsom)
lattice(::QuadSpaceWithIsom, ::MatElem{ <:RationalUnion})
lattice_in_same_ambient_space(::ZZLatWithIsom, ::MatElem)
```

## Attributes and first operations

Given a lattice with isometry $Lf := (L, f)$, one can have access most of the
attributes of $L$ and $f$ by calling the similar function for the pair. For
instance, in order to know the genus of $L$, one can simply call `genus(Lf)`.
Here is a list of what are the current accessible attributes:

```@docs
basis_matrix(::ZZLatWithIsom)
characteristic_polynomial(::ZZLatWithIsom)
degree(::ZZLatWithIsom)
det(::ZZLatWithIsom)
discriminant(::ZZLatWithIsom)
genus(::ZZLatWithIsom)
gram_matrix(::ZZLatWithIsom)
is_definite(::ZZLatWithIsom)
is_even(::ZZLatWithIsom)
is_elementary(::ZZLatWithIsom, ::IntegerUnion)
is_elementary_with_prime(::ZZLatWithIsom)
is_integral(::ZZLatWithIsom)
is_positive_definite(::ZZLatWithIsom)
is_primary(::ZZLatWithIsom, ::IntegerUnion)
is_primary_with_prime(::ZZLatWithIsom)
is_negative_definite(::ZZLatWithIsom)
is_unimodular(::ZZLatWithIsom)
minimum(::ZZLatWithIsom)
minimal_polynomial(::ZZLatWithIsom)
norm(::ZZLatWithIsom)
rank(::ZZLatWithIsom)
rational_span(::ZZLatWithIsom)
scale(::ZZLatWithIsom)
signature_tuple(::ZZLatWithIsom)
```

Similarly, some basic operations on $\mathbb Z$-lattices and matrices are
available for integer lattices with isometry.

```@docs
Base.:^(::ZZLatWithIsom, ::Int)
direct_sum(::Vector{ZZLatWithIsom})
dual(::ZZLatWithIsom)
lll(::ZZLatWithIsom)
orthogonal_submodule(::ZZLatWithIsom, ::QQMatrix)
rescale(::ZZLatWithIsom, ::RationalUnion)
```

## Type for finite order isometries

Given a lattice with isometry $Lf := (L, f)$ where $f$ is of finite order $n$,
one can compute the *type* of $Lf$.

```@docs
type(::ZZLatWithIsom)
```

Since determining whether two pairs of lattices with isometry are isomorphic is
a challenging task, one can perform a coarser comparison by looking at the
type. This set of data keeps track of some local and global invariants of the
pair $(L, f)$ with respect to the action of $f$ on $L$.

```@docs
is_of_type(::ZZLatWithIsom, t::Dict)
is_of_same_type(::ZZLatWithIsom, ::ZZLatWithIsom)
```

Finally, if the minimal polynomial of $f$ is irreducible, then we say that the
pair $(L, f)$ is of *hermitian type*. The type of a lattice with isometry of
hermitian type is called *hermitian* (note that the type is only defined for
finite order isometries).

These names follow from the fact that, by the trace equivalence, one can
associate to the pair $(L, f)$ a hermitian lattice over the equation order of
$f$, if it is maximal in the associated number field $\mathbb{Q}[f]$.

```@docs
is_of_hermitian_type(::ZZLatWithIsom)
is_hermitian(::Dict)
```

## Hermitian structures and trace equivalence

As mentioned in the previous section, to a lattice with isometry $Lf := (L, f)$
such that the minimal polynomial of $f$ is irreducible, one can associate a
hermitian lattice $\mathfrak{L}$ over the equation order of $f$, if it is
maximal, for which $Lf$ is the associated trace lattice. Hecke provides the
tools to perform the trace equivalence for lattices with isometry of hermitian
type.

```@docs
hermitian_structure(::ZZLatWithIsom)
```

## Discriminant groups

Given an integral lattice with isometry $Lf := (L, f)$, if one denotes by $D_L$
the discriminant group of $L$, there exists a natural map
$\pi\colon O(L) \to O(D_L)$ sending any isometry to its induced action on the
discriminant group of $L$. In general, this map is neither injective nor
surjective. If we denote by $D_f := \pi(f)$ then $\pi$ induces a map between
centralizers $O(L, f)\to O(D_L, D_f)$. Again, this induced map is in general
neither injective nor surjective, and we denote its image by $G_{L,f}$.

```@docs
discriminant_group(::ZZLatWithIsom)
```

For simple cases as for definite lattices, $f$ being plus-or-minus the identity
or if the rank of $L$ is equal to the Euler totient of the order of $f$ (in the
finite case), $G_{L,f}$ can be easily computed. For the remaining cases, we use
the hermitian version of *Miranda--Morrison theory* as presented in
[BH23](@cite). The general computation of $G_{L, f}$ has been implemented in
this project and it can be indirectly used through the general following
method:

```@docs
image_centralizer_in_Oq(::ZZLatWithIsom)
```

Important note: hermitian Miranda-Morrison is only available for even lattices.

For an implementation of the regular Miranda-Morrison theory, we refer to the
function `image_in_Oq` which actually computes the image of
$\pi$ in both the definite and the indefinite case.

More generally, for a finitely generated subgroup $G$ of $O(L)$, we have
implemented a function which computes the representation of $G$ on $D_L$:

```@docs
discriminant_representation(::ZZLat, ::MatrixGroup)
```

We will see later in the section about enumeration of lattices with isometry
that one can compute $G_{L,f}$ in some particular cases arising from
equivariant primitive embeddings of lattices with isometries.

## Kernel sublattices

As for single integer lattices, it is possible to compute kernel sublattices of
some $\mathbb{Z}$-module homomorphisms. We provide here the possibility to
compute $\ker(p(f))$ as a sublattice of $L$ equipped with the induced action of
$f$, where $p$ is a polynomial with rational coefficients.

```@docs
kernel_lattice(::ZZLatWithIsom, ::Union{ZZPolyRingElem, QQPolyRingElem})
kernel_lattice(::ZZLatWithIsom, ::Integer)
```

Note that such sublattices are by definition primitive in $L$ since $L$ is
non-degenerate. As particular kernel sublattices of $L$, one can also compute
the so-called *invariant* and *coinvariant* lattices of $(L, f)$:

```@docs
coinvariant_lattice(::ZZLatWithIsom)
invariant_lattice(::ZZLatWithIsom)
invariant_coinvariant_pair(::ZZLatWithIsom)
```

Similarly, we provide the possibility to compute invariant and coinvariant
sublattices given an orthogonal representation `G` in matrix form of a finite
group on a given lattice `L`:

```@docs
coinvariant_lattice(::ZZLat, ::MatrixGroup)
invariant_lattice(::ZZLat, ::MatrixGroup)
invariant_coinvariant_pair(::ZZLat, ::MatrixGroup)
```

## Signatures

We conclude this introduction about standard functionalities for lattices with
isometry by introducing a last invariant for lattices with finite isometry of
hermitian type $(L, f)$, called the *signatures*. These signatures are
intrinsically connected to the local archimedean invariants of the
hermitian structure associated to $(L, f)$ via the trace equivalence.

```@docs
signatures(::ZZLatWithIsom)
```

## Spinor norm

Given an integer lattice with isometry $(L, f)$, one often would like to know
the *spinor norm* of $f$ seen as an isometry of the rational quadratic space
$L\times \mathbb{Q}$. See [`rational_spinor_norm(::QuadSpaceWithIsom)`](@ref)
for a definition.

```@docs
rational_spinor_norm(::ZZLatWithIsom)
```

## Equality

We choose as a convention that two pairs $(L, f)$ and $(L', f')$ of integer
lattices with isometries are *equal* if their ambient quadratic space with
isometry of type [`QuadSpaceWithIsom`](@ref) are equal, and if the underlying
lattices $L$ and $L'$ are equal as $\mathbb Z$-modules in the common ambient
quadratic space.
