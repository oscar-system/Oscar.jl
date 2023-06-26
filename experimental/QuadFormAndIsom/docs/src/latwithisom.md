```@meta
CurrentModule = Oscar
```

# Lattice with isometry

We call *lattice with isometry* any pair $(L, f)$ consisting of an integer
lattice $L$ together with an isometry $f \in O(L)$. We refer to the section
about integer lattices of the documentation for new users.

On Oscar, such a pair is contained into a type called `ZZLatWithIsom`:

```@docs
ZZLatWithIsom
```

and it is seen as a quadruple $(Vf, L, f, n)$ where $Vf = (V, f_a)$ consists on
the ambient rational quadratic space $V$ of $L$ and an isometry $f_a$ of $V$
preversing $L$ and inducing $f$ on $L$. $n$ is the order of $f$, which is a
divisor of the order of the isometry $f_a\in O(V)$.

Given a lattice with isometry $(L, f)$, we provide the following accessors to the
elements of the previously described quadruple:

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

## Constructor

We provide two ways to construct a pair $Lf = (L,f)$ consisting on an integer
lattice endowed with an isometry. One way to construct an object of type
`ZZLatWithIsom` is through the methods `integer_lattice_with_isometry`. These
two methods does not require as input an ambient quadratic space with isometry.

```@docs
integer_lattice_with_isometry(::ZZLat, ::QQMatrix)
integer_lattice_with_isometry(::ZZLat)
```

By default, the first constructor will always check whether the entry matrix
defines an isometry of the lattice, or its ambient space. We recommend not to
disable this parameter to avoid any further issues. Note that as in the case of
quadratic space with isometries, both isometries of integer lattices of *finite
order* and *infinite order* are supported.

Another way of constructing such lattices with isometry is by fixing an ambient
quadratic space with isometry, of type `QuadSpaceWithIsom`, and specify a basis
for an integral lattice in that space. If this lattice is preserved by the fixed
isometry of the quadratic space considered, then we endow it with the induced
action.

```@docs
lattice(::QuadSpaceWithIsom)
lattice(::QuadSpaceWithIsom, ::MatElem{ <:RationalUnion})
lattice_in_same_ambient_space(::ZZLatWithIsom, ::MatElem)
```
### Examples

```@repl 2
using Oscar # hide
L = root_lattice(:E, 6);
f = matrix(QQ, 6, 6, [ 1  2  3  2  1  1;
                      -1 -2 -2 -2 -1 -1;
                       0  1  0  0  0  0;
                       1  0  0  0  0  0;
                      -1 -1 -1  0  0 -1;
                       0  0  1  1  0  1]);
Lf = integer_lattice_with_isometry(L, f)
```

## Attributes and first operations

Given a lattice with isometry $Lf := (L, f)$, one can have access most of the
attributes of $L$ and $f$ by calling the similar function to the pair. For
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
is_integral(::ZZLatWithIsom)
is_positive_definite(::ZZLatWithIsom)
is_negative_definite(::ZZLatWithIsom)
minimum(::ZZLatWithIsom)
minimal_polynomial(::ZZLatWithIsom)
norm(::ZZLatWithIsom)
rank(::ZZLatWithIsom)
rational_span(::ZZLatWithIsom)
scale(::ZZLatWithIsom)
signature_tuple(::ZZLatWithIsom)
```

Similarly, some basic operations on $\mathbb Z$-lattices are available for
lattices with isometry.

```@docs
biproduct(::Vector{ZZLatWithIsom})
direct_product(::Vector{ZZLatWithIsom})
direct_sum(::Vector{ZZLatWithIsom})
dual(::ZZLatWithIsom)
lll(::ZZLatWithIsom)
rescale(::ZZLatWithIsom, ::RationalUnion)
```

## Type for finite order isometries

Given a lattice with isometry $Lf := (L, f)$ where $f$ is of finite order $n$,
one can compute the *type* of $Lf$, which can be seen as an equivalent of the
*genus* used to classified single lattices.

```@docs
type(::Lf)
```

Since determining whether two pairs of lattices with isometry are isomorphic is
a challenging task, one can perform a coarser comparison by looking at the type.
This set of data keeps track of some local and global invariants of the pair $(L,
f)$ with respect to the action of $f$ on $L$.

```@docs
is_of_type(::ZZLatWithIsom, t:Dict)
is_of_same_type(::ZZLatWithIsom, ::ZZLatWithIsom)
```

### Examples

```@repl 2
using Oscar # hide
L = root_lattice(:E, 6);
f = matrix(QQ, 6, 6, [ 1  2  3  2  1  1;
                      -1 -2 -2 -2 -1 -1;
                       0  1  0  0  0  0;
                       1  0  0  0  0  0;
                      -1 -1 -1  0  0 -1;
                       0  0  1  1  0  1]);
Lf = integer_lattice_with_isometry(L, f)
type(Lf)
```

Finally, if the minimal polynomial of $f$ is cyclotomic, i.e. the $n$-th
cyclotomic polynomial, then we say that the pair $(L, f)$ is *of hermitian
type*. The type of a lattice with isometry of hermitian type is called
*hermitian*.

These namings follow from the fact that, by the trace equivalence, one can
associate to the pair $(L, f)$ a hermitian lattice over the ring of integers of
the $n$-th cyclotomic field

```@docs
is_of_hermitian_type(::ZZLatWithIsom)
is_hermitian()
```

## Hermitian structure and trace equivalence

As mentioned in the previous section, to a lattice with isometry $Lf := (L, f)$
such that the minimal polynomial of $f$ is the $n$-th cyclotomic polynomial, one
can associate a hermitian lattice $\mathfrak{L}$ over the ring of integers of
the $n$-th cyclotomic field for which $Lf$ is the associated trace lattice (see
[`trace_lattice_with_isometry(::AbstractLat)`](@ref)). Hecke provides the tools
to perform the trace equivalence for lattices with isometry of hermitian type.

```@docs
hermitian_structure(::ZZLatWithIsom)
```

## Discriminant group

Given an integral lattice with isometry $Lf := (L, f)$, if one denotes $D_L$ the
discriminant group of $L$, there exists a natural map $\pi\colon O(L) \to O(D_L)$
sending any isometry to its induced action on the discriminant form of $L$. In
general, this map is neither injective nor surjective. If we denote $D_f :=
\pi(f)$ then $\pi$ induces a map between centralizers $O(L, f)\to O(D_L, D_f)$.
Again, this induced map is in general neither injective nor surjective, and we
denote its image $G_{L,f}$.

```@docs
discriminant_group(::ZZLatWithIsom)
```

For simple cases as for definite lattices, $f$ being plus-or-minus the identity
or if the rank of $L$ is equal to the totient of the order of $f$ (in the
finite case), $G_{L,f}$ can be easily computed. The only other case which can
be currently handled is for lattices with isometry of hermitian type following
the *hermitian Miranda-Morisson theory* from [BH23]. This has been implemented
in this project and it can be indirectly used through the general following method:

```@docs
image_centralizer_in_Oq(::ZZLatWithIsom)
```

For an implementation of the regular Miranda-Morisson theory, we refer to the
function [`image_in_Oq(::ZZLat)`](@ref) which actually compute the image of
$\pi$ in both the definite and the indefinite case.

We will see later in the section about enumeration of lattices with isometry
that one can compute $G_{L,f}$ in some particular cases arising from equivariant
primitive embeddings of lattices with isometries.

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
isometry by introducing a last invariant for lattices with isometry of
hermitain type $(L, f)$, called the *signatures*. These signatures are
are intrinsequely connected to the local archimedean invariants of the
hermitian structure associated to $(L, f)$ via the trace equivalence.

```@docs
signatures(::ZZLatWithIsom)
```

## Equality

We choose as a convention that two pairs $(L, f)$ and $(L', f')$ of integer
lattices with isometries are *equal* if their ambient quadratic space with
isometry of type `QuadSpaceWithIsom` are equal, and if the underlying lattices
$L$ and $L'$ are equal as $\mathbb Z$-modules in the common ambient quadratic
space.

## Tips for users

### Report an issue

If you are working with some lattices with isometry, of type `ZZLatWithIsom`, and
you need to report an issue, you can produce directly some lines of codes
helping to reconstruct your "non-working" example. We have implemented a method
`to_oscar` which prints 5 lines for reconstructing your example.

```@repl 2
using Oscar # hide
```

### Make the code more talkative

Within the code, there are more hidden messages and testing which are disabled
by default. If you plan to experiment with the codes with your favourite
examples, you may want to be able to detect some issues to be reported, as well
as knowing what the code is doing. Indeed, some functions might take time in
term of compilation but also computations. For this, you can enable these extra
tests and printings by setting:

```julia
Oscar.set_lwi_level(2)
```
