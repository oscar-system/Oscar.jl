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

and it is seen as a quadruple $(L, f, f_a, n)$ where $n$ is the order of $f$ and
$f_a$ is an isometry of the ambient quadratic space of $L$ inducing $f$ on $L$.
Note that $f_a$ might not always be of order $n$.

Given a lattice with isometry $(L, f)$, we provide the following accessors to the
elements of the previously described quadruple:

```@docs
lattice(::ZZLatWithIsom)
isometry(::ZZLatWithIsom)
ambient_isometry(::ZZLatWithIsom)
order_of_isometry(::ZZLatWithIsom)
```

Note that for some computations, it is more convenient to work either with the
isometry of the lattice itself, or with an isometry of the ambient quadratic
space inducing it on the lattice.

## Constructor

For simplicity, we have gathered the main constructors under the same name
`lattice_with_isometry`. The user has then the choice on the parameters
depending on what they attend to do:

```@docs
lattice_with_isometry(::ZZLat, ::QQMatrix)
lattice_with_isometry(::ZZLat)
```

By default, the first constructor will always check whether the entry matrix
defines an isometry of the lattice, or its ambient space. We recommend not to
disable this parameter to avoid any further issues. Both isometries of *finite
order* and *infinite order* are supported.

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
Lf = lattice_with_isometry(L, f)
```

## Attributes and first operations

Given a lattice with isometry $Lf := (L, f)$, one can have access most of the
attributes of $L$ and $f$ by calling the similar function to the pair. For
instance, in order to know the genus of $L$, one can simply call `genus(Lf)`.
Here is a list of what are the current accessible attributes:

```@docs
ambient_space(::ZZLatWithIsom)
basis_matrix(::ZZLatWithIsom)
charpoly(::ZZLatWithIsom)
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
minpoly(::ZZLatWithIsom)
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
Lf = lattice_with_isometry(L, f)
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

## Signatures

We conclude this introduction about standard functionalities for lattices with
isometry by introducing a last invariant for lattices with isometry of
hermitain type $(L, f)$, called the *signatures*. These signatures are
are intrinsequely connected to the local archimedean invariants of the
hermitian structure associated to $(L, f)$ via the trace equivalence.

```@docs
signatures(::ZZLatWithIsom)
```

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
