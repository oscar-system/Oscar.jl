```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Quadratic spaces with isometry

We call *quadratic space with isometry* any pair $(V, f)$ consisting of a
non-degenerate quadratic space $V$ together with an isometry $f\in O(V)$.
We refer to the section about [Spaces](@ref Spaces2) of the documentation for
new users.

Note that currently, we support only rational quadratic forms, i.e.
quadratic spaces defined over $\mathbb{Q}$.

In Oscar, such a pair is encoded by the type called `QuadSpaceWithIsom`:

```@docs
QuadSpaceWithIsom
```

It is seen as a triple $(V, f, n)$ where $n$ is the order of $f$. We
actually support isometries of finite and infinite order. In the case where
$f$ is of infinite order, then `n = PosInf`. If $V$ has rank 0, then any
isometry $f$ of $V$ is trivial and we set by default `n = -1`.

Given a quadratic space with isometry $(V, f)$, we provide the following
accessors to the elements of the previously described triple:

```@docs
isometry(::QuadSpaceWithIsom)
order_of_isometry(::QuadSpaceWithIsom)
space(::QuadSpaceWithIsom)
```

The main purpose of the definition of such objects is to define a
contextual ambient space for quadratic lattices endowed with an isometry.
Indeed, as we will see in the next section, *lattices with isometry* are
attached to an ambient quadratic space with an isometry inducing the one on the
lattice.

## Constructors

For simplicity, we have gathered the main constructors for objects of type
`QuadSpaceWithIsom` under the same name `quadratic_space_with_isometry`. The
user has then the choice on the parameters depending on what they intend to do:

```@docs
quadratic_space_with_isometry(::Hecke.QuadSpace, ::QQMatrix)
quadratic_space_with_isometry(::Hecke.QuadSpace)
```

By default, the first constructor always checks whether the matrix defines
an isometry of the quadratic space. We recommend not to disable this parameter
to avoid any complications. Note however that in the rank 0 case, the checks are
avoided since all isometries are necessarily trivial.

## Attributes and first operations

Given a quadratic space with isometry $Vf := (V, f)$, one has access to
most of the attributes of $V$ and $f$ by calling the similar functions on the
pair $(V, f)$ itself. For instance, in order to know the rank of $V$, one can
simply call `rank(Vf)`. Here is a list of what are the current accessible
attributes:

```@docs
characteristic_polynomial(::QuadSpaceWithIsom)
det(::QuadSpaceWithIsom)
diagonal(::QuadSpaceWithIsom)
dim(::QuadSpaceWithIsom)
discriminant(::QuadSpaceWithIsom)
gram_matrix(::QuadSpaceWithIsom)
is_definite(::QuadSpaceWithIsom)
is_positive_definite(::QuadSpaceWithIsom)
is_negative_definite(::QuadSpaceWithIsom)
minimal_polynomial(::QuadSpaceWithIsom)
rank(::QuadSpaceWithIsom)
signature_tuple(::QuadSpaceWithIsom)
```

Similarly, some basic operations on quadratic spaces and matrices are available
for quadratic spaces with isometry.

```@docs
Base.:^(::QuadSpaceWithIsom, ::Int)
biproduct(::Vector{QuadSpaceWithIsom})
direct_product(::Vector{QuadSpaceWithIsom})
direct_sum(::Vector{QuadSpaceWithIsom})
rescale(::QuadSpaceWithIsom, ::RationalUnion)
```

## Spinor norm

Given a rational quadratic space $(V, \Phi)$, and given an integer $b\in\mathbb{Q}$,
we define the *rational spinor norm* $\sigma$ on $(V, b\Phi)$ to be the group
homomorphism

$\sigma\colon O(V, b\Phi) = O(V, \Phi)\to \mathbb{Q}^\ast/(\mathbb{Q}^\ast)^2$

defined as follows. For $f\in O(V, b\Phi)$, there exist elements $v_1,\ldots,
v_r\in V$ where $1\leq r\leq \text{rank}(V)$ such that $f =
\tau_{v_1}\circ\cdots\circ \tau_{v_r}$ is equal to the product of the associated
reflections. We define

$\sigma(f) := (-\frac{b\Phi(v_1, v_1)}{2})\cdots(-\frac{b\Phi(v_r,v_r)}{2}) \mod (\mathbb{Q}^{\ast})^2.$

```@docs
rational_spinor_norm(::QuadSpaceWithIsom)
```

## Equality

We choose as a convention that two pairs $(V, f)$ and $(V', f')$ of quadratic
spaces with isometries are *equal* if $V$ and $V'$ are the same space, and $f$
and $f'$ are represented by the same matrix with respect to the standard basis
of $V = V'$.
