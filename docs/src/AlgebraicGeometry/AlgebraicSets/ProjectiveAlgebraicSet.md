# Projective Algebraic Sets
For finitely many homogeneous polynomials $f_1,\dots f_r \in k[x_0,\dots x_n]$,
and $I=(f_1,\dots , f_n) \leq k[x_0,\dots x_n]$ the homogeneous ideal they generate, we denote by $X = V(I) \subseteq \mathbb{P}^n$ the
projective algebraic set defined by $I$ and call $k$ its base field.

Let $\mathbb{P}^n(k)=(k^{n+1}\setminus\{0\})/k^*$ be the set of $k$-points of projective space of dimension $n$.
If $k \subseteq K$ is any field extension, we denote the set of $K$-points of $X$ by

$$\begin{aligned}X(K) &= \{ P \in \mathbb{P}^n(K) \mid f_1(P)=\dots = f_n(P)=0\}\\
&=\{P \in \mathbb{P}^n(K) \mid \forall f\in I : f(P)=0\}.\end{aligned}$$

Most properties of the projective variety $X$ refer to $X(K)$ where $K$ is an
algebraically closed field.
Just like for affine schemes there are a few exceptions to this rule,
for instance, whether $X$ is irreducible or not depends on its base field.
See [`is_irreducible(X::AbsProjectiveScheme)`](@ref) for details.
Further exceptions are documented in the individual methods.

## Relation to schemes

One can view a projective algebraic set as a scheme.
See [Projective schemes](@ref).

More formally we define a projective algebraic set as follows:
```@docs
AbsProjectiveAlgebraicSet
```

## Constructors
Projective algebraic sets can be created from homogeneous polynomials and
homogeneous ideals in standard graded rings.
```@docs
algebraic_set(I::MPolyIdeal{<:MPolyDecRingElem}; check::Bool=true)
algebraic_set(p::MPolyDecRingElem; check::Bool=true)
```
Algebraic sets can also be constructed from projective schemes.
```@docs
algebraic_set(X::AbsProjectiveScheme; check::Bool=true)
```

```@docs
set_theoretic_intersection(X::AbsProjectiveAlgebraicSet, Y::AbsProjectiveAlgebraicSet)
irreducible_components(X::AbsProjectiveAlgebraicSet)
geometric_irreducible_components(X::AbsProjectiveAlgebraicSet)
```

## Attributes
In addition to the attributes inherited from [Projective schemes](@ref)
the following are available.
```@docs
vanishing_ideal(X::AbsProjectiveAlgebraicSet)
fat_ideal(X::AbsProjectiveAlgebraicSet)
```

## Methods
Inherited from [Projective schemes](@ref)
## Properties
Inherited from [Projective schemes](@ref)
