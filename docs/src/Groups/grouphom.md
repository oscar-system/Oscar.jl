```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Group homomorphisms

A *group homomorphism* from a group $G$ to a group $H$ is a map $f$ with
domain $G$ and codomain $H$ that respects the group structure,
that is, $(g*h)^f = g^f * h^f$ and $(g^{-1})^f = (g^f)^{-1}$ hold for
all $g, h \in G$.

## Creation of group homomorphisms

A homomorphism from $G$ to $H$ can be defined
by prescribing the images of some generators of $G$
(see [`hom(G::GAPGroup, H::GAPGroup, gensG::Vector, imgs::Vector)`](@ref)
or by prescribing a Julia function that maps the elements as required
(see [`hom(G::GAPGroup, H::GAPGroup, img::Function)`](@ref).

For homomorphisms to permutation groups that are defined by the action
of $G$ on a set, see [`action_homomorphism`](@ref).

For the special cases of an identity mapping between groups
or the mapping whose image is a trivial group,
use [`id_hom(G::GAPGroup)`](@ref)
and [`trivial_morphism(G::GAPGroup, H::GAPGroup = G)`](@ref),
respectively.

The restriction of a homomorphism to a subgroup of its domain is given by
[`restrict_homomorphism(f::GAPGroupHomomorphism, H::GAPGroup)`](@ref),

```@docs
hom(G::GAPGroup, H::GAPGroup, img::Function)
hom(G::GAPGroup, H::GAPGroup, gensG::Vector, imgs::Vector)
id_hom(G::GAPGroup)
trivial_morphism(G::GAPGroup, H::GAPGroup = G)
restrict_homomorphism(f::GAPGroupHomomorphism, H::GAPGroup)
```

## Applying group homomorphisms

Group homomorphisms can be used to compute images and preimages of
elements and subgroups.

For a group homomorphism `f` and a group element `x`,
the image of `x` under `f` can be computed as
[`image(f::GAPGroupHomomorphism, x::GAPGroupElem)`](@ref)
or `f(x)` or `x^f`,
and the preimage of an element `y` in the image can be computed as
[`preimage(f::GAPGroupHomomorphism, x::GAPGroupElem)`](@ref).

```@docs
image(f::GAPGroupHomomorphism, x::GAPGroupElem)
preimage(f::Union{GAPGroupHomomorphism, GAPGroupEmbedding}, x::GAPGroupElem)
has_preimage_with_preimage(f::GAPGroupHomomorphism, x::GAPGroupElem; check::Bool = true)
```

## Operations on homomorphisms

OSCAR supports the following operations on homomorphisms.

* `inv(f)` is the inverse of `f`.
  An exception is thrown if `f` is not bijective.
* `f^n` is the homomorphism `f` composed `n` times with itself.
  An exception is thrown if the domain and the codomain of `f` do not coincide
  (unless `n=1`).
  If `n` is negative, the result is the inverse of `f` composed `n` times
  with itself.
* `compose(f, g)` is the composition of `f` and `g`.
  This works only if the codomain of `f` coincides with the domain of `g`.
  A shorter equivalent expressions is `f*g`.

### Examples
```jldoctest
julia> S = symmetric_group(4)
Symmetric group of degree 4

julia> f = hom(S, S, x->x^S[1])
Group homomorphism
  from symmetric group of degree 4
  to symmetric group of degree 4

julia> g = hom(S, S, x->x^S[2])
Group homomorphism
  from symmetric group of degree 4
  to symmetric group of degree 4

julia> f*g == hom(S, S, x->x^(S[1]*S[2]))
true

julia> f == f^-3
true
```

!!! note
    The composition operation `*` has to be read from the right to the left. So, `(f*g)(x)` is equivalent to `g(f(x))`.

## Properties of homomorphisms

OSCAR implements the following attributes of homomorphisms,
in addition to the usual `domain` and `codomain`.

```@docs
is_injective(f::GAPGroupHomomorphism)
is_surjective(f::GAPGroupHomomorphism)
is_bijective(f::GAPGroupHomomorphism)
is_invertible(f::GAPGroupHomomorphism)
is_invariant(f::GAPGroupHomomorphism, H::GAPGroup)
```

## Subgroups described by homomorphisms

The following functions compute subgroups or quotients of either the domain or the codomain. Analogously to the functions described in Sections [Subgroups](@ref subgroups) and [Quotients](@ref quotient), the output consists of a pair (`H`, `g`), where `H` is a subgroup (resp. quotient) and `g` is its embedding (resp. projection) homomorphism.

```@docs
kernel(f::GAPGroupHomomorphism)
image(f::GAPGroupHomomorphism)
image(f::GAPGroupHomomorphism{<: GAPGroup, <: GAPGroup}, H::GAPGroup)
cokernel(f::GAPGroupHomomorphism)
preimage(f::GAPGroupHomomorphism{<: GAPGroup, <: GAPGroup}, H::GAPGroup)
```

## Group isomorphisms

For all functions that return group isomorphisms,
we have the following rule about the direction of the result.

If two groups are given as inputs then the domain of the returned isomorphism
is the first given group and the codomain is the second.

If one group is given then the domain of the result is this group,
and the codomain is some new group constructed by the function.

```@docs
is_isomorphic(G::GAPGroup, H::GAPGroup)
is_isomorphic_with_map(G::GAPGroup, H::GAPGroup)
isomorphism(G::GAPGroup, H::GAPGroup)
isomorphic_subgroups(H::GAPGroup, G::GAPGroup)
```

```@docs
isomorphism(::Type{T}, G::Group) where T <: Group
isomorphism(::Type{FinGenAbGroup}, G::GAPGroup)
```

## Technicalities

```@docs
GAPGroupHomomorphism
GAPGroupEmbedding
```
