```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["grouphom.md"]
```

# Group homomorphisms

In Oscar, a group homomorphism from `G` to `H` is an object of parametric type `GAPGroupHomomorphism{S,T}`, where `S` and `T` are the types of `G` and `H` respectively.

A homomorphism from `G` to `H` can be defined in two ways.

* Writing explicitly the images of the generators of `G`:
```julia
f = hom(G,H,[x1,x2,...],[y1,y2,...])
```
Here, `[x1,x2,...]` must be a generating set for `G` (not necessarily minimal) and `[y1,y2,...]` is a vector of elements of `H` of the same length of `[x1,x2,...]`. This assigns to `f` the value of the group homomorphism sending `x_i` into `y_i`.

An exception is thrown if such a homomorphism does not exist.

* Taking an existing function `g` satisfying the group homomorphism properties:
```julia
f = hom(G,H,g)
```
An exception is thrown if the function `g` does not define a group homomorphism.

  **Example:**
The following procedures define the same homomorphism (conjugation by `x`) in the two ways explained above.
```jldoctest
julia> S=symmetric_group(4);

julia> x=S[1];

julia> f=hom(S,S,gens(S),[S[1]^x,S[2]^x]);

julia> g=hom(S,S,y->y^x);

julia> f==g
true
```

```@docs
hom(G::GAPGroup, H::GAPGroup, img::Function)
hom(G::GAPGroup, H::GAPGroup, gensG::Vector, imgs::Vector)
image(f::GAPGroupHomomorphism, x::GAPGroupElem)
preimage(f::GAPGroupHomomorphism, x::GAPGroupElem)
restrict_homomorphism(f::GAPGroupHomomorphism, H::GAPGroup)
```

Oscar has also the following standard homomorphism.
```@docs
id_hom
trivial_morphism
```

To evaluate the homomorphism `f` in the element `x` of `G`, it is possible to use the instruction
```julia
image(f,x)
```
or the more compact notations `f(x)` and `x^f`.

  **Example:**
```jldoctest
julia> S=symmetric_group(4);

julia> f=hom(S,S,x->x^S[1]);

julia> x=cperm(S,[1,2]);

julia> image(f,x)
(2,3)

julia> f(x)
(2,3)

julia> x^f
(2,3)
```

A sort of "inverse" of the evaluation is the following
```@docs
haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem; check::Bool = true)
```
  **Example:**
```jldoctest
julia> S=symmetric_group(4);

julia> f=hom(S,S,x->x^S[1]);

julia> x=cperm(S,[1,2]);

julia> haspreimage(f,x)
(true, (1,4))
```

!!! warning
    Do not confuse `haspreimage` with the function `has_preimage`, which works on variable of type `GrpGenToGrpGenMor`.

## Operations on homomorphisms

Oscar supports the following operations on homomorphisms.

* `inv(f)` = the inverse of `f`.
  An exception is thrown if `f` is not bijective.
* `f^n` = the homomorphism `f` composed `n` times with itself.
  An exception is thrown if the domain and the codomain of `f` do not coincide
  (unless `n=1`). If `n` is negative, the result is the inverse of `f` composed `n` times with itself.
* `compose(f, g)` = composition of `f` and `g`. This works only if the codomain of `f` coincides with the domain of `g`. Shorter equivalent expressions are `f*g` and `g(f)`.

  **Example:**
```jldoctest
julia> S=symmetric_group(4);

julia> f=hom(S,S,x->x^S[1]);

julia> g=hom(S,S,x->x^S[2]);

julia> f*g==hom(S,S,x->x^(S[1]*S[2]))
true

julia> f==f^-3
true
```

!!! note
    The composition operation `*` has to be read from the right to the left. So, `(f*g)(x)` is equivalent to `g(f(x))`.

## Properties of homomorphisms

Oscar implements the following attributes of homomorphisms,
in addition to the usual `domain` and `codomain`.

```@docs
isinjective(f::GAPGroupHomomorphism)
issurjective(f::GAPGroupHomomorphism)
isbijective(f::GAPGroupHomomorphism)
isinvertible(f::GAPGroupHomomorphism)
isinvariant(f::GAPGroupHomomorphism, H::GAPGroup)
```

## Subgroups described by homomorphisms

The following functions compute subgroups or quotients of either the domain or the codomain. Analogously to the functions described in Sections [Subgroups](@ref subgroups) and [Quotients](@ref quotient), the output consists of a pair (`H`, `g`), where `H` is a subgroup (resp. quotient) and `g` is its embedding (resp. projection) homomorphism.

```@docs
kernel(f::GAPGroupHomomorphism)
image(f::GAPGroupHomomorphism)
image(f::GAPGroupHomomorphism{S, T}, H::S) where S <: GAPGroup where T <: GAPGroup
cokernel(f::GAPGroupHomomorphism)
preimage(f::GAPGroupHomomorphism{S, T}, H::T) where S <: GAPGroup where T <: GAPGroup
```

## Group isomorphisms

```@docs
isisomorphic(G::GAPGroup, H::GAPGroup)
isisomorphic_with_map(G::GAPGroup, H::GAPGroup)
isomorphism(G::GAPGroup, H::GAPGroup)
```

```@docs
isomorphism(::Type{T}, G::GAPGroup) where T <: Union{FPGroup, PcGroup, PermGroup}
isomorphism(::Type{GrpAbFinGen}, G::GAPGroup)
simplified_fp_group(G::FPGroup)
```
