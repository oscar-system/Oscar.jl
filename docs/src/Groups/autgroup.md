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
Pages = ["autgroup.md"]
```

# Groups of automorphisms

Groups of automorphisms over a group `G` have parametric type `AutomorphismGroup{T}`, where `T` is the type of `G`. The group of automorphisms over a group `G` is defined by the following instruction:
```julia
AutomorphismGroup{T}
automorphism_group(G)
```
The evaluation of the automorphism `f` in the element `x` is analogous to the homomorphism evaluation: it can be obtained by typing either `f(x)` or `x^f`.

It is possible to turn an automorphism `f` into a homomorphism by typing `hom(f)`. The converse is also possible: if `g` is a bijective homomorphism from the group `G` to itself and `A` is the automorphism group of `G`, then the instruction `A(g)` returns `g` as automorphism of `G`. This is the standard way to build explicitly an automorphism (another way, available for inner automorphisms, is shown in Section [Inner_automorphisms](@ref inner_automorphisms)).

```@docs
automorphism_group(G::GAPGroup)
```

  **Example:**
```jldoctest
julia> S=symmetric_group(4);

julia> A=automorphism_group(S);

julia> g=hom(S,S,x->x^S[1]);

julia> g in A
false

julia> au=A(g);

julia> au in A
true

julia> g==hom(au)
true

julia> x=cperm(S,[1,2,3]);

julia> au(x)
(2,3,4)

julia> g(x)==au(x)
true
```
In Oscar it is possible to multiply homomorphisms and automorphisms (whenever it makes sense); in such cases, the output is always a variable of type `GAPGroupHomomorphism{S,T}`.
```@repl oscar
S=symmetric_group(4);
A=automorphism_group(S);
g=hom(S,S,x->x^S[1]);
f=A(g);
typeof(g*f)
```

The following functions are available for automorphisms, some of them similar to the corresponding functions for homomorphisms of groups.
```@docs
is_invariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T<:GAPGroup
restrict_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T, A=automorphism_group(H)) where T <: GAPGroup
induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupHomomorphism)
hom(x::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup
```

## [Inner automorphisms](@id inner_automorphisms)

Oscar disposes of the following functions to handle inner automorphisms of a group.
```@docs
inner_automorphism(g::GAPGroupElem)
is_inner_automorphism(f::GAPGroupHomomorphism)
inner_automorphisms_group(A::AutomorphismGroup{T}) where T <: GAPGroup
```
