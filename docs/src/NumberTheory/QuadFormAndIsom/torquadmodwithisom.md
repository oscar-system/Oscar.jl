```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Torsion quadratic modules with isometry

We call *torsion quadratic module with isometry* any pair $(T, f)$ consisting
of a torsion quadratic module $T$, of type `TorQuadModule`, together with an
isometry $f\in O(T)$, of type `TorQuadModuleMap`. We refer to the section
about [Torsion Quadratic Modules](@ref) of the manual for new users.

```@docs
TorQuadModuleWithIsom
```

Given a torsion quadratic module with isometry $(T, f)$, we provide the
following accessors:

```@docs
underlying_module(::TorQuadModuleWithIsom)
torsion_quadratic_module(::TorQuadModuleWithIsom)
isometry(::TorQuadModuleWithIsom)
```

Note that the isometry ``f`` is of finite order, but this order is not computed
by default. One can use the following function, which stores this order of
``f`` after being computed once:

```@docs
order_of_isometry(::TorQuadModuleWithIsom)
```

## Constructors

There are two standard ways to describe a torsion quadratic module:
  * As a quotient $L/S$ where $L$ is $\mathbb{Z}$-lattice and $S\subseteq L$
    is a finite index sublattice;
  * As an finite abelian group $A$ equipped with a torsion quadratic/bilinear
    form $q$.

These different viewpoints bring different ways to construct objects of types
`TorQuadModuleWithIsom`, especially when it comes to describing the associated
isometry. We provide the following constructors:

```@docs
torsion_quadratic_module_with_isometry(::TorQuadModule, ::TorQuadModuleMap)
torsion_quadratic_module_with_isometry(::QQMatrix, ::ZZMatrix)
```

Given a fixed torsion quadratic module with isometry $(T, f)$, there isanother
way to construct objects of type `TorQuadModuleWithIsom` as $f$-stable
submodules of $T$.

```@docs
sub(::TorQuadModuleWithIsom, ::Vector{TorQuadModuleElem})
primary_part(::TorQuadModuleWithIsom, ::IntegerUnion)
orthogonal_submodule(::TorQuadModuleWithIsom, ::TorQuadModule)
submodules(::TorQuadModuleWithIsom)
```

## (Anti-)Isomorphism

Given two torsion quadratic module $(T, f)$ and $(S, g)$, we call an abelian
group isomorphism $\psi\colon T\to S$ an *(anti-)isomorphism* between $(T, f)$ and
$(S, g)$ if $\psi$ defines an (anti-)isometry of $T$ and $S$ and if satisfies
$\psi\circ f = g\circ \psi$. We denote by $O(T, f)$ the *automorphism group* of
$(T, f)$, i.e. the group of isomorphisms from $(T, f)$ to $(T, f)$.

```@docs
automorphism_group_with_inclusion(::TorQuadModuleWithIsom)
automorphism_group(::TorQuadModuleWithIsom)
is_isomorphic_with_map(::TorQuadModuleWithIsom, ::TorQuadModuleWithIsom)
is_anti_isomorphic_with_map(::TorQuadModuleWithIsom, ::TorQuadModuleWithIsom)
```

## Equality

We choose as a convention that two pairs $(T, f)$ and $(S, g)$ of torsion
quadratic modules with isometry are *equal* if the underlying modules $T$ and
$S$ of type[`TorQuadModule`](@ref) are the same julia object, and if the
associated isometries have the same matrix representation.
