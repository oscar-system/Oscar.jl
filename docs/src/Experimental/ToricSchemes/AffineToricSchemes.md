```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["AffineToricSchemes.md"]
```

# Affine Toric Schemes

## Constructors

We can construct an affine toric scheme as follows:

```julia
julia> C = positive_hull([1 0; 0 1])
A polyhedral cone in ambient dimension 2

julia> antv = AffineNormalToricVariety(C)
A normal, affine toric variety

julia> affine_toric_scheme = ToricSpec(antv)
Spec of an affine toric variety with cone spanned by RayVector{fmpq}[[1, 0], [0, 1]]
```


## Properties

We currently support the following properties:
```@docs
is_simplicial(X::ToricSpec)
is_smooth(X::ToricSpec)
```


## Attributes

We currently support the following attributes:
```@docs
underlying_scheme(X::ToricSpec)
affine_normal_toric_variety(X::ToricSpec)
cone(X::ToricSpec)
dual_cone(X::ToricSpec)
hilbert_basis(X::ToricSpec)
```


## A note on the torus inclusion and the torus action

Note that the `torus_inclusions(X::ToricSpec)` is not yet supported.
For an affine toric scheme ``X``, we envision that this function should
return a list `l` containing the inclusions ``Tʳ⁽ⁱ⁾ ↪ X`` of the different tori. 

Similarly, `torus_action(X::ToricSpec)` is not yet supported.
For an affine toric scheme ``X`` with a dense open torus ``T``
this method should returns a quintuple of morphisms `(pT, pX, incX, mult)` 
consisting of the following:
 * the projection ``T × X → T`` of the product with the torus ``T`` to ``T``,
 * the projection ``T × X → X``,
 * the inclusion ``X ↪ T × X`` taking ``x`` to ``(1, x)``,
 * the group action ``T × X → X``.
