```@meta
CurrentModule = Oscar
```

# Affine Toric Schemes

## Constructors

We provide the following constructors for affine toric schemes:
```@docs
toric_spec(antv::AffineNormalToricVariety)
```


## Attributes

An affine toric scheme has all attributes of a normal toric scheme.
In addition, there are the following special attributes, that we
overload from the corresponding affine toric variety:
* ``cone(X::ToricSpec)``,
* ``dual_cone(X::ToricSpec)``,
* ``hilbert_basis(X::ToricSpec)``,
* ``toric_ideal(R::MPolyRing, X::ToricSpec)``,
* ``toric_ideal(X::ToricSpec)``.


## Properties

For an affine toric scheme, all properties of normal toric
schemes are supported.


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
