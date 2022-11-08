```@meta
CurrentModule = FTheoryTools
```

```@contents
Pages = ["Weierstrass.md"]
```

# Global Weierstrass models

A global Weierstrass model describes a particular form of an elliptic fibration.
We focus on elliptic fibrations over base 3-folds ``B3``. Consider
the weighted projective space ``\mathbb{P}^{2,3,1}`` with coordinates
``x, y, z``. In addition, consider
- ``f \in H^0( B_3, \overline{K}_{B_3}^{\otimes 4} )``,
- ``g \in H^0( B_3, \overline{K}_{B_3}^{\otimes 6} )``,
Then form a ``\mathbb{P}^{2,3,1}``-bundle over ``B3`` such that
- ``x`` transforms as a section of ``2 \overline{K}_{B_3}``,
- ``y`` transforms as a section of ``3 \overline{K}_{B_3}``,
- ``z`` transforms as a section of ``0 \overline{K}_{B_3} = \mathcal{O}_{B_3}``.
In this 5-fold ambient space, a global Weierstrass model is the hypersurface defined
by the vanishing of the Weierstrass polynomial ``P_W = x^3 - y^2 + f x z^4 + g z^6``.
For such geometries, we support functionality.

## Constructors

We support the following constructors:
```@docs
GenericGlobalWeierstrassModel(base::Oscar.AbstractNormalToricVariety)
GenericGlobalWeierstrassModelOverProjectiveSpace()
SpecificGlobalWeierstrassModel(f::MPolyElem_dec{fmpq, fmpq_mpoly}, g::MPolyElem_dec{fmpq, fmpq_mpoly}, base::Oscar.AbstractNormalToricVariety)
```

## Attributes

```@docs
poly_f(w::GlobalWeierstrassModel)
poly_g(w::GlobalWeierstrassModel)
pw(w::GlobalWeierstrassModel)
toric_base_space(w::GlobalWeierstrassModel)
toric_ambient_space(w::GlobalWeierstrassModel)
```
