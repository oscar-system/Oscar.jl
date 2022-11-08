```@meta
CurrentModule = FTheoryTools
```

```@contents
Pages = ["Tate.md"]
```

# Global Tate models

A global Tate model describes a particular form of an elliptic fibration.
We focus on elliptic fibraitons over base 3-folds ``B3``. Consider
the weighted projective space ``\mathbb{P}_{2,3,1}`` with coordinates
``x, y, z``. In addition, consider
- ``a_1 \in H^0( B_3, \overline{K}_{B_3} )``,
- ``a_2 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 2} )``,
- ``a_3 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 3} )``,
- ``a_4 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 4} )``,
- ``a_6 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 6} )``.
Then form a ``\mathbb{P}_{2,3,1}``-bundle over ``B3`` such that
- ``x`` transforms as sections of ``2 \overline{K}_{B_3}``,
- ``y`` transforms as sections of ``3 \overline{K}_{B_3}``,
- ``z`` transforms as sections of ``0 \overline{K}_{B_3} = \mathcal{O}_{B_3}``.
In this 5-fold ambient space, a global Tate model is the hypersurface defined
by the vanishing of the Tate polynomial
``P_T = x^3 - y^2 + x y z a_1 + x^2 z^2 a_2 + y z^3 a_3 + x z^4 a_4 + z^6 a_6``.
For such geometries, we support functionality.


## Constructors

```@docs
GenericGlobalTateModel(base::Oscar.AbstractNormalToricVariety)
GenericGlobalTateModelOverProjectiveSpace()
SpecificGlobalTateModel(ais::Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}, base::Oscar.AbstractNormalToricVariety)
```


## Attributes

```@docs
a1(t::GlobalTateModel)
a2(t::GlobalTateModel)
a3(t::GlobalTateModel)
a4(t::GlobalTateModel)
a6(t::GlobalTateModel)
pt(t::GlobalTateModel)
toric_base_space(t::GlobalTateModel)
toric_ambient_space(t::GlobalTateModel)
```
