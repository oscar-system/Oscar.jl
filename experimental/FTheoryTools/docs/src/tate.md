```@meta
CurrentModule = Oscar
```

# Global Tate models

A global Tate model describes a particular form of an elliptic fibration.
We focus on elliptic fibrations over base 3-folds ``B3``. Consider
the weighted projective space ``\mathbb{P}^{2,3,1}`` with coordinates
``x, y, z``. In addition, consider
* ``a_1 \in H^0( B_3, \overline{K}_{B_3} )``,
* ``a_2 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 2} )``,
* ``a_3 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 3} )``,
* ``a_4 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 4} )``,
* ``a_6 \in H^0( B_3, \overline{K}_{B_3}^{\otimes 6} )``.
Then form a ``\mathbb{P}^{2,3,1}``-bundle over ``B3`` such that
* ``x`` transforms as a section of ``2 \overline{K}_{B_3}``,
* ``y`` transforms as a section of ``3 \overline{K}_{B_3}``,
* ``z`` transforms as a section of ``0 \overline{K}_{B_3} = \mathcal{O}_{B_3}``.
In this 5-fold ambient space, a global Tate model is the hypersurface defined
by the vanishing of the Tate polynomial
``P_T = x^3 - y^2 - x y z a_1 + x^2 z^2 a_2 - y z^3 a_3 + x z^4 a_4 + z^6 a_6``.

Crucially, for non-trivial F-theory settings, the elliptic fibration in question must
be singular. In fact, by construction, one usually engineers certain singularities.
For this, vanishing orders of the sections ``a_i`` above need to specified. The
following table---often referred to as the Tate table and taken from
[Wei10](@cite)---summarizes the singularities introduced by certain vanishing orders:

| sing. type | ``\mathrm{ord}(\Delta)`` | singularity | group ``G`` | ``a_1`` | ``a_2`` | ``a_3`` | ``a_4`` | ``a_6`` |
| :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- |
| ``I_0`` | ``0`` |  |  | ``0`` | ``0`` | ``0`` | ``0`` | ``0`` |
| ``I_1`` | ``1`` |  |  | ``0`` | ``0`` | ``1`` | ``1`` | ``1`` |
| ``I_2`` | ``2`` | ``A_1`` | ``SU(2)`` | ``0`` | ``0`` | ``1`` | ``1`` | ``2`` |
| ``I_{2k}^{ns}`` | ``2k`` | ``C_k`` | ``Sp(k)`` | ``0`` | ``0`` | ``k`` | ``k`` | ``2k`` |
| ``I_{2k}^s`` | ``2k`` | ``A_{2k-1}`` | ``SU(2k)`` | ``0`` | ``1`` | ``k`` | ``k`` | ``2k`` |
| ``I_{2k+1}^{ns}`` | ``2k+1`` | | ``Sp(k)`` | ``0`` | ``0`` | ``k+1`` | ``k+1`` | ``2k+1`` |
| ``I_{2k+1}^{s}`` | ``2k+1`` | ``A_{2k}`` | ``SU(2k+1)`` | ``0`` | ``1`` | ``k`` | ``k+1`` | ``2k+1`` |
| ``II`` | ``2`` | | | ``1`` | ``1`` | ``1`` | ``1`` | ``1`` |
| ``III`` | ``3`` | ``A_1`` | ``SU(2)`` | ``1`` | ``1`` | ``1`` | ``1`` | ``2`` |
| ``IV^{ns}`` | ``4`` |  | ``Sp(1)`` | ``1`` | ``1`` | ``1`` | ``2`` | ``2`` |
| ``IV^s`` | ``4`` | ``A_2`` | ``SU(3)`` | ``1`` | ``1`` | ``1`` | ``2`` | ``3`` |
| ``I_0^{*ns}`` | ``6`` | ``G_2`` | ``G_2`` | ``1`` | ``1`` | ``2`` | ``2`` | ``3`` |
| ``I_0^{*ss}`` | ``6`` | ``B_3`` | ``SO(7)`` | ``1`` | ``1`` | ``2`` | ``2`` | ``4`` |
| ``I_0^{*s}`` | ``6`` | ``D_4`` | ``SO(8)`` | ``1`` | ``1`` | ``2`` | ``2`` | ``4`` |
| ``I_1^{*ns}`` | ``7`` | ``B_4`` | ``SO(9)`` | ``1`` | ``1`` | ``2`` | ``3`` | ``4`` |
| ``I_1^{*s}`` | ``7`` | ``D_5`` | ``SO(10)`` | ``1`` | ``1`` | ``2`` | ``3`` | ``5`` |
| ``I_2^{*ns}`` | ``8`` | ``B_5`` | ``SO(11)`` | ``1`` | ``1`` | ``3`` | ``3`` | ``5`` |
| ``I_2^{*s}`` | ``8`` | ``D_6`` | ``SO(12)`` | ``1`` | ``1`` | ``3`` | ``3`` | ``5`` |
| ``I_{2k-3}^{*ns}`` | ``2k+3`` | ``B_{2k}`` | ``SO(4k+1)`` | ``1`` | ``1`` | ``k`` | ``k+1`` | ``2k`` |
| ``I_{2k-3}^{*s}`` | ``2k+3`` | ``D_{2k+1}`` | ``SO(4k+2)`` | ``1`` | ``1`` | ``k`` | ``k+1`` | ``2k+1`` |
| ``I_{2k-2}^{*ns}`` | ``2k+4`` | ``B_{2k+1}`` | ``SO(4k+3)`` | ``1`` | ``1`` | ``k+1`` | ``k+1`` | ``2k+1`` |
| ``I_{2k-2}^{*s}`` | ``2k+4`` | ``D_{2k+2}`` | ``SO(4k+4)`` | ``1`` | ``1`` | ``k+1`` | ``k+1`` | ``2k+1`` |
| ``IV^{*ns}`` | ``8`` | ``F_4`` | ``F_4`` | ``1`` | ``2`` | ``2`` | ``3`` | ``4`` |
| ``IV^{*s}`` | ``8`` | ``E_6`` | ``E_6`` | ``1`` | ``2`` | ``2`` | ``3`` | ``5`` |
| ``III^*`` | ``9`` | ``E_7`` | ``E_7`` | ``1`` | ``2`` | ``3`` | ``3`` | ``5`` |
| ``II^*`` | ``10`` | ``E_8`` | ``E_8`` | ``1`` | ``2`` | ``3`` | ``4`` | ``5`` |
| non-min. | ``12`` |  |  | ``1`` | ``2`` | ``3`` | ``4`` | ``6`` |


## ... over concrete bases

For some F-theory constructions, one fixes a particular
 base space. We support the following constructors:
```@docs
global_tate_model(base::AbstractNormalToricVariety)
global_tate_model_over_projective_space()
global_tate_model(ais::Vector{T}, base::AbstractNormalToricVariety) where {T<:MPolyRingElem{QQFieldElem}}
```

## ... over not fully specified bases

This method constructs a global Tate model over a base space that is not
fully specified. Rather, it assumes that a base space exists such that
the Tate sections ``a_i`` are well-defined so that the Tate model in
question is well-defined.

For many practical applications, one wishes to assume a further factorization
of the Tate sections ``a_i``. This has the advantage that one can engineer
singularity loci or even the singularity type over a specific locus. This is
the backbone of many F-theory constructions.

To this end, this method accepts a polynomial ring whose variables are the sections
used in the desired factorization of the Tate sections ``a_i``. For example, if we
desired a factorization:
* ``a_1 = a_{10} w^0``,
* ``a_2 = a_{21} w^1``,
* ``a_3 = a_{32} w^2``,
* ``a_4 = a_{43} w^3``,
* ``a_6 = a_{65} w^5``,
then the polynomial ring in question is the ring with indeterminates
``a_{10}``, ``a_{21}``, ``a_{32}``, ``a_{43}``, ``a_{65}`` and ``w``.

In theory, one can consider these indeterminates as local coordinate of the base space.
For the computer implementation, this ring will therefore serve as the coordinate
ring of an auxiliary toric base space, namely an affine space with those coordinates.

Such geometries can be constructed with the following constructor:
```@docs
global_tate_model(ais::Vector{T}, auxiliary_base_ring::MPolyRing, d::Int) where {T<:MPolyRingElem{QQFieldElem}}
```


## Standard constructions

Certain Tate models have been studied in the physics literature over and over again. Thereby,
these constructions became famous and some were given special names. We aim to provide
support for such standard constructions. Currently, we have support for the following:
```@docs
su5_tate_model_over_arbitrary_3d_base()
```


## Attributes

For all Tate models, we support the following attributes:
```@docs
tate_section_a1(t::GlobalTateModel)
tate_section_a2(t::GlobalTateModel)
tate_section_a3(t::GlobalTateModel)
tate_section_a4(t::GlobalTateModel)
tate_section_a6(t::GlobalTateModel)
tate_polynomial(t::GlobalTateModel)
calabi_yau_hypersurface(t::GlobalTateModel)
global_weierstrass_model(t::GlobalTateModel)
discriminant(t::GlobalTateModel)
```
In case the Tate model is constructed over a not fully specified base, it is
nonetheless possible to construct an auxiliary base space as well as an
auxiliary ambient space. The (auxiliary) base and ambient space can
be accessed with the following functions:
```@docs
toric_base_space(t::GlobalTateModel)
toric_ambient_space(t::GlobalTateModel)
```
To tell if they are auxiliary or not, one can use `base_fully_specified(t)`,
which returns `true` in case the Tate model is defined over a concrete base and
`false` otherwise.

The user can decide to get an information whenever an auxiliary base space,
auxiliary ambient space or auxiliary hypersurface have been computed.
To this end, one invokes `set_verbosity_level(:GlobalTateModel, 1)`.
More background information is available
[here](http://www.thofma.com/Hecke.jl/dev/features/macros/).


## Properties

```@docs
base_fully_specified(t::GlobalTateModel)
```

## Singular loci

For applications in F-theory, singular elliptic fibrations are key
(cf. [Wei18](@cite) and references therein). The general approach is
to not work with the singular space directly. Rather, one resolves
the singularities in order to obtain a smooth space instead.
Subsequently, one performs computations on this smooth space.

In this sense, knowledge of the singular loci is a first step to
a successful analysis of such a geometry. For this, we provide
the following functionality.
```@docs
singular_loci(t::GlobalTateModel)
```
