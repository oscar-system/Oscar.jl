```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["Weierstrass.md"]
```

# Global Weierstrass models

A global Weierstrass model describes a particular form of an elliptic fibration.
We focus on elliptic fibrations over base 3-folds ``B3``. Consider
the weighted projective space ``\mathbb{P}^{2,3,1}`` with coordinates
``x, y, z``. In addition, consider
* ``f \in H^0( B_3, \overline{K}_{B_3}^{\otimes 4} )``,
* ``g \in H^0( B_3, \overline{K}_{B_3}^{\otimes 6} )``,
Then form a ``\mathbb{P}^{2,3,1}``-bundle over ``B3`` such that
* ``x`` transforms as a section of ``2 \overline{K}_{B_3}``,
* ``y`` transforms as a section of ``3 \overline{K}_{B_3}``,
* ``z`` transforms as a section of ``0 \overline{K}_{B_3} = \mathcal{O}_{B_3}``.
In this 5-fold ambient space, a global Weierstrass model is the hypersurface defined
by the vanishing of the Weierstrass polynomial ``P_W = x^3 - y^2 + f x z^4 + g z^6``.

Crucially, for non-trivial F-theory settings, the elliptic fibration in question must
be singular. In fact, by construction, one usually engineers certain singularities.
This can be read-off from the Weierstrass table, which we have reproduced from
[Wei18](@cite):

| type | ``\mathrm{ord}(f)`` | ``\mathrm{ord}(g)`` | ``\mathrm{ord}(\Delta)`` | sing. | monodromy cover | algebra ``\mathfrak{g}`` | comp. |
| :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- |
| ``I_0`` | ``\geq 0`` | ``\geq 0`` | 0 | 
| ``I_1`` | ``0`` | ``0`` | ``1`` |
| ``II`` | ``\geq 1`` | ``1`` | ``2`` |
| ``III`` | ``1`` | ``\geq 2`` | ``3`` | ``A_1`` | | ``\mathfrak{su}(2)`` |
| ``IV^{ns}`` | ``\geq 2`` | ``2`` | ``4`` | ``A_2`` | ``\left. \psi^2 - \frac{g}{w^2} \right\|_{w = 0}`` | ``\mathfrak{sp}(1)`` | ``1`` |
| ``IV^{s}`` | ``\geq 2`` | ``2`` | ``4`` | ``A_2`` | ``\left. \psi^2 - \frac{g}{w^2} \right\|_{w = 0}`` | ``\mathfrak{su}(3)`` | ``2`` |
| ``I_m^{ns}`` | ``0`` | ``0`` | ``m`` | ``A_{m - 1}`` | ``\left. \psi^2 + \frac{9g}{2f} \right\|_{w = 0}`` | ``\mathfrak{sp}(\lfloor \frac{m}{2} \rfloor])`` | ``1`` |
| ``I_m^{s}`` | ``0`` | ``0`` | ``m`` | ``A_{m - 1}`` | ``\left. \psi^2 + \frac{9g}{2f} \right\|_{w = 0}`` | ``\mathfrak{su}(m)`` | ``2`` |
| ``I_0^{*ns}`` | ``\geq 2`` | ``\geq 3`` | ``6`` | ``D_4`` | ``\left. \psi^3 + \psi \cdot \frac{f}{w^2} + \frac{g}{w^3} \right\|_{w = 0}`` | ``\mathfrak{g}_2`` | ``1`` |
| ``I_0^{*ss}`` | ``\geq 2`` | ``\geq 3`` | ``6`` | ``D_4`` | ``\left. \psi^3 + \psi \cdot \frac{f}{w^2} + \frac{g}{w^3} \right\|_{w = 0}`` | ``\mathfrak{so}(7)`` | ``2`` |
| ``I_0^{*s}`` | ``\geq 2`` | ``\geq 3`` | ``6`` | ``D_4`` | ``\left. \psi^3 + \psi \cdot \frac{f}{w^2} + \frac{g}{w^3} \right\|_{w = 0}`` | ``\mathfrak{so}(8)`` | ``3`` |
| ``I_{2n-5}^{*ns}`` (``n \geq 3``) | ``2`` | ``3`` | ``2n+1`` | ``D_{2n-1}`` | ``\left. \psi^2 + \frac{1}{4} \left( \frac{\Delta}{w^{2n+1}} \right) \left( \frac{2wf}{9g} \right)^3 \right\|_{w = 0}`` | ``\mathfrak{so}(4n-3)`` | ``1`` |
| ``I_{2n-5}^{*s}`` (``n \geq 3``) | ``2`` | ``3`` | ``2n+1`` | ``D_{2n-1}`` | ``\left. \psi^2 + \frac{1}{4} \left( \frac{\Delta}{w^{2n+1}} \right) \left( \frac{2wf}{9g} \right)^3 \right\|_{w = 0}`` | ``\mathfrak{so}(4n-2)`` | ``2`` |
| ``I_{2n-4}^{*ns}`` (``n \geq 3``) | ``2`` | ``3`` | ``2n+2`` | ``D_{2n}`` | ``\left. \psi^2 + \left( \frac{\Delta}{w^{2n+2}} \right) \left( \frac{2wf}{9g} \right)^2 \right\|_{w = 0}`` | ``\mathfrak{so}(4n-1)`` | 1 |
| ``I_{2n-4}^{*s}`` (``n \geq 3``) | ``2`` | ``3`` | ``2n+2`` | ``D_{2n}`` | ``\left. \psi^2 + \left( \frac{\Delta}{w^{2n+2}} \right) \left( \frac{2wf}{9g} \right)^2 \right\|_{w = 0}`` | ``\mathfrak{so}(4n)`` | 2 |
| ``IV^{*ns}`` | ``\geq 3`` | ``4`` | ``8`` | ``E_6`` | ``\left. \psi^2 - \frac{g}{w^4} \right\|_{w = 0}`` | ``\mathfrak{f}_4`` | ``1`` |
| ``IV^{*s}`` | ``\geq 3`` | ``4`` | ``8`` | ``E_6`` | ``\left. \psi^2 - \frac{g}{w^4} \right\|_{w = 0}`` | ``\mathfrak{e}_6`` | ``2`` |
| ``III^*`` | ``3`` | ``\geq 5`` | ``9`` | ``E_7`` | | ``\mathfrak{e}_7`` | |
| ``II^*`` | ``\geq 4`` | ``5`` | ``10`` | ``E_8`` | | ``\mathfrak{e}_8`` | |
| non-min. | ``\geq 4`` | ``\geq 6`` | ``\geq 12`` | non-can. | | | |


## ... over concrete bases

To construct a Weierstrass model as a hypersurface in an ambient space,
we first need to construct the ambient space in question. For a toric base,
one way to achieve this is to first focus on the Cox ring of the toric
ambient space. This ring must be graded such that the Weierstrass polynomial
is homogeneous and cuts out a Calabi-Yau hypersurface. Given this grading,
one can perform a triangulation. This triangulation will typically
take a long time to complete and yield a large number of candidate ambient
spaces. Typically, one wishes to focus on those spaces which contain the
toric base space in a manifest way. But even this criterion usually
allows for many ambient spaces to be used.

To circumvent this obstacle, all our constructors operate in the opposite direction.
That is, they begins by extracting the rays and maximal cones of the chosen
toric base space. Subsequently, those rays and cones are extended
to form one of the many toric ambient spaces. Note that the so-constructed
ambient space need not be smooth.

We support the following constructors:
```@docs
global_weierstrass_model(base::AbstractNormalToricVariety)
global_weierstrass_model_over_projective_space()
global_weierstrass_model(f::MPolyRingElem{QQFieldElem}, g::MPolyRingElem{QQFieldElem}, base::AbstractNormalToricVariety)
```

## ... over not fully specified bases

A global Weierstrass model can also be constructed over a base space that
is not fully specified. Rather, it assumes that a base space exists such that
the Weierstrass sections ``f`` and ``g`` are well-defined, so that the
Weierstrass model in question is well-defined.

For many practical applications, one wishes to assume a further specialize
the Weierstrass sections ``f`` and ``g``. This has the advantage that one can
engineer singularity loci or even the singularity type over a specific locus.
To some extend, this is the backbone of many F-theory constructions.

Consequently, the construction of such models accepts a polynomial ring whose
variables are the sections used in the desired factorization of the Weierstrass
sections ``f`` and ``g``.

In theory, one can consider the indeterminates of this auxiliary ring as a
(possilby redundant) set of local coordinate of the base space. For the computer
implementation, this ring will therefore serve as the coordinate ring of an auxiliary
toric base space, namely an affine space with those coordinates.

For such geometries, we support the following functionality.
```@docs
global_weierstrass_model(poly_f::MPolyRingElem{QQFieldElem}, poly_g::MPolyRingElem{QQFieldElem}, auxiliary_base_ring::MPolyRing, d::Int)
```


## Attributes

```@docs
weierstrass_section_f(w::GlobalWeierstrassModel)
weierstrass_section_g(w::GlobalWeierstrassModel)
weierstrass_polynomial(w::GlobalWeierstrassModel)
cy_hypersurface(w::GlobalWeierstrassModel)
discriminant(w::GlobalWeierstrassModel)
```
In case the Weierstrass model is constructed over a not fully specified base,
it is nonetheless possibly to construct an auxiliary base space as well as an
auxiliary ambient space. The (auxiliary) base and ambient space can
be accessed with the following functions:
```@docs
toric_base_space(w::GlobalWeierstrassModel)
toric_ambient_space(t::GlobalWeierstrassModel)
```
To tell if they are auxiliary or not, one can use `base_fully_specified(w)`,
which returns `true` in case the Tate model is defined over a concrete base and
`false` otherwise.

The user can decide to get an information whenever an auxiliary base space,
auxiliary ambient space or auxiliary hypersurface have been computed.
To this end, one invokes `set_verbosity_level(:GlobalWeierstrassModel, 1)`.
More background information is available
[here](http://www.thofma.com/Hecke.jl/dev/features/macros/).


## Properties

```@docs
base_fully_specified(t::GlobalWeierstrassModel)
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
singular_loci(t::GlobalWeierstrassModel)
```
