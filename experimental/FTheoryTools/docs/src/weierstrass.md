```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Weierstrass models

A Weierstrass model describes a particular form of an elliptic fibration.
We focus on an elliptic fibration over a complete base ``B``. Consider
the weighted projective space ``\mathbb{P}^{2,3,1}`` with coordinates
``x, y, z``. In addition, consider
* ``f \in H^0( B, \overline{K}_{B}^{\otimes 4} )``,
* ``g \in H^0( B, \overline{K}_{B}^{\otimes 6} )``,
Then form a ``\mathbb{P}^{2,3,1}``-bundle over ``B`` such that
* ``x`` transforms as a section of ``2 \overline{K}_{B}``,
* ``y`` transforms as a section of ``3 \overline{K}_{B}``,
* ``z`` transforms as a section of ``0 \overline{K}_{B} = \mathcal{O}_{B}``.
In this 5-fold ambient space, a Weierstrass model is the hypersurface defined
by the vanishing of the Weierstrass polynomial ``P_W = x^3 - y^2 + f x z^4 + g z^6``.

Crucially, for non-trivial F-theory settings, the elliptic fibration in question must
be singular. In fact, by construction, one usually engineers certain singularities.
This can be read-off from the Weierstrass table, which we have reproduced from
[Wei18](@cite) with small corrections:

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
| ``I_{2n-4}^{*ns}`` (``n \geq 3``) | ``2`` | ``3`` | ``2n+2`` | ``D_{2n}`` | ``\left. \psi^2 + \left( \frac{\Delta}{w^{2n+2}} \right) \left( \frac{2wf}{9g} \right)^2 \right\|_{w = 0}`` | ``\mathfrak{so}(4n-1)`` | ``1`` |
| ``I_{2n-4}^{*s}`` (``n \geq 3``) | ``2`` | ``3`` | ``2n+2`` | ``D_{2n}`` | ``\left. \psi^2 + \left( \frac{\Delta}{w^{2n+2}} \right) \left( \frac{2wf}{9g} \right)^2 \right\|_{w = 0}`` | ``\mathfrak{so}(4n)`` | ``2`` |
| ``IV^{*ns}`` | ``\geq 3`` | ``4`` | ``8`` | ``E_6`` | ``\left. \psi^2 - \frac{g}{w^4} \right\|_{w = 0}`` | ``\mathfrak{f}_4`` | ``1`` |
| ``IV^{*s}`` | ``\geq 3`` | ``4`` | ``8`` | ``E_6`` | ``\left. \psi^2 - \frac{g}{w^4} \right\|_{w = 0}`` | ``\mathfrak{e}_6`` | ``2`` |
| ``III^*`` | ``3`` | ``\geq 5`` | ``9`` | ``E_7`` | | ``\mathfrak{e}_7`` | |
| ``II^*`` | ``\geq 4`` | ``5`` | ``10`` | ``E_8`` | | ``\mathfrak{e}_8`` | |
| non-min. | ``\geq 4`` | ``\geq 6`` | ``\geq 12`` | non-can. | | | |


## Constructors

We aim to provide support for Weierstrass models over the following bases:
* a toric variety,
* a toric scheme,
* a (covered) scheme.
Often, one also wishes to obtain information about a Weierstrass model without
explicitly specifying the base space. Also for this application, we provide support.
Finally, we provide support for some standard constructions.

Before we detail these constructors, we must comment on the constructors over toric base
spaces. Namely, in order to construct a Weierstrass model as a hypersurface in an ambient
space, we first wish to construct the ambient space in question. For a toric base, one way to
achieve this is to first focus on the Cox ring of the toric ambient space. This ring must be
graded such that the Weierstrass polynomial is homogeneous and cuts out a Calabi-Yau
hypersurface. Given this grading, one can perform a triangulation task. Typically, this
combinatorial task is very demanding, consumes a lot of computational power and takes a
long time to complete. Even more, it will yield a large, often huge, number of candidate
ambient spaces of which the typical user will only pick one. For instance, a common and
often appropriate choice is a toric ambient space which contains the toric base space in a
manifest way.

To circumvent this very demanding computation, our toric constructors operate in the
opposite direction. That is, they begin by extracting the rays and maximal cones of the chosen
toric base space. Subsequently, those rays and cones are extended to form one of the many
toric ambient spaces. This proves hugely superior in performance than going through the
triangulation task of enumerating all possible toric ambient spaces. One downside of this strategy
is that the so-constructed ambient space need not be smooth.

### A toric variety as base space

We require that the provided toric base space is complete. This is a technical limitation as of now.
The functionality of OSCAR only allows us to compute a section basis (or a finite subset thereof)
for complete toric varieties. In the future, this could be extended.

However, completeness is an expensive check. Therefore, we provide an optional argument which
 one can use to disable this check if desired. To this end, one passes the optional argument
 `completeness_check = false` as last argument to the constructor. The following examples
 demonstrate this:
```@docs
weierstrass_model(base::NormalToricVariety; completeness_check::Bool = true)
weierstrass_model(base::NormalToricVariety, f::MPolyRingElem, g::MPolyRingElem; completeness_check::Bool = true)
```


### A (covered) scheme as base space

This functionality does not yet exist.

### Base space not specified

A Weierstrass model can also be constructed over a base space that
is not fully specified. Rather, it assumes that a base space exists such that
the Weierstrass sections ``f`` and ``g`` are well-defined, so that the
Weierstrass model in question is well-defined.

For many practical applications, one wishes to assume a further specialize
the Weierstrass sections ``f`` and ``g``. This has the advantage that one can
engineer singularity loci or even the singularity type over a specific locus.
To some extend, this is the backbone of many F-theory constructions. It is
useful to consider a polynomial ring whose variables are the sections
used in the desired factorization of the Weierstrass sections ``f`` and ``g``.
In theory, one can consider the indeterminates of this polynomial ring as local
coordinate of an auxiliary base space. Indeed, for our computer implementation
the polynomial ring with these indeterminates serves as the coordinate ring of
a family of base spaces. We support the following constructor:
```@docs
weierstrass_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, weierstrass_f::MPolyRingElem, weierstrass_g::MPolyRingElem)
```

### Standard constructions

We provide convenient constructions of Weierstrass models
over famous base spaces. Currently, we support the following:
```@docs
weierstrass_model_over_projective_space(d::Int)
weierstrass_model_over_hirzebruch_surface(r::Int)
weierstrass_model_over_del_pezzo_surface(b::Int)
```


### Literature models

Certain Weierstrass models have been studied in the physics literature over and over again.
Thereby, these constructions became famous and some were given special names. We aim
to provide support for such standard constructions. Currently, we provide support for the following:
```@docs
su5_weierstrass_model_over_arbitrary_3d_base()
```


## Attributes

### Basic attributes

For all Weierstrass models -- irrespective over whether the base is toric or not -- we support
the following attributes:
```@docs
weierstrass_section_f(w::WeierstrassModel)
weierstrass_section_g(w::WeierstrassModel)
weierstrass_polynomial(w::WeierstrassModel)
weierstrass_ideal_sheaf(w::WeierstrassModel)
```
The base space can be obtained with `base_space`, the ambient space with `ambient_space` and the
fiber ambient space with `fiber_ambient_space`. Recall that `is_base_space_fully_specified` will
tell if the model has been constructed over a concrete space (in which case the function returns
`true`) or a family of spaces (returning `false`).


### Advanced attributes

The following attributes are currently only supported in a toric setting:
```@docs
calabi_yau_hypersurface(w::WeierstrassModel)
```
Note that for applications in F-theory, *singular* elliptic fibrations are key
(cf. [Wei18](@cite) and references therein). Consequently the discriminant
locus as well as the singular loci of the fibration in question are of ample
importance:
```@docs
discriminant(w::WeierstrassModel)
singular_loci(w::WeierstrassModel)
```

## Methods

### Blowup

We can blow up a Weierstrass model with the `blow_up` function. The resulting model
will thereafter be partially resolved. No checks are currently implemented to test
if a model is completely resolved. However, `is_partially_resolved` will return `true`
if a blowup has been applied to the model in question.


### Tuning

Often, one wishes to tune an existing model, e.g. in an attempt to engineer a
larger gauge group. We support the following functionality:
```@docs
tune(w::WeierstrassModel, special_section_choices::Dict{String, <:MPolyRingElem}; completeness_check::Bool = true)
```
See also the `tune` function described in [Functionality for all F-theory models](@ref).
