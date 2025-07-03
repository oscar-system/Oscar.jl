```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Weierstrass models

Weierstrass models are central to many constructions in F-theory. Such a model describes a particular form of an elliptic fibration. Our primary focus is on Weierstrass models built over **toric base spaces**. These models receive the most comprehensive support. Weierstrass models over more general bases, such as arbitrary schemes or unspecified base spaces, are under development.

---

## Background: What Is a Weierstrass Model?

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


---

## Constructing Weierstrass Models

Weierstrass models are most robustly supported over **toric base spaces**. While there are plans to extend support to toric schemes, covered schemes, and unspecified base families, these directions remain experimental and are reserved for future development. This section briefly discusses constructing Weierstrass models over arbitrary bases but primarily focuses on working with a **concrete toric base**.

As mentioned, a Weierstrass model is defined as a hypersurface within an ambient space. Thus, constructing the ambient space is the first step. When working with a family of base spaces, this ambient space can be treated symbolically---flexible enough to represent an entire family while still allowing algebraic manipulations. In particular, it is often useful to construct a Weierstrass model **without explicitly specifying the base space**. This is especially common when engineering singularities through specific factorizations of the Weierstrass sections `f` and `g`, as seen in the Weierstrass table. This method is prevalent in F-theory model building, where singularities are engineered independently of the full base geometry.

To facilitate such workflows, we allow users to define the Weierstrass sections `f` and `g` as elements in a multivariate polynomial ring. The variables in this polynomial ring are interpreted as coordinates of an **auxiliary base space**---or as sections of line bundles over it.

The mathematical meaning and implementation details of such symbolic base families are subtle and under active discussion. We do not elaborate on them further here. Rather, let us mention that users can create such models with the following constructor, but support for these models may be limited.

```@docs
weierstrass_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, weierstrass_f::MPolyRingElem, weierstrass_g::MPolyRingElem)
```

Full support exists for constructing Weierstrass models over **concrete, complete toric base spaces**. Completeness is a technical assumption. It ensure that the set of global sections of a line bundle are finite-dimensional vector spaces, enabling OSCAR to handle these sets efficiently. In the future, this restriction may be relaxed. For now, completeness checks---while sometimes slow---are executed in many methods involving Weierstrass models. To gain performance, you can skip these checks using the optional keyword argument:
```julia
completeness_check = false
```

We proceed under the assumption that the base space is a fixed, complete toric variety.

In this setting, constructing a suitable ambient space becomes significantly more efficient. First, we focus on toric ambient spaces. Second, rather than enumerating a large number of such ambient spaces from the triangulations of polytopes---a resource-intensive process that generates many ambient spaces---we take a different, **performance-optimized** route. Starting from the rays and maximal cones of the toric base, we **extend** them to construct a compatible polyhedral fan that defines the toric ambient space. This avoids triangulation entirely and produces a single suitable ambient toric variety. The resulting space may not be smooth, and it may differ from choices made in the F-theory literature, but it is guaranteed to be compatible with the elliptic fibration structure.

Users can construct Weierstrass models over such concrete toric bases with the following constructors:

```@docs
weierstrass_model(base::NormalToricVariety; completeness_check::Bool = true)
weierstrass_model(base::NormalToricVariety, f::MPolyRingElem, g::MPolyRingElem; completeness_check::Bool = true)
```

For convenience ---ideal for quick experiments and educational purposes--- we also support constructors for Weierstrass models over commonly chosen base spaces:

```@docs
weierstrass_model_over_projective_space(d::Int)
weierstrass_model_over_hirzebruch_surface(r::Int)
weierstrass_model_over_del_pezzo_surface(b::Int)
```

A number of Weierstrass models have gained popularity within the F-theory community. Those models are often associated with specific publications, maybe informally referred to by the author names or particular buzz words. For these established constrctions, we provide support through the specialized `literature_model` interface, which is discussed on a dedicated page within this documentation.


---

## Attributes and Properties of Weierstrass Models

Weierstrass, global Tate and hypersurface models are very similar: They all encode an elliptic fibration as a hypersurface in an ambient space. They differ in details, but in brought strokes are similar. As such, they share attributes and properties. The shared properties and attributes are collected on the page [Functionality for all F-theory models](@ref). Among others, they include `base_space`, `ambient_space` and `fiber_ambient_space`.

What follows, is therefore "just" a list of those attributes that are special to Weierstrass models and thus do not apply to global Tate and/or hypersurface models.
```@docs
weierstrass_section_f(w::WeierstrassModel)
weierstrass_section_g(w::WeierstrassModel)
weierstrass_polynomial(w::WeierstrassModel)
weierstrass_ideal_sheaf(w::WeierstrassModel)
calabi_yau_hypersurface(w::WeierstrassModel)
```


---

## Studying and Resolving Singularities of Weierstrass Models

Let us emphasize again that in F-theory, *singular* elliptic fibrations are of central importance (cf. [Wei18](@cite) and references therein): singularities resemble non-trivial physics. One crucial quantity in this regard is the discriminant locus of the fibration, i.e. the locus of the base space over which the elliptic fibers become singular:
```@docs
discriminant(w::WeierstrassModel)
```
Even more critical than the discriminant itself is its decomposition: the discriminant locus may split into several irreducible components. Over each such component, the type of singularity in the elliptic fiber is determined. The following function singular_loci returns this list of irreducible base loci along with the corresponding singularity types:
```@docs
singular_loci(w::WeierstrassModel)
```
!!! warning
    The classification of singularities is performed using a Monte Carlo algorithm, i.e. it involves random sampling. While the implementation has been extensively tested and meets high standards, its probabilistic nature may occasionally yield non-deterministic results.

Sometimes, one might wish to take an existing model and then change some of its parameters, more specifically its tunable sections, in an attempt to create a different singularity structure. We provide support for this through the following method:
```@docs
tune(w::WeierstrassModel, special_section_choices::Dict{String, <:MPolyRingElem}; completeness_check::Bool = true)
```
More details ---in particular what tunable sections are--- are discussed in the exposition of the `tune` function described in [Functionality for all F-theory models](@ref).

While the singularities are critical for non-trivial physics, our understanding of how we can read-off the physics directly from quantities of the singular geometry is limited. Instead, it is standard in the F-theory model building to employ a crepant resolution and then to investigate the resulting resolved space. Of course, many resolution techniques do exist. Among them, we focus on blowups, which can be facilitated by the `blow_up` function. The resulting model will thereafter be partially resolved --- crepant resolutions may not exist/resolve only certain singularities. Currently, no checks are executed to test if a model is smooth following a blowup or a sequence therefore. Instead, the function `is_partially_resolved` will always return `true` if a blowup has been applied to the model in question. More details can be found in [Functionality for all F-theory models](@ref).
