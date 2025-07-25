```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Weierstrass Models](@id weierstrass_models)

Weierstrass models are central to many constructions in F-theory. Such a model describes
a particular form of an elliptic fibration. Our primary focus is on Weierstrass models
built over **toric base spaces**. These models receive the most comprehensive support.
Weierstrass models over more general bases, such as arbitrary schemes or unspecified
base spaces, are under development.

---

## What Is a Weierstrass Model?

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

Weierstrass models are most robustly supported over **toric base spaces**. While there are
plans to extend support to toric schemes, covered schemes, and unspecified base families,
these directions remain experimental and are reserved for future development. This section
briefly discusses constructing Weierstrass models over arbitrary bases, but primarily focuses
on working with a **concrete toric base**.

### Unspecified Base Spaces

As mentioned, a Weierstrass model is defined as a hypersurface within an ambient space. Thus,
constructing the ambient space is the first step. When working with a family of base spaces,
this ambient space can be treated symbolically—flexible enough to represent an entire family
while still allowing algebraic manipulations. In particular, it is often useful to construct
a Weierstrass model **without explicitly specifying the base space**. This is especially common
when engineering singularities through specific factorizations of the Weierstrass sections `f`
and `g`, as seen in the Weierstrass table. This method is prevalent in F-theory model building,
where singularities are engineered independently of the full base geometry.

To facilitate such workflows, users can define the Weierstrass sections `f` and `g` as elements
in a multivariate polynomial ring. The variables in this ring are interpreted as sections of
line bundles on the **auxiliary base space**.

The mathematical meaning and implementation details of symbolic base families are subtle and
under active discussion. We do not elaborate on them further here. Let us only mention that
users can create such models with the following constructor, though support for these models
may be limited:

```@docs
weierstrass_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, weierstrass_f::MPolyRingElem, weierstrass_g::MPolyRingElem)
```

To check whether a Weierstrass model has been constructed over a concrete base space or over a
symbolic family, use `is_base_space_fully_specified`. This returns `true` if the base is fully
specified, and `false` otherwise.

### Concrete Toric Base Spaces

Full support exists for constructing Weierstrass models over **concrete, complete toric base spaces**.
Completeness is a technical assumption: it ensures that the set of global sections of a line
bundle forms a finite-dimensional vector space, enabling OSCAR to handle these sets efficiently.
In the future, this restriction may be relaxed. For now, completeness checks—though sometimes
slow—are performed in many methods involving Weierstrass models. To skip them and improve performance,
use the optional keyword:

```julia
completeness_check = false
```

We proceed under the assumption that the base space is a fixed, complete toric variety.

In this setting, constructing a suitable ambient space becomes significantly more efficient. First,
we focus on toric ambient spaces. Second, rather than enumerating large numbers of ambient spaces from
polytope triangulations—a resource-intensive process—we adopt a **performance-optimized** approach.
Starting from the rays and maximal cones of the toric base, we **extend** them to construct a compatible
polyhedral fan defining the toric ambient space. This avoids triangulation entirely and yields a single
suitable ambient toric variety. The resulting space may not be smooth and may differ from those commonly
used in the F-theory literature, but it is guaranteed to be compatible with the elliptic fibration structure.

Users can construct Weierstrass models over such concrete toric bases with the following constructors:

```@docs
weierstrass_model(base::NormalToricVariety; completeness_check::Bool = true)
weierstrass_model(base::NormalToricVariety, f::MPolyRingElem, g::MPolyRingElem; completeness_check::Bool = true)
```

For convenience—ideal for quick experiments and educational use—we also support constructors for Weierstrass
models over commonly used base spaces:

```@docs
weierstrass_model_over_projective_space(d::Int)
weierstrass_model_over_hirzebruch_surface(r::Int)
weierstrass_model_over_del_pezzo_surface(b::Int)
```

### Famous Weierstrass Models

Several Weierstrass models have gained popularity in the F-theory community. These models are often
associated with specific publications and may be informally referred to by author names or recognizable
keywords. For these established constructions, we provide support through the specialized `literature_model`
interface, which is discussed on the page [Literature Models](@ref literature_models).

---

## Attributes of Weierstrass Models

Weierstrass models are one way to represent an elliptic fibration as a hypersurface in an ambient space.
While different representations may vary in implementation details, they share a common structure in broad
strokes. As such, many attributes and properties are shared across model representations. These shared
components—such as `base_space`, `ambient_space`, and `fiber_ambient_space`—are documented on the page
[Common Model Ops](@ref common_model_ops).

Below, we list the attributes that are **specific to Weierstrass models** and do not generally apply to
other representations (such as [Global Tate Models](@ref global_tate_models) or
[Hypersurface Models](@ref hypersurface_models)):

```@docs
weierstrass_section_f(w::WeierstrassModel)
weierstrass_section_g(w::WeierstrassModel)
weierstrass_polynomial(w::WeierstrassModel)
hypersurface_equation(w::WeierstrassModel)
weierstrass_ideal_sheaf(w::WeierstrassModel)
calabi_yau_hypersurface(w::WeierstrassModel)
```

---

## Singularities in Weierstrass Models

Let us emphasize again that in F-theory, *singular* elliptic fibrations are of central importance
(cf. [Wei18](@cite) and references therein): singularities signal non-trivial physics. 

A key quantity in this context is the discriminant locus of the fibration—i.e., the subset of the
base space over which the elliptic fibers degenerate:

```@docs
discriminant(w::WeierstrassModel)
```

Even more critical than the discriminant itself is its **decomposition**: the discriminant locus can
split into multiple irreducible components. Over each such component, one can determine the singularity
type of the fiber. The function `singular_loci` returns a list of these irreducible base loci, together
with the corresponding singularity types:

```@docs
singular_loci(w::WeierstrassModel)
```

We discuss singularities in greater depth—including how to deform models to achieve a desired singularity
structure and how to resolve them—in [Resolving F-Theory Models](@ref resolving_f_theory_models).
