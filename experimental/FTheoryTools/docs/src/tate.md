```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Global Tate Models

Global Tate models provide a powerful framework to systematically engineer elliptic fibrations
with prescribed fiber singularities. They are widely used in F-theory model building because
they make it easier to realize certain gauge symmetries and matter spectra geometrically. Their
structure is governed by Tate’s algorithm, which relates the vanishing orders of specific polynomial
coefficients to fiber singularities—classified via the refined Tate table.

---

## What Is a Global Tate Model?

A global Tate model describes a particular form of an elliptic fibration. We focus on an elliptic
fibration over a base ``B``. Consider the weighted projective space ``\mathbb{P}^{2,3,1}`` with
coordinates ``x, y, z``. In addition, consider the Tate sections

- ``a_1 \in H^0(B, \overline{K}_{B})``,
- ``a_2 \in H^0(B, \overline{K}_{B}^{\otimes 2})``,
- ``a_3 \in H^0(B, \overline{K}_{B}^{\otimes 3})``,
- ``a_4 \in H^0(B, \overline{K}_{B}^{\otimes 4})``,
- ``a_6 \in H^0(B, \overline{K}_{B}^{\otimes 6})``.

Then form a ``\mathbb{P}^{2,3,1}``-bundle over ``B`` such that:

- ``x`` transforms as a section of ``2\overline{K}_{B}``,
- ``y`` transforms as a section of ``3\overline{K}_{B}``,
- ``z`` transforms as a section of ``\mathcal{O}_B``.

In this 5-fold ambient space, a global Tate model is the hypersurface defined by the vanishing of
the Tate polynomial:

$P_T = x^3 - y^2 - x y z a_1 + x^2 z^2 a_2 - y z^3 a_3 + x z^4 a_4 + z^6 a_6\,.$

In a sufficiently small open neighborhood ``U \subset B`` of a point ``p \in B``, the polynomial
locally takes the form:

$y^2 + a_1(q) x y z + a_3(q) y z^3 = x^3 + a_2(q) x^2 z^2 + a_4(q) x z^4 + a_6(q) z^6\,,$

where the ``a_i(q)`` are the values of the Tate sections at a local coordinate ``q \in U``. This
formulation is known as a **Tate model**. One may define the elliptic fibration globally by specifying
a single Tate polynomial as above. Such constructions, known as **global Tate models**, are strictly less
general than [Weierstrass Models](@ref weierstrass_models) but have proven especially useful for model
building in F-theory: Engineering the desired fiber singularities—crucial for F-theory model building—is
often simpler for global Tate models than for [Weierstrass Models](@ref weierstrass_models).

Much like the Kodaira classification, in a global Tate model the singularities of the elliptic fiber are
characterized by the vanishing orders of the Tate sections ``a_i``, although no additional monodromy needs
to be taken into account (with a few notable exceptions). The complete classification is summarized in the
following table—commonly referred to as the **Tate table**—which is taken from [Wei10](@cite):

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


---

## Constructing Global Tate Models

Global Tate models are most robustly supported over **toric base spaces**. While there are plans to extend
support to toric schemes, covered schemes, and unspecified base families, these directions remain experimental
and are reserved for future development. This section briefly outlines constructing global Tate models over
arbitrary bases, but primarily focuses on working with a **concrete toric base**.

### Unspecified Base Spaces

As with [Weierstrass Models](@ref weierstrass_models), constructing a global Tate model begins with building
a suitable **ambient space**. When the base is not fully specified, the ambient space can be treated symbolically.
This symbolic flexibility is particularly useful when **engineering singular fibers** via specified vanishing orders
or factorizations of the Tate sections ``a_i``. This method is foundational to F-theory model building.

To support such workflows, the Tate sections ``a_i`` are defined as indeterminates of a multivariate polynomial
ring. The indeterminates are interpreted as sections of line bundles on an **auxiliary base space**. This allows
engineering desired vanishing profiles and exploring families of models over a space of base geometries. A typical
use case involves factoring the Tate sections ``a_i`` to control the singularity type over a specific divisor. For
instance, to realize an ``I_5`` fiber over the divisor ``\{w = 0\}`` in the unspecified base space, one might set:

- ``a_1 = a_{10} \times w^0``,
- ``a_2 = a_{21} \times w^1``,
- ``a_3 = a_{32} \times w^2``,
- ``a_4 = a_{43} \times w^3``,
- ``a_6 = a_{65} \times w^5``.

This formalism underlies many symbolic constructions in F-theory model building. The mathematical interpretation of
symbolic bases is subtle and under active development. Users can construct such models with the following constructor:

```@docs
global_tate_model(auxiliary_base_ring::MPolyRing, auxiliary_base_grading::Matrix{Int64}, d::Int, ais::Vector{T}) where {T<:MPolyRingElem}
```

To check whether a global Tate model has been constructed over a concrete base space or over a symbolic family, use
`is_base_space_fully_specified`. This returns `true` if the base is fully specified, and `false` otherwise.

### Concrete Toric Base Spaces

Full support exists for constructing global Tate models over **concrete, complete toric base varieties**. Completeness
is a technical assumption: it ensures that the set of global sections of a line bundle forms a finite-dimensional vector
space, enabling OSCAR to handle these sets efficiently. In the future, this restriction may be relaxed. For now,
completeness checks—though sometimes slow—are performed in many methods involving global Tate models. To skip them and
improve performance, use the optional keyword:

```julia
completeness_check = false
```

We proceed under the assumption that the base space is a fixed, complete toric variety.

In this setting, constructing a suitable ambient space becomes significantly more efficient. First, we focus on toric
ambient spaces. Second, rather than enumerating large numbers of ambient spaces from polytope triangulations—a
resource-intensive process—we adopt a **performance-optimized** approach. Starting from the rays and maximal cones of
the toric base, we **extend** them to construct a compatible polyhedral fan defining the toric ambient space. This
avoids triangulation entirely and yields a single suitable ambient toric variety. The resulting space may not be smooth
and may differ from those commonly used in the F-theory literature, but it is guaranteed to be compatible with the
elliptic fibration structure.

Users can construct global Tate models over such concrete toric bases with the following constructors:

```@docs
global_tate_model(base::NormalToricVariety; completeness_check::Bool = true)
global_tate_model(base::NormalToricVariety, ais::Vector{T}; completeness_check::Bool = true) where {T<:MPolyRingElem}
```

For convenience—ideal for quick experiments and educational use—we also support constructors for global Tate models
over commonly used base spaces:

```@docs
global_tate_model_over_projective_space(d::Int)
global_tate_model_over_hirzebruch_surface(r::Int)
global_tate_model_over_del_pezzo_surface(b::Int)
```

### Famous Global Tate Models

Several global Tate models have gained popularity in the F-theory community. These models are often
associated with specific publications and may be informally referred to by author names or recognizable
keywords. For these established constructions, we provide support through the specialized `literature_model`
interface, which is discussed on the page [Literature Models](@ref literature_models).

---

## Attributes of Global Tate Models

Global Tate models are one way to represent an elliptic fibration as a hypersurface in an ambient space. While
different representations may vary in implementation details, they share a common structure in broad strokes. As
such, many attributes and properties are shared across model representations. These shared components—such as
`base_space`, `ambient_space`, and `fiber_ambient_space`—are documented on the page
[Functionality for all F-theory models](@ref "Functionality for all F-theory models").

Below, we list the attributes that are **specific to global Tate models** and do not generally apply to other
representations (such as [Weierstrass Models](@ref weierstrass_models) or [Hypersurface Models](@ref "Hypersurface Models")):

```@docs
tate_section_a1(t::GlobalTateModel)
tate_section_a2(t::GlobalTateModel)
tate_section_a3(t::GlobalTateModel)
tate_section_a4(t::GlobalTateModel)
tate_section_a6(t::GlobalTateModel)
tate_polynomial(t::GlobalTateModel)
hypersurface_equation(t::GlobalTateModel)
tate_ideal_sheaf(t::GlobalTateModel)
weierstrass_model(t::GlobalTateModel)
calabi_yau_hypersurface(t::GlobalTateModel)
```

---

## Singularities in Global Tate Models

Let us emphasize again that in F-theory, *singular* elliptic fibrations are of central importance (cf. [Wei18](@cite)
and references therein): singularities signal non-trivial physics.

A key step in analyzing an elliptic fibration is identifying its singular fibers—those whose structure degenerates
over certain loci in the base. The **discriminant locus** is the subset of the base space over which the fibers degenerate:

```@docs
discriminant(t::GlobalTateModel)
```

More informative than the discriminant itself is its **decomposition** into irreducible components. Each component
corresponds to a locus where the fiber exhibits a distinct singularity structure. These can be classified using:

```@docs
singular_loci(t::GlobalTateModel)
```

We discuss singularities in greater depth—including how to deform models to achieve a desired singularity
structure and how to resolve them—in [Resolving F-Theory Models](@ref "Resolving F-Theory Models").
