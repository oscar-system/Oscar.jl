```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Hypersurface Models

A **hypersurface model** is a description of an elliptic fibration whose
total space is defined as the vanishing locus of a single polynomial in a
suitable ambient space. Prominent examples include [Weierstrass Models](@ref)
and [Global Tate Models](@ref). Our algorithmic framework is rooted in the
constructions presented in [KM-POPR15](@cite).

---

## What is a Hypersurface Model?

Every elliptic fibration is birationally equivalent to an elliptic fibration
represented as [Weierstrass model](@ref) (in characteristics other than 2 or 3).
However, for practical purposes it is often more convenient to work with alternative
descriptions. This is the main reason for working with [Global Tate Models](@ref),
and extends more generally to the constructions presented in [KM-POPR15](@cite).
Our implementation is focused on exactly these constructions.

The setup of a hypersurface model following [KM-POPR15](@cite) consists of the
following ingredients:

- A **base space** ``B``, over which the elliptic fibration is defined.
- A **fiber ambient space** ``F``, in which the elliptic fiber appears as a hypersurface.
- Two divisor classes ``D_1`` and ``D_2`` in ``\text{Cl}(B)``, and a choice of two homogeneous coordinates of the fiber ambient space ``F``. These two coordinates transform over the base ``B`` as sections of the line bundles associated to ``D_1`` and ``D_2``, respectively. All remaining homogeneous coordinates of ``F`` transform as sections of the trivial line bundle over ``B``.
- A hypersurface equation defining the total space of the elliptic fibration as a section of the anti-canonical bundle ``\overline{K}_A`` of the full ambient space ``A``, which combines both the fiber ambient space ``F`` and the base space ``B``. This ensures that the hypersurface equation is Calabi–Yau.

It is worth noting that any elliptic fibration, for which the fiber ambient space is toric,
can be cast into this form [KM-POPR15](@cite). Consequently, this approach allows for a
uniform interface for constructing and manipulating hypersurface models in a way that
generalizes [Global Tate Models](@ref) and [Weierstrass Models](@ref) naturally. Our standing
assumption is therefore that the fiber ambient space ``F`` is toric.

---

## Constructing Hypersurface Models

### Unspecified Base Spaces

We support the construction of hypersurface models over unspecified base spaces, though
computational capabilities are significantly limited in this setting. These models can
be constructed using the following interface:

```@docs
hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, p::MPolyRingElem)
hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, indices::Vector{Int}, p::MPolyRingElem)
```

### Concrete Toric Base Spaces

Full support exists for constructing hypersurface models over **concrete, complete toric base spaces**.
Completeness is a technical assumption: it ensures that the set of global sections of a line bundle
forms a finite-dimensional vector space, enabling OSCAR to handle these sets efficiently. In the
future, this restriction may be relaxed. For now, completeness checks—though sometimes slow—are
performed in many methods involving hypersurface models. To skip them and improve performance, use
the optional keyword:

```julia
completeness_check = false
```

We proceed under the assumption that the base space is a fixed, complete toric variety.

Under this assumption, a toric ambient space ``A`` can be constructed algorithmically. Similar to
our approach for [Weierstrass Models](@ref) and [Global Tate Models](@ref), our approach is optimized
for performance. In general, several such ambient spaces ``A`` may exist. Instead of enumerating a
large number of such ambient spaces, we merely compute a single one. As such, the ambient space ``A``
computed by our methods may differ from explicit choices in the literature. However, the obtained
space ``A`` is guaranteed to be consistent with the fibration structure.

Users can construct hypersurface models over such concrete toric bases with the following constructor:

```@docs
hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::MPolyRingElem; completeness_check::Bool = true)
hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, indices::Vector{Int}, p::MPolyRingElem; completeness_check::Bool = true)
```

### Famous Hypersurface Models

Several hypersurface models have gained popularity in the F-theory community. These models are often
associated with specific publications and may be informally referred to by author names or recognizable
keywords. For these established constructions, we provide support through the specialized `literature_model`
interface, which is discussed on the page [Literature Models](@ref).

---

## Attributes of Hypersurface Models

Hypersurface models represent an elliptic fibration as a hypersurface in an ambient space. While different types
of models may vary in implementation, they share a broadly similar structure. Common attributes—such as `base_space`,
`ambient_space`, and `fiber_ambient_space`—are documented on the page [Functionality for all F-theory models](@ref).

The following attributes are **specific to hypersurface models** and do not generally apply to other representations
(such as [Weierstrass Models](@ref) or [Global Tate Models](@ref)):

```@docs
hypersurface_equation(h::HypersurfaceModel)
calabi_yau_hypersurface(h::HypersurfaceModel)
weierstrass_model(h::HypersurfaceModel)
global_tate_model(h::HypersurfaceModel)
```

Currently, we do not provide automatic functionality to convert a hypersurface model into a [Weierstrass Model](@ref) or
a [Global Tate Model](@ref). However, such relations may be known or derived in the literature. If desired, users can
manually establish the connection using the functions below:

```@docs
set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)
set_global_tate_model(h::HypersurfaceModel, w::GlobalTateModel)
```

---

## Singularities in Hypersurface Models

Let us emphasize again that in F-theory, *singular* elliptic fibrations are of central importance (cf. [Wei18](@cite)
and references therein): singularities signal non-trivial physics.

### Detecting Singularities

A key step in analyzing an elliptic fibration is identifying its singular fibers—those whose structure degenerates
over certain loci in the base. The **discriminant locus** is the subset of the base space over which the fibers degenerate.

For hypersurface models, we provide this functionality only if corresponding [Weierstrass Models](@ref) or
[Global Tate Models](@ref) are known.

```@docs
discriminant(h::HypersurfaceModel)
```

More informative than the discriminant itself is its **decomposition** into irreducible components. Each component
corresponds to a locus where the fiber exhibits a distinct singularity structure. These can be classified using:

```@docs
singular_loci(h::HypersurfaceModel)
```

We discuss singularities in greater depth—including how to deform models to achieve a desired singularity
structure and how to resolve them—in [Resolving F-Theory Models](@ref).
