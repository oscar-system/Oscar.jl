```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Hypersurface Models

A **hypersurface model** is a description of an elliptic fibration whose total space is defined as the vanishing locus of a single polynomial in a suitable ambient space. Prominent examples include [Weierstrass Models](@ref) and [Global Tate Models](@ref). Our algorithmic framework is rooted in the constructions presented in [KM-POPR15](@cite).

---

## What is a Hypersurface Model?

Every elliptic fibration is birationally equivalent to an elliptic fibration represented as [Weierstrass model](@ref) (in characteristics other than 2 or 3). However, for practical purposes it is often more convenient to work with alternative descriptions. This is the main reason for working with [Global Tate Models](@ref), and extends more generally to the constructions presented in [KM-POPR15](@cite). Our implementation is focused on exactly these constructions.

The setup of a hypersurface model following [KM-POPR15](@cite) consists of the following ingredients:

- A **base space** ``B``, over which the elliptic fibration is defined.
- A **fiber ambient space** ``F``, in which the elliptic fiber appears as a hypersurface.
- Two divisor classes ``D_1`` and ``D_2`` in ``\text{Cl}(B)``, and a choice of two homogeneous coordinates of the fiber ambient space ``F``. These two coordinates transform over the base ``B`` as sections of the line bundles associated to ``D_1`` and ``D_2``, respectively. All remaining homogeneous coordinates of ``F`` transform as sections of the trivial line bundle over ``B``.
- A hypersurface equation defining the total space of the elliptic fibration as a section of the anti-canonical bundle ``\overline{K}_A`` of the full ambient space ``A``, which combines both the fiber ambient space ``F`` and the base space ``B``. This ensures that the hypersurface equation is Calabi–Yau.

It is worth noting that any elliptic fibration, for which the fiber ambient space is toric, can be cast into this form [KM-POPR15](@cite). Consequently, this approach allows for a uniform interface for constructing and manipulating hypersurface models in a way that generalizes [Global Tate Models](@ref) and [Weierstrass Models](@ref) naturally. Our standing assumption is therefore that the fiber ambient space ``F`` is toric.

---

## Constructing Hypersurface Models

### Unspecified Base Spaces

We support the construction of hypersurface models over unspecified base spaces, though computational capabilities are significantly limited in this setting. These models can be constructed using the following interface:

```@docs
hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, p::MPolyRingElem)
```

THIS METHOD MUST ALSO SUPPORT TWO DIVISOR CLASSES, WHICH ARE THEN SET FOR THE FIRST TWO COORDINATES AND ZEROS OTHERWISE. AND WE MUST ALSO HAVE A METHOD TO SET ANY TWO COORDINATES TO TRANSFORM...

### Concrete Toric Base Spaces

Full support exists for constructing hypersurface models over **concrete, complete toric base spaces**. Completeness is a technical assumption: it ensures that the set of global sections of a line bundle forms a finite-dimensional vector space, enabling OSCAR to handle these sets efficiently. In the future, this restriction may be relaxed. For now, completeness checks—though sometimes slow—are performed in many methods involving hypersurface models. To skip them and improve performance, use the optional keyword:

```julia
completeness_check = false
```

We proceed under the assumption that the base space is a fixed, complete toric variety.

Under this assumption, a toric ambient space ``A`` can be constructed algorithmically. Similar to our approach for [Weierstrass Models](@ref) and [Global Tate Models](@ref), our approach is optimized for performance. In general, several such ambient spaces ``A`` may exist. Instead of enumerating a large number of such ambient spaces, we merely compute a single one. As such, the ambient space ``A`` computed by our methods may differ from explicit choices in the literature. However, the obtained space ``A`` is guaranteed to be consistent with the fibration structure.

Users can construct hypersurface models over such concrete toric bases with the following constructor:

```@docs
hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::MPolyRingElem; completeness_check::Bool = true)
hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, indices::Vector{Int}, p::MPolyRingElem; completeness_check::Bool = true)
```

---

## Attributes

### Basic attributes

All hypersurface models come with a hypersurface equation:
```@docs
hypersurface_equation(h::HypersurfaceModel)
```
Note that this equation is a polynomial in the coordinate ring of a suitable
ambient space. This ambient space is computed by our constructors. Consequently,
the coordinate ring, in which the hypersurface equation is an element, is only
available after the model has been constructed.

Some hypersurface models parametrize the hypersurface equation with sections
of line bundles of the base space. One can obtain those sections and their
explicit polynomial expressions over the base with `explicit_model_sections`.
The parametrization of the hypersurface equation by these sections is found as
follows:
```@docs
hypersurface_equation_parametrization(h::HypersurfaceModel)
```

To specify a custom value for the hypersurface equation, first use one of the above
constructors. Once completed, employ the `tune` function described in the methods
section to set the hypersurface equation to your desired value.

The base space can be obtained with `base_space`, the ambient space with `ambient_space` and the
fiber ambient space with `fiber_ambient_space`. Recall that `is_base_space_fully_specified` will
tell if the model has been constructed over a concrete space (in which case the function returns
`true`) or a family of spaces (returning `false`).


### Attributes in toric settings

If the base space of the hypersurface model is a toric space, then we
also provide a special type for the Calabi-Yau hypersurface:
```@docs
calabi_yau_hypersurface(h::HypersurfaceModel)
```

### Attributes based on the corresponding global Tate and Weierstrass models

Currently, we do not provide functionality to convert a hypersurface model
into a Weierstrass or global Tate model. Still, for some constructions this might
be known or detailed in the literature. If the user wishes, one can then associate
a corresponding Weierstrass or global Tate model as follows:
```@docs
set_weierstrass_model(h::HypersurfaceModel, w::WeierstrassModel)
set_global_tate_model(h::HypersurfaceModel, w::GlobalTateModel)
```
These models can then be accessed with the following functions:
```@docs
weierstrass_model(h::HypersurfaceModel)
global_tate_model(h::HypersurfaceModel)
```
Provided that the corresponding Weierstrass model is known for a hypersurface model, we can compute the discriminant locus of this corresponding Weierstrass model:
```@docs
discriminant(h::HypersurfaceModel)
```
Even more critical than the discriminant itself is its decomposition into irreducible components. Over each of these components, the singularities of the elliptic fiber are classified. The function below returns the list of irreducible base loci -- arising from the discriminant of the associated Weierstrass model -- along with the corresponding singularity types:
```@docs
singular_loci(h::HypersurfaceModel)
```
!!! warning
    The classification of singularities is performed using a Monte Carlo algorithm, i.e. it involves random sampling. While the implementation has been extensively tested and meets high standards, its probabilistic nature may occasionally yield non-deterministic results.



## Methods

### Tuning

Tuning is possible with the `tune` function, cf. [Functionality for all F-theory models](@ref).
