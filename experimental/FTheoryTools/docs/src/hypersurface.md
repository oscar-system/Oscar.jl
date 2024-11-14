```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Hypersurface models

## Introduction

A hypersurface model describes a particular form of an elliptic fibration.
We make the following assumptions:
* The generic fiber is a hypersurface in a 2-dimensional toric fiber ambient space ``F``.
* The first two (homogeneous) coordinates in the coordinate ring of ``F`` transform in the
divisor classes ``D_1`` and ``D_2`` over the base ``B`` of the elliptic fibration. The
remaining coordinates of ``F`` are assumed to transform in the trivial bundle over ``B``.
See [KM-POPR15](@cite) for more details.

Our tools are most powerful if the base space ``B`` is a toric variety. Then, based on the
above information, it is possible to compute a toric ambient space ``A`` for the elliptic
fibration. The elliptic fibration is then a hypersurface in this toric space ``A``. Furthermore,
since we assume that this fibration is Calabi-Yau, it is clear that the hypersurface equation is
a (potentially very special) section of the ``\overline{K}_A``. This hypersurface equation
completes the information required about a hypersurface model.

Certainly, one often wishes to extend beyond this setting:
* Oftentimes, a special hypersurface equation is chosen. In the toric setting above, this will
be a polynomial in the Cox ring of the toric ambient space. Since said ambient space, and therefore
also its Cox ring, are computed in the cause of the above construction, we do not support one
constructor, which immediately accepts a special hypersurface equation. Rather, in this case,
the user must first use one of the general constructors described in the next section. Subsequently,
the `tune` function, described in the methods section below, can be employed.
* Bases other than toric spaces matter. Often, the F-theory literature will not even assume
one particular base space but rather an entire family of base spaces. We extend our machinery to
these cases, but such more general settings are typically significantly more limited than the
toric setting. This limitation originates from the nature of the matter -- the more general the geometry,
the less is known.



## Constructors

We aim to provide support for hypersurface models over the following bases:
* a toric variety,
* a toric scheme,
* a (covered) scheme.

[Often, one also wishes to obtain information about a hypersurface model without
explicitly specifying the base space. For this, please use our `tune` function,
described in the methods section below.]

Finally, we provide support for some standard constructions.

Before we detail these constructors, we must comment on the constructors over toric base
spaces. Namely, in order to construct a hypersurface model, we first have to construct
the ambient space in question. For a toric base, one way to achieve this is by means of
triangulations. However, this is a rather time consuming and computationally challenging
task, which leads to a huge number of ambient spaces. Even more, typically one wishes to
only pick one of thees many ambient spaces. For instance, a common and often appropriate
choice is a toric ambient space which contains the toric base space in a manifest way.

To circumvent this very demanding computation, our constructors operate in the opposite direction.
That is, they begin by extracting the rays and maximal cones of the chosen toric base space.
Subsequently, those rays and cones are extended to form one of the many toric ambient spaces.
This proves hugely superior in performance than going through the triangulation task of enumerating
all possible toric ambient spaces. One downside of this strategy is that the so-constructed ambient
space need not be smooth.

### A toric variety as base space

We require that the provided toric base space is complete. This is a technical limitation as of now.
The functionality of OSCAR only allows us to compute a section basis (or a finite subset thereof)
for complete toric varieties. In the future, this could be extended.

Completeness is an expensive check. Therefore, we provide an optional argument which
one can use to disable this check if desired. To this end, one passes the optional argument
`completeness_check = false` as last argument to the constructor. Here is how we can construct
a hypersurface model in OSCAR:
```@docs
hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass}, p::MPolyRingElem; completeness_check::Bool = true)
```
For convenience, it also possible to provide the hypersurface polynomial `p` as a string.


### A (covered) scheme as base space

This functionality does not yet exist.

### Base space not specified

This method constructs a hypersurface model over a base space, where
this base space is not (fully) specified. We currently provide the following constructors:
```@docs
hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{Vector{Int64}}, p::MPolyRingElem)
```
For convenience, the fiber_twist_divisor_classes can also be provided as `ZZMatrix`.


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
Provided that the corresponding Weierstrass model is known for a hypersurface
model, the following functionality is available. It returns the attribute in question
of the corresponding Weierstrass model.
```@docs
discriminant(h::HypersurfaceModel)
singular_loci(h::HypersurfaceModel)
```


## Methods

### Tuning

Tuning is possible with the `tune` function, cf. [Functionality for all F-theory models](@ref).
