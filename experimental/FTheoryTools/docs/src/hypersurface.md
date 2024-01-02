```@meta
CurrentModule = Oscar
```

# Hypersurface models

## Introduction

A hypersurface model describes a particular form of an elliptic fibration.
For now, we consider such models in toric settings only, that is we restrict
to a toric base space ``B`` and assume that the generic fiber is a
hypersurface in a 2-dimensional toric fiber ambient space ``F``.

In addition, we shall assume that the first two homogeneous coordinates of
``F`` transform in the divisor classes ``D_1`` and ``D_2`` over the base ``B``
of the elliptic fibration. The remaining coordinates of ``F`` are assumed to
transform in the trivial bundle over ``B``. See [KM-POPR15](@cite) for more details.

Given this set of information, it is possible to compute a toric ambient space ``A``
for the elliptic fibration. The elliptic fibration is then a hypersurface in this toric
space ``A``. Furthermore, since we assume that this fibration is Calabi-Yau,
it is clear that the hypersurface equation is a (potentially very special) section
of the ``\overline{K}_A``. This hypersurface equation completes the information
required about a hypersurface model.


## Constructors

We aim to provide support for hypersurface models over the following bases:
* a toric variety,
* a toric scheme,
* a (covered) scheme.

[Often, one also wishes to obtain information about a hypersurface model without
explicitly specifying the base space. Also for this application, we provide support.]

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
 `completeness_check = false` as last argument to the constructor. The following examples
 demonstrate this:
```@docs
hypersurface_model(base::NormalToricVariety; completeness_check::Bool = true)
hypersurface_model(base::NormalToricVariety, fiber_ambient_space::NormalToricVariety, D1::ToricDivisorClass, D2::ToricDivisorClass; completeness_check::Bool = true)
```

### A (covered) scheme as base space

This functionality does not yet exist.

### Base space not specified

This method constructs a hypersurface model over a base space, where
this base space is not (fully) specified. We currently provide the following constructors:
```@docs
hypersurface_model(auxiliary_base_vars::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int, fiber_ambient_space::NormalToricVariety, D1::Vector{Int64}, D2::Vector{Int64}, p::MPolyRingElem)
```

### Standard constructions

We provide convenient constructions of hypersurface models over
famous base spaces. Currently, we support the following:
```@docs
hypersurface_model_over_projective_space(d::Int)
hypersurface_model_over_hirzebruch_surface(r::Int)
hypersurface_model_over_del_pezzo_surface(b::Int)
```


## Attributes

### Basic attributes

For all hypersurface models -- irrespective over whether the base is toric or not -- we support
the following attributes:
```@docs
hypersurface_equation(h::HypersurfaceModel)
```
One can also decide to specify a custom hypersurface equation:
```@docs
set_hypersurface_equation(h::HypersurfaceModel, p::MPolyRingElem)
```
The fiber ambient space can be accessed via
```@docs
fiber_ambient_space(h::HypersurfaceModel)
```
In case the hypersurface model is constructed over a not fully specified base,
recall that we construct an auxiliary (toric) base space as well as an
auxiliary (toric) ambient space. The (auxiliary) base and ambient space can
be accessed with the following functions:
```@docs
base_space(h::HypersurfaceModel)
ambient_space(h::HypersurfaceModel)
```
The following method allows to tell if the base/ambient space is auxiliary or not:
```@docs
is_base_space_fully_specified(h::HypersurfaceModel)
```
The user can decide to get an information whenever an auxiliary base space,
auxiliary ambient space or auxiliary hypersurface have been computed.
To this end, one invokes `set_verbosity_level(:HypersurfaceModel, 1)`.
More background information is available
[here](http://www.thofma.com/Hecke.jl/dev/features/macros/).

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

## Properties

```@docs
is_partially_resolved(h::HypersurfaceModel)
```
