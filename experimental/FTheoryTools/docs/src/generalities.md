```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Common Model Ops](@id common_model_ops)

All F-theory models describe elliptic (or genus-one) fibrations. While implementation details vary by
model—e.g. [Weierstrass Model](@ref weierstrass_model), [Global Tate Model](@ref global_tate_model),
[Hypersurface Model](@ref hypersurface_model), or [Literature Model](@ref literature_model)—a core set
of functionality is shared across them. This page documents that common interface.

---

## Printouts and Verbosity Control

To receive diagnostic information (e.g. for debugging or introspection), set the verbosity level as follows:

```julia
set_verbosity_level(:FTheoryModelPrinter, 1)
```

More details on verbosity settings are available [here](https://nemocas.github.io/AbstractAlgebra.jl/dev/assertions/#AbstractAlgebra.@vprint).

---

## Model Geometry and Base Specification

These attributes describe the geometric building blocks of the F-theory model—its ambient space, fiber
ambient space, base space and the Calabi–Yau hypersurface.

```@docs
ambient_space(m::AbstractFTheoryModel)
base_space(m::AbstractFTheoryModel)
fiber_ambient_space(m::AbstractFTheoryModel)
calabi_yau_hypersurface(m::AbstractFTheoryModel)
```

You can also check whether the base space is fully specified as a concrete geometry, or whether it is
a [Families of Spaces](@ref family_of_spaces):

```@docs
is_base_space_fully_specified(m::AbstractFTheoryModel)
```

---

## Chern Classes and Euler Characteristic

The Chern classes of the variety can be computed—or retrieved if precomputed (e.g. in
[Literature Models](@ref literature_models))—using:

```@docs
chern_class(m::AbstractFTheoryModel, k::Int; check::Bool = true)
chern_classes(m::AbstractFTheoryModel; check::Bool = true)
```

These classes allow for a consistency check on the Calabi–Yau condition:

```@docs
is_calabi_yau(m::AbstractFTheoryModel; check::Bool = true)
```

The Euler characteristic is obtained by integrating the top Chern class:

```@docs
euler_characteristic(m::AbstractFTheoryModel; check::Bool = true)
```

---

## [Hodge Numbers](@id non_yet_algorithmic_advanced_attributes)

Hodge numbers are essential topological invariants of Calabi–Yau spaces. While not yet (July 2025)
computed algorithmically, they are stored for certain [Literature Models](@ref literature_models)
and can be accessed via:

```@docs
hodge_h11(m::AbstractFTheoryModel)
hodge_h12(m::AbstractFTheoryModel)
hodge_h13(m::AbstractFTheoryModel)
hodge_h22(m::AbstractFTheoryModel)
```

Once algorithmic computation is implemented, these same functions will trigger it automatically.

If Hodge numbers are available, they can be used to verify the Euler characteristic independently:

```@docs
verify_euler_characteristic_from_hodge_numbers(m::AbstractFTheoryModel; check::Bool = true)
```

---

## [Base Intersection Numbers](@id base_top_data)

The following method returns the triple self-intersection number of the anticanonical class
``\overline{K}_{B_3}`` of the 3-dimensional base:

```@docs
kbar3(m::AbstractFTheoryModel)
```

This number is of ample importance to the F-theory QSMs introduced in [CHLLT19](@cite).

---

## Modifying or Instantiating Models

These methods allow the user to modify tunable sections or instantiate a model over a concrete base.

```@docs
put_over_concrete_base(m::AbstractFTheoryModel, concrete_data::Dict{String, <:Any})
tune(w::WeierstrassModel, special_section_choices::Dict{String, <:MPolyRingElem})
tune(t::GlobalTateModel, special_ai_choices::Dict{String, <:Any})
tune(h::HypersurfaceModel, input_sections::Dict{String, <:Any})
```

---

## [Zero Section](@id zero_section_data)

These following attributes specify the zero section both as a multivariate polynomial tuple and, when
applicable, as a divisor class or coordinate ring variable index.

Return the zero section of the given model. If no zero section is known, an error is raised. This information
is not typically stored as an attribute for Weierstrass and global Tate models, whose zero sections are known.

```@docs
zero_section(::AbstractFTheoryModel)
```

Return the zero section class of a model as a cohomology class in the toric ambient space. If no zero section
class is known, an error is raised. This information is always available for Weierstrass and global Tate models,
whose zero section classes are known.

```@docs
zero_section_class(::AbstractFTheoryModel)
```

Return the index of the generator of the Cox ring of the ambient space, whose corresponding vanishing locus defines
the zero section of a model. If no zero section class is known, an error is raised. This attribute is always set
simultaneously with zero_section_class. This information is always available for Weierstrass and global Tate models,
whose zero section classes are known.

```@docs
zero_section_index(::AbstractFTheoryModel)
```

---

## [Mordell–Weil Group](@id mordell_weil_group_data)

Returns the known generators and torsion sections of the Mordell–Weil group. This includes both sections
visible in the original model and those that may arise after resolution.

```@docs
generating_sections(::AbstractFTheoryModel)
torsion_sections(::AbstractFTheoryModel)
```

You can also add new instances of these two with the following functionality.
```@docs
add_generating_section!(m::AbstractFTheoryModel, addition::Vector{String})
add_torsion_section!(m::AbstractFTheoryModel, addition::Vector{String})
```

---

## [Gauge Group](@id gauge_group_data)

Returns the (possibly reducible) gauge algebra and global gauge group structure. The gauge group is
determined as a central quotient of the algebra, specified via known discrete identifications.

```@docs
gauge_algebra(::AbstractFTheoryModel)
```

Return list of lists of matrices, where each list of matrices corresponds to a gauge factor of the same
index given by `gauge_algebra(m)`. These matrices are elements of the center of the corresponding gauge factor
and quotienting by them replicates the action of some discrete group on the center of the lie algebra. This
list combined with `gauge_algebra(m)` completely determines the gauge group of the model. If no gauge quotients
are known, an error is raised.

```@docs
global_gauge_group_quotient(::AbstractFTheoryModel)
```
