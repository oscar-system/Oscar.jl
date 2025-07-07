```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Functionality for all F-theory models

All F-theory models focus on elliptic (or genus-one) fibrations. Details depend on the specific
way in which the fibration is constructed/described. Still, some functionality is
common among all or at least the majority of all supported models. We will document
such common functionality here.



## Tuning Singularities

Often, one may wish to start with an existing model and modify some of its parameters—specifically its
*tunable sections*—to generate a different singularity structure. This is supported through the following
method:

```@docs
tune(w::WeierstrassModel, special_section_choices::Dict{String, <:MPolyRingElem}; completeness_check::Bool = true)
tune(t::GlobalTateModel, special_ai_choices::Dict{String, <:Any}; completeness_check::Bool = true)
tune(h::HypersurfaceModel, input_sections::Dict{String, <:Any}; completeness_check::Bool = true)
```


## Printouts

The user can decide to get information whenever a family of spaces is being used.
To this end, one invokes `set_verbosity_level(:FTheoryModelPrinter, 1)`.
More information is available [here](http://www.thofma.com/Hecke.jl/dev/features/macros/).



## Attributes of all (or most) F-theory models

```@docs
ambient_space(m::AbstractFTheoryModel)
base_space(m::AbstractFTheoryModel)
fiber_ambient_space(m::AbstractFTheoryModel)
gauge_algebra(m::AbstractFTheoryModel)
global_gauge_group_quotient(m::AbstractFTheoryModel)
chern_class(m::AbstractFTheoryModel, k::Int; check::Bool = true)
chern_classes(m::AbstractFTheoryModel; check::Bool = true)
euler_characteristic(m::AbstractFTheoryModel; check::Bool = true)
```


## Properties of all (or most) F-theory models

```@docs
is_base_space_fully_specified(m::AbstractFTheoryModel)
is_calabi_yau(m::AbstractFTheoryModel; check::Bool = true)
verify_euler_characteristic_from_hodge_numbers(m::AbstractFTheoryModel; check::Bool = true)
```


## Methods for all (or most) F-theory models

```@docs
put_over_concrete_base(m::AbstractFTheoryModel, concrete_data::Dict{String, <:Any}; completeness_check::Bool = true)
```
