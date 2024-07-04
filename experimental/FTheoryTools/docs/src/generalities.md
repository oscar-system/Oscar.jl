```@meta
CurrentModule = Oscar
```

# Functionality for all F-theory models

All F-theory models focus on elliptic (or genus-one) fibrations. Details depend on the specific
way in which the fibration is constructed/described. Still, some functionality is
common among all or at least the majority of all supported models. We will document
such common functionality here.


## Family of Spaces

Many F-theory constructions (e.g. in the literature) work without fully specifying
the base space of the elliptic fibrations. Put differently, those works consider
an entire family of base spaces. We aim to support this data structure.

Note of caution: This data structure is subject to discussion. The exact implementation
details may change drastically in the future. Use with care (as all experimental
code, of course).


### Constructors

We currently support the following constructor:
```@docs
family_of_spaces(coordinate_ring::MPolyRing, grading::Matrix{Int64}, dim::Int)
```

### Attributes

We currently support the following attributes:
```@docs
coordinate_ring(f::FamilyOfSpaces)
weights(f::FamilyOfSpaces)
dim(f::FamilyOfSpaces)
stanley_reisner_ideal(f::FamilyOfSpaces)
irrelevant_ideal(f::FamilyOfSpaces)
ideal_of_linear_relations(f::FamilyOfSpaces)
```

### Printouts

The user can decide to get information whenever a family of spaces is being used.
To this end, one invokes `set_verbosity_level(:FTheoryModelPrinter, 1)`.
More information is available [here](http://www.thofma.com/Hecke.jl/dev/features/macros/).



## Attributes of all (or most) F-theory models

```@docs
ambient_space(m::AbstractFTheoryModel)
base_space(m::AbstractFTheoryModel)
fiber_ambient_space(m::AbstractFTheoryModel)
explicit_model_sections(m::AbstractFTheoryModel)
defining_section_parametrization(m::AbstractFTheoryModel)
```


## Properties of all (or most) F-theory models

```@docs
is_base_space_fully_specified(m::AbstractFTheoryModel)
is_partially_resolved(m::AbstractFTheoryModel)
```


## Methods for all (or most) F-theory models

```@docs
blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e")
tune(m::AbstractFTheoryModel, p::MPolyRingElem; completeness_check::Bool = true)
put_over_concrete_base(m::AbstractFTheoryModel, concrete_data::Dict{String, <:Any}; completeness_check::Bool = true)
```
