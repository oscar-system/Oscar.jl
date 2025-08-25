```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Resolving F-Theory Models](@id resolving_f_theory_models)

In F-theory, the standard approach to handling singular geometries is to replace them with **smooth** ones
via **crepant resolutions**. This process preserves the Calabi–Yau condition and ensures the correct encoding
of physical data. However, several important caveats apply:

- Not all singularities admit crepant resolutions; some singularities are obstructed from being resolved without violating the Calabi–Yau condition.
- No general algorithm is known for computing a crepant resolution of a given singular geometry. In practice, one applies all known resolution techniques, guided by mathematical structure and physical expectations. A particularly prominent strategy is a sequence of **blowups**, which we discuss below.

After applying a resolution strategy, one obtains a **partially resolved** model. For the reasons stated above,
we do currently no verify whether the model has been fully resolved — i.e., whether all resolvable singularities
have been removed crepantly to the extent possible. Instead, the following method returns `true` if resolution
techniques were applied and `false` otherwise.

```@docs
is_partially_resolved(m::AbstractFTheoryModel)
```

---

## Blowups

### Manually Applying Individual Blowups

You can execute individual blowups, whether toric or not, using the following methods:

```@docs
blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e")
blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; coordinate_name::String = "e")
blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; coordinate_name::String = "e")
```

### [Registering And Extracting Known Resolution Sequences](@id working_with_resolution_sequences)

Typically, a full resolution requires a sequence of blowups. The known (weighted) blowup resolution
sequences, can be accessed with the following methods.

```@docs
resolutions(::AbstractFTheoryModel)
weighted_resolutions(::AbstractFTheoryModel)
```

If beyond this, a resolution sequence is known, it is advantageous to
register it with the model:

```@docs
add_resolution!(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
add_weighted_resolution!(m::AbstractFTheoryModel, centers::Vector{Tuple{Vector{String}, Vector{Int64}}}, exceptionals::Vector{String})
```

Once registered, the following method applies the complete sequence of blowups in a single step:

```@docs
resolve(m::AbstractFTheoryModel, index::Int)
```

---

## [Resolution Metadata Functions](@id resolution_meta_data)

These methods retrieve known resolution associated section data.

```@docs
resolution_zero_sections(::AbstractFTheoryModel)
resolution_generating_sections(::AbstractFTheoryModel)
weighted_resolution_zero_sections(::AbstractFTheoryModel)
weighted_resolution_generating_sections(::AbstractFTheoryModel)
```

You can also add to this information, if more generating sections or zero sections are known.

```@docs
add_resolution_zero_section!(m::AbstractFTheoryModel, addition::Vector{Vector{String}})
add_resolution_generating_section!(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
add_weighted_resolution_zero_section!(m::AbstractFTheoryModel, addition::Vector{Vector{String}})
add_weighted_resolution_generating_section!(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
```

---

## [Exceptional Divisors](@id exceptional_divisors)

The following methods return information on exceptional divisors introduced by toric blowups.
Therefore, these attributes are only available for models build over a concrete toric base space.

```@docs
exceptional_classes(::AbstractFTheoryModel)
exceptional_divisor_indices(::AbstractFTheoryModel)
```

---

## Analyzing the Resolved Fiber Structure

> **Note**: This functionality is currently only supported for [Global Tate Models](@ref global_tate_models).

After resolution, one typically studies the structure of the resolved fibers to extract intersection numbers
and representation-theoretic data. The method below computes the fiber components and their intersection graph
in codimension one of the base space:

```@docs
analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})
```
