```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Resolving F-Theory Models](@id resolving_f_theory_models)

In F-theory, the standard approach to handling singular geometries is to replace them with **smooth** ones
via **crepant resolutions**. This process preserves the Calabi–Yau condition and ensures the correct encoding
of physical data. However, several important caveats apply:

- Not all singularities admit crepant resolutions, rather some singularities are obstructed from being resolved without violating the Calabi–Yau condition.
- No general algorithm is known for computing a crepant resolution of a given singular geometry. In practice, one applies all known resolution techniques, guided by mathematical structure and physical expectations. A particularly prominent strategy is a sequence of **blowups**, which we will discuss below.

After applying a resolution strategy, one obtains a **partially resolved** model. For the reasons stated above,
we do not currently verify whether the model has been fully resolved—i.e., whether all resolvable
singularities have been removed crepantly to the extent possible:

```@docs
is_partially_resolved(m::AbstractFTheoryModel)
```

## Blowups

We can execute individual blowups, be they toric or not, with the following functions:

```@docs
blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e")
blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; coordinate_name::String = "e")
blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; coordinate_name::String = "e")
```

Typically, an entire sequence of blowups is needed to resolve an F-Theory model. If known, it
can therefore be advantageous to inform the model of the resolution:

```@docs
add_resolution(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
```

Subsequently, the following method can be used to execute the entire sequence of blowups all
at once:

```@docs
resolve(m::AbstractFTheoryModel, index::Int)
```

---

## Resolutions and Section Data

The following functions describe known crepant resolutions of singularities and their associated sections. The
resolution data is represented as blow-up centers and corresponding exceptional divisors, either in standard or
weighted format. Generating and zero sections are tracked for each resolution when known.

Return the list of all known resolutions for the given model. If no resolutions are known, an error is raised.

```@docs
resolutions(::AbstractFTheoryModel)
```

Return a list of known Mordell–Weil zero sections for the given model after each blowup of each known resolution.
Each element of the list corresponds to a known resolution (in the same order). If no resolution zero sections are
known, an error is raised.

```@docs
resolution_zero_sections(::AbstractFTheoryModel)
```

Return a list of lists of known Mordell–Weil generating sections for the given model after each known resolution.
Each element of the outer list corresponds to a known resolution (in the same order), and each element of the list
associated to a given resolution corresponds to a known generating section (in the same order). If no resolution
generating sections are known, an error is raised.

```@docs
resolution_generating_sections(::AbstractFTheoryModel)
```

Return the list of all known weighted resolutions for the given model. If no weighted resolutions are known, an error
is raised.

```@docs
weighted_resolutions(::AbstractFTheoryModel)
```

Return a list of known Mordell–Weil zero sections for the given model after each known weighted resolution. Each element of
the list corresponds to a known weighted resolution (in the same order). If no weighted resolution zero sections are known,
an error is raised.

```@docs
weighted_resolution_zero_sections(::AbstractFTheoryModel)
```

Return a list of lists of known Mordell–Weil generating sections for the given model after each known weighted resolution.
Each element of the outer list corresponds to a known weighted resolution (in the same order), and each element of the list
associated to a given weighted resolution corresponds to a known generating section (in the same order). If no weighted
resolution generating sections are known, an error is raised.

```@docs
weighted_resolution_generating_sections(::AbstractFTheoryModel)
```

---

## Exceptional Divisors

```@docs
exceptional_classes(::AbstractFTheoryModel)
exceptional_divisor_indices(::AbstractFTheoryModel)
```

These functions return information on exceptional divisors introduced by toric blowups. Available only for models over a concrete
base that is a `NormalToricVariety`. The cohomology classes are expressed in the toric ambient space.

---

## Analyzing the resolved Fiber Structure

The following is currently only supported for [Global Tate Models](@ref global_tate_models).

After resolution, one typically studies the structure of the (resolved) fibers to extract intersection numbers and
representation theory information. The following method computes the fiber components and their intersection graph:

```@docs
analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})
```
