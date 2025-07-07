```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Resolving F-Theory Models

In F-theory, the standard approach to handling singular geometries is to replace them with **smooth** ones
via **crepant resolutions**. This process preserves the Calabi–Yau condition and ensures the correct encoding
of physical data. However, several important caveats apply:

- Not all singularities admit crepant resolutions, rather some singularities are obstructed from being resolved without violating the Calabi–Yau condition. No algorithm is known to the authors that determines whether a given singularity admits a crepant resolution.
- Likewise, no general algorithm is known for computing a crepant resolution of a given singular geometry. In practice, one applies all known resolution techniques, guided by mathematical structure and physical expectations. A particularly prominent strategy is a sequence of **blowups**, which we will discuss below.

After applying a resolution strategy, one obtains a **partially resolved** model. For the reasons stated above,
we do not currently verify whether the model has been fully resolved—i.e., whether all resolvable
singularities have been removed crepantly to the extent possible:

```@docs
is_partially_resolved(m::AbstractFTheoryModel)
```

We can execute individual blowups with the following functions:

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


### Analyzing the resolved Fiber Structure

The following is currently only supported for [Global Tate Models](@ref).

After resolution, one typically studies the structure of the (resolved) fibers to extract intersection numbers and
representation theory information. The following method computes the fiber components and their intersection graph:

```@docs
analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})
```
