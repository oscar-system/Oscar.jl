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

A **blowup** is a birational modification of an algebraic variety or scheme that replaces a chosen subvariety (or subscheme), called
the **center** of the blowup, with an exceptional divisor. This process provides a tool for resolving singularities in a controlled way.

In the **weighted** case, the exceptional divisor is introduced with a grading by positive integers. This allows for finer control
over the resolution process.

### Manually Applying Individual Blowups

You can execute individual blowups, whether toric or not, using the following methods:

```@docs
blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e")
blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; coordinate_name::String = "e")
blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; coordinate_name::String = "e")
```

### Data Format for Resolutions

Typically, a resolution requires a sequence of blowups.

The resolution metadata functions store information about sections in terms of their *homogeneous coordinates* after a given sequence of blowups.  

Our framework is tailored towards toric blowups.

- A **resolution entry** consists of two parts:
  1. The sequence of blowup centers (each a list of coordinate names).
  2. The list of new exceptional coordinate names introduced at each step.
- The zero and generating sections are then tracked through these blowups and recorded in terms of the homogeneous factors that describe their proper transforms.

For example, a sequence of blowups might look like

```julia
[
  [["x","y","w"], ["y","e1"], ["x","e4"], ["y","e2"], ["x","y"]],
  ["e1","e4","e2","e3","s"]
]
```

This represents the blowup sequence:

```text
(x, y, w | e1)
(y, e1   | e4)
(x, e4   | e2)
(y, e2   | e3)
(x, y    | s)
```

Each tuple ``(g_1, \dotsc, g_n | e)`` indicates that we blow up the locus ``g_1 = \dotsb = g_n = 0`` by replacing it with a new exceptional locus ``e = 0``. This is done by replacing ``g_i \mapsto e g_i`` and introducing a new homogeneous factor ``[g_1 : \dotsb : g_n]``. Thus, in the case of the above example blowup sequence, the map from the coordinates of the original to the resolved space is
```math
x \mapsto e_1 e_2^2 e_3^2 e_4 s x \\
y \mapsto e_1 e_2^2 e_3^3 e_4^2 s y \\
w \mapsto e_1 e_2 e_3 e_4 w
```
with all other coordinates being unchanged, and the new homogeneous factors (along with the original ambient weighted projective factor) are
```math
[e_1 e_2^2 e_3^2 e_4 s x : e_1 e_2^2 e_3^3 e_4^2 s y : z] [e_2 e_3 s x : e_2 e_3^2 e_4 s y : w] [e_3 s y : e_1] [s x : e_4] [s y : e_2] [x : y]\,.
```
Each blowup ``(g_1, \dotsc, g_n | e)`` has an associated rescaling ``(g_1, \dotsc, g_n, e) \sim (\lambda g_1, \dotsc, \lambda g_n, \lambda^{-1} e)``, under which the products ``e g_i`` are invariant.

For a weighted blowup, along with a center and exceptional coordinate ``(g_1, \dotsc, g_n | e)``, a grading vector ``\mu = (\mu_1, \dotsc, \mu_n)`` is specified, and the blowup is carried out by the replacement ``g_i \mapsto e^{\mu_i} g_i``. The associated rescaling is then ``(g_1, \dotsc, g_n, e) \sim (\lambda^{\mu_1} g_1, \dotsc, \lambda^{\mu_n} g_n, \lambda^{-1} e)``, under which the products ``e^{\mu_i} g_i`` are invariant.

### [Registering And Extracting Known Resolution Sequences](@id working_with_resolution_sequences)

The known (weighted) blowup resolution sequences, can be accessed with the following methods.

```@docs
resolutions(::AbstractFTheoryModel)
weighted_resolutions(::AbstractFTheoryModel)
```

If beyond this, a resolution sequence is known, it is advantageous to register it with the model:

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

The following methods retrieve known resolution associated section data.

```@docs
resolution_zero_sections(::AbstractFTheoryModel)
resolution_generating_sections(::AbstractFTheoryModel)
weighted_resolution_zero_sections(::AbstractFTheoryModel)
weighted_resolution_generating_sections(::AbstractFTheoryModel)
```

You can also add to this information, if more resolution generating sections or zero sections are known.

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
