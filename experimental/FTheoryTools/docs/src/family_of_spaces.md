```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Families of Spaces](@id family_of_spaces)

Many F-theory constructions in the literature work without fully specifying the base space of
the elliptic fibration. Put differently, they describe an **entire family of base spaces**
rather than a single one. OSCAR supports such families via a dedicated data structure.

> **Note of caution**: This data structure is still experimental and under discussion. The implementation details may change significantly in the future. Use with care.

---

## Constructors

The following constructors are currently supported for creating a `FamilyOfSpaces` object:

```@docs
family_of_spaces(coordinate_ring::MPolyRing, grading::Matrix{Int64}, dim::Int)
family_of_spaces(coordinate_ring::MPolyDecRing{QQFieldElem, QQMPolyRing}, dim::Int)
```

---

## Attributes

The following methods provide access to core attributes of a family of spaces:

```@docs
coordinate_ring(f::FamilyOfSpaces)
weights(f::FamilyOfSpaces)
dim(f::FamilyOfSpaces)
irrelevant_ideal(f::FamilyOfSpaces)
ideal_of_linear_relations(f::FamilyOfSpaces)
```

---

## Printouts and Verbosity Control

To receive diagnostic information whenever a `FamilyOfSpaces` object is used (e.g. for debugging or introspection),
set the verbosity level as follows:

```julia
set_verbosity_level(:FTheoryModelPrinter, 1)
```

More details on verbosity settings are available in the [Hecke documentation](http://www.thofma.com/Hecke.jl/dev/features/macros/).
