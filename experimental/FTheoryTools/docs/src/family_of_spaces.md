```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Families of Spaces

Many F-theory constructions (e.g. in the literature) work without fully specifying
the base space of the elliptic fibrations. Put differently, those works consider
an entire family of base spaces. We aim to support this data structure.

Note of caution: This data structure is subject to discussion. The exact implementation
details may change drastically in the future. Use with care (as all experimental
code, of course).


## Constructors

We currently support the following constructor:
```@docs
family_of_spaces(coordinate_ring::MPolyRing, grading::Matrix{Int64}, dim::Int)
```

## Attributes

We currently support the following attributes:
```@docs
coordinate_ring(f::FamilyOfSpaces)
weights(f::FamilyOfSpaces)
dim(f::FamilyOfSpaces)
irrelevant_ideal(f::FamilyOfSpaces)
ideal_of_linear_relations(f::FamilyOfSpaces)
```
