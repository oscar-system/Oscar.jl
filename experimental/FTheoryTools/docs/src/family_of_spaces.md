```@meta
CurrentModule = Oscar
```

# Family of Spaces

Literature constructions in F-theory often work without fully specifying the
base space of the elliptic fibrations. Put differently, those works consider
an entire family of base spaces. We aim to support this data structure.

Note of caution: This is subject to discussion and the exact implementation
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
stanley_reisner_ideal(f::FamilyOfSpaces)
irrelevant_ideal(f::FamilyOfSpaces)
ideal_of_linear_relations(f::FamilyOfSpaces)
```
