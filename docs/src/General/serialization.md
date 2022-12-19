```@contents
Pages = ["serialization.md"]
```

# Serialization

## Introduction

For some of our datatypes we provide a way to save them in and load them from
JSON format. This is still experimental and it will take some time until all
corners of OSCAR are covered by this effort. The goal of this effort is
threefold:
  - Avoid recomputation by providing an easy way to store data.
  - Increase portability by giving a convenient possibility to transport data.
  - Increase overall software quality by testing against existing data and
    tracking errors through data computed by different versions.

```@docs
save
load
```

## Objects that can be serialized

In this section we will list objects that may be (de-)serialized. This list may
be incomplete.

Many low level objects may be stored and these in turn allow serializing higher
level objects. Such low level objects are various types of matrices, vectors
and sets.

### Combinatorics
```julia
Graph
SimplicialComplex
```

### Commutative Algebra
```julia
Ideal
Polynomial
PolynomialRing
```

### Polyhedral Geometry
```julia
Cone
LinearProgram
PolyhedralFan
PolyhedralComplex
Polyhedron
SubdivisionOfPoints
```

### Toric Geometry
```julia
NormalToricVariety
ToricDivisor
```

### Tropical Geometry
```julia
TropicalCurve
TropicalHypersurface
```
