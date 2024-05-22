# Saving and loading files

## Introduction

For some of our datatypes we provide a way to save them in and load them from
JSON format. This is still experimental and it will take some time until all
corners of OSCAR are covered by this effort. The goal of this effort is
threefold:
  - Avoid recomputation by providing an easy way to store data.
  - Increase portability by giving a convenient possibility to transport data.
  - Increase overall software quality by testing against existing data and
    tracking errors through data computed by different versions.

Work on serialization is supported by the MaRDI project. You can find out more about
Mardi [here.](https://www.mardi4nfdi.de/about/mission)

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
polynomial_ring
```

### Groups
```julia
FPGroup
FinGenAbGroup
PcGroup
PermGroup
SubPcGroup
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

## Reading from and writting to external software

We write to files following the MaRDI file format which is an extension of JSON
defined [here](https://www.oscar-system.org/schemas/mrdi.json). The MaRDI file format
aims to be language agnostic. If you would like to understand more about how our files
are written or how you could implement a way to load files serialized with OSCAR
into your computer algebra setup see the developer documentation on [serialization](@ref dev_serialization).
