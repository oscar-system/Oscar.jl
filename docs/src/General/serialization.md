```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Saving and loading files](@id serialization)

## Introduction

For some of our datatypes we provide a way to save them in and load them from
JSON format.  The most common OSCAR types are supported, but it will take some time until all
corners of OSCAR are covered by this effort. Our overall goal is threefold:
  - Avoid recomputation by providing an easy way to store data.
  - Increase portability by giving a convenient possibility to transport data.
  - Increase overall software quality by testing against existing data and
    tracking errors through data computed by different versions of OSCAR (or other computer algebra systems).

For more details read the [developer documentation](@ref dev_serialization).
Work on serialization is supported by the [MaRDI project](https://www.mardi4nfdi.de/about/mission). You can find out more about its Task Area 1 (Computer Algebra) [here.](https://portal.mardi4nfdi.de/wiki/Portal/TA1)

```@docs
save
load
```

## Objects that can be serialized

In this section we will list a selection of objects that may be (de-)serialized. 

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
PolyRing
PolyRingElem
MPolyRing
MPolyRingElem
```

### Groups
```julia
FPGroup
FinGenAbGroup
PcGroup
PermGroup
SubFPGroup
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

## Listing all serializable types of the current session

If you are curious about whether your type can already be serialized given your version of Oscar
you can run the following command in your active session.

```@setup oscar
using Oscar
```

```@repl oscar
Oscar.Serialization.reverse_type_map
```
