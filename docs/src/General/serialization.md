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

## Serializers

OSCAR provides several serializers for different use cases. The serializer is
passed via the `serializer` keyword argument to `save` and `load`.

### JSONSerializer (default)

The default serializer. Writes a single `.mrdi` JSON file. Referenced objects
(e.g. the parent ring of a polynomial) are stored inline under the `_refs` key.

```julia
R, x = QQ[:x]
p = x^2 + 1

# default: serializer=JSONSerializer() is implicit
save("poly.mrdi", p)
load("poly.mrdi"; params=R)

# disable inline refs (useful when the parent ring is already serialized elsewhere)
save("poly.mrdi", p; serializer=Oscar.Serialization.JSONSerializer(serialize_refs=false))
```

### MultiFileRefSerializer

Splits refs into separate files alongside the main file, using the path as a
prefix. Each referenced object is written to its own `<prefix>_<UUID>.mrdi`
file. The main file is `<prefix>.mrdi`. Use this when objects share large
sub-objects.

The `filename` argument is used as a prefix directly — no extension is added or
stripped. On overwrite, stale ref files from the previous save are removed
automatically.

```julia
Qx, x = QQ[:x]
F, a = number_field(x^2 + 2)
R, (y, z) = F[:y, :z]
p = a * y - z

save("polydata", p; serializer=Oscar.Serialization.MultiFileRefSerializer())
# creates: polydata.mrdi  polydata_<UUID>.mrdi  ...

Oscar.reset_global_serializer_state()
load("polydata"; serializer=Oscar.Serialization.MultiFileRefSerializer(), params=R)
```

Gzip compression applies to main and all ref files:

```julia
save("polydata", p; serializer=Oscar.Serialization.MultiFileRefSerializer(), compression=:gzip)
# creates: polydata.mrdi.gz  polydata_<UUID>.mrdi.gz  ...

Oscar.reset_global_serializer_state()
load("polydata"; serializer=Oscar.Serialization.MultiFileRefSerializer(), params=R)
```

### LPSerializer

Specialized for `LinearProgram{QQFieldElem}`. Writes the LP data to an external
`.lp` file (the standard LP file format) and stores only a reference to it in
the `.mrdi` file. The `basepath` argument is used as the filename prefix for the
`.lp` file.

```julia
P = cube(3)
LP = linear_program(P, [3, -2, 4]; k=2, convention=:min)

serializer = Oscar.Serialization.LPSerializer("/tmp/mydata/lp")
save("/tmp/mydata/lp.mrdi", LP; serializer=serializer)

Oscar.reset_global_serializer_state()
load("/tmp/mydata/lp.mrdi"; serializer=serializer)
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
