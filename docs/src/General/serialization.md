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

```jldoctest; setup=:(current=pwd(); cd(mktempdir())), teardown=:(cd(current))
julia> R, x = QQ[:x];

julia> p = x^2 + 1
x^2 + 1

julia> save("poly.mrdi", p)

julia> save("poly.mrdi", p; serializer=Oscar.Serialization.JSONSerializer(serialize_refs=false))
```

### MultiFileRefSerializer

Splits refs into separate files alongside the main file, using the path as a
prefix. Each referenced object is written to its own `<prefix>_<UUID>.mrdi`
file. The main file is `<prefix>.mrdi`. Use this when objects share large
sub-objects.

The `filename` argument is used as a prefix directly — no extension is added or
stripped. On overwrite, stale ref files from the previous save are removed
automatically.

Without compression:

```jldoctest; setup=:(current=pwd(); cd(mktempdir())), teardown=:(cd(current)), filter = r"_[-0-9a-f]*\.mrdi"
julia> Qx, x = QQ[:x];

julia> F, a = number_field(x^2 + 2);

julia> R, (y, z) = F[:y, :z];

julia> p = a * y - z;

julia> poly_dir = mkdir("poly_dir");

julia> save(joinpath(poly_dir, "polydata"), p; serializer=Oscar.Serialization.MultiFileRefSerializer())

julia> readdir(poly_dir)
4-element Vector{String}:
 "polydata.mrdi"
 "polydata_52e4c5e2-ff94-4b1c-9832-a61d8a54331a.mrdi"
 "polydata_c8be128f-a28c-4d08-a67f-5e1c2ad9409d.mrdi"
 "polydata_f2bd8b4b-e6a7-4961-943b-1d68306889a2.mrdi"
```

With `gzip` compression:

```jldoctest; setup=:(current=pwd(); cd(mktempdir())), teardown=:(cd(current)), filter = r"_[-0-9a-f]*\.mrdi"
julia> Qx, x = QQ[:x];

julia> F, a = number_field(x^2 + 2);

julia> R, (y, z) = F[:y, :z];

julia> p = a * y - z;

julia> poly_dir = mkdir("poly_dir");

julia> save(joinpath(poly_dir, "polydata"), p; serializer=Oscar.Serialization.MultiFileRefSerializer(), compression=:gzip)

julia> readdir(poly_dir)
4-element Vector{String}:
 "polydata.mrdi.gz"
 "polydata_52e4c5e2-ff94-4b1c-9832-a61d8a54331a.mrdi.gz"
 "polydata_c8be128f-a28c-4d08-a67f-5e1c2ad9409d.mrdi.gz"
 "polydata_f2bd8b4b-e6a7-4961-943b-1d68306889a2.mrdi.gz"
```

### LPSerializer

Specialized for `LinearProgram{QQFieldElem}`. Writes the LP data to an external
`.lp` file (the standard LP file format) and stores only a reference to it in
the `.mrdi` file. The `basepath` argument is used as the filename prefix for the
`.lp` file.

```jldoctest; setup=:(current=pwd(); cd(mktempdir())), teardown=:(cd(current)), filter = r"_[-0-9a-f]*\.lp"
julia> P = cube(3);

julia> LP = linear_program(P, [3, -2, 4]; k=2, convention=:min);

julia> lp_dir = mkdir("lp_dir");

julia> serializer = Oscar.Serialization.LPSerializer(joinpath(lp_dir, "lp"));

julia> save(joinpath(lp_dir, "lp.mrdi"), LP; serializer=serializer);

julia> readdir(lp_dir)
2-element Vector{String}:
 "lp.mrdi"
 "lp_10612199771096508645.lp"
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
