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

## Container Types
The container types used for storing and organizing data of the types listed above can also be serialized. 
There are some caveats, for example when constructing a `Vector` with entries of different types, julia will use `typejoin` to find the closest common ancestor among the types.
This may result in a type that can no longer be serialized as a `Vector`, this is due to the fact that the type parameters for serialization are not necessarily the type parameters for the data. 
To serialize a `Vector` all entries must have the output of [`type_params`](@ref) be equal (`==` returns true), this allows us to provide a size efficient storage of homogeneous data. The `Tuple` and `NamedTuple` types will store type information (type names and parameters) for each entry. 
Storing `Dict` types depends on the key and value types, for example if the keys are either `String`, `Symbol` or `Int` then the type parameters for the value may vary by key. 
In the other cases the type parameters for the keys and the values must be the same for all keys and values.

```
julia> Qxy, (x, y) = QQ[:x, :y]
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> Qt, t = QQ[:t]
(Univariate polynomial ring in t over QQ, t)

julia> typeof([t, x, y])
Vector{RingElem} (alias for Array{RingElem, 1})

julia> save("/path/to/file.mrdi", [t, x, y])
ERROR: ArgumentError: Not all type parameters of Vector or Matrix entries are the same, consider using a Tuple for serialization

...
```

## Type Parameters

The `TypeParams` type is used internally, and is a container struct for the type and the parameters need to serialize an object. Since it can be useful for users to create them when overriding parmaeters on load we include a small section here.

```@docs
type_params
```



## Listing all serializable types of the current session

If you are curious about whether your type can already be serialized given your version of Oscar
you can run the following command in your active session.

```@setup oscar
using Oscar
```

```@repl oscar
Oscar.reverse_type_map
```

