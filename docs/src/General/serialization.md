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


## Type Parameters

The `TypeParams` type is used internally, and is a struct for the julia type and the parameters needed to serialize an object. 
The parameters needed to serialize OSCAR types are not necessarily the same as the parameters of the type itself.
For example, take `p` with julia type `MPolyRingElem{QQField}`, this julia type is parametrized by the type of its coefficients, however
the parameter needed for serializing is its parent, that is `parent(p)`. The specific parameters needed vary depending on the type,
for more details about how `TypeParams` are used during (de)serialization see the [developer documentation](@ref dev_serialization).
Since it can be useful for users to create them when overriding parameters on load we include a small section here.

```@docs
type_params
```

## Container Types
We use the term "container types" to refer to the types `Vector`, `Tuple`, `NamedTuple`, `Set`, `Dict` and `Matrix`. 
Container types that contain serializable OSCAR types can themselves be serialized.
However, there are some caveats. Although your data type might be a `Vector` of an OSCAR type say `Vector{MPolyRingElem}`, this type might not be serializable.
To serialize a `Vector` all entries must have the output of [`type_params`](@ref) be equal (`==` returns true), this allows us to provide a size efficient storage of homogeneous data. The `Tuple` and `NamedTuple` types will store type information (type names and parameters) for each entry separately, this is the recommended way of serializing collections of types
with different `TypeParams`. 
Storing `Dict` types depends on the key and value types, for example if the keys are either `String`, `Symbol` or `Int` then the type parameters for the value may vary by key. 
In the other cases the type parameters for all keys must be equal, and likewise all values must have equal type parameters.

```jldoctest; setup=:(current=pwd(); cd(mktempdir())), teardown=:(cd(current))
julia> Qxy, (x, y) = QQ[:x, :y]
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> Qt, t = QQ[:t]
(Univariate polynomial ring in t over QQ, t)

julia> typeof([t, x, y])
Vector{RingElem} (alias for Array{RingElem, 1})

julia> save("bad-example.mrdi", [t, x, y])
ERROR: ArgumentError: Not all type parameters of Vector or Matrix entries are the same, consider using a Tuple for serialization
[...]
```
