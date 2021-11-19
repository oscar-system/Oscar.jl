```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```

# Introduction

The polyhedral geometry part of OSCAR provides functionality for handling
- convex polytopes, unbounded polyhedra and cones
- polyhedral fans
- linear programs

General textbooks offering details on theory and algorithms include:
- [JT13](@cite)
- [Zie95](@cite)


## Serialization

Most objects from the polyhedral geometry section can be saved through the
polymake interface in the background. These functions are documented in the
subsections on the different objects. The format of the files is JSON and you
can find details of the specification
[here](https://polymake.org/schemas/data.json).

More details on the serialization, albeit concerning the older XML format, can be
found in [GHJ16](@cite). Even though the underlying format changed to JSON, the
abstract mathematical structure of the data files is still the same.
