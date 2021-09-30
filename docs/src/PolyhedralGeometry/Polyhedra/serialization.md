```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["serialization.md"]
```


# Saving and loading

Objects of type `Polyhedron` can be saved to a file and loaded from a file in the
following way:
```@repl oscar
square = cube(2)
save_polyhedron(square, "square.poly")
s = load_polyhedron("square.poly")
s == square
```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and Oscar.

```@docs
save_polyhedron(P::Polyhedron, filename::String)
load_polyhedron(filename::String)
```
