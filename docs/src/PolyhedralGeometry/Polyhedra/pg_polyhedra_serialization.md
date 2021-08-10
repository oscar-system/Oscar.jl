```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_polyhedra_serialization.md"]
```


# Saving and loading

Objects of type Polyhedron can be saved to a file and loaded from a file in the
following way:
```@repl oscar
square = cube(2)
save_polyhedron(square, "square.poly")
s = load_polyhedron("square.poly")
s == square
```
The file is in json format and contains all the underlying polymake object. In
particular, this file can now be read by both polymake and Oscar.


