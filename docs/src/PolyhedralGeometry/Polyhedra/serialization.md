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
save("square.poly", square)
s = load("square.poly")
s == square
```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and OSCAR.

