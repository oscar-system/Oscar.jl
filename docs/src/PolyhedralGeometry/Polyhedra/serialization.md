```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
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
```jldoctest
julia> square = cube(2)

julia> save("square.poly", square)

julia> s = load("square.poly")

julia> s == square

```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and OSCAR.

