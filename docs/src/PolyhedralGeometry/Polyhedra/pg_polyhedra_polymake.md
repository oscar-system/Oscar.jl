```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_polyhedra_polymake.md"]
```


# `Polyhedron` and `polymake`'s `Polytope`

Many polyhedral computations are done through `polymake`.
This is visible in the structure `Polyhedron` via a pointer  `pm_polytope` to the corresponding `polymake` object.
This allows to apply all functionality of `polymake` by calling a suitable `Polymake.jl` function on this pointer.

To allow both `Oscar`'s and `polymake`'s functionality to be applicable to a polyhedron object, it can be converted back and forth:

```@docs
Polyhedron(::Polymake.BigObject)
pm_polytope(::Polyhedron)
```

There are several design differences between `polymake` and `Oscar`.

!!! warning
    Polyhedra in `polymake` and `Polymake.jl` use homogeneous coordinates. The polyhedra in `Oscar` use affine coordinates.
