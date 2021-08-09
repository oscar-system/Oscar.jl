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
This is visible in the structure `Polyhedron` via a pointer `pm_polytope` to the corresponding `polymake` object.
This allows to apply all functionality of `polymake` by calling a suitable `Polymake.jl` function on this pointer.

To allow both `Oscar`'s and `polymake`'s functionality to be applicable to a polyhedron object, it can be converted back and forth:

```@docs
Polyhedron(::Polymake.BigObject)
pm_polytope(::Polyhedron)
```

The following shows all the data "known to" a `Polyhedron`.

```@repl oscar
c = cube(3)
c.pm_polytope
```

The following example shows a purely combinatorial construction of a polytope.
Inspite of being given no coordinates, we can see that this is a simple polytope; i.e., each vertex is contained in dimension many facets.

```@repl oscar
X = Polymake.@pm polytope.Polytope(VERTICES_IN_FACETS=c.pm_polytope.VERTICES_IN_FACETS);
X.SIMPLE
```

However, without coordinates, e.g., computing the volume cannot work:
```@repl oscar
X.VOLUME
```

There are several design differences between `polymake` and `Oscar`.

!!! warning
    Polyhedra in `polymake` and `Polymake.jl` use homogeneous coordinates. The polyhedra in `Oscar` use affine coordinates.
	Indices in `polymake` are zero-based, whereas in `Oscar` they are one-based.

