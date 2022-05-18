```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["polymake.md"]
```


# `Polyhedron` and `polymake`'s `Polytope`

Many polyhedral computations are done through `polymake`. `polymake`
([GJ00](@cite), [polymake.org](https://polymake.org)) is open source software
for research in polyhedral geometry and is attached to Julia via `Polymake.jl`
([KLT20](@cite), [Polymake.jl](https://github.com/oscar-system/Polymake.jl)).
This is visible in the structure `Polyhedron` via a pointer `pm_polytope` to
the corresponding `polymake` object.
Using `Polymake.jl` one can apply all functionality of `polymake` to the
`polymake` object hidden behind this pointer.

Sometimes it can be necessary to directly invoke some `polymake` functions on
an Oscar `Polyhedron` object (e.g. because some functionality has not yet been
made available via Oscar's interface). In that case, the following two
functions allow extracting the underlying `Polymake.jl` object from a
`Polyhedron`, respectively wrapping a `Polymake.jl` object representing a
polyhedron into a high-level `Polyhedron` object.

```@docs
Polyhedron{T}(::Polymake.BigObject) where T<:scalar_types
```

The following shows all the data currently known for a `Polyhedron`.

```@repl oscar
C = cube(3)
C.pm_polytope
```

`polymake` allows for an interactive visualization of 3-dimensional polytopes in the browser: `Polymake.visual(C.pm_polytope)`.


!!! warning
    There are several design differences between `polymake` and `Oscar`.
    Polyhedra in `polymake` and `Polymake.jl` use homogeneous coordinates. The polyhedra in `Oscar` use affine coordinates.
	 Indices in `polymake` are zero-based, whereas in `Oscar` they are one-based.

The next example shows a purely combinatorial construction of a polytope (here: a square).
In spite of being given no coordinates, `polymake` can check for us that this is a simple polytope; i.e., each vertex is contained in dimension many facets.

```@repl oscar
Q = Polymake.polytope.Polytope(VERTICES_IN_FACETS=[[0,2],[1,3],[0,1],[2,3]]);
Q.SIMPLE
```

However, without coordinates, some operations such as computing the volume cannot work:
```@repl oscar
Q.VOLUME
```

