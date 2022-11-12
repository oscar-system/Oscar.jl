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
an OSCAR `Polyhedron` object (e.g. because some functionality has not yet been
made available via OSCAR's interface). In that case, the following two
functions allow extracting the underlying `Polymake.jl` object from a
`Polyhedron`, respectively wrapping a `Polymake.jl` object representing a
polyhedron into a high-level `Polyhedron` object.

```@docs
Polyhedron{T}(::Polymake.BigObject) where T<:scalar_types
```

The following shows all the data currently known for a `Polyhedron`.

```jldoctest
julia> C = cube(3)
A polyhedron in ambient dimension 3

julia> C.pm_polytope
type: Polytope<Rational>
description: cube of dimension 3

AFFINE_HULL


BOUNDED
	true

CONE_AMBIENT_DIM
	4

CONE_DIM
	4

FACETS
  1   1   0   0
  1  -1   0   0
  1   0   1   0
  1   0  -1   0
  1   0   0   1
  1   0   0  -1

VERTICES_IN_FACETS
	{0 2 4 6}
	{1 3 5 7}
	{0 1 4 5}
	{2 3 6 7}
	{0 1 2 3}
	{4 5 6 7}

```

`polymake` allows for an interactive visualization of 3-dimensional polytopes in the browser: `Polymake.visual(C.pm_polytope)`.


!!! warning
    There are several design differences between `polymake` and `OSCAR`.
    Polyhedra in `polymake` and `Polymake.jl` use homogeneous coordinates. The polyhedra in `OSCAR` use affine coordinates.
	 Indices in `polymake` are zero-based, whereas in `OSCAR` they are one-based.

The next example shows a purely combinatorial construction of a polytope (here: a square).
In spite of being given no coordinates, `polymake` can check for us that this is a simple polytope; i.e., each vertex is contained in dimension many facets.

```jldoctest
julia> Q = Polymake.polytope.Polytope(VERTICES_IN_FACETS=[[0,2],[1,3],[0,1],[2,3]]);

julia> Q.SIMPLE
true

```

However, without coordinates, some operations such as computing the volume cannot work:
```jldoctest
julia> Q.VOLUME
ERROR: UndefVarError: Q not defined
Stacktrace:
 [1] top-level scope
   @ none:1

```

