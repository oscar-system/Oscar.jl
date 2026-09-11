# polymake Integration

This section explains how OSCAR interacts with polymake.

## The Julia package [Polymake.jl](https://github.com/oscar-system/Polymake.jl)

This package provides an interface between polymake and Julia.
Its [documentation](https://oscar-system.github.io/Polymake.jl/stable/)
describes how polymake objects and their properties are accessed from
Julia, and how polymake's data types are converted.
Inside OSCAR, the package is available as `Polymake`.

## Using polymake from OSCAR

!!! warning "A route of last resort"
    Calling polymake directly is an escape hatch for functionality that
    OSCAR does not provide yet.
    Prefer OSCAR's own functions whenever they exist, and please
    [open an issue](https://github.com/oscar-system/Oscar.jl/issues)
    for what is missing, so that it can be added to OSCAR itself.

[Polyhedra](@ref "`Polyhedron` and `polymake`'s `Polytope`"),
[cones](@ref "Cones"), [polyhedral fans](@ref "Polyhedral Fans"),
[matroids](@ref "Matroids"),
[simplicial complexes](@ref "Simplicial Complexes"), and
[graphs](@ref "Graphs") in OSCAR wrap polymake objects, and
`Oscar.pm_object` returns the wrapped object.
Properties of polymake objects are accessed with the dot syntax,
and the functions of a polymake application `app` are available as
`Polymake.app.<name>`.
The results are polymake objects; `polyhedron` turns a polymake polytope
back into an OSCAR polyhedron, and `matrix` converts polymake matrices.
```jldoctest
julia> P = cube(3);

julia> pm = Oscar.pm_object(P);

julia> pm.F_VECTOR
pm::Vector<pm::Integer>
8 12 6

julia> Polymake.polytope.ambient_dim(pm)
3

julia> matrix(QQ, pm.VERTICES)
[1   -1   -1   -1]
[1    1   -1   -1]
[1   -1    1   -1]
[1    1    1   -1]
[1   -1   -1    1]
[1    1   -1    1]
[1   -1    1    1]
[1    1    1    1]

julia> polyhedron(Polymake.polytope.hypersimplex(2, 4))
Polytope in ambient dimension 4
```
Note that `Oscar.pm_object` is not part of OSCAR's official interface
and may change.
Prefer the OSCAR functions whenever they exist
(`f_vector(P)`, `vertices(P)`, and `hypersimplex(2, 4)` in the example
above), they take care of the differences described below, and please
[open an issue](https://github.com/oscar-system/Oscar.jl/issues)
if you need polymake functionality that OSCAR does not provide.
More details can be found in the
[Polymake.jl manual](https://oscar-system.github.io/Polymake.jl/stable/).

### Accessing properties and user methods

`user_method`s cannot be accessed via Julia's dot syntax, i.e. something like

```julia
c = Polymake.polytope.cube(3)
c.AMBIENT_DIM
```

will not work. Instead `user_method`s are attached as Julia functions in
their respective application. They are always written in lowercase. In the
example the following works:

```julia
c = Polymake.polytope.cube(3)
Polymake.polytope.ambient_dim(c)
```

