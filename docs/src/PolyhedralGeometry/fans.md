```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["fans.md"]
```

# Polyhedral Fans

## Introduction

Let $\mathbb{F}$ be an ordered field; the default is that
$\mathbb{F}=\mathbb{Q}$ is the field of rational numbers and other fields are
not yet supported everywhere in the implementation.

A nonempty finite collection $\mathcal{F}$ of (polyhedral) cones in
$\mathbb{F}^n$, for $n$ fixed, is a *(polyhedral) fan* if

- the set $\mathcal{F}$ is closed with respect to taking faces and
- if $C,D\in\mathcal{F}$ then $C\cap D$ is a face of both, $C$ and $D$.

## Construction

To construct a polyhedral fan, you must pass the rays of each cone in the fan,
along with an `IncidenceMatrix` encoding which rays generate which cones.

```@docs
PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
```

```@docs
normal_fan(P::Polyhedron{T}) where T<:scalar_types
face_fan(P::Polyhedron{T}) where T<:scalar_types
```

## Saving and loading

Objects of type `PolyhedralFan` can be saved to a file and loaded from a file
in the following way:
```@repl oscar
square = cube(2)
fan = normal_fan(square)
save_polyhedralfan(fan, "F.fan")
f = load_polyhedralfan("F.fan")
collect(rays(f))
```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and Oscar.

```@docs
save_polyhedralfan(PF::PolyhedralFan, filename::String)
load_polyhedralfan(filename::String)
```

## Auxiliary functions
```@docs
ambient_dim(PF::PolyhedralFan)
dim(PF::PolyhedralFan)
f_vector(PF::PolyhedralFan)
iscomplete(PF::PolyhedralFan)
ispointed(PF::PolyhedralFan)
isregular(PF::PolyhedralFan)
issimplicial(PF::PolyhedralFan)
issmooth(PF::PolyhedralFan{fmpq})
lineality_dim(PF::PolyhedralFan)
lineality_space(PF::PolyhedralFan{T}) where T<:scalar_types
maximal_cones(PF::PolyhedralFan{T}) where T<:scalar_types
cones(PF::PolyhedralFan{T}, cone_dim::Int) where T<:scalar_types
n_maximal_cones(PF::PolyhedralFan)
nrays(PF::PolyhedralFan)
rays(PF::PolyhedralFan)
primitive_collections(PF::PolyhedralFan)
starsubdivision(PF::PolyhedralFan{T}, n::Int) where T<:scalar_types
*(PF1::PolyhedralFan, PF2::PolyhedralFan)
```

### Visualization
```@docs
visualize(PF::PolyhedralFan)
```
