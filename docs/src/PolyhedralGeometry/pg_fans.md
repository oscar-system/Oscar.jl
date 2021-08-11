```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_fans.md"]
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
along with an `IncidenceMatrix` encoding which rays are extremal with respect
to which cones.

```@docs
PolyhedralFan(Rays::Union{Oscar.MatElem,AbstractMatrix}, Incidence::IncidenceMatrix)
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
The file is in json format and contains all the underlying polymake object. In
particular, this file can now be read by both polymake and Oscar.
