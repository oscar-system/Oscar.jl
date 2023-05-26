```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
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
polyhedral_fan
polyhedral_fan_from_rays_action
```

```@docs
normal_fan(P::Polyhedron{T}) where T<:scalar_types
face_fan(P::Polyhedron{T}) where T<:scalar_types
```


## Auxiliary functions
```@docs
ambient_dim(PF::PolyhedralFan)
dim(PF::PolyhedralFan)
f_vector(PF::PolyhedralFan)
is_complete(PF::PolyhedralFan)
is_pointed(PF::PolyhedralFan)
is_regular(PF::PolyhedralFan)
is_simplicial(PF::PolyhedralFan)
is_smooth(PF::PolyhedralFan{QQFieldElem})
lineality_dim(PF::PolyhedralFan)
lineality_space(PF::PolyhedralFan{T}) where T<:scalar_types
maximal_cones(PF::PolyhedralFan{T}) where T<:scalar_types
cones(PF::PolyhedralFan{T}, cone_dim::Int) where T<:scalar_types
cones(PF::PolyhedralFan)
n_maximal_cones(PF::PolyhedralFan)
n_cones(PF::PolyhedralFan)
nrays(PF::PolyhedralFan)
rays(PF::PolyhedralFan{T}) where T<:scalar_types
rays_modulo_lineality(PF::PolyhedralFan{T}) where T<:scalar_types
primitive_collections(PF::PolyhedralFan)
star_subdivision(PF::PolyhedralFan{T}, n::Int) where T<:scalar_types
*(PF1::PolyhedralFan, PF2::PolyhedralFan)
```
