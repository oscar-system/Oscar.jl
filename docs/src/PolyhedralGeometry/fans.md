```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
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
arrangement_polynomial
dim(PF::PolyhedralFan)
f_vector(PF::PolyhedralFan)
is_complete(PF::PolyhedralFan)
is_pointed(PF::PolyhedralFan)
is_regular(PF::PolyhedralFan)
is_simplicial(PF::PolyhedralFan)
is_smooth(PF::PolyhedralFan{QQFieldElem})
lineality_dim(PF::PolyhedralFan)
lineality_space(PF::PolyhedralFan)
maximal_cones(PF::PolyhedralFan)
cones(PF::PolyhedralFan, cone_dim::Int)
cones(PF::PolyhedralFan)
minimal_supercone_coordinates(PF::PolyhedralFan, v::AbstractVector{<:RationalUnion})
minimal_supercone_indices(PF::PolyhedralFan, v::AbstractVector{<:RationalUnion})
is_minimal_supercone_coordinate_vector(PF::PolyhedralFan, v::AbstractVector{<:RationalUnion})
standard_coordinates(PF::PolyhedralFan, coords::AbstractVector{<:RationalUnion})
n_maximal_cones(PF::PolyhedralFan)
n_cones(PF::PolyhedralFan)
n_rays(PF::PolyhedralFan)
rays(PF::PolyhedralFan)
rays_modulo_lineality(PF::PolyhedralFan)
primitive_collections(PF::PolyhedralFan)
star_subdivision(PF::PolyhedralFan, n::Int)
star_subdivision(PF::PolyhedralFan, exceptional_ray::AbstractVector{<:IntegerUnion})
*(PF1::PolyhedralFan{QQFieldElem}, PF2::PolyhedralFan{QQFieldElem})
```
