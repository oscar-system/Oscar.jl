# Tropical varieties

## Introduction
Tropial varieties (in OSCAR) are weighted polyhedral complexes.  They may arise as tropicalizations of polynomial ideals or from operations on the more specialized types of tropical varieties, such as the stable intersection of tropical hypersurfaces.  For more on the first, see
- Chapter 3 in [MS15](@cite)
- Chapter 2.8 in [Jos21](@cite)

#### Note:
- Objects of type `TropicalVariety` need to be embedded, abstract tropical varieties are currently not supported.
- The type `TropicalVariety` can be thought of as supertype of `TropicalHypersurface`, `TropicalCurve`, and `TropicalLinearSpace` in the sense that the latter three should have all properties and features of the former.
- Embedded tropical varieties are polyhedral complexes with multiplicities and should have all properties of polyhedral complexes

## Constructor
Objects of type `TropicalVariety` can be constructed as follows:
```@docs
tropical_variety(Sigma::PolyhedralComplex, mult::Vector{ZZRingElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
```

## Properties
Objects of type `TropicalVariety` (and `TropicalHypersurface`, `TropicalCurve`, `TropicalLinearSpace`) have the following properties:
```@docs
polyhedral_complex(TropV::TropicalVariety)
ambient_dim(TropV::TropicalVariety)
codim(TropV::TropicalVariety)
dim(TropV::TropicalVariety)
f_vector(TropV::TropicalVariety)
lineality_dim(TropV::TropicalVariety)
lineality_space(TropV::TropicalVariety)
maximal_polyhedra(TropV::TropicalVariety)
maximal_polyhedra_and_multiplicities(TropV::TropicalVariety)
minimal_faces(TropV::TropicalVariety)
multiplicities(TropV::TropicalVariety)
n_maximal_polyhedra(TropV::TropicalVariety)
n_polyhedra(TropV::TropicalVariety)
n_vertices(TropV::TropicalVariety)
is_pure(TropV::TropicalVariety)
is_simplicial(TropV::TropicalVariety)
rays(TropV::TropicalVariety)
rays_modulo_lineality(TropV::TropicalVariety)
stable_intersection(::Union{Tuple{minOrMax}, Tuple{Oscar.TropicalVarietySupertype{minOrMax, true}, Oscar.TropicalVarietySupertype{minOrMax, true}}, Tuple{Oscar.TropicalVarietySupertype{minOrMax, true}, Oscar.TropicalVarietySupertype{minOrMax, true}, Union{Nothing, Vector{Int64}}}} where minOrMax)
tropical_prevariety
vertices_and_rays(TropV::TropicalVariety)
vertices(TropV::TropicalVariety)
visualize(TropV::TropicalVariety)
```
