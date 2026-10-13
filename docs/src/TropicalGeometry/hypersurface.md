# Tropical hypersurfaces

## Introduction
A tropical hypersurface is a balanced polyhedral complex of codimension one.  It is dual to a regular subdivision of a Newton polytope. For more on tropical hypersurfaces, see
- Chapter 3.1 in [MS15](@cite)
- Chapter 1 in [Jos21](@cite)

#### Note:
- Objects of type `TropicalHypersurface` need to be embedded, abstract tropical hypersurfaces are currently not supported.
- The type `TropicalHypersurface` can be thought of as subtype of `TropicalVariety` in the sense that it should have all properties and features of the latter.


## Constructors
In addition to converting from `TropicalVariety`, objects of type `TropicalHypersurface` can be constructed from:
- polynomials over a tropical semiring,
- polynomials over a field and a tropical semiring map,
- subdivision of points and a choice of min- or max-convention.
```@docs
tropical_hypersurface
```

## Properties
In addition to the properties inherited from `TropicalVariety`, objects of type `TropicalHypersurface` have the following exclusive properties:
```@docs
algebraic_polynomial(TropH::TropicalHypersurface)
tropical_polynomial(TropH::TropicalHypersurface)
dual_subdivision(TropH::TropicalHypersurface)
```

## Example
The following code sets up an example and prints the vertices and rays of the tropical hypersurface (in no particular order).
```jldoctest exampleHypersurface
julia> T = tropical_semiring();

julia> Tx,(x1,x2) = polynomial_ring(T, 2);

julia> g = 1 + 2*x1^2 + 2*x2^2 + 1*x1*x2;

julia> TropH = tropical_hypersurface(g);

julia> vertRays = vertices_and_rays(TropH)
5-element SubObjectIterator{Union{PointVector{QQFieldElem}, RayVector{QQFieldElem}}}:
 [-1, -1]
 [1, 0]
 [0, 1]
 [-1//2, 1//2]
 [1//2, -1//2]
```
By broadcasting the `typeof()` command, we can see, which are vertices, and which are rays.
```jldoctest exampleHypersurface
julia> typeof.(vertRays)
5-element Vector{DataType}:
 RayVector{QQFieldElem}
 RayVector{QQFieldElem}
 RayVector{QQFieldElem}
 PointVector{QQFieldElem}
 PointVector{QQFieldElem}
```
The maximal polyhedra of our tropical hypersurface is simply the edges (both bounded and unbounded). The command `maximal_polyhedra()` gives us a list of these edges (in no particular order).
```jldoctest exampleHypersurface
julia> maxPol = maximal_polyhedra(TropH)
5-element SubObjectIterator{Polyhedron{QQFieldElem}}:
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2
 Polyhedron in ambient dimension 2
 Polytope in ambient dimension 2
 Polyhedron in ambient dimension 2
```
The polyhedrons are the unbounded edges, and the polytopes are the bounded edges.

The incidence matrix of the maximal polyhedra is a list of the relations between the elements of `vertices_and_rays(TropH)`. 
From these relations, we can draw the hypersurface. However, one should be careful, as there is no distinction between vertices and rays in the incidence matrix.
```jldoctest exampleHypersurface
julia> IncidenceMatrix(maxPol)
5×5 IncidenceMatrix
 [1, 4]
 [1, 5]
 [3, 4]
 [4, 5]
 [2, 5]
```
This is made clearer if we ask for the vertices of each of the maximal polyhedra (the bounded edges have a vertex at each end, while the unbounded only have one vertex).
```jldoctest exampleHypersurface
julia> vertices.(maxPol)
5-element Vector{SubObjectIterator{PointVector{QQFieldElem}}}:
 [[-1//2, 1//2]]
 [[1//2, -1//2]]
 [[-1//2, 1//2]]
 [[-1//2, 1//2], [1//2, -1//2]]
 [[1//2, -1//2]]
```
Instead of being between two vertices, the unbounded edges are defined by a vertex and a ray. These rays can be seen in the following way.
```jldoctest exampleHypersurface
julia> rays.(maxPol)
5-element Vector{SubObjectIterator{RayVector{QQFieldElem}}}:
 [[-1, -1]]
 [[-1, -1]]
 [[0, 1]]
 0-element SubObjectIterator{RayVector{QQFieldElem}}
 [[1, 0]]
```
