#TODO: inward/outward options? via polymake changes?

"""
    normal_fan(P::Polyhedron)

Return the normal fan of `P`. The maximal cones of the normal fan of `P` are
dual to the edge cones at the vertices of `P`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C)
A polyhedral fan in ambient dimension 3

julia> rays(NF)
6-element SubObjectIterator{RayVector{fmpq}}:
 [1, 0, 0]
 [-1, 0, 0]
 [0, 1, 0]
 [0, -1, 0]
 [0, 0, 1]
 [0, 0, -1]
```
"""
function normal_fan(P::Polyhedron{T}) where T<:scalar_types
   pmp = pm_object(P)
   pmnf = Polymake.fan.normal_fan(pmp)
   return PolyhedralFan{T}(pmnf)
end

"""
    face_fan(P::Polyhedron)

Return the face fan of `P`. The polytope `P` has to contain the origin, then
the maximal cones of the face fan of `P` are the cones over the facets of `P`.

# Examples
By definition, this bounded polyhedron's number of facets equals the amount of
maximal cones of its face fan.
```jldoctest
julia> C = cross(3);

julia> FF = face_fan(C)
A polyhedral fan in ambient dimension 3

julia> n_maximal_cones(FF) == nfacets(C)
true
```
"""
function face_fan(P::Polyhedron{T}) where T<:scalar_types
   pmp = pm_object(P)
   pmff = Polymake.fan.face_fan(pmp)
   return PolyhedralFan{T}(pmff)
end
