#TODO: inward/outward options? via polymake changes?

"""
    normal_fan(P::Polyhedron)

Return the normal fan of `P`.

# Examples
The rays of a normal fan of a cube point in every positive and negative unit
direction.
```jldoctest
julia> C = cube(3);

julia> NF = normal_fan(C)
A polyhedral fan in ambient dimension 3

julia> collect(rays(NF))
6-element Vector{Polymake.Vector{Polymake.Rational}}:
 pm::Vector<pm::Rational>
1 0 0
 pm::Vector<pm::Rational>
-1 0 0
 pm::Vector<pm::Rational>
0 1 0
 pm::Vector<pm::Rational>
0 -1 0
 pm::Vector<pm::Rational>
0 0 1
 pm::Vector<pm::Rational>
0 0 -1
```
"""
function normal_fan(P::Polyhedron)
   pmp = pm_polytope(P)
   pmnf = Polymake.fan.normal_fan(pmp)
   return PolyhedralFan(pmnf)
end

"""
    face_fan(P::Polyhedron)

Return the face fan of `P`.

# Examples
By definition, this bounded polyhedron's number of facets equals the amount of
maximal cones of its face fan.
```jldoctest
julia> C = cross(3);

julia> FF = face_fan(C)
A polyhedral fan in ambient dimension 3

julia> nmaximal_cones(FF) == nfacets(C)
true
```
"""
function face_fan(P::Polyhedron)
   pmp = pm_polytope(P)
   pmff = Polymake.fan.face_fan(pmp)
   return PolyhedralFan(pmff)
end
