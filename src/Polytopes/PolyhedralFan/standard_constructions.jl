#TODO: inward/outward options? via polymake changes?

"""
   normal_fan(P::Polyhedron)

Returns the normal fan of a polyhedron.
"""
function normal_fan(P::Polyhedron)
   pmp = pm_polytope(P)
   pmnf = Polymake.fan.normal_fan(pmp)
   return PolyhedralFan(pmnf)
end

"""
   face_fan(P::Polyhedron)

Returns the face fan of a polyhedron.
"""
function face_fan(P::Polyhedron)
   pmp = pm_polytope(P)
   pmff = Polymake.fan.face_fan(pmp)
   return PolyhedralFan(pmff)
end

