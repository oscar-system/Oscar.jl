##############################################################################
"""
   save_cone(Cone, String)

Save a cone to a file in JSON format. The first argument is the cone, the
second argument is the filename.
"""
function save_cone(C::Cone, filename::String)
   bigobject = pm_cone(C)
   Polymake.save_bigobject(bigobject, filename)
end

"""
   load_cone(String)

Load a cone stored in JSON format, given the filename as input.
"""
function load_cone(filename::String)
   bigobject = Polymake.load_bigobject(filename)
   if typename[1:4] != "Cone"
      throw(ArgumentError("Loaded object is not of polymake type Cone, it has type " * typename))
   end
   return Cone(bigobject)
end

##############################################################################
"""
   save_polyhedron(Cone, String)

Save a polyhedron to a file in JSON format. The first argument is the
polyhedron, the second argument is the filename.
"""
function save_polyhedron(P::Polyhedron, filename::String)
   bigobject = pm_polytope(P)
   Polymake.save_bigobject(bigobject, filename)
end

"""
   load_polyhedron(String)

Load a polyhedron stored in JSON format, given the filename as input.
"""
function load_polyhedron(filename::String)
   bigobject = Polymake.load_bigobject(filename)
   typename = Polymake.type_name(bigobject)
   if typename[1:8] != "Polytope"
      throw(ArgumentError("Loaded object is not of polymake type Polytope, it has type " * typename))
   end
   return Polyhedron(bigobject)
end

##############################################################################
"""
   save_polyhedralfan(PolyhedralFan, String)

Save a polyhedral fan to a file in JSON format. The first argument is the
polyhedral fan, the second argument is the filename.
"""
function save_polyhedralfan(PF::PolyhedralFan, filename::String)
   bigobject = pm_ralfan(PF)
   Polymake.save_bigobject(bigobject, filename)
end

"""
   load_polyhedralfan(String)

Load a polyhedral fan stored in JSON format, given the filename as input.
"""
function load_polyhedralfan(filename::String)
   bigobject = Polymake.load_bigobject(filename)
   if typename[1:13] != "PolyhedralFan"
      throw(ArgumentError("Loaded object is not of polymake type PolyhedralFan, it has type " * typename))
   end
   return PolyhedralFan(bigobject)
end
