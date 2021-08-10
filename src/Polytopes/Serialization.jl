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
   typename = Polymake.type_name(bigobject)
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
   bigobject = pm_fan(PF)
   Polymake.save_bigobject(bigobject, filename)
end

"""
    load_polyhedralfan(String)

Load a polyhedral fan stored in JSON format, given the filename as input.
"""
function load_polyhedralfan(filename::String)
   bigobject = Polymake.load_bigobject(filename)
   typename = Polymake.type_name(bigobject)
   if typename[1:13] != "PolyhedralFan"
      throw(ArgumentError("Loaded object is not of polymake type PolyhedralFan, it has type " * typename))
   end
   return PolyhedralFan(bigobject)
end


##############################################################################
"""
    save_linearprogram(LinearProgram, String)

Save a cone to a file in JSON format. The first argument is the cone, the
second argument is the filename.
"""
function save_linearprogram(LP::LinearProgram, filename::String)
   bigobject = pm_polytope(feasible_region(LP))
   Polymake.save_bigobject(bigobject, filename)
end

"""
    load_linearprogram(String)

Load a cone stored in JSON format, given the filename as input.
"""
function load_linearprogram(filename::String)
   fr = Polymake.load_bigobject(filename)
   typename = Polymake.type_name(fr)
   if typename[1:8] != "Polytope"
      throw(ArgumentError("Loaded object is not of polymake type LinearProgram."))
   end
   if !Polymake.exists(fr, "LP")
      throw(ArgumentError("Loaded object is not of polymake type LinearProgram."))
   end
   lp = fr.LP
   conv = Polymake.get_attachment(lp, "convention")
   if isnothing(conv)
      throw(ArgumentError("Loaded object is not of polymake type LinearProgram."))
   end
   if conv == "max"
       return LinearProgram(Polyhedron(fr), lp, :max)
   elseif conv == "min"
       return LinearProgram(Polyhedron(fr), lp, :min)
   end
end


##############################################################################
"""
    save_subdivisionofpoints(SubdivisionOfPoints, String)

Save a subdivision of points to a file in JSON format. The first argument is
the subdivision of points, the second argument is the filename.
"""
function save_subdivisionofpoints(SOP::SubdivisionOfPoints, filename::String)
   bigobject = pm_subdivision(SOP)
   Polymake.save_bigobject(bigobject, filename)
end

"""
    load_subdivisionofpoints(String)

Load a subdivision of points stored in JSON format, given the filename as input.
"""
function load_subdivisionofpoints(filename::String)
   bigobject = Polymake.load_bigobject(filename)
   typename = Polymake.type_name(bigobject)
   if typename[1:19] != "SubdivisionOfPoints"
      throw(ArgumentError("Loaded object is not of polymake type SubdivisionOfPoints, it has type " * typename))
   end
   return SubdivisionOfPoints(bigobject)
end
