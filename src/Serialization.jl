export load_linearprogram,
   load,
   save



##############################################################################
"""
    save(Cone, String)

Save a cone to a file in JSON format. The first argument is the cone, the
second argument is the filename.
"""
function save(C::Cone, filename::String)
   bigobject = pm_cone(C)
   Polymake.save_bigobject(bigobject, filename)
end

##############################################################################
"""
    save(Cone, String)

Save a polyhedron to a file in JSON format. The first argument is the
polyhedron, the second argument is the filename.
"""
function save(P::Polyhedron, filename::String)
   bigobject = pm_polytope(P)
   Polymake.save_bigobject(bigobject, filename)
end

##############################################################################
"""
    save(PolyhedralFan, String)

Save a polyhedral fan to a file in JSON format. The first argument is the
polyhedral fan, the second argument is the filename.
"""
function save(PF::PolyhedralFan, filename::String)
   bigobject = pm_fan(PF)
   Polymake.save_bigobject(bigobject, filename)
end

##############################################################################
"""
    save(LinearProgram, String)

Save a linear program to a file in JSON format. The first argument is the
linear program, the second argument is the filename.
"""
function save(LP::LinearProgram, filename::String)
   bigobject = pm_polytope(feasible_region(LP))
   Polymake.attach(bigobject, "oscar_type", "LinearProgram")
   Polymake.save_bigobject(bigobject, filename)
end

function _load(::Type{LinearProgram}, fr::Polymake.BigObject)
   Polymake.remove_attachment(p_obj, "oscar_type")
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
    save(SubdivisionOfPoints, String)

Save a subdivision of points to a file in JSON format. The first argument is
the subdivision of points, the second argument is the filename.
"""
function save(SOP::SubdivisionOfPoints, filename::String)
   bigobject = pm_subdivision(SOP)
   Polymake.save_bigobject(bigobject, filename)
end

##############################################################################

# Unify loading

function _load(p_obj::Polymake.BigObject)
  typename = Polymake.type_name(p_obj)
  if typename[1:4] == "Cone"
    return Cone(p_obj)
  elseif typename[1:8] == "Polytope"
    return Polyhedron(p_obj)
  elseif typename[1:13] == "PolyhedralFan"
    return PolyhedralFan(p_obj)
  elseif typename[1:19] == "SubdivisionOfPoints"
    return SubdivisionOfPoints(p_obj)
  end
  throw(ArgumentError("Polymake type not embedded in Oscar."))
end

function load(filename::String)
  p_obj = Polymake.load_bigobject(filename)
  o_type = Polymake.get_attachment(p_obj, "oscar_type")
  if isnothing(o_type)
    return _load(p_obj)
  else
    return _load(eval(Meta.parse(o_type)), p_obj)
  end
end