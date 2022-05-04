# Loading polymake objects into the appropriate types.

##############################################################################
"""
    load_from_polymake(::Type{Cone}, jsondict::Dict{Symbol, Any})

Load a cone stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{Cone}, jsondict::Dict{Symbol, Any})
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    typename = Polymake.type_name(bigobject)
    if typename != "Cone<Rational>"
        throw(ArgumentError("Loaded object is not of polymake type Cone, it has type " * typename))
    end
    return Cone(bigobject)
end


##############################################################################
"""
    load_from_polymake(::Type{Polyhedron}, jsondict::Dict{Symbol, Any})

Load a polyhedron stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{Polyhedron}, jsondict::Dict{Symbol, Any})
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    typename = Polymake.type_name(bigobject)
    if typename != "Polytope<Rational>"
        throw(ArgumentError("Loaded object is not of polymake type Polytope, it has type " * typename))
    end
    return Polyhedron(bigobject)
end



##############################################################################

"""
    load_from_polymake(::Type{PolyhedralFan}, jsondict::Dict{Symbol, Any})

Load a polyhedral fan stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{PolyhedralFan}, jsondict::Dict{Symbol, Any})
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    typename = Polymake.type_name(bigobject)
    if typename != "PolyhedralFan<Rational>"
        throw(ArgumentError("Loaded object is not of polymake type PolyhedralFan, it has type " * typename))
    end
    return PolyhedralFan(bigobject)
end




##############################################################################

"""
    load_from_polymake(::Type{LinearProgram}, jsondict::Dict{Symbol, Any})

Load a linear program stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{LinearProgram}, jsondict::Dict{Symbol, Any})
    fr = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    typename = Polymake.type_name(fr)
    if typename != "Polytope<Rational>"
        throw(ArgumentError("Loaded object is not of polymake type LinearProgram."))
    end
    if !Polymake.exists(fr, "LP")
        throw(ArgumentError("Incomplete or corrupt data file."))
    end
    lp = fr.LP
    conv = Polymake.get_attachment(lp, "convention")
    if isnothing(conv)
        throw(ArgumentError("Incomplete or corrupt data file."))
    end
    if conv == "max"
        return LinearProgram(Polyhedron(fr), lp, :max)
    elseif conv == "min"
        return LinearProgram(Polyhedron(fr), lp, :min)
    end
end






##############################################################################
"""
    load_from_polymake(::Type{SubdivisionOfPoints}, jsondict::Dict{Symbol, Any})

Load a subdivision of points stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{SubdivisionOfPoints}, jsondict::Dict{Symbol, Any})
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    typename = Polymake.type_name(bigobject)
    if typename != "SubdivisionOfPoints<Rational>"
        throw(ArgumentError("Loaded object is not of polymake type SubdivisionOfPoints, it has type " * typename))
    end
    return SubdivisionOfPoints(bigobject)
end



##############################################################################

"""
    load_from_polymake(::Type{PolyhedralComplex}, jsondict::Dict{Symbol, Any})

Load a polyhedral complex stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{PolyhedralComplex}, jsondict::Dict{Symbol, Any})
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    typename = Polymake.type_name(bigobject)
    if typename != "PolyhedralComplex<Rational>"
        throw(ArgumentError("Loaded object is not of polymake type PolyhedralComplex, it has type " * typename))
    end
    return PolyhedralComplex(bigobject)
end


##############################################################################
"""
    load_from_polymake(::Type{Graphs.Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Graphs.Directed, Graphs.Undirected}}

Load a graph stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{Graphs.Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Graphs.Directed, Graphs.Undirected}}
    polymake_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    return Graphs.Graph{T}(polymake_object)
end




# Distinguish between the various polymake datatypes.
function load_from_polymake(jsondict::Dict{Symbol, Any})
    typename = jsondict[:_type]
    if typename == "polytope::Cone<Rational>"
        return load_from_polymake(Cone, jsondict)
    elseif typename == "polytope::Polytope<Rational>"
        if haskey(jsondict, "is_OSCAR_LP")
            load_from_polymake(LinearProgram, jsondict)
        else
            load_from_polymake(Polyhedron, jsondict)
        end
    elseif typename == "fan::PolyhedralFan<Rational>"
        return load_from_polymake(PolyhedralFan, jsondict)
    elseif typename == "fan::PolyhedralComplex<Rational>"
        return load_from_polymake(PolyhedralComplex, jsondict)
    elseif typename == "fan::SubdivisionOfPoints<Rational>"
        return load_from_polymake(SubdivisionOfPoints, jsondict)
    elseif typename == "topaz::SimplicialComplex"
        return load_from_polymake(SimplicialComplex, jsondict)
    elseif typename == "common::GraphAdjacency<Undirected>"
        return load_from_polymake(Graphs.Graph{Graphs.Undirected}, jsondict)
    elseif typename == "common::GraphAdjacency<Directed>"
        return load_from_polymake(Graphs.Graph{Graphs.Directed}, jsondict)
    else 
        # We just try to default to something from Polymake.jl
        return Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    end
end

