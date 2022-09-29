# Loading polymake objects into the appropriate types.


##############################################################################
"""
    load_from_polymake(::Type{Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Directed, Undirected}}

Load a graph stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Directed, Undirected}}
    polymake_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    return Graph{T}(polymake_object)
end


const polymake2OscarTypes = Dict{String, Type}([
    "polytope::Cone<Rational>" => Cone{fmpq},
    "polytope::Polytope<Rational>" => Polyhedron{fmpq},
    "fan::PolyhedralFan<Rational>" => PolyhedralFan{fmpq},
    "fan::PolyhedralComplex<Rational>" => PolyhedralComplex{fmpq},
    "fan::SubdivisionOfPoints<Rational>" => SubdivisionOfPoints{fmpq},
    "topaz::SimplicialComplex" => SimplicialComplex,
    "common::GraphAdjacency<Undirected>" => Graph{Undirected},
    "common::GraphAdjacency<Directed>" => Graph{Directed},
])

@registerSerializationType(Polymake.BigObjectAllocated, "Polymake.BigObject")

function load_from_polymake(::Type{T}, jsondict::Dict{Symbol, Any}) where {
        S<:scalar_types,
        T<:Union{Cone{S}, Polyhedron{S}, PolyhedralFan{S}, 
            PolyhedralComplex{S}, SubdivisionOfPoints{S}, SimplicialComplex}}
    
    inner_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    return T(inner_object)
end


# Distinguish between the various polymake datatypes.
function load_from_polymake(jsondict::Dict{Symbol, Any})
    typename = jsondict[:_type]
    if haskey(polymake2OscarTypes, typename)
        oscar_type = polymake2OscarTypes[typename]
        return load_from_polymake(oscar_type, jsondict)
    else 
        # We just try to default to something from Polymake.jl
        deserialized = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
        if !isa(deserialized, Polymake.BigObject)
            @warn "No function for converting the deserialized Polymake type to Oscar type: $(typeof(deserialized))"
            return deserialized
        end
        
        if Polymake.type_name(deserialized) == "Ideal"
            return convert(MPolyIdeal{fmpq_mpoly}, deserialized)
        else
            @warn "No function for converting the deserialized Polymake type to Oscar type: $(typeof(deserialized))"
            return deserialized
        end
    end
end

