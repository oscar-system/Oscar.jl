using JSON

###############################################################################
## Graphs
###############################################################################
@registerSerializationType(Graphs.Graph{Graphs.Directed}, "Oscar.Graphs.Graph{Oscar.Graphs.Directed}")
@registerSerializationType(Graphs.Graph{Graphs.Undirected}, "Oscar.Graphs.Graph{Oscar.Graphs.Undirected}")

function save_internal(s::SerializerState, g::Graphs.Graph{T}) where {T <: Union{Graphs.Directed, Graphs.Undirected}}
    smallobject = pm_object(g)
    serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, smallobject)
    jsonstr = Polymake.call_function(:common, :encode_json, serialized)
    return JSON.parse(jsonstr)
end


function load_internal(s::DeserializerState, g::Type{Graphs.Graph{T}}, dict::Dict) where {T <: Union{Graphs.Directed, Graphs.Undirected}}
    smallobj = Polymake.call_function(:common, :deserialize_json_string, json(dict))
    return g(smallobj)
end

###############################################################################
## IncidenceMatrix
###############################################################################
@registerSerializationType(Polymake.IncidenceMatrixAllocated{Polymake.NonSymmetric})

function save_internal(s::SerializerState, IM::IncidenceMatrix)
    serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, IM)
    jsonstr = Polymake.call_function(:common, :encode_json, serialized)

    return JSON.parse(jsonstr)
end

function load_internal(s::DeserializerState, ::Type{<: IncidenceMatrix}, dict::Dict)
    IM = Polymake.call_function(:common, :deserialize_json_string, json(dict))
    return IM
end


###############################################################################
## SimplicialComplex
###############################################################################
@registerSerializationType(SimplicialComplex)

function save_internal(s::SerializerState, K::SimplicialComplex)
    bo = pm_object(K)
    return bigobject_to_dict(bo)
end

function load_internal(s::DeserializerState, K::Type{SimplicialComplex}, dict::Dict)
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(dict))
    return K(bigobject)
end
