###############################################################################
## Graphs
###############################################################################

function save_intern(s::SerializerState, g::Graphs.Graph{T}) where {T <: Union{Graphs.Directed, Graphs.Undirected}}
    smallobject = pm_object(g)
    serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, smallobject)
    jsonstr = Polymake.call_function(:common, :encode_json, serialized)
    return JSON.parse(jsonstr)
end


function load_intern(s::DeserializerState, g::Type{Graphs.Graph{T}}, dict::Dict) where {T <: Union{Graphs.Directed, Graphs.Undirected}}
    smallobj = Polymake.call_function(:common, :deserialize_json_string, json(dict))
    return g(smallobj)
end


###############################################################################
## SimplicialComplex
###############################################################################
function save_intern(s::SerializerState, K::SimplicialComplex)
    bo = pm_object(K)
    return bigobject_to_dict(bo)
end

function load_intern(s::DeserializerState, K::Type{SimplicialComplex}, dict::Dict)
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(dict))
    return K(bigobject)
end
