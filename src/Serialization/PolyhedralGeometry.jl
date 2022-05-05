using JSON

function bigobject_to_dict(bo::Polymake.BigObject)
    serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, bo)
    jsonstr = Polymake.call_function(:common, :encode_json, serialized)
    return JSON.parse(jsonstr)
end

function save_internal(s::SerializerState, p::Polymake.BigObject)
    return bigobject_to_dict(p)
end

function load_internal(s::DeserializerState, ::Type{Polymake.BigObject}, dict::Dict)
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(dict))
    return bigobject
end


##############################################################################
function save_internal(s::SerializerState, lp::LinearProgram)
    lpcoeffs = lp.polymake_lp.LINEAR_OBJECTIVE
    serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, lpcoeffs)
    jsonstr = Polymake.call_function(:common, :encode_json, serialized)
    return Dict(
        :feasible_region => save_type_dispatch(s, lp.feasible_region),
        :convention => lp.convention,
        :lpcoeffs => JSON.parse(jsonstr)
    )
end

function load_internal(s::DeserializerState, ::Type{LinearProgram{T}}, dict::Dict) where T
    fr = load_type_dispatch(s, Polyhedron{T}, dict[:feasible_region])
    conv = dict[:convention]
    lpcoeffs = Polymake.call_function(:common, :deserialize_json_string, json(dict[:lpcoeffs]))
    all = Polymake._lookup_multi(pm_object(fr), "LP")
    index = 0
    for i in 1:length(all)
        if all[i].LINEAR_OBJECTIVE == lpcoeffs
            index = i
            break
        end
    end
    lp = Polymake._lookup_multi(pm_object(fr), "LP", index-1)
    return LinearProgram{T}(fr, lp, Symbol(conv))
end


