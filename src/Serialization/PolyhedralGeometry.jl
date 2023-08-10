
using JSON

function bigobject_to_jsonstr(bo::Polymake.BigObject)
    serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, bo)
    return Polymake.call_function(:common, :encode_json, serialized)
end

function bigobject_to_dict(bo::Polymake.BigObject)
    jsonstr = bigobject_to_jsonstr(bo)
    return JSON.parse(jsonstr)
end

function save_object(s::SerializerState, p::Polymake.BigObject)
    data_json(s, bigobject_to_jsonstr(p))
end

function load_object(s::DeserializerState, ::Type{Polymake.BigObjectAllocated}, dict::Dict)
    return load_from_polymake(dict)
end

function load_object(s::DeserializerState, ::Type{Polymake.BigObject}, dict::Dict)
    bigobject = Polymake.call_function(:common, :deserialize_json_string, json(dict))
    return bigobject
end

function load_object(s::DeserializerState, ::Type{Polymake.BigObject}, str::String)
    return load_ref(s, str)
end

##############################################################################
# Abstract Polyhedral Object
type_needs_params(::Type{<:PolyhedralObject}) = true

function save_type_params(s::SerializerState, obj::T, key::Symbol) where T <: PolyhedralObject
    s.key = key
    data_dict(s) do
        save_object(s, encode_type(T), :name)
        save_typed_object(s, coefficient_field(obj), :params)
    end
end

function save_object(s::SerializerState, obj::PolyhedralObject)
    save_object(s, pm_object(obj))
end

function load_type_params(s::DeserializerState, ::Type{<:PolyhedralObject}, dict)
    return load_typed_object(s, dict)
end

function load_object_with_params(s::DeserializerState, T::Type{<:PolyhedralObject},
                                 dict::Dict, params::Any)
    pm_obj = load_from_polymake(dict)
    return T(pm_obj, params)
end
##############################################################################
@registerSerializationType(LinearProgram{QQFieldElem})

function save_internal(s::SerializerState, lp::LinearProgram)
    lpcoeffs = lp.polymake_lp.LINEAR_OBJECTIVE
    serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, lpcoeffs)
    jsonstr = Polymake.call_function(:common, :encode_json, serialized)
    save_type_dispatch(s, lp.feasible_region, :feasible_region)
    save_type_dispatch(s, lp.convention, :convention)
    save_type_dispatch(s, JSON.parse(jsonstr), :lpcoeffs)
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

##############################################################################
@registerSerializationType(MixedIntegerLinearProgram{QQFieldElem})

function save_internal(s::SerializerState, milp::MixedIntegerLinearProgram)
    milp_coeffs = milp.polymake_milp.LINEAR_OBJECTIVE
    int_vars = milp.polymake_milp.INTEGER_VARIABLES
    coeffs_serialized = Polymake.call_function(
        Symbol("Core::Serializer"), :serialize, milp_coeffs)
    int_vars_serialized = Polymake.call_function(
        Symbol("Core::Serializer"), :serialize, int_vars)
    coeffs_jsonstr = Polymake.call_function(:common, :encode_json, coeffs_serialized)
    int_vars_jsonstr = Polymake.call_function(:common, :encode_json, int_vars_serialized)

    save_type_dispatch(s, milp.feasible_region, :feasible_region)
    save_type_dispatch(s, milp.convention, :convention)
    save_type_dispatch(s, JSON.parse(coeffs_jsonstr), :milp_coeffs)
    save_type_dispatch(s, JSON.parse(int_vars_jsonstr), :int_vars)
end

function load_internal(s::DeserializerState, ::Type{MixedIntegerLinearProgram{T}}, dict::Dict) where T
    fr = load_type_dispatch(s, Polyhedron{T}, dict[:feasible_region])
    conv = dict[:convention]
    milp_coeffs = Polymake.call_function(
        :common,
        :deserialize_json_string,
        json(dict[:milp_coeffs])
    )
    int_vars = Polymake.call_function(
        :common,
        :deserialize_json_string,
        json(dict[:int_vars])
    )

    all = Polymake._lookup_multi(pm_object(fr), "MILP")
    index = 0
    for i in 1:length(all)
        if all[i].LINEAR_OBJECTIVE == milp_coeffs && all[i].INTEGER_VARIABLES == int_vars
            index = i
            break
        end
    end
    lp = Polymake._lookup_multi(pm_object(fr), "MILP", index-1)
    return MixedIntegerLinearProgram{T}(fr, lp, Symbol(conv))
end

# use generic serialization for the other types:
@registerSerializationType(Cone)



@registerSerializationType(PolyhedralComplex{QQFieldElem})
save_internal(s::SerializerState, obj::PolyhedralComplex) = save_internal_generic(s, obj)
load_internal(s::DeserializerState, ::Type{T}, dict::Dict) where T <: PolyhedralComplex = load_internal_generic(s, T, dict)

@registerSerializationType(Polyhedron{QQFieldElem})
save_internal(s::SerializerState, obj::Polyhedron) = save_internal_generic(s, obj)
load_internal(s::DeserializerState, ::Type{T}, dict::Dict) where T <: Polyhedron = load_internal_generic(s, T, dict)

@registerSerializationType(PolyhedralFan{QQFieldElem})
save_internal(s::SerializerState, obj::PolyhedralFan) = save_internal_generic(s, obj)
load_internal(s::DeserializerState, ::Type{T}, dict::Dict) where T <: PolyhedralFan = load_internal_generic(s, T, dict)

@registerSerializationType(SubdivisionOfPoints{QQFieldElem})
save_internal(s::SerializerState, obj::SubdivisionOfPoints) = save_internal_generic(s, obj)
load_internal(s::DeserializerState, ::Type{T}, dict::Dict) where T <: SubdivisionOfPoints = load_internal_generic(s, T, dict)
