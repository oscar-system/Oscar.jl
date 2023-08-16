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
  save_json(s, bigobject_to_jsonstr(p))
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

function save_type_params(s::SerializerState, obj::T) where T <: PolyhedralObject
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

function load_object(s::DeserializerState, T::Type{<:PolyhedralObject},
                     dict::Dict, field::Field) 
  return  load_from_polymake(T{elem_type(field)}, dict)
end

function load_object(s::DeserializerState, T::Type{<:PolyhedralObject{S}},
                     dict::Dict, field::Field) where S <: FieldElem
  return load_from_polymake(T, dict)
end

##############################################################################
@registerSerializationType(LinearProgram)

function save_object(s::SerializerState, lp::LinearProgram)
  lpcoeffs = lp.polymake_lp.LINEAR_OBJECTIVE
  serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, lpcoeffs)
  jsonstr = Polymake.call_function(:common, :encode_json, serialized)
  data_dict(s) do 
    save_object(s, lp.feasible_region, :feasible_region)
    save_object(s, lp.convention, :convention)
    save_json(s, jsonstr, :lpcoeffs)
  end
end

function load_object(s::DeserializerState, ::Type{<:LinearProgram},
                                 dict::Dict, field::Field)
  coeff_type = elem_type(field)
  fr = load_object(s, Polyhedron, dict[:feasible_region], field)
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
  return LinearProgram{coeff_type}(fr, lp, Symbol(conv))
end

##############################################################################
@registerSerializationType(MixedIntegerLinearProgram)

function save_object(s::SerializerState, milp::MixedIntegerLinearProgram)
  milp_coeffs = milp.polymake_milp.LINEAR_OBJECTIVE
  int_vars = milp.polymake_milp.INTEGER_VARIABLES
  coeffs_serialized = Polymake.call_function(
    Symbol("Core::Serializer"), :serialize, milp_coeffs)
  int_vars_serialized = Polymake.call_function(
    Symbol("Core::Serializer"), :serialize, int_vars)
  coeffs_jsonstr = Polymake.call_function(:common, :encode_json, coeffs_serialized)
  int_vars_jsonstr = Polymake.call_function(:common, :encode_json, int_vars_serialized)
  data_dict(s) do
    save_object(s, milp.feasible_region, :feasible_region)
    save_object(s, milp.convention, :convention)
    save_json(s, coeffs_jsonstr, :milp_coeffs)
    save_json(s, int_vars_jsonstr, :int_vars)
  end
end

function load_object(s::DeserializerState, ::Type{<: MixedIntegerLinearProgram},
                                 dict::Dict, field::Field) 
  fr = load_object(s, Polyhedron, dict[:feasible_region], field)
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
  T = elem_type(field)
  return MixedIntegerLinearProgram{T}(fr, lp, Symbol(conv))
end

# use generic serialization for the other types:
@registerSerializationType(Cone)
@registerSerializationType(PolyhedralComplex)
@registerSerializationType(Polyhedron)
@registerSerializationType(PolyhedralFan)
@registerSerializationType(SubdivisionOfPoints)
