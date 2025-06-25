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

function load_object(s::DeserializerState, ::Type{Polymake.BigObjectAllocated})
  dict = Dict{Symbol, Any}(s.obj)
  return load_from_polymake(dict)
end

function load_object(s::DeserializerState, ::Type{Polymake.BigObject})
  dict = Dict{Symbol, Any}(s.obj)
  bigobject = Polymake.call_function(:common, :deserialize_json_string, JSON3.write(dict))
  return bigobject
end

function load_object(s::DeserializerState, ::Type{Polymake.BigObject}, str::String)
  return load_ref(s, str)
end

##############################################################################
# Abstract Polyhedral Object

type_params(obj::T) where {S <: Union{QQFieldElem, Float64}, T <: PolyhedralObject{S}} = TypeParams(T, coefficient_field(obj))

function type_params(obj::T) where {S, T <: PolyhedralObject{S}}
  p_dict = _polyhedral_object_as_dict(obj)
  field = p_dict[:_coeff]
  delete!(p_dict, :_coeff)
  return TypeParams(
    T,
    :field => field,
    :pm_params => type_params(p_dict))
end

function save_object(s::SerializerState, obj::PolyhedralObject{S}) where S <: Union{QQFieldElem, Float64}
  save_object(s, pm_object(obj))
end

function save_object(s::SerializerState, obj::PolyhedralObject{<:FieldElem})
  if typeof(obj) <: Union{MixedIntegerLinearProgram, LinearProgram}
    T = typeof(obj)
    error("Unsupported type $T for serialization")
  end
  p_dict = _polyhedral_object_as_dict(obj)
  delete!(p_dict, :_coeff)
  save_data_dict(s) do
    for (k, v) in p_dict
      if !Base.issingletontype(typeof(v))
        save_object(s, v, k)
      end
    end
  end
end

function load_object(s::DeserializerState, T::Type{<:PolyhedralObject},
                     field::U) where {U <: Union{QQField, AbstractAlgebra.Floats}}
  return load_from_polymake(T{elem_type(field)}, Dict{Symbol, Any}(s.obj))
end

function load_object(s::DeserializerState, T::Type{<:PolyhedralObject{S}},
    field::U) where {S <: Union{QQFieldElem, Float64}, U <: Union{QQField, AbstractAlgebra.Floats}}
  return load_from_polymake(T, Dict{Symbol, Any}(s.obj))
end

function load_object(s::DeserializerState, T::Type{<:PolyhedralObject}, dict::Dict)
  field = dict[:field]
  polymake_dict = load_object(s, Dict{Symbol, Any}, dict[:pm_params])
  bigobject = _dict_to_bigobject(polymake_dict)

  return T{elem_type(field)}(bigobject, field)
end

function load_object(s::DeserializerState, T::Type{<:PolyhedralObject{S}},
                     dict::Dict) where S <: FieldElem
  field = dict[:field]
  polymake_dict = load_object(s, Dict{Symbol, Any}, dict[:pm_params])
  bigobject = _dict_to_bigobject(polymake_dict)

  return T(bigobject, field)
end

##############################################################################
@register_serialization_type LinearProgram

function save_object(s::SerializerState, lp::LinearProgram{QQFieldElem})
  lpcoeffs = lp.polymake_lp.LINEAR_OBJECTIVE
  serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, lpcoeffs)
  jsonstr = Polymake.call_function(:common, :encode_json, serialized)
  save_data_dict(s) do
    save_object(s, lp.feasible_region, :feasible_region)
    save_object(s, lp.convention, :convention)
    save_json(s, jsonstr, :lpcoeffs)
  end
end

function save_object(s::SerializerState{<: LPSerializer}, lp::LinearProgram{QQFieldElem})
  lp_filename = basepath(s.serializer) * "-$(objectid(lp)).lp"
  save_lp(lp_filename, lp)

  save_object(s, basename(lp_filename))
end

function load_object(s::DeserializerState, ::Type{<:LinearProgram}, field::QQField)
  if s.obj isa String
    error("Loading this file requires using the LPSerializer")
  end
  coeff_type = elem_type(field)
  fr = load_object(s, Polyhedron, field, :feasible_region)
  conv = load_object(s, String, :convention)
  lpcoeffs = load_node(s, :lpcoeffs) do lpcoeffs
    Polymake.call_function(:common, :deserialize_json_string, json(lpcoeffs))
  end
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

function load_object(s::DeserializerState{LPSerializer},
                     ::Type{<:LinearProgram}, field::QQField)
  load_node(s) do _
    lp_filename = dirname(basepath(s.serializer)) * "/$(s.obj)"
    pm_lp = load_lp(lp_filename)
  end
end

##############################################################################
@register_serialization_type MixedIntegerLinearProgram

function save_object(s::SerializerState, milp::MixedIntegerLinearProgram{QQFieldElem})
  milp_coeffs = milp.polymake_milp.LINEAR_OBJECTIVE
  int_vars = milp.polymake_milp.INTEGER_VARIABLES
  coeffs_serialized = Polymake.call_function(
    Symbol("Core::Serializer"), :serialize, milp_coeffs)
  int_vars_serialized = Polymake.call_function(
    Symbol("Core::Serializer"), :serialize, int_vars)
  coeffs_jsonstr = Polymake.call_function(:common, :encode_json, coeffs_serialized)
  int_vars_jsonstr = Polymake.call_function(:common, :encode_json, int_vars_serialized)
  save_data_dict(s) do
    save_object(s, milp.feasible_region, :feasible_region)
    save_object(s, milp.convention, :convention)
    save_json(s, coeffs_jsonstr, :milp_coeffs)
    save_json(s, int_vars_jsonstr, :int_vars)
  end
end

function load_object(s::DeserializerState, ::Type{<: MixedIntegerLinearProgram}, field::QQField)
  fr = load_object(s, Polyhedron, field, :feasible_region)
  conv = load_object(s, String, :convention)
  milp_coeffs = load_node(s, :milp_coeffs) do coeffs
    Polymake.call_function(
      :common,
      :deserialize_json_string,
      json(coeffs)
    )
  end
  int_vars = load_node(s, :int_vars) do vars
    Polymake.call_function(
      :common,
      :deserialize_json_string,
      json(vars)
    )
  end

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
  return MixedIntegerLinearProgram{T}(fr, lp, Symbol(conv), field)
end

# use generic serialization for the other types:
@register_serialization_type Cone
@register_serialization_type PolyhedralComplex
@register_serialization_type Polyhedron
@register_serialization_type PolyhedralFan
@register_serialization_type SubdivisionOfPoints
