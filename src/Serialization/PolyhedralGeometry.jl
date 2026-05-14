import JSON

function bigobject_to_jsonstr(bo::Polymake.BigObject)
  serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, bo)
  return Polymake.call_function(:common, :encode_json, serialized)
end

# unused, is this needed?
function bigobject_to_dict(bo::Polymake.BigObject)
  jsonstr = bigobject_to_jsonstr(bo)
  return JSON.parse(jsonstr)
end

function save_object(s::SerializerState, p::Polymake.BigObject)
  save_json(s, bigobject_to_jsonstr(p))
end

function load_object(s::DeserializerState, ::Type{Polymake.BigObjectAllocated})
  dict = load_json(s, Dict{String, Any})
  return load_from_polymake(dict)
end

function load_object(s::DeserializerState, ::Type{Polymake.BigObject})
  dict = load_json(s, Dict{String, Any})
  bigobject = Polymake.call_function(:common, :deserialize_json_string, JSON.json(dict))
  return bigobject
end

##############################################################################
# Abstract Polyhedral Object

type_params(obj::T) where {T <: PolyhedralObject} = TypeParams(T, coefficient_field(obj))
type_params(obj::T) where {T <: Union{LinearProgram, MixedIntegerLinearProgram}} = TypeParams(T, coefficient_field(obj))


function save_object(s::SerializerState, obj::PolyhedralObject{S}) where S <: Union{QQFieldElem, Float64}
  save_object(s, pm_object(obj))
end

function save_object(s::SerializerState, obj::PolyhedralObject{<:FieldElem})
  p_dict = _polyhedral_object_as_dict(obj)
  delete!(p_dict, "_coeff")
  save_data_dict(s) do
    for (k, v) in p_dict
      if !Base.issingletontype(typeof(v))
        save_typed_object(s, v, Symbol(k))
      end
    end
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:PolyhedralObject, <:Union{QQField, AbstractAlgebra.Floats}})
  field = Oscar.params(tp)
  return load_from_polymake(tp.type{elem_type(field)}, load_json(s, Dict{String, Any}))
end

function load_object(s::DeserializerState, tp::TypeParams{<:PolyhedralObject{S}, <:Union{QQField, AbstractAlgebra.Floats}}) where S <: Union{QQFieldElem, Float64}
  return load_from_polymake(tp.type, load_json(s, Dict{String, Any}))
end

function load_object(s::DeserializerState, tp::TypeParams{T, <:Field}) where T <: PolyhedralObject
  field = Oscar.params(tp)
  polymake_dict = Dict{Symbol, Any}()
  for k in propertynames(s.obj)
    polymake_dict[k] = load_typed_object(s, k)
  end
  bigobject = _dict_to_bigobject(polymake_dict)
  T_concrete = T isa UnionAll ? T{elem_type(field)} : T
  return T_concrete(bigobject, field)
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

function save_object(s::SerializerState, lp::LinearProgram{<:FieldElem})
  lpcoeffs = lp.polymake_lp.LINEAR_OBJECTIVE
  save_data_dict(s) do
    save_object(s, lp.feasible_region, :feasible_region)
    save_object(s, lp.convention, :convention)
    save_object(s, _pmdata_for_oscar(lpcoeffs, coefficient_field(lp)), :lpcoeffs)
  end
end

function save_object(s::SerializerState{<: LPSerializer}, lp::LinearProgram{QQFieldElem})
  lp_filename = basepath(s.serializer) * "-$(objectid(lp)).lp"
  save_lp(lp_filename, lp)

  save_object(s, basename(lp_filename))
end

function load_object(s::DeserializerState, tp::TypeParams{<:LinearProgram, QQField})
  if is_string(s)
    error("Loading this file requires using the LPSerializer")
  end
  field = Oscar.params(tp)
  coeff_type = elem_type(field)
  fr = load_object(s, TypeParams(Polyhedron, field), :feasible_region)
  conv = load_object(s, String, :convention)
  lpcoeffs = load_node(s, :lpcoeffs) do _
    Polymake.call_function(:common, :deserialize_json_string, JSON.json(load_json(s, Dict{String, Any})))
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

function load_object(s::DeserializerState, tp::TypeParams{<:LinearProgram, <:Field})
  if is_string(s)
    error("Loading this file requires using the LPSerializer")
  end
  field = Oscar.params(tp)
  coeff_type = elem_type(field)
  fr = load_object(s, TypeParams(Polyhedron, field), :feasible_region)
  conv = load_object(s, String, :convention)
  lpcoeffs = load_object(s, TypeParams(Vector{coeff_type}, field), :lpcoeffs)
  all = Polymake._lookup_multi(pm_object(fr), "LP")
  lp = nothing
  for i in 1:length(all)
    lo = _pmdata_for_oscar(all[i].LINEAR_OBJECTIVE, field)
    if lpcoeffs == lo
      lp = all[i]
      break
    end
  end
  @req lp !== nothing "could not identify LP subobject"
  return LinearProgram{coeff_type}(fr, lp, Symbol(conv), field)
end

function load_object(s::DeserializerState{LPSerializer},
                     tp::TypeParams{<:LinearProgram, QQField})
  load_node(s) do _
    lp_filename = dirname(basepath(s.serializer)) * "/$(load_json(s, String))"
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
    # we can probably changes this line to just store a vector{Int} but will need upgrade
    save_json(s, int_vars_jsonstr, :int_vars) 
  end
end

function save_object(s::SerializerState, milp::MixedIntegerLinearProgram{<:FieldElem})
  milp_coeffs = milp.polymake_milp.LINEAR_OBJECTIVE
  int_vars = milp.polymake_milp.INTEGER_VARIABLES
  save_data_dict(s) do
    save_object(s, milp.feasible_region, :feasible_region)
    save_object(s, milp.convention, :convention)
    save_object(s, _pmdata_for_oscar(milp_coeffs, coefficient_field(milp)), :milp_coeffs)
    save_object(s, _pmdata_for_oscar(int_vars, QQ), :int_vars)
  end
end

function load_object(s::DeserializerState, tp::TypeParams{<:MixedIntegerLinearProgram, QQField})
  field = Oscar.params(tp)
  fr = load_object(s, TypeParams(Polyhedron, field), :feasible_region)
  conv = load_object(s, String, :convention)
  milp_coeffs = load_node(s, :milp_coeffs) do _
    Polymake.call_function(
      :common,
      :deserialize_json_string,
      JSON.json(load_json(s, Dict{String, Any}))
    )
  end
  int_vars = load_node(s, :int_vars) do _
    Polymake.call_function(
      :common,
      :deserialize_json_string,
      JSON.json(load_json(s, Dict{String, Any}))
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

function load_object(s::DeserializerState, tp::TypeParams{<:MixedIntegerLinearProgram, <:Field})
  field = Oscar.params(tp)
  coeff_type = elem_type(field)
  conv = load_object(s, String, :convention)
  fr = load_object(s, TypeParams(Polyhedron, field), :feasible_region)
  milp_coeffs = load_object(s, TypeParams(Vector{coeff_type}, field), :milp_coeffs)
  int_vars = load_object(s, Vector{Int}, :int_vars)

  all = Polymake._lookup_multi(pm_object(fr), "MILP")
  index = 0
  for i in 1:length(all)
    if _pmdata_for_oscar(all[i].LINEAR_OBJECTIVE, field) == milp_coeffs && all[i].INTEGER_VARIABLES == Set(int_vars)
      index = i
      break
    end
  end
  milp = Polymake._lookup_multi(pm_object(fr), "MILP", index-1)
  return MixedIntegerLinearProgram{coeff_type}(fr, milp, Symbol(conv), field)
end

# use generic serialization for the other types:
@register_serialization_type Cone
@register_serialization_type PolyhedralComplex
@register_serialization_type Polyhedron
@register_serialization_type PolyhedralFan
@register_serialization_type SubdivisionOfPoints
