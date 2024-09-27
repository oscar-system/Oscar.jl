@register_serialization_type MPolyAnyMap uses_params

function save_type_params(s::SerializerState, phi::T) where T <: MPolyAnyMap
  save_data_dict(s) do
    save_object(s, encode_type(T), :name)

    save_data_dict(s, :params) do
      save_typed_object(s, domain(phi), :domain)
      save_typed_object(s, codomain(phi), :codomain)
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:MPolyAnyMap})
  d = load_typed_object(s, :domain)
  c = load_typed_object(s, :codomain)
  return (d, c)
end

function save_object(s::SerializerState{T}, phi::MPolyAnyMap) where {T}
  save_data_dict(s) do
    save_object(s, _images(phi), :images)

    cm = coefficient_map(phi)

    @req _check_whether_this_is_admissible_serialization(typeof(cm), T) "the type of the coefficient map is not supported for long term storage"
    if !isnothing(cm)
      save_object(s, cm, :coeff_map)
    end
  end
end

_check_whether_this_is_admissible_serialization(::Type{CMType}, ::Type{SerializeParam}) where {CMType, SerializeParam <: OscarSerializer} = false
_check_whether_this_is_admissible_serialization(::Type{Nothing}, ::Type{SerializeParam}) where {SerializeParam <: OscarSerializer} = true
# TODO This list can be extended for more types of coefficient maps if needed
_check_whether_this_is_admissible_serialization(::Type{CMType}, ::Type{SerializeParam}) where {CMType, SerializeParam <: IPCSerializer} = true
_check_whether_this_is_admissible_serialization(::Type{Nothing}, ::Type{SerializeParam}) where {SerializeParam <: IPCSerializer} = true

function load_object(s::DeserializerState, ::Type{<:MPolyAnyMap},
    params::Tuple{<:Union{<:MPolyRing, <:MPolyQuoRing}, <:Ring}
  )
  d = params[1]
  c = params[2]
  imgs = load_object(s, Vector, c, :images)

  if haskey(s, :coeff_map)
    return hom(d, c, load_object(s, :coeff_map), img_gens)
  end
  return MPolyAnyMap(d, c, nothing, imgs)
end
