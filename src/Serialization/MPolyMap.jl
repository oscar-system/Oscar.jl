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

function save_object(s::SerializerState, phi::MPolyAnyMap)
  save_data_dict(s) do
    save_object(s, _images(phi), :images)

    cm = coefficient_map(phi)
    if !isnothing(cm)
      save_object(s, cm, :coeff_map)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:MPolyAnyMap},
                     params::Tuple{MPolyRing, MPolyRing})
  d = params[1]
  c = params[2]
  imgs = load_object(s, Vector, c, :images)

  if haskey(s, :coeff_map)
    throw("MPolyAnyMap with coefficient map serialization unimplemented")
  end
  return MPolyAnyMap(d, c, nothing, imgs)
end
