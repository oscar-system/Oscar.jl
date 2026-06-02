using Oscar: _images, coefficient_map

@register_serialization_type MPolyAnyMap

function type_params(phi::MPolyAnyMap{S, T, U, V}) where {S, T, U, V}
  return TypeParams(
    MPolyAnyMap,
    :domain => domain(phi),
    :codomain => codomain(phi)
  )
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

function load_object(s::DeserializerState,
                     tp::TypeParams{<:MPolyAnyMap, <:Tuple{Vararg{Pair}}})
  d = tp[:domain]
  c = tp[:codomain]
  T = elem_type(c)
  imgs = load_object(s, TypeParams(Vector{T}, c), :images)

  if haskey(s, :coeff_map)
    coeff_map = load_object(s, TypeParams(MPolyAnyMap, :domain => base_ring(d), :codomain => base_ring(c)), :coeff_map)
    return MPolyAnyMap(d, c, coeff_map, imgs)
  end
  return MPolyAnyMap(d, c, nothing, imgs)
end
