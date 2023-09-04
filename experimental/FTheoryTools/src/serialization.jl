@register_serialization_type GlobalTateModel uses_params

function save_type_params(s::SerializerState, gtm::GlobalTateModel)
  save_data_dict(s) do
    save_object(s, encode_type(GlobalTateModel), :name)
    base = base_space(gtm)
    ambient = ambient_space(gtm)
    global_poly_parent = parent(gtm.tate_polynomial)
    sections_poly_parent = parent(gtm.tate_a1)
    save_data_dict(s, :params) do
      if serialize_with_id(base)
        parent_ref = save_as_ref(s, base)
        save_object(s, parent_ref, :base_space)
      else
        save_typed_object(s, base, :base_space)
      end

      if serialize_with_id(ambient)
        parent_ref = save_as_ref(s, ambient)
        save_object(s, parent_ref, :ambient_space)
      else
        save_typed_object(s, ambient, :base_space)
      end

      if serialize_with_id(global_poly_parent)
        parent_ref = save_as_ref(s, global_poly_parent)
        save_object(s, parent_ref, :global_poly_parent)
      else
        save_typed_object(s, global_poly_parent, :global_poly_parent)
      end

      if serialize_with_id(sections_poly_parent)
        parent_ref = save_as_ref(s, sections_poly_parent)
        save_object(s, parent_ref, :sections_poly_parent)
      else
        save_typed_object(s, sections_poly_parent, :sections_poly_parent)
      end
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<: GlobalTateModel}, dict::Dict)
  return (
    load_typed_object(s, dict[:base_space]),
    load_typed_object(s, dict[:ambient_space]),
    load_typed_object(s, dict[:global_poly_parent]),
    load_typed_object(s, dict[:sections_poly_parent])
  )
end

function save_object(s::SerializerState, gtm::GlobalTateModel)
  save_data_dict(s) do
    save_data_array(s, :section_polys) do
      save_object(s, tate_section_a1(gtm))
      save_object(s, tate_section_a2(gtm))
      save_object(s, tate_section_a3(gtm))
      save_object(s, tate_section_a4(gtm))
      save_object(s, tate_section_a6(gtm))
    end
    save_object(s, tate_polynomial(gtm), :tate_polynomial)
  end
end

function load_object(s::DeserializerState, ::Type{<: GlobalTateModel}, dict::Dict,
                     params::Tuple{NormalToricVariety,
                                   NormalToricVariety,
                                   MPolyDecRing,
                                   MPolyDecRing})
  section_polys = load_object(s, Vector, dict[:section_polys], params[3])
  p = load_object(s, MPolyDecRingElem, dict[:tate_polynomial], params[4])

  model = GlobalTateModel(section_polys..., p, params[1], params[2])

  # not 100% sure about this line?
  set_attribute!(model, :base_fully_specified, true)
  return model
end
