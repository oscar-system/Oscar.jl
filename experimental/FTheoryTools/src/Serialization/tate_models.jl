@register_serialization_type GlobalTateModel uses_params



###########################################################################
# This function saves the types of the data that define a global Tate model
###########################################################################

function save_type_params(s::SerializerState, gtm::GlobalTateModel)
  save_data_dict(s) do
    save_object(s, encode_type(GlobalTateModel), :name)
    base = base_space(gtm)
    ambient = ambient_space(gtm)
    tate_polynomial_ring = parent(gtm.tate_polynomial)
    tate_section_ring = parent(gtm.tate_a1)
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
        save_typed_object(s, ambient, :ambient_space)
      end

      if serialize_with_id(tate_polynomial_ring)
        parent_ref = save_as_ref(s, tate_polynomial_ring)
        save_object(s, parent_ref, :tate_polynomial_ring)
      else
        save_typed_object(s, tate_polynomial_ring, :tate_polynomial_ring)
      end

      if serialize_with_id(tate_section_ring)
        parent_ref = save_as_ref(s, tate_section_ring)
        save_object(s, parent_ref, :tate_section_ring)
      else
        save_typed_object(s, tate_section_ring, :tate_section_ring)
      end
    end
  end
end



###########################################################################
# This function loads the types of the data that define a global Tate model
###########################################################################

function load_type_params(s::DeserializerState, ::Type{<: GlobalTateModel}, dict::Dict)
  return (
    load_typed_object(s, dict[:base_space]),
    load_typed_object(s, dict[:ambient_space]),
    load_typed_object(s, dict[:tate_polynomial_ring]),
    load_typed_object(s, dict[:tate_section_ring])
  )
end



#########################################
# This function saves a global Tate model
#########################################

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



#########################################
# This function loads a global Tate model
#########################################

function load_object(s::DeserializerState, ::Type{<: GlobalTateModel}, dict::Dict,
                     params::Tuple{NormalToricVariety,
                                   NormalToricVariety,
                                   MPolyDecRing,
                                   MPolyDecRing})
  p = load_object(s, MPolyDecRingElem, dict[:tate_polynomial], params[3])
  section_polys = load_object(s, Vector, dict[:section_polys], params[4])

  model = GlobalTateModel(section_polys..., p, params[1], params[2])

  # not 100% sure about this line?
  set_attribute!(model, :base_fully_specified, true)
  return model
end
