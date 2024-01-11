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

function load_type_params(s::DeserializerState, ::Type{<: GlobalTateModel})
  return (
    load_typed_object(s, :base_space),
    load_typed_object(s, :ambient_space),
    load_typed_object(s, :tate_polynomial_ring),
    load_typed_object(s, :tate_section_ring)
  )
end

#########################################
# This function saves a global Tate model
#########################################

function save_object(s::SerializerState, gtm::GlobalTateModel)
  # Currently, only serialize Tate models with toric defining data
  @req typeof(base_space(gtm)) == NormalToricVariety "Currently, we only serialize Tate models defined over a toric base space"
  @req typeof(ambient_space(gtm)) == NormalToricVariety "Currently, we only serialize Tate models defined within a toric ambient space"

  # Save information
  save_data_dict(s) do
    # Tate sections
    save_data_array(s, :section_polys) do
      save_object(s, tate_section_a1(gtm))
      save_object(s, tate_section_a2(gtm))
      save_object(s, tate_section_a3(gtm))
      save_object(s, tate_section_a4(gtm))
      save_object(s, tate_section_a6(gtm))
    end

    # Tate polynomial
    save_object(s, tate_polynomial(gtm), :tate_polynomial)

    # A couple of boolean values, that are always known for Tate models
    save_data_array(s, :boolean_data) do
      save_object(s, is_partially_resolved(gtm))
      save_object(s, base_fully_specified(gtm))
    end
  end
end

#########################################
# This function loads a global Tate model
#########################################

function load_object(s::DeserializerState, ::Type{<: GlobalTateModel}, params::Tuple{NormalToricVariety, NormalToricVariety, MPolyDecRing, MPolyDecRing})
  # Extract base and ambient space
  base_space = params[1]
  ambient_space = params[2]
  pt = load_object(s, MPolyDecRingElem, params[3], :tate_polynomial)
  tate_sections = load_object(s, Vector, params[4], :section_polys)
  model = GlobalTateModel(tate_sections..., pt, base_space, ambient_space)

  # Set boolean attributes
  bools = load_object(s, Vector, Bool, :boolean_data)
  set_attribute!(model, :partially_resolved, bools[1])
  set_attribute!(model, :base_fully_specified, bools[2])
  return model
end
