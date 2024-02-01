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
    explicit_model_section_ring = parent(gtm.explicit_model_sections["a1"])
    parametrizing_sections = collect(keys(gtm.defining_section_parametrization))
    if length(parametrizing_sections) > 0
      defining_section_parametrization_ring = parent(gtm.defining_section_parametrization[parametrizing_sections[1]])
    else
      defining_section_parametrization_ring = explicit_model_section_ring
    end

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

      if serialize_with_id(explicit_model_section_ring)
        parent_ref = save_as_ref(s, explicit_model_section_ring)
        save_object(s, parent_ref, :explicit_model_section_ring)
      else
        save_typed_object(s, explicit_model_section_ring, :explicit_model_section_ring)
      end

      if serialize_with_id(defining_section_parametrization_ring)
        parent_ref = save_as_ref(s, defining_section_parametrization_ring)
        save_object(s, parent_ref, :defining_section_parametrization_ring)
      else
        save_typed_object(s, defining_section_parametrization_ring, :defining_section_parametrization_ring)
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
    load_typed_object(s, :explicit_model_section_ring),
    load_typed_object(s, :defining_section_parametrization_ring)
  )
end

#########################################
# This function saves a global Tate model
#########################################

function save_object(s::SerializerState, gtm::GlobalTateModel)
  # Currently, only serialize Tate models with toric defining data
  @req base_space(gtm) isa NormalToricVariety "Currently, we only serialize Tate models defined over a toric base space"
  @req ambient_space(gtm) isa NormalToricVariety "Currently, we only serialize Tate models defined within a toric ambient space"

  # Save information
  save_data_dict(s) do

    # Save keys of explicit_model_sections
    save_data_array(s, :explicit_model_section_keys) do
      for (key, value) in explicit_model_sections(gtm)
        save_object(s, key)
      end
    end

    # Save values of explicit_model_sections
    save_data_array(s, :explicit_model_section_values) do
      for (key, value) in explicit_model_sections(gtm)
        save_object(s, value)
      end
    end

    # Save keys of defining_section_parametrization
    save_data_array(s, :defining_section_parametrization_keys) do
      for (key, value) in defining_section_parametrization(gtm)
        save_object(s, key)
      end
    end

    # Save keys of defining_section_parametrization
    save_data_array(s, :defining_section_parametrization_values) do
      for (key, value) in defining_section_parametrization(gtm)
        save_object(s, value)
      end
    end

    # Tate polynomial
    save_object(s, tate_polynomial(gtm), :tate_polynomial)

    # Boolean values, that are always known for Tate models
    save_data_array(s, :boolean_data) do
      save_object(s, is_partially_resolved(gtm))
    end
  end
end

#########################################
# This function loads a global Tate model
#########################################

function load_object(s::DeserializerState, ::Type{<: GlobalTateModel}, params::Tuple{NormalToricVariety, NormalToricVariety, MPolyDecRing, MPolyDecRing, <:MPolyRing})

  # Extract base and ambient space
  base_space = params[1]
  ambient_space = params[2]

  # Extract Tate polynomial
  pt = load_object(s, MPolyDecRingElem, params[3], :tate_polynomial)

  # Extract explicit_model_sections
  values = load_object(s, Vector, params[4], :explicit_model_section_values)
  keys = load_object(s, Vector, String, :explicit_model_section_keys)
  explicit_model_sections = Dict{String, MPolyRingElem}()
  for i in 1:length(keys)
    explicit_model_sections[keys[i]] = values[i]
  end

  # Extract defining_section_parametrization
  values = load_object(s, Vector, params[5], :defining_section_parametrization_values)
  keys = load_object(s, Vector, String, :defining_section_parametrization_keys)
  defining_section_parametrization = Dict{String, MPolyRingElem}()
  for i in 1:length(keys)
    defining_section_parametrization[keys[i]] = values[i]
  end

  # Construct the model
  model = GlobalTateModel(explicit_model_sections, defining_section_parametrization, pt, base_space, ambient_space)

  # Set boolean attributes
  bools = load_object(s, Vector, Bool, :boolean_data)
  set_attribute!(model, :partially_resolved, bools[1])

  # Return the loaded model
  return model
end
