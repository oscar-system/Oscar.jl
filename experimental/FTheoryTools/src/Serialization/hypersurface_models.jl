@register_serialization_type HypersurfaceModel uses_params



###########################################################################
# This function saves the types of the data that define a hypersurface model
###########################################################################

function save_type_params(s::SerializerState, h::HypersurfaceModel)
  save_data_dict(s) do
    save_object(s, encode_type(HypersurfaceModel), :name)
    base = base_space(h)
    ambient = ambient_space(h)
    fiber_amb_space = fiber_ambient_space(h)
    hypersurface_equation_ring = parent(hypersurface_equation(h))
    hypersurface_equation_parametrization_ring = parent(h.hypersurface_equation_parametrization)
    explicit_model_section_ring = parent(hypersurface_equation(h))
    explicit_model_section_keys = collect(keys(h.explicit_model_sections))
    if length(explicit_model_section_keys) > 0
      explicit_model_section_ring = parent(explicit_model_sections(h)[explicit_model_section_keys[1]])
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

      if serialize_with_id(fiber_amb_space)
        parent_ref = save_as_ref(s, fiber_amb_space)
        save_object(s, parent_ref, :fiber_ambient_space)
      else
        save_typed_object(s, fiber_amb_space, :fiber_ambient_space)
      end

      if serialize_with_id(hypersurface_equation_ring)
        parent_ref = save_as_ref(s, hypersurface_equation_ring)
        save_object(s, parent_ref, :hypersurface_equation_ring)
      else
        save_typed_object(s, hypersurface_equation_ring, :hypersurface_equation_ring)
      end

      if serialize_with_id(hypersurface_equation_parametrization_ring)
        parent_ref = save_as_ref(s, hypersurface_equation_parametrization_ring)
        save_object(s, parent_ref, :hypersurface_equation_parametrization_ring)
      else
        save_typed_object(s, hypersurface_equation_parametrization_ring, :hypersurface_equation_parametrization_ring)
      end

      if serialize_with_id(explicit_model_section_ring)
        parent_ref = save_as_ref(s, explicit_model_section_ring)
        save_object(s, parent_ref, :explicit_model_section_ring)
      else
        save_typed_object(s, explicit_model_section_ring, :explicit_model_section_ring)
      end

    end
  end
end

############################################################################
# This function loads the types of the data that define a hypersurface model
############################################################################

function load_type_params(s::DeserializerState, ::Type{<: HypersurfaceModel})
  return (
    load_typed_object(s, :base_space),
    load_typed_object(s, :ambient_space),
    load_typed_object(s, :fiber_ambient_space),
    load_typed_object(s, :hypersurface_equation_ring),
    load_typed_object(s, :hypersurface_equation_parametrization_ring),
    load_typed_object(s, :explicit_model_section_ring)
  )
end

##########################################
# This function saves a hypersurface model
##########################################

function save_object(s::SerializerState, h::HypersurfaceModel)
  # Currently, only serialize hypersurface models with toric defining data
  @req base_space(h) isa NormalToricVariety "Currently, we only serialize hypersurface models defined over a toric base space"
  @req ambient_space(h) isa NormalToricVariety "Currently, we only serialize hypersurface models defined within a toric ambient space"

  # Save information
  save_data_dict(s) do

    # hypersurace equation and parametrization
    save_object(s, hypersurface_equation(h), :hypersurface_equation)
    save_object(s, hypersurface_equation_parametrization(h), :hypersurface_equation_parametrization)

    # Save keys of explicit_model_sections
    save_data_array(s, :explicit_model_section_keys) do
      for (key, value) in explicit_model_sections(h)
        save_object(s, key)
      end
    end

    # Save values of explicit_model_sections
    save_data_array(s, :explicit_model_section_values) do
      for (key, value) in explicit_model_sections(h)
        save_object(s, value)
      end
    end

    # Boolean values, that are always known for Tate models
    save_data_array(s, :boolean_data) do
      save_object(s, is_partially_resolved(h))
    end
  end
end

##########################################
# This function loads a hypersurface model
##########################################

function load_object(s::DeserializerState, ::Type{<: HypersurfaceModel}, params::Tuple{NormalToricVariety, NormalToricVariety, NormalToricVariety, <:MPolyRing, <:MPolyRing, <:MPolyRing})
  # load basic data
  base_space = params[1]
  ambient_space = params[2]
  fiber_ambient_space = params[3]
  defining_equation = load_object(s, MPolyRingElem, params[4], :hypersurface_equation)
  defining_equation_parametrization = load_object(s, MPolyRingElem, params[5], :hypersurface_equation_parametrization)

  # Extract explicit_model_sections
  values = load_object(s, Vector, params[6], :explicit_model_section_values)
  keys = load_object(s, Vector, String, :explicit_model_section_keys)
  explicit_model_sections = Dict{String, MPolyRingElem}()
  for i in 1:length(keys)
    explicit_model_sections[keys[i]] = values[i]
  end

  # create and return model
  model = HypersurfaceModel(explicit_model_sections, defining_equation_parametrization, defining_equation, base_space, ambient_space, fiber_ambient_space)
  bools = load_object(s, Vector, Bool, :boolean_data)
  set_attribute!(model, :partially_resolved, bools[1])
  return model
end
