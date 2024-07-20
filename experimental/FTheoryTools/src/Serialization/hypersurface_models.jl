@register_serialization_type HypersurfaceModel uses_params

############################################################################
# This function saves the types of the data that define a hypersurface model
############################################################################

function save_type_params(s::SerializerState, h::HypersurfaceModel)
  save_data_dict(s) do
    save_object(s, encode_type(HypersurfaceModel), :name)
    save_data_dict(s, :params) do
      save_typed_object(s, base_space(h), :base_space)
      save_typed_object(s, ambient_space(h), :ambient_space)
      save_typed_object(s, fiber_ambient_space(h), :fiber_ambient_space)
      save_typed_object(s, parent(hypersurface_equation(h)), :hypersurface_equation_ring)
      save_typed_object(s, parent(hypersurface_equation_parametrization(h)), :hypersurface_equation_parametrization_ring)
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
    load_typed_object(s, :hypersurface_equation_parametrization_ring)
  )
end


##########################################
# This function saves a hypersurface model
##########################################

function save_object(s::SerializerState, h::HypersurfaceModel)
  @req base_space(h) isa NormalToricVariety "Currently, we only serialize hypersurface models defined over a toric base space"
  @req ambient_space(h) isa NormalToricVariety "Currently, we only serialize hypersurface models defined within a toric ambient space"

  save_data_dict(s) do
    for (data, key) in [
        (explicit_model_sections(h), :explicit_model_sections),
        (defining_classes(h), :defining_classes)
        ]
      !isempty(data) && save_typed_object(s, data, key)
    end
    save_object(s, hypersurface_equation(h), :hypersurface_equation)
    save_object(s, hypersurface_equation_parametrization(h), :hypersurface_equation_parametrization)
    save_data_array(s, :boolean_data) do
      save_object(s, is_partially_resolved(h))
    end
  end
end


##########################################
# This function loads a hypersurface model
##########################################

function load_object(s::DeserializerState, ::Type{<: HypersurfaceModel}, params::Tuple{NormalToricVariety, NormalToricVariety, NormalToricVariety, <:MPolyRing, <:MPolyRing})
  base_space, ambient_space, fiber_ambient_space, R1, R2 = params
  defining_equation = load_object(s, MPolyRingElem, R1, :hypersurface_equation)
  defining_equation_parametrization = load_object(s, MPolyRingElem, R2, :hypersurface_equation_parametrization)
  explicit_model_sections = haskey(s, :explicit_model_sections) ? load_typed_object(s, :explicit_model_sections) : Dict{String, MPolyRingElem}()
  defining_classes = haskey(s, :defining_classes) ? load_typed_object(s, :defining_classes) : Dict{String, ToricDivisorClass}()
  model = HypersurfaceModel(explicit_model_sections, defining_equation_parametrization, defining_equation, base_space, ambient_space, fiber_ambient_space)
  model.defining_classes = defining_classes
  set_attribute!(model, :partially_resolved, load_object(s, Vector{Bool}, :boolean_data)[1])
  return model
end
