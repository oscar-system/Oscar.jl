@register_serialization_type GlobalTateModel uses_params

###########################################################################
# This function saves the types of the data that define a global Tate model
###########################################################################

function save_type_params(s::SerializerState, gtm::GlobalTateModel)
  save_data_dict(s) do
    save_object(s, encode_type(GlobalTateModel), :name)
    base, ambient, tp_ring = base_space(gtm), ambient_space(gtm), parent(tate_polynomial(gtm))

    save_data_dict(s, :params) do
      for (obj, key) in [(base, :base_space), (ambient, :ambient_space), (tp_ring, :tate_polynomial_ring)]
        if serialize_with_id(obj)
          save_object(s, save_as_ref(s, obj), key)
        else
          save_typed_object(s, obj, key)
        end
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
    load_typed_object(s, :tate_polynomial_ring)
  )
end


#########################################
# This function saves a global Tate model
#########################################

function save_object(s::SerializerState, gtm::GlobalTateModel)
  @req base_space(gtm) isa NormalToricVariety "Currently, we only serialize Tate models defined over a toric base space"
  @req ambient_space(gtm) isa NormalToricVariety "Currently, we only serialize Tate models defined within a toric ambient space"

  save_data_dict(s) do
    for (data, key) in [
        (explicit_model_sections(gtm), :explicit_model_sections),
        (defining_section_parametrization(gtm), :defining_section_parametrization),
        (defining_classes(gtm), :defining_classes)
        ]
      !isempty(data) && save_typed_object(s, data, key)
    end
    save_object(s, tate_polynomial(gtm), :tate_polynomial)
    save_data_array(s, :boolean_data) do
      save_object(s, is_partially_resolved(gtm))
    end
  end
end


#########################################
# This function loads a global Tate model
#########################################

function load_object(s::DeserializerState, ::Type{<: GlobalTateModel}, params::Tuple{NormalToricVariety, NormalToricVariety, MPolyDecRing})
  base_space, ambient_space, tp_ring = params
  pt = load_object(s, MPolyDecRingElem, tp_ring, :tate_polynomial)
  explicit_model_sections = haskey(s, :explicit_model_sections) ? load_typed_object(s, :explicit_model_sections) : Dict{String, MPolyRingElem}()
  defining_section_parametrization = haskey(s, :defining_section_parametrization) ? load_typed_object(s, :defining_section_parametrization) : Dict{String, MPolyRingElem}()
  defining_classes = haskey(s, :defining_classes) ? load_typed_object(s, :defining_classes) : Dict{String, ToricDivisorClass}()
  model = GlobalTateModel(explicit_model_sections, defining_section_parametrization, pt, base_space, ambient_space)
  model.defining_classes = defining_classes
  set_attribute!(model, :partially_resolved, load_object(s, Vector{Bool}, :boolean_data)[1])
  return model
end
