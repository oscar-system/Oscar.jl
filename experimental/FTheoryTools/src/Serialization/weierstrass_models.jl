@register_serialization_type WeierstrassModel uses_params

###########################################################################
# This function saves the types of the data that define a Weierstrass model
###########################################################################

function save_type_params(s::SerializerState, w::WeierstrassModel)
  save_data_dict(s) do
    save_object(s, encode_type(WeierstrassModel), :name)
    base, ambient, wp_ring = base_space(w), ambient_space(w), parent(weierstrass_polynomial(w))

    save_data_dict(s, :params) do
      for (obj, key) in [(base, :base_space), (ambient, :ambient_space), (wp_ring, :weierstrass_polynomial_ring)]
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
# This function loads the types of the data that define a Weierstrass model
###########################################################################

function load_type_params(s::DeserializerState, ::Type{<: WeierstrassModel})
  return (
    load_typed_object(s, :base_space),
    load_typed_object(s, :ambient_space),
    load_typed_object(s, :weierstrass_polynomial_ring)
  )
end


#########################################
# This function saves a Weierstrass model
#########################################

function save_object(s::SerializerState, w::WeierstrassModel)
  @req base_space(w) isa NormalToricVariety "We only serialize Weierstrass models defined over a toric base space"
  @req ambient_space(w) isa NormalToricVariety "We only serialize Weierstrass models defined within a toric ambient space"

  save_data_dict(s) do
    for (data, key) in [
        (explicit_model_sections(w), :explicit_model_sections),
        (defining_section_parametrization(w), :defining_section_parametrization),
        (defining_classes(w), :defining_classes)
        ]
      !isempty(data) && save_typed_object(s, data, key)
    end
    save_object(s, weierstrass_polynomial(w), :weierstrass_polynomial)
    save_data_array(s, :boolean_data) do
      save_object(s, is_partially_resolved(w))
    end
  end
end


#########################################
# This function loads a Weierstrass model
#########################################

function load_object(s::DeserializerState, ::Type{<: WeierstrassModel}, params::Tuple{NormalToricVariety, NormalToricVariety, MPolyDecRing})
  base_space, ambient_space, wp_ring = params
  pw = load_object(s, MPolyDecRingElem, wp_ring, :weierstrass_polynomial)
  explicit_model_sections = haskey(s, :explicit_model_sections) ? load_typed_object(s, :explicit_model_sections) : Dict{String, MPolyRingElem}()
  defining_section_parametrization = haskey(s, :defining_section_parametrization) ? load_typed_object(s, :defining_section_parametrization) : Dict{String, MPolyRingElem}()
  defining_classes = haskey(s, :defining_classes) ? load_typed_object(s, :defining_classes) : Dict{String, ToricDivisorClass}()
  model = WeierstrassModel(explicit_model_sections, defining_section_parametrization, pw, base_space, ambient_space)
  model.defining_classes = defining_classes
  set_attribute!(model, :partially_resolved, load_object(s, Vector{Bool}, :boolean_data)[1])
  return model
end
