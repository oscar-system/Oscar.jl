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
    attrs_dict = Dict{Symbol, Any}()
    for (key, value) in w.__attrs
      if value isa String || value isa Vector{String} || value isa Bool
        attrs_dict[key] = value
      end
    end
    if has_resolutions(w)
      res = resolutions(w)
      resolution_loci = [k[1] for k in res]
      exceptional_divisors = [k[2] for k in res]
      attrs_dict[:resolution_loci] = resolution_loci
      attrs_dict[:exceptional_divisors] = exceptional_divisors
    end
    !isempty(attrs_dict) && save_typed_object(s, attrs_dict, :__attrs)
  end
end


#########################################
# This function loads a Weierstrass model
#########################################

function load_object(s::DeserializerState, ::Type{<: WeierstrassModel}, params::Tuple{NormalToricVariety, NormalToricVariety, MPolyDecRing})
  base_space, amb_space, wp_ring = params
  pw = load_object(s, MPolyDecRingElem, wp_ring, :weierstrass_polynomial)
  explicit_model_sections = haskey(s, :explicit_model_sections) ? load_typed_object(s, :explicit_model_sections) : Dict{String, MPolyRingElem}()
  defining_section_parametrization = haskey(s, :defining_section_parametrization) ? load_typed_object(s, :defining_section_parametrization) : Dict{String, MPolyRingElem}()
  defining_classes = haskey(s, :defining_classes) ? load_typed_object(s, :defining_classes) : Dict{String, ToricDivisorClass}()
  model = WeierstrassModel(explicit_model_sections, defining_section_parametrization, pw, base_space, amb_space)
  model.defining_classes = defining_classes
  attrs_data = haskey(s, :__attrs) ? load_typed_object(s, :__attrs) : Dict{Symbol, Any}()
  for (key, value) in attrs_data
    if (key != :resolution_loci) && (key != :exceptional_divisors)
      set_attribute!(model, Symbol(key), value)
    end
  end
  if haskey(attrs_data, :resolution_loci)
    resolution_loci = attrs_data[:resolution_loci]
    exceptional_divisors = attrs_data[:exceptional_divisors]
    @req length(exceptional_divisors) == length(exceptional_divisors) "Inconsistency upon loading resolutions"
    set_attribute!(model, :resolutions, [[resolution_loci[i], exceptional_divisors[i]] for i in 1:length(resolution_loci)])
  end
  @req cox_ring(ambient_space(model)) == parent(weierstrass_polynomial(model)) "Weierstrass polynomial not in Cox ring of toric ambient space"
  return model
end
