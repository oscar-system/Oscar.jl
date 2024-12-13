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

    # Save resolutions, if they are known.
    if has_resolutions(w)
      res = resolutions(w)
      resolution_loci = [k[1] for k in res]
      exceptional_divisors = [k[2] for k in res]
      attrs_dict[:resolution_loci] = resolution_loci
      attrs_dict[:exceptional_divisors] = exceptional_divisors
    end

    # Have intersection numbers been computed?
    if has_attribute(w, :inter_dict)
      attrs_dict[:inter_dict] = get_attribute(w, :inter_dict)
    end

    # Have special intersections been remembered?
    if has_attribute(w, :s_inter_dict)
      attrs_dict[:s_inter_dict] = get_attribute(w, :s_inter_dict)
    end

    # Do we know which pairs of ambient space divisors intersect non-trivially with the hypersurface?
    if has_attribute(w, :_ambient_space_divisor_pairs_to_be_considered)
      attrs_dict[:_ambient_space_divisor_pairs_to_be_considered] = _ambient_space_divisor_pairs_to_be_considered(w)
    end

    # Do we know which pairs of ambient space base divisors intersect non-trivially with the hypersurface?
    if has_attribute(w, :_ambient_space_base_divisor_pairs_to_be_considered)
      attrs_dict[:_ambient_space_base_divisor_pairs_to_be_considered] = _ambient_space_base_divisor_pairs_to_be_considered(w)
    end

    # Do we know ambient space models for g4-fluxes?
    if has_attribute(w, :ambient_space_models_of_g4_fluxes_indices)
      attrs_dict[:g4_flux_tuple_list] = get_attribute(w, :ambient_space_models_of_g4_fluxes_indices)
    end

    # Do we know the well-quantized G4-fluxes
    if has_attribute(w, :well_quantized_ambient_space_models_of_g4_fluxes)
      res = well_quantized_ambient_space_models_of_g4_fluxes(w, check = false)
      attrs_dict[:well_quantized_integral] = res[1]
      attrs_dict[:well_quantized_rational] = res[2]
    end

    # Do we know the well-quantized and vertical G4-fluxes
    if has_attribute(w, :well_quantized_and_vertical_ambient_space_models_of_g4_fluxes)
      res = well_quantized_and_vertical_ambient_space_models_of_g4_fluxes(w, check = false)
      attrs_dict[:well_quantized_and_vertical_integral] = res[1]
      attrs_dict[:well_quantized_and_vertical_rational] = res[2]
    end

    # Save all of the above data...
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
  @req cox_ring(ambient_space(model)) == parent(weierstrass_polynomial(model)) "Weierstrass polynomial not in Cox ring of toric ambient space"

  # Set "generic" attributes
  attrs_data = haskey(s, :__attrs) ? load_typed_object(s, :__attrs) : Dict{Symbol, Any}()
  for (key, value) in attrs_data
    if (key != :resolution_loci) && (key != :exceptional_divisors)
      set_attribute!(model, Symbol(key), value)
    end
  end

  # Resolution loci known? If so, set them.
  if haskey(attrs_data, :resolution_loci)
    resolution_loci = attrs_data[:resolution_loci]
    exceptional_divisors = attrs_data[:exceptional_divisors]
    @req length(exceptional_divisors) == length(exceptional_divisors) "Inconsistency upon loading resolutions"
    set_attribute!(model, :resolutions, [[resolution_loci[i], exceptional_divisors[i]] for i in 1:length(resolution_loci)])
  end

  # Some intersection numbers known? That is, is it known that the intersection of the toric divisors
  # (i1, i2, i3, i4) with the hypersurface consists of k points (counted with multiplicity)? If so, set it.
  # That is, we have a dictionary which assigns the tuple (i1, i2, i3, i4) - sorted in ascending order - the
  # intersection number k. In general, k could be a rational number! Cf. orbifold singularities.
  if haskey(attrs_data, :inter_dict)
    set_attribute!(model, :inter_dict, attrs_data[:inter_dict])
  end

  # Some special intersection numbers known? If we intersect the toric divisors
  # (i1, i2, i3, i4) with the hypersurface p, then we can set a number of variables to 1
  # by the scaling relations, which simplifies the hypersurface. The remaining variables
  # are subject to a number of remaining SR-ideal generators and some remaining scaling
  # relations. The special intersection algorithm remembers the intersection points of such
  # a locus in a dictionary. The key in this dictionary is the string formed from
  # (simplified_hypersurface, remaining variables, remaining scaling relations, remaining SR-generators)
  # and its value is the number of intersection points. If such information is known,
  # then set it.
  if haskey(attrs_data, :s_inter_dict)
    set_attribute!(model, :s_inter_dict, attrs_data[:s_inter_dict])
  end

  # For the computation off all well-quantized G4-fluxes we compute intersection numbers in the toric
  # ambient space. Specifically, we intersect the ambient space G4-flux candidates (a product of algebraic cycles,
  # each of which is associated to a toric divisor) with another pair of (algebraic cycles associated to) toric
  # divisors and the hypersurface. Of course, only certain pairs of toric divisors restrict non-trivially
  # to the CY-hypersurface in question. We make a selection of divisor pairs that (likely - we only make a naive
  # test regarding a non-trivial restriction) restrict non-trivially to the hypersurface. If the list of those
  # toric divisor pairs is known, set it.
  if haskey(attrs_data, :_ambient_space_divisor_pairs_to_be_considered)
    set_attribute!(model, :_ambient_space_divisor_pairs_to_be_considered, attrs_data[:_ambient_space_divisor_pairs_to_be_considered])
  end
  if haskey(attrs_data, :_ambient_space_base_divisor_pairs_to_be_considered)
    set_attribute!(model, :_ambient_space_base_divisor_pairs_to_be_considered, attrs_data[:_ambient_space_base_divisor_pairs_to_be_considered])
  end

  # Likewise, there are only so many ambient space G4-flux candidates. If those are known, set them.
  # In the serialization, we remember only the integer tuples that tell us which toric divisors to consider,
  # as this saves memory. If known, set this data and recompute the corresponding cohomology classes and
  # set those classes as well.
  if haskey(attrs_data, :g4_flux_tuple_list)
    tuple_list = attrs_data[:g4_flux_tuple_list]
    S = cohomology_ring(amb_space, check = false)
    c_ds = [k.f for k in gens(S)]
    ambient_space_models_of_g4_fluxes = [cohomology_class(amb_space, MPolyQuoRingElem(c_ds[my_tuple[1]]*c_ds[my_tuple[2]], S)) for my_tuple in tuple_list]
    set_attribute!(model, :ambient_space_models_of_g4_fluxes, ambient_space_models_of_g4_fluxes)
    set_attribute!(model, :ambient_space_models_of_g4_fluxes_indices, tuple_list)
  end
    
  # Are the well-quanitzed G4-fluxes known?
  if haskey(attrs_data, :well_quantized_integral) && haskey(attrs_data, :well_quantized_rational)
    quant_tuple = (attrs_data[:well_quantized_integral], attrs_data[:well_quantized_rational])
    set_attribute!(model, :well_quantized_ambient_space_models_of_g4_fluxes, quant_tuple)
  end

  # Are the well-quantized and vertical G4-fluxes known?
  if haskey(attrs_data, :well_quantized_and_vertical_integral) && haskey(attrs_data, :well_quantized_and_vertical_rational)
    quant_tuple = (attrs_data[:well_quantized_and_vertical_integral], attrs_data[:well_quantized_and_vertical_rational])
    set_attribute!(model, :well_quantized_and_vertical_ambient_space_models_of_g4_fluxes, quant_tuple)
  end

  return model
end
