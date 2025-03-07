@register_serialization_type WeierstrassModel uses_id [
  # Intersections,
  :inter_dict,
  # special intersections,
  :s_inter_dict,
  # ambient space divisors intersecting non-trivially with the hypersurface
  :_ambient_space_divisor_pairs_to_be_considered,
  # ambient space base divisors intersecting non-trivially with the hypersurface,
  :_ambient_space_base_divisor_pairs_to_be_considered,
  # ambient space models for g4-fluxes
  :ambient_space_models_of_g4_fluxes_indices
]

###########################################################################
# This function saves the types of the data that define a Weierstrass model
###########################################################################

type_params(w::WeierstrassModel) = TypeParams(
  WeierstrassModel,
  :base_space => base_space(w),
  :ambient_space => ambient_space(w),
  :wp_ring => parent(weierstrass_polynomial(w))
)

#########################################
# This function saves a Weierstrass model
#########################################

function save_object(s::SerializerState, w::WeierstrassModel)
  @req base_space(w) isa NormalToricVariety "We only serialize Weierstrass models defined over a toric base space"
  @req ambient_space(w) isa NormalToricVariety "We only serialize Weierstrass models defined within a toric ambient space"

  save_data_dict(s) do
    for (data, key) in [
        (explicit_model_sections(w), :explicit_model_sections),
        (model_section_parametrization(w), :model_section_parametrization),
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

    # Do we know the well-quantized G4-fluxes
    if has_attribute(w, :well_quantized_ambient_space_models_of_g4_fluxes)
      res = well_quantized_ambient_space_models_of_g4_fluxes(w, check = false)
      attrs_dict[:well_quantized_integral] = matrix_integral(res)
      attrs_dict[:well_quantized_rational] = matrix_rational(res)
      if has_attribute(res, :d3_tadpole_constraint)
        attrs_dict[:well_quantized_d3_tadpole] = d3_tadpole_constraint(res)
      end
    end

    # Do we know the well-quantized and vertical G4-fluxes
    if has_attribute(w, :well_quantized_and_vertical_ambient_space_models_of_g4_fluxes)
      res = well_quantized_and_vertical_ambient_space_models_of_g4_fluxes(w, check = false)
      attrs_dict[:well_quantized_and_vertical_integral] = matrix_integral(res)
      attrs_dict[:well_quantized_and_vertical_rational] = matrix_rational(res)
      if has_attribute(res, :d3_tadpole_constraint)
        attrs_dict[:well_quantized_and_vertical_d3_tadpole] = d3_tadpole_constraint(res)
      end
    end

    # Do we know the well-quantized, vertical G4-fluxes that do not break the non-Abelian gauge group?
    if has_attribute(w, :well_quantized_and_vertical_and_no_non_abelian_gauge_group_breaking_ambient_space_models_of_g4_fluxes)
      res = well_quantized_and_vertical_and_no_non_abelian_gauge_group_breaking_ambient_space_models_of_g4_fluxes(w, check = false)
      attrs_dict[:well_quantized_vertical_no_break_integral] = matrix_integral(res)
      attrs_dict[:well_quantized_vertical_no_break_rational] = matrix_rational(res)
      if has_attribute(res, :d3_tadpole_constraint)
        attrs_dict[:well_quantized_vertical_and_no_breaking_d3_tadpole] = d3_tadpole_constraint(res)
      end
    end

    # Save all of the above data...
    !isempty(attrs_dict) && save_typed_object(s, attrs_dict, :__attrs)
  end
end


#########################################
# This function loads a Weierstrass model
#########################################

function load_object(s::DeserializerState, ::Type{<: WeierstrassModel}, params::Dict)
  base_space = params[:base_space]
  amb_space = params[:ambient_space]
  wp_ring = params[:wp_ring]
  pw = load_object(s, MPolyDecRingElem, wp_ring, :weierstrass_polynomial)
  explicit_model_sections = haskey(s, :explicit_model_sections) ? load_typed_object(s, :explicit_model_sections) : Dict{String, MPolyRingElem}()
  model_section_parametrization = haskey(s, :model_section_parametrization) ? load_typed_object(s, :model_section_parametrization) : Dict{String, MPolyRingElem}()
  defining_classes = haskey(s, :defining_classes) ? load_typed_object(s, :defining_classes) : Dict{String, ToricDivisorClass}()
  model = WeierstrassModel(explicit_model_sections, model_section_parametrization, pw, base_space, amb_space)
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
    # We want this inter_dict to be of type Dict{NTuple{4, Int64}, ZZRingElem}().
    # Sadly, serializing and loading turns NTuple{4, Int64} into Tuple.
    # So we need to massage this... Not at all good, as it doubles memory usage!
    original_dict = attrs_data[:inter_dict]
    new_dict = Dict{NTuple{4, Int64}, ZZRingElem}()
    for (key, value) in original_dict
      new_key = NTuple{4, Int64}(key)
      new_dict[new_key] = value
    end
    set_attribute!(model, :inter_dict, new_dict)
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
    fgs = family_of_g4_fluxes(model, quant_tuple[1], quant_tuple[2])
    set_attribute!(fgs, :is_well_quantized, true)
    set_attribute!(fgs, :passes_transversality_checks, false)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
    if haskey(attrs_data, :well_quantized_d3_tadpole)
      set_attribute!(fgs, :d3_tadpole_constraint, attrs_dict[:well_quantized_d3_tadpole])
    end
    set_attribute!(model, :well_quantized_ambient_space_models_of_g4_fluxes, fgs)
  end

  # Are the well-quantized and vertical G4-fluxes known?
  if haskey(attrs_data, :well_quantized_and_vertical_integral) && haskey(attrs_data, :well_quantized_and_vertical_rational)
    quant_tuple = (attrs_data[:well_quantized_and_vertical_integral], attrs_data[:well_quantized_and_vertical_rational])
    fgs = family_of_g4_fluxes(model, quant_tuple[1], quant_tuple[2])
    set_attribute!(fgs, :is_well_quantized, true)
    set_attribute!(fgs, :passes_transversality_checks, true)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, true)
    if haskey(attrs_data, :well_quantized_and_vertical_d3_tadpole)
      set_attribute!(fgs, :d3_tadpole_constraint, attrs_dict[:well_quantized_and_vertical_d3_tadpole])
    end
    set_attribute!(model, :well_quantized_and_transversal_ambient_space_models_of_g4_fluxes, fgs)
  end

  # Are the well-quantized, vertical G4-fluxes which do not break a non-abelian gauge group known?
  if haskey(attrs_data, :well_quantized_vertical_no_break_integral) && haskey(attrs_data, :well_quantized_vertical_no_break_rational)
    quant_tuple = (attrs_data[:well_quantized_vertical_no_break_integral], attrs_data[:well_quantized_vertical_no_break_rational])
    fgs = family_of_g4_fluxes(model, quant_tuple[1], quant_tuple[2])
    set_attribute!(fgs, :is_well_quantized, true)
    set_attribute!(fgs, :passes_transversality_checks, true)
    set_attribute!(fgs, :breaks_non_abelian_gauge_group, false)
    if haskey(attrs_data, :well_quantized_vertical_and_no_breaking_d3_tadpole)
      set_attribute!(fgs, :d3_tadpole_constraint, attrs_dict[:well_quantized_vertical_and_no_breaking_d3_tadpole])
    end
    set_attribute!(model, :well_quantized_and_transversal_and_no_non_abelian_gauge_group_breaking_ambient_space_models_of_g4_fluxes, fgs)
  end

  return model
end
