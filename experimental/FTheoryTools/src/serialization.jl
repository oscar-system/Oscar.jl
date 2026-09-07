###########################################################################
# (1) Register the FTheory model types with their fields for serialization
###########################################################################

const COMMON_FIELDS = [
  :kbar3,
  :_ambient_space_base_divisor_pairs_to_be_considered,
  :_ambient_space_divisor_pairs_to_be_considered,
  :ambient_space_models_of_g4_fluxes,
  :ambient_space_models_of_g4_fluxes_indices,
  :associated_literature_models,
  :basis_of_h22_ambient,
  :basis_of_h22_ambient_indices,
  :gens_of_h22_hypersurface,
  :gens_of_h22_hypersurface_indices,
  :birational_literature_models,
  :chern_classes,
  :classes_of_model_sections,
  :classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes,
  :components_of_dual_graph,
  :components_of_simplified_dual_graph,
  :converter_dict_h22_ambient,
  :converter_dict_h22_hypersurface,
  :degrees_of_kbar_restrictions_to_ci_curves,
  :degrees_of_kbar_restrictions_to_components_of_dual_graph,
  :degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph,
  :dual_graph,
  :estimated_number_of_triangulations,
  :euler_characteristic,
  :exceptional_classes,
  :exceptional_divisor_indices,
  :g4_flux_tuple_list,
  :gauge_algebra,
  :generating_sections,
  :genera_of_ci_curves,
  :genera_of_components_of_dual_graph,
  :genera_of_components_of_simplified_dual_graph,
  :global_gauge_group_quotient,
  :global_tate_model,
  :has_quick_triangulation,
  :hodge_h11,
  :hodge_h12,
  :hodge_h13,
  :hodge_h22,
  :indices_of_trivial_ci_curves,
  :inter_dict,
  :topological_intersection_numbers_among_ci_curves,
  :topological_intersection_numbers_among_nontrivial_ci_curves,
  :is_calabi_yau,
  :literature_identifier,
  :matrix_integral_quant_transverse,
  :matrix_rational_quant_transverse,
  :matrix_integral_quant_transverse_nobreak,
  :matrix_rational_quant_transverse_nobreak,
  :max_lattice_pts_in_facet,
  :model_parameters,
  :model_sections,
  :offset_quant_transverse,
  :offset_quant_transverse_nobreak,
  :partially_resolved,
  :polytope_index,
  :resolutions,
  :resolution_generating_sections,
  :resolution_zero_sections,
  :s_inter_dict,
  :simplified_dual_graph,
  :torsion_sections,
  :tunable_sections,
  :vertices,
  :weierstrass_model,
  :weighted_resolutions,
  :weighted_resolution_generating_sections,
  :weighted_resolution_zero_sections,
  :zero_section,
  :zero_section_class,
  :zero_section_index,
  # the following attributes should be moved into a meta data framework eventually
  :arxiv_doi,
  :arxiv_id,
  :arxiv_link,
  :arxiv_model_equation_number,
  :arxiv_model_page,
  :arxiv_model_section,
  :arxiv_version,
  :journal_doi,
  :journal_link,
  :journal_model_equation_number,
  :journal_model_page,
  :journal_model_section,
  :journal_name,
  :journal_pages,
  :journal_report_numbers,
  :journal_volume,
  :journal_year,
  :model_description,
  :model_index,
  :paper_authors,
  :paper_buzzwords,
  :paper_description,
  :paper_title,
]
@register_serialization_type WeierstrassModel uses_id COMMON_FIELDS
@register_serialization_type HypersurfaceModel uses_id COMMON_FIELDS
@register_serialization_type GlobalTateModel uses_id COMMON_FIELDS

const LEGACY_FTHEORY_MODEL_ATTRIBUTE_NAMES = Dict(
  :Kbar3 => :kbar3,
  :degree_of_Kbar_of_tv_restricted_to_ci => :degrees_of_kbar_restrictions_to_ci_curves,
  :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph =>
    :degrees_of_kbar_restrictions_to_components_of_dual_graph,
  :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph =>
    :degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph,
  :genus_ci => :genera_of_ci_curves,
  :genus_of_components_of_dual_graph => :genera_of_components_of_dual_graph,
  :genus_of_components_of_simplified_dual_graph =>
    :genera_of_components_of_simplified_dual_graph,
  :global_gauge_group_quotients => :global_gauge_group_quotient,
  :h11 => :hodge_h11,
  :h12 => :hodge_h12,
  :h13 => :hodge_h13,
  :h22 => :hodge_h22,
  :index_facet_interior_divisors => :indices_of_trivial_ci_curves,
  :intersection_number_among_ci_cj => :topological_intersection_numbers_among_ci_curves,
  :intersection_number_among_nontrivial_ci_cj =>
    :topological_intersection_numbers_among_nontrivial_ci_curves,
  :poly_index => :polytope_index,
  :triang_quick => :has_quick_triangulation,
)

function Serialization.load_attrs(
  s::DeserializerState, m::Union{WeierstrassModel,GlobalTateModel,HypersurfaceModel}
)
  !with_attrs(s) && return nothing

  haskey(s, :attrs) && load_node(s, :attrs) do d
    serialized_attributes = keys(d)
    for old_name in serialized_attributes
      new_name = get(LEGACY_FTHEORY_MODEL_ATTRIBUTE_NAMES, old_name, old_name)
      old_name != new_name && new_name in serialized_attributes && continue
      set_attribute!(m, new_name, Serialization.load_typed_object(s, old_name))
    end
  end
end

###########################################################################
# (2) Type parameters
###########################################################################

function _fmodel_params(m::Union{WeierstrassModel,GlobalTateModel,HypersurfaceModel})
  params = [:base_space => base_space(m),
    :ambient_space => ambient_space(m),
    :fiber_ambient_space => fiber_ambient_space(m),
    :hypersurface_equation_ring => parent(hypersurface_equation(m)),
    if m isa HypersurfaceModel
      (
        :hypersurface_equation_parametrization_ring =>
          parent(hypersurface_equation_parametrization(m))
      )
    else
      nothing
    end,
    if !isempty(explicit_model_sections(m))
      (:explicit_model_sections => type_and_params(explicit_model_sections(m)))
    else
      nothing
    end,
    if !isempty(defining_classes(m))
      (:defining_classes => type_and_params(defining_classes(m)))
    else
      nothing
    end,
    if !isempty(model_section_parametrization(m))
      (:model_section_parametrization => type_and_params(model_section_parametrization(m)))
    else
      nothing
    end]
  return filter(!isnothing, params)
end

type_and_params(m::T) where {T<:Union{WeierstrassModel,GlobalTateModel,HypersurfaceModel}} =
  TypeAndParams(
    T, _fmodel_params(m)...
  )

###########################################################################
# (3) Saving
###########################################################################

function _check_serializable(m::Union{WeierstrassModel,GlobalTateModel,HypersurfaceModel})
  @req base_space(m) isa NormalToricVariety "For serialization, base space must be toric"
  @req ambient_space(m) isa NormalToricVariety "For serialization, ambient space must be toric"
  @req hypersurface_equation(m) !== nothing "For serialization, hypersurface equation must be known (not only the ideal sheaf)"
end

function _save_common_data(
  s::SerializerState, m::Union{WeierstrassModel,GlobalTateModel,HypersurfaceModel}
)
  save_object(s, hypersurface_equation(m), :hypersurface_equation)
  m isa HypersurfaceModel && save_object(
    s, hypersurface_equation_parametrization(m), :hypersurface_equation_parametrization
  )
  for (data, key) in [(explicit_model_sections(m), :explicit_model_sections),
    (defining_classes(m), :defining_classes),
    (model_section_parametrization(m), :model_section_parametrization)]
    !isempty(data) && save_object(s, data, key)
  end
end

function save_object(
  s::SerializerState, m::Union{WeierstrassModel,GlobalTateModel,HypersurfaceModel}
)
  _check_serializable(m)
  save_data_dict(() -> _save_common_data(s, m), s)
end

###########################################################################
# (4) Loading
###########################################################################

function _maybe_load(
  s::DeserializerState, ::Type{T}, key::Symbol, tp::TypeAndParams
) where {S,T<:Dict{String,S}}
  if haskey(s, key)
    val = tp[key]
    if val isa TypeAndParams
      return load_object(s, val, key)
    end
    dict_params = val isa Dict ? val : Dict(:key_params => nothing, :value_params => val)
    return load_object(s, TypeAndParams(T, dict_params), key)
  end
  return T()
end

function _load_common_parts(s::DeserializerState, tp::TypeAndParams)
  # TODO: The following should eventually be handled by a serializatino upgrade script.
  # TODO: Discussed with Antony on slack today (August 22, 2025), and we decided to postpone.
  # TODO: This can/should/may be fixed after the OSCAR 1.5.0 release and the corresponding upgrade script.
  ring_key, poly_key = first((
    (rk, pk) for (rk, pk) in (
      (:hypersurface_equation_ring, :hypersurface_equation),
      (:weierstrass_polynomial_ring, :weierstrass_polynomial),
      (:tate_polynomial_ring, :tate_polynomial),
    ) if haskey(tp, rk)
  ))
  def_poly = load_object(s, TypeAndParams(MPolyDecRingElem, tp[ring_key]), poly_key)
  @req coordinate_ring(tp[:ambient_space]) == parent(def_poly) "Hypersurface equation not in Cox ring of toric ambient space"

  explicit_model_sections = _maybe_load(
    s,
    Dict{String,MPolyDecRingElem{QQFieldElem,QQMPolyRingElem}},
    :explicit_model_sections,
    tp,
  )
  model_section_parametrization = _maybe_load(
    s, Dict{String,MPolyRingElem}, :model_section_parametrization, tp
  )
  defining_classes = _maybe_load(
    s, Dict{String,ToricDivisorClass}, :defining_classes, tp
  )
  return def_poly, explicit_model_sections, model_section_parametrization, defining_classes
end

function load_object(s::DeserializerState, tp::TypeAndParams{<:WeierstrassModel, <:Tuple{Vararg{Pair}}})
  def_poly, explicit_model_sections, model_section_parametrization, defining_classes = _load_common_parts(s, tp)
  model = WeierstrassModel(
    explicit_model_sections,
    model_section_parametrization,
    def_poly,
    tp[:base_space],
    tp[:ambient_space],
  )
  model.defining_classes = defining_classes
  return model
end

function load_object(s::DeserializerState, tp::TypeAndParams{<:GlobalTateModel, <:Tuple{Vararg{Pair}}})
  def_poly, explicit_model_sections, model_section_parametrization, defining_classes = _load_common_parts(s, tp)
  model = GlobalTateModel(
    explicit_model_sections,
    model_section_parametrization,
    def_poly,
    tp[:base_space],
    tp[:ambient_space],
  )
  model.defining_classes = defining_classes
  return model
end

function load_object(s::DeserializerState, tp::TypeAndParams{<:HypersurfaceModel, <:Tuple{Vararg{Pair}}})
  def_poly, explicit_model_sections, model_section_parametrization, defining_classes = _load_common_parts(s, tp)
  defining_equation_parametrization = load_object(
    s,
    TypeAndParams(MPolyRingElem, tp[:hypersurface_equation_parametrization_ring]),
    :hypersurface_equation_parametrization,
  )
  model = HypersurfaceModel(
    explicit_model_sections,
    defining_equation_parametrization,
    def_poly,
    tp[:base_space],
    tp[:ambient_space],
    tp[:fiber_ambient_space],
  )
  model.model_section_parametrization = model_section_parametrization
  model.defining_classes = defining_classes
  return model
end
