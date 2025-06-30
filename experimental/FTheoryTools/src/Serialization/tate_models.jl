@register_serialization_type GlobalTateModel uses_id [
  :Kbar3,
  :_ambient_space_base_divisor_pairs_to_be_considered,
  :_ambient_space_divisor_pairs_to_be_considered,
  :ambient_space_models_of_g4_fluxes,
  :ambient_space_models_of_g4_fluxes_indices,
  :associated_literature_models,
  :basis_of_h22_ambient,
  :basis_of_h22_ambient_indices,
  :basis_of_h22_hypersurface,
  :basis_of_h22_hypersurface_indices,
  :birational_literature_models,
  :chern_classes,
  :classes_of_model_sections,
  :classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes,
  :components_of_dual_graph,
  :components_of_simplified_dual_graph,
  :converter_dict_h22_ambient,
  :converter_dict_h22_hypersurface,
  :degree_of_Kbar_of_tv_restricted_to_ci,
  :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph,
  :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph,
  :dual_graph,
  :estimated_number_of_triangulations,
  :euler_characteristic,
  :exceptional_classes,
  :exceptional_divisor_indices,
  :g4_flux_tuple_list,
  :gauge_algebra,
  :generating_sections,
  :genus_ci,
  :genus_of_components_of_dual_graph,
  :genus_of_components_of_simplified_dual_graph,
  :global_gauge_quotients,
  :h11,
  :h12,
  :h13,
  :h22,
  :index_facet_interior_divisors,
  :inter_dict,
  :intersection_number_among_ci_cj,
  :intersection_number_among_nontrivial_ci_cj,
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
  :poly_index,
  :resolutions,
  :resolution_generating_sections,
  :resolution_zero_sections,
  :s_inter_dict,
  :simplified_dual_graph,
  :torsion_sections,
  :triang_quick,
  :tunable_sections,
  :vertices,
  :weierstrass_model,
  :weighted_resolutions,
  :weighted_resolution_generating_sections,
  :weighted_resolution_zero_sections,
  :zero_section,
  :zero_section_class,
  :zero_section_index,

  # these attributes should be moved into some form of meta data
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



###########################################################################
# This function saves the types of the data that define a global Tate model
###########################################################################

function type_params(m::GlobalTateModel)
  extra_params = [data[2] => type_params(data[1]) for data in 
    [(explicit_model_sections(m), :explicit_model_sections), 
     (defining_classes(m), :defining_classes), 
     (model_section_parametrization(m), :model_section_parametrization)] 
    if !isempty(data[1])]

  TypeParams(
    GlobalTateModel,
    :base_space => base_space(m),
    :ambient_space => ambient_space(m),
    :fiber_ambient_space => fiber_ambient_space(m),
    :tate_polynomial_ring => parent(tate_polynomial(m)),
    extra_params...
  )
end



#########################################
# This function saves a global Tate model
#########################################

function save_object(s::SerializerState, m::GlobalTateModel)
  @req base_space(m) isa NormalToricVariety "Currently, we only serialize global Tate models defined over a toric base space"
  @req ambient_space(m) isa NormalToricVariety "Currently, we only serialize global Tate models defined within a toric ambient space"
  @req m.tate_polynomial !== nothing "Currently, we only serialize global Tate models for which the Tate polynomial (and not only the Tate ideal sheaf) is known"
  save_data_dict(s) do
    for (data, key) in [
        (explicit_model_sections(m), :explicit_model_sections),
        (defining_classes(m), :defining_classes),
        (model_section_parametrization(m), :model_section_parametrization)
        ]
      !isempty(data) && save_object(s, data, key)
    end
    save_object(s, tate_polynomial(m), :tate_polynomial)
  end
end



#########################################
# This function loads a global Tate model
#########################################

function load_object(s::DeserializerState, ::Type{<: GlobalTateModel}, params::Dict)
  base_space = params[:base_space]
  amb_space = params[:ambient_space]
  tp_ring = params[:tate_polynomial_ring]
  pt = load_object(s, MPolyDecRingElem, tp_ring, :tate_polynomial)
  explicit_model_sections = Dict{String, MPolyRingElem}()
  if haskey(s, :explicit_model_sections)
    explicit_model_sections = load_object(s, Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}, params[:explicit_model_sections], :explicit_model_sections)
  end
  model_section_parametrization = Dict{String, MPolyRingElem}()
  if haskey(s, :model_section_parametrization)
    model_section_parametrization = load_object(s, Dict{String, MPolyRingElem}, params[:model_section_parametrization], :model_section_parametrization)
  end
  defining_classes = Dict{String, ToricDivisorClass}()
  if haskey(s, :defining_classes)
    defining_classes = load_object(s, Dict{String, ToricDivisorClass}, params[:defining_classes], :defining_classes)
  end
  model = GlobalTateModel(explicit_model_sections, model_section_parametrization, pt, base_space, amb_space)
  model.defining_classes = defining_classes
  @req cox_ring(ambient_space(model)) == parent(tate_polynomial(model)) "Tate polynomial not in Cox ring of toric ambient space"
  return model
end
