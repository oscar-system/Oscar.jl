@register_serialization_type QSMModel

function save_object(s::SerializerState, qsm::QSMModel)
  save_data_dict(s) do
    save_object(s, qsm.vertices, :vertices)
    save_object(s, qsm.poly_index, :poly_inx)

    save_object(s, qsm.triang_quick, :triang_quick)
    save_object(s, qsm.max_lattice_pts_in_facet, :max_lattice_pts_in_facet)
    save_object(s, qsm.estimated_number_of_triangulations, :estimated_number_of_triangulations)
    save_typed_object(s, qsm.hs_model, :hs_model)

    save_object(s, qsm.Kbar3, :Kbar3)
    save_object(s, qsm.h11, :h11)
    save_object(s, qsm.h12, :h12)
    save_object(s, qsm.h13, :h13)
    save_object(s, qsm.h22, :h22)

    save_typed_object(s, qsm.genus_ci, :genus_ci)
    save_typed_object(s, qsm.degree_of_Kbar_of_tv_restricted_to_ci, :degree_of_Kbar_of_tv_restricted_to_ci)
    save_object(s, qsm.intersection_number_among_ci_cj, :intersection_number_among_ci_cj)
    save_object(s, qsm.index_facet_interior_divisors, :index_facet_interior_divisors)
    save_object(s, qsm.intersection_number_among_nontrivial_ci_cj, :intersection_number_among_nontrivial_ci_cj)

    save_object(s, qsm.dual_graph, :dual_graph)
    save_object(s, qsm.components_of_dual_graph, :components_of_dual_graph)
    save_object(s, qsm.degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph, :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph)
    save_object(s, qsm.genus_of_components_of_dual_graph, :genus_of_components_of_dual_graph)

    save_object(s, qsm.simplified_dual_graph, :simplified_dual_graph)
    save_object(s, qsm.components_of_simplified_dual_graph, :components_of_simplified_dual_graph)
    save_object(s, qsm.degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph, :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph)
    save_object(s, qsm.genus_of_components_of_simplified_dual_graph, :genus_of_components_of_simplified_dual_graph)
  end
end

function load_object(s::DeserializerState, ::Type{QSMModel})
  vertices = load_object(s, Vector{Vector{QQFieldElem}}, :vertices)
  poly_index = load_object(s, Int, :poly_inx)
  triang_quick = load_object(s, Bool, :triang_quick)
  max_lattice_pts_in_facet = load_object(s, Int, :max_lattice_pts_in_facet)
  estimated_number_of_triangulations = load_object(s, Int, :estimated_number_of_triangulations)
  hs_model = load_typed_object(s, :hs_model)

  Kbar3 = load_object(s, Int, :Kbar3)
  h11 = load_object(s, Int, :h11)
  h12 = load_object(s, Int, :h12)
  h13 = load_object(s, Int, :h13)
  h22 = load_object(s, Int, :h22)

  genus_ci = load_typed_object(s, :genus_ci)
  degree_of_Kbar_of_tv_restricted_to_ci = load_typed_object(s, :degree_of_Kbar_of_tv_restricted_to_ci)
  intersection_number_among_ci_cj = load_object(s, Matrix, Int, :intersection_number_among_ci_cj)
  index_facet_interior_divisors = load_object(s, Vector{Int}, :index_facet_interior_divisors)
  intersection_number_among_nontrivial_ci_cj = load_object(s, Matrix, Int, :intersection_number_among_nontrivial_ci_cj)

  dual_graph = load_object(s, Graph{Undirected}, :dual_graph)
  components_of_dual_graph = load_object(s, Vector{String}, :components_of_dual_graph)
  degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph = load_object(
    s, Dict{String, Int64}, :degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph)
  genus_of_components_of_dual_graph = load_object(s, Dict{String, Int64}, :genus_of_components_of_dual_graph)

  simplified_dual_graph = load_object(s, Graph{Undirected}, :simplified_dual_graph)
  components_of_simplified_dual_graph = load_object(s, Vector{String}, :components_of_simplified_dual_graph)

  degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph = load_object(
    s, Dict{String, Int},
    :degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph)

  genus_of_components_of_simplified_dual_graph = load_object(
    s, Dict{String, Int},
    :genus_of_components_of_simplified_dual_graph
  )

  h = hypersurface_equation(hs_model)
  S = cox_ring(ambient_space(hs_model))
  var_names = symbols(parent(h))
  S.R.S = var_names
  hs_model.hypersurface_equation = Oscar.eval_poly(string(h), S)

  return QSMModel(vertices,
                  poly_index,
                  triang_quick,
                  max_lattice_pts_in_facet,
                  estimated_number_of_triangulations,
                  hs_model,
                  Kbar3,
                  h11,
                  h12,
                  h13,
                  h22,
                  genus_ci,
                  degree_of_Kbar_of_tv_restricted_to_ci,
                  intersection_number_among_ci_cj,
                  index_facet_interior_divisors,
                  intersection_number_among_nontrivial_ci_cj,
                  dual_graph,
                  components_of_dual_graph,
                  degree_of_Kbar_of_tv_restricted_to_components_of_dual_graph,
                  genus_of_components_of_dual_graph,
                  simplified_dual_graph,
                  components_of_simplified_dual_graph,
                  degree_of_Kbar_of_tv_restricted_to_components_of_simplified_dual_graph,
                  genus_of_components_of_simplified_dual_graph)
end
