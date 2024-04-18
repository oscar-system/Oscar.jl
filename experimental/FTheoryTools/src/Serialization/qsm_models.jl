@register_serialization_type QSMModel

function save_object(s::SerializerState, qsm::QSMModel)
  save_data_dict(s) do
    save_object(s, qsm.vertices, :vertices)
    save_object(s, qsm.poly_index, :poly_inx)

    save_object(s, qsm.triang_quick, :triang_quick)
    save_object(s, qsm.max_lattice_pts_in_facet, :max_lattice_pts_in_facet)
    save_object(s, qsm.estimated_number_of_triangulations, :estimated_number_of_triangulations)
    save_object(s, qsm.hs_model, :hs_model)

    save_object(s, qsm.Kbar3, :Kbar3)
    save_object(s, qsm.h11, :h11)
    save_object(s, qsm.h12, :h12)
    save_object(s, qsm.h13, :h13)
    save_object(s, qsm.h22, :h22)

    save_object(s, qsm.genus_ci, :genus_ci)
    save_object(s, qsm.degree_of_Kbar_of_tv_restricted_to_ci, :degree_of_Kbar_of_tv_restricted_to_ci)
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
