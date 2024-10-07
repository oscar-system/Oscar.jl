# Updates to the QSM Database might be needed from time to time.
# Currently, this data is added as tar.gz folder to the following release tag of OSCAR:
# https://github.com/oscar-system/Oscar.jl/releases/tag/archive-tag-1
# When updating, note that you must update Artifacts.toml in the OSCAR root folder with
# the latest sha256 and git-tree-sha1, which can be computed as outlined here:
# https://pkgdocs.julialang.org/v1/artifacts/
# For recreating/updating the data, you could execute code like the following on the latest master branch of OSCAR:

#=
for k in 1:5000
  try
    qsm_model = literature_model(arxiv_id = "1903.00009", model_parameters = Dict("k" => k))
    equ = hypersurface_equation(qsm_model)
    S = cox_ring(ambient_space(qsm_model))
    my_ring_map = hom(parent(equ), S, gens(S))
    println("Working on k = $k...")
    var_names = symbols(parent(equ))
    S.R.S = var_names
    qsm_model.hypersurface_equation = my_ring_map(equ)
    if is_empty(explicit_model_sections(qsm_model)) == false 
      println("Problem with explicit model sections!")
    end
    if my_ring_map(equ) != hypersurface_equation(qsm_model)
      println("Problem with mapping of hypersurface equation!")
    end
    if parent(degree(hypersurface_equation(qsm_model))) != class_group(ambient_space(qsm_model))
      println("Problem with class group!")
    end
    if cox_ring(ambient_space(qsm_model)).D != class_group(ambient_space(qsm_model))
      println("Inconsistency among grading group of cox ring and class group!")
    end
    if parent(hypersurface_equation(qsm_model)) != cox_ring(ambient_space(qsm_model))
      println("Problem with parent of hypersurface equation!")
    end
    if is_calabi_yau(qsm_model, check = false) == false
      println("Problem with is_calabi_yau!")
    end
    corresponding_qsm_model = QSMModel(
      vertices(qsm_model),
      polytope_index(qsm_model),
      has_quick_triangulation(qsm_model),
      max_lattice_pts_in_facet(qsm_model),
      estimated_number_of_triangulations(qsm_model),
      qsm_model,
      kbar3(qsm_model),
      hodge_h11(qsm_model),
      hodge_h12(qsm_model),
      hodge_h13(qsm_model),
      hodge_h22(qsm_model),
      genera_of_ci_curves(qsm_model),
      degrees_of_kbar_restrictions_to_ci_curves(qsm_model),
      topological_intersection_numbers_among_ci_curves(qsm_model),
      indices_of_trivial_ci_curves(qsm_model),
      topological_intersection_numbers_among_nontrivial_ci_curves(qsm_model),
      dual_graph(qsm_model),
      components_of_dual_graph(qsm_model),
      degrees_of_kbar_restrictions_to_components_of_dual_graph(qsm_model),
      genera_of_components_of_dual_graph(qsm_model),
      simplified_dual_graph(qsm_model),
      components_of_simplified_dual_graph(qsm_model),
      degrees_of_kbar_restrictions_to_components_of_simplified_dual_graph(qsm_model),
      genera_of_components_of_simplified_dual_graph(qsm_model)
    )
    save("/YOUR_CHOSEN_PATH/$(k).mrdi", corresponding_qsm_model)
    println("Saved...")
    println("")
  catch e
    #println("Error for k = $k: $e. Skipping...")
    continue
  end
end
=#
  
# To verify, that the new data can be read, switch to your development branch and execute something like the following:
  
#=
for k in 1:500
  try
    qsm_model = load("/YOUR_CHOSEN_PATH/$(k).mrdi")
    println("Working on k = $k...")
    hs_model = qsm_model.hs_model
    set_attribute!(ambient_space(hs_model), :class_group, cox_ring(ambient_space(hs_model)).D)
    if parent(degree(hypersurface_equation(hs_model))) != class_group(ambient_space(hs_model))
      println("Problem with class group!")
    end
    if parent(hypersurface_equation(hs_model)) != cox_ring(ambient_space(hs_model))
      println("Problem with parent of hypersurface equation!")
    end
    if is_calabi_yau(hs_model, check = false) == false
      println("Problem with is_calabi_yau!")
    end
    println("Success for k = $k...")
    println("")
  catch e
    #println("Error for k = $k: $e. Skipping...")
    continue
  end
end
=#

# Do not forget to zip the new data and update https://pkgdocs.julialang.org/v1/artifacts/.

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

  #h = hypersurface_equation(hs_model)
  #S = cox_ring(ambient_space(hs_model))
  #var_names = symbols(parent(h))
  #S.R.S = var_names
  #hs_model.hypersurface_equation = Oscar.eval_poly(string(h), S)

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
