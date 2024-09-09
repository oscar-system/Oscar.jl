@register_serialization_type PhylogeneticModel 
@register_serialization_type GroupBasedPhylogeneticModel 

function save_object(s::SerializerState, pm::PhylogeneticModel)
  save_data_dict(s) do 
    save_typed_object(s, graph(pm), :tree) 
    save_object(s, number_states(pm), :states) 
    save_typed_object(s, probability_ring(pm), :probability_ring) 
    save_typed_object(s, root_distribution(pm), :root_distribution) 
    save_typed_object(s, transition_matrices(pm), :transition_matrices) 
  end
end

function save_object(s::SerializerState, pm::GroupBasedPhylogeneticModel)
  save_data_dict(s) do
    save_typed_object(s, phylogenetic_model(pm), :phylomodel) 
    save_typed_object(s, fourier_ring(pm), :fourier_ring) 
    save_typed_object(s, fourier_parameters(pm), :fourier_params) 
    save_typed_object(s, group_of_model(pm), :group) 
  end
end

function load_object(s::DeserializerState, g::Type{PhylogeneticModel})
  graph = load_typed_object(s, :tree)
  n_states = load_object(s, Int64, :states)
  prob_ring = load_typed_object(s, :probability_ring)
  root_distr = load_typed_object(s, :root_distribution)
  trans_matrices = load_typed_object(s, :transition_matrices)
  
  return PhylogeneticModel(graph, n_states, prob_ring, root_distr, trans_matrices)
end 

function load_object(s::DeserializerState, g::Type{GroupBasedPhylogeneticModel})
  phylomodel = load_typed_object(s, :phylomodel)
  fourier_ring = load_typed_object(s, :fourier_ring)
  fourier_params = load_typed_object(s, :fourier_params)
  group = load_typed_object(s, :group)
    
  return GroupBasedPhylogeneticModel(phylomodel, fourier_ring, fourier_params, group)
end
