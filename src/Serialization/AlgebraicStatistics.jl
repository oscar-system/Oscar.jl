@register_serialization_type PhylogeneticModel uses_params 
@register_serialization_type GroupBasedPhylogeneticModel uses_params 

function save_object(s::SerializerState, pm::PhylogeneticModel)
    save_data_dict(s) do 
        save_object(s, pm.graph, :tree) # already implemented 
        save_object(s, pm.n_states, :states) # already implemented 
        save_object(s, pm.prob_ring, :probability_ring) # already implemented
        save_object(s, pm.root_distr, :root_distribution) # needs to be implemented, or we change the entry type to QQFieldElem
        save_object(s, pm.trans_matrices, :transition_matrices) # needs to be implemented
    end
end

function save_object(s::SerializerState, pm::GroupBasedPhylogeneticModel)
    save_data_dict(s) do
        save_object(s,pm.phylo_model,:phylomodel) # details see above
        save_object(s, pm.fourier_ring, :fourier_ring) # already implemented
        save_object(s,pm.fourier_params,:fourier_params) # needs to be implemented
        save_object(s,pm.group,:group) # already implemented
    end
end


function load_object(s::DeserializerState, ::Type{<: PhylogeneticModel})
    graph = load_object(s, Graph{Directed}, :tree)
    n_states = load_object(s, Int64, :states)
    prob_ring = load_object(s, QQMPolyRing, :probability_ring)
    root_distr = load_object(s, Vector{Any}, :root_distribution)
    trans_matrices = load_object(s, Dict{Edge, MatElem{QQMPolyRingElem}}, :trans_matrices)
    
    return PhylogeneticModel(graph, n_states, prob_ring, root_distr, trans_matrices)
end 

function load_object(s::DeserializerState, ::Type{<: GroupBasedPhylogeneticModel})
    phylomodel = load_object(s, PhylogeneticModel, :phylomodel)
    fourier_ring = load_object(s, QQMPolyRing, :fourier_ring)
    fourier_params = load_object(s, Dict{Edge, Vector{QQMPolyRingElem}} , :fourier_params)
    group = load_object(s, Vector{FinGenAbGroupElem}, :group)
    
    return GroupBasedPhylogeneticModel(phylomodel, fourier_ring, fourier_params, group)
end

### Structure of both Types: ###

#= struct PhylogeneticModel
  graph::Graph{Directed}
  n_states::Int
  prob_ring::MPolyRing{QQFieldElem}
  root_distr::Vector{Any}
  trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}
end

struct GroupBasedPhylogeneticModel
  phylo_model::PhylogeneticModel
  fourier_ring::MPolyRing{QQFieldElem}
  fourier_params::Dict{Edge, Vector{QQMPolyRingElem}}
  group::Vector{FinGenAbGroupElem}
end =#
