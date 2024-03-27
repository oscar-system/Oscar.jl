#### PARAMETRISATION IN PROBABLITY COORDINATES #### 


function monomial_parametrization(pm::PhylogeneticModel, states::Dict{Int, Int})
    gr = graph(pm)
    tr_mat = transition_matrices(pm)
  
    monomial = 1
    for edge in edges(gr)
        stateParent = states[src(edge)]
        stateChild = states[dst(edge)]
        monomial = monomial * tr_mat[edge][stateParent, stateChild]
    end
  
    return(monomial)
end
  
function probability_parametrization(pm::PhylogeneticModel, leaves_states::Vector{Int})
    gr = graph(pm)
    int_nodes = interior_nodes(gr)
    lvs_nodes = leaves(gr)
    n_states = number_states(pm)
  
    interior_indices = collect.(Iterators.product([collect(1:n_states) for _ in int_nodes]...))  
    nodes_states = Dict(lvs_nodes[i] => leaves_states[i] for i in 1:length(lvs_nodes))
  
    poly = 0
    # Might be useful in the future to use a polynomial ring context
    for labels in interior_indices
      for (int_node, label) in zip(int_nodes, labels)
        nodes_states[int_node] = label
      end
      poly = poly + monomial_parametrization(pm, nodes_states)
    end 
    return poly
end 
  
function probability_map(pm::PhylogeneticModel)
    lvs_nodes = leaves(graph(pm))
    n_states = number_states(pm)
  
    leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
    probability_coordinates = Dict(leaves_states => probability_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
    return(probability_coordinates)
end


#### FOURIER PARAMETRISATION #### 
  

function monomial_fourier(pm::PhylogeneticModel, leaves_states::Vector{Int})
    gr = graph(pm)
    param = fourier_parameters(pm)
    monomial = 1
    for edge in edges(gr)
      dsc = vertex_descendants(dst(edge), gr, [])
      elem = group_sum(pm, leaves_states[dsc])
      monomial = monomial * param[edge][which_group_element(pm, elem)]
    end
  
    return(monomial)
end
  
function fourier_parametrization(pm::PhylogeneticModel, leaves_states::Vector{Int})
    S = fourier_ring(pm)
    if group_sum(pm, leaves_states) == [0,0]
      poly = monomial_fourier(pm, leaves_states)
    else 
      poly = S(0)
    end
  
    return poly
end 
  
function fourier_map(pm::PhylogeneticModel)
    lvs_nodes = leaves(graph(pm))
    n_states = number_states(pm)
  
    leaves_indices = collect.(Iterators.product([collect(1:n_states) for _ in lvs_nodes]...))
    fourier_coordinates = Dict(leaves_states => fourier_parametrization(pm, leaves_states) for leaves_states in leaves_indices)
    return(fourier_coordinates)
end

