################################################
#### PHYLOGENETIC DATA STRUCTURES & METHODS ####
################################################

# -------------------- #
## Phylogenetic Model ##
# -------------------- #

@attributes mutable struct PhylogeneticModel{L} <: GraphicalModel{Directed, L}
  # do need to for T to be directed here? YES!
  base_field::Field
  graph::Graph{Directed}
  labelings::L
  trans_mat_structure::Matrix{<: VarName}
  root_distribution::Vector
  n_states::Int
  model_parameter_name::VarName

  
  function PhylogeneticModel(F::Field,
                             G::Graph{Directed},
                             trans_mat_structure::Matrix{<: VarName},
                             root_distribution::Union{Nothing, Vector} = nothing,
                             n_states::Union{Nothing, Int} = nothing,
                             varname::VarName="p")
    if isnothing(n_states)
      n_states = size(trans_mat_structure)[1]
    end
    if isnothing(root_distribution)
      root_distribution = repeat([1//n_states], outer = n_states)
    end
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{typeof(graph_maps)}(F, G,
                                   graph_maps,
                                   trans_mat_structure,
                                   root_distribution,
                                   n_states,
                                   varname)
  end

  function PhylogeneticModel(G::Graph{Directed},
                             trans_mat_structure::Matrix{<: VarName},
                             root_distribution::Union{Nothing, Vector} = nothing,
                             n_states::Union{Nothing, Int} = nothing,
                             varname::VarName="p")
    return PhylogeneticModel(QQ, G, trans_mat_structure, root_distribution, n_states, varname)
  end
end

n_states(PM::PhylogeneticModel) = PM.n_states
transition_matrix(PM::PhylogeneticModel) = PM.trans_mat_structure
base_field(PM::PhylogeneticModel) = PM.base_field
varname(PM::PhylogeneticModel) = PM.model_parameter_name
root_distribution(PM::PhylogeneticModel) = PM.root_distribution

@attr Tuple{MPolyRing, GenDict} function parameter_ring(PM::PhylogeneticModel; cached=false)
  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), (j,e) in enumerate(sort_edges(graph(PM)))
  )
end

@attr Tuple{MPolyRing, Array} function model_ring(PM::PhylogeneticModel; cached=false)
  leave_indices = leaves_indices(PM)

  return polynomial_ring(base_field(PM),
                         ["$(varname(PM))[$(join(x, ", "))]" for x in leave_indices];
                         cached=cached)
end

function entry_transition_matrix(PM::PhylogeneticModel, i::Int, j::Int, e::Edge)
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][tr_mat[i,j], e]
end

function entry_transition_matrix(PM::PhylogeneticModel, i::Int, j::Int, u::Int, v::Int,)
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][tr_mat[i,j], Edge(u,v)]
end

@attr MPolyAnyMap function parametrization(PM::PhylogeneticModel)
  gr = graph(PM)

  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = leaves_indices(PM::PhylogeneticModel)

  map = [leaves_probability(PM, Dict(i => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))

end

function Base.show(io::IO, pm::PhylogeneticModel)
  gr = graph(pm)
  ns = n_states(pm)
  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  #root_dist = join(Oscar.root_distribution(pm), ", " )
  #M = collect(values(transition_matrix(pm)))[1]

  # commenting out root distribution for now
  print(io, "Phylogenetic model on a tree with $(nl) leaves and $(ne) edges") # \n with distribution at the root [$(root_dist)] \n")
  print(io, " and transition matrices of the form \n ")

  # printing this matrix can probably be improved
  print(io, "$(pm.trans_mat_structure)")
end


# -------------------------------- #
## Group-based Phylogenetic Model ##
# -------------------------------- #

@attributes mutable struct GroupBasedPhylogeneticModel{L} <: GraphicalModel{Directed, L}
  # do need to for T to be directed here? YES!
  phylo_model::PhylogeneticModel 
  fourier_param_structure::Vector{<: VarName}
  group::Vector{FinGenAbGroupElem}
  model_parameter_name::VarName
  
  function GroupBasedPhylogeneticModel(F::Field, 
                                       G::Graph{Directed},
                                       trans_mat_structure::Matrix{<: VarName},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       root_distribution::Union{Nothing, Vector} = nothing,
                                       n_states::Union{Nothing, Int} = nothing,
                                       varname_phylo_model::VarName="p",
                                       varname_group_based::VarName="q")
    if isnothing(group)
      group = collect(abelian_group(2,2))
      group = [group[1],group[3],group[2],group[4]]
    end

    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{typeof(graph_maps)}(PhylogeneticModel(F, G, trans_mat_structure, 
                                                     root_distribution, n_states, 
                                                     varname_phylo_model),
                                   fourier_param_structure,
                                   group, varname_group_based)
  end

  # F, G,
  # trans_mat_structure, fourier_param_structure,
  # n_states, root_distribution, group,
  # varname_phylo_model, varname_group_based

  function GroupBasedPhylogeneticModel(G::Graph{Directed},
                                       trans_mat_structure::Matrix{<: VarName},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       root_distribution::Union{Nothing, Vector} = nothing,
                                       n_states::Union{Nothing, Int} = nothing,
                                       varname_phylo_model::VarName="p",
                                       varname_group_based::VarName="q")
    return GroupBasedPhylogeneticModel(QQ, G, trans_mat_structure, fourier_param_structure,
                                       group, root_distribution, n_states, 
                                       varname_phylo_model, varname_group_based)
  end

  # # ? Antony
  # function GroupBasedPhylogeneticModel(G::Graph{Directed},
  #                                      trans_mat_structure::Matrix{<: VarName},
  #                                      fourier_param_structure::Vector{<: VarName},
  #                                      group::Vector{FinGenAbGroupElem},
  #                                      root_distribution::Union{Nothing, Vector} = nothing,
  #                                      n_states::Union{Nothing, Int} = nothing,
  #                                      varname_phylo_model::VarName="p",
  #                                      varname_group_based::VarName="q")
  #   return GroupBasedPhylogeneticModel(QQ, G, trans_mat_structure, fourier_param_structure,
  #                                      root_distribution, n_states, group, 
  #                                      varname_phylo_model, varname_group_based)
  # end


  # can I do smth like this?
  # function GroupBasedPhylogeneticModel(G::Graph{Directed},
  #                                      trans_mat_structure::Matrix{<: VarName},
  #                                      fourier_param_structure::Vector{<: VarName},
  #                                      group::Vector{FinGenAbGroupElem} = nothing,
  #                                      varname_group_based::VarName="q")
  #   return GroupBasedPhylogeneticModel(QQ, G, trans_mat_structure, nothing, nothing, "p",
  #                                      fourier_param_structure, group, varname_group_based)
  # end

  #Is this necessary? Or it's the same as the struct constructor?
  # function GroupBasedPhylogeneticModel(PM::PhylogeneticModel,
  #                                      fourier_param_structure::Vector{<: VarName},
  #                                      group::Vector{FinGenAbGroupElem},
  #                                      varname_group_based::VarName="q")
  #   return GroupBasedPhylogeneticModel(base_field(PM), graph(PM), transition_matrix(PM),
  #                                      n_states(PM), root_distribution(PM), varname(PM),
  #                                      fourier_param_structure, group, varname_group_based)
  # end

end

graph(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.graph # ? Antony

n_states(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.n_states
transition_matrix(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.trans_mat_structure
base_field(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.base_field
root_distribution(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.root_distribution
varname_probabilities(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.model_parameter_name # ? Antony

phylogenetic_model(PM::GroupBasedPhylogeneticModel) = PM.phylo_model
group(PM::GroupBasedPhylogeneticModel) = PM.group
fourier_parameters(PM::GroupBasedPhylogeneticModel) = PM.fourier_param_structure
varname_fourier(PM::GroupBasedPhylogeneticModel) = PM.model_parameter_name # ? Antony


@attr Tuple{MPolyRing, GenDict} function parameter_ring(PM::GroupBasedPhylogeneticModel; cached=false)
  vars = unique(fourier_parameters(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), (j,e) in enumerate(sort_edges(graph(PM)))
  )
end

@attr Tuple{
  MPolyRing, 
  Array
  } function model_ring(PM::GroupBasedPhylogeneticModel; cached=false)
  leave_indices = leaves_indices(PM) # ? Antony

  return polynomial_ring(base_field(PM),
                         ["$(varname_fourier(PM))[$(join(x, ", "))]" for x in leave_indices];
                         cached=cached)
end

function entry_fourier_parameter(PM::GroupBasedPhylogeneticModel, i::Int, e::Edge,)
  x = fourier_parameters(PM)
  parameter_ring(PM)[2][x[i], e]
end

@attr MPolyAnyMap function parametrization(PM::GroupBasedPhylogeneticModel)
  gr = graph(PM)

  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = leaves_indices(PM)

  map = [leaves_fourier(PM, Dict(i => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))

end

function Base.show(io::IO, pm::GroupBasedPhylogeneticModel)
  gr = graph(pm)
  ns = n_states(pm)
  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))

  # commenting out root distribution for now
  print(io, "Group-based phylogenetic model on a tree with $(nl) leaves and $(ne) edges") # \n with distribution at the root [$(root_dist)] \n")
  
  # printing this matrix can probably be improved
  print(io, " with transition matrices of the form \n ")
  print(io, "$(pm.phylo_model.trans_mat_structure) \n")

  print(io, " and fourier parameters of the form \n ")
  print(io, "$(pm.fourier_param_structure)")
end

## Equivalent classes ##
# -------------------- #

const EquivDict = Dict{T, Vector{T}} where T <: MPolyRingElem
@attr EquivDict function equivalent_classes(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}) 
  
  _, ps = model_ring(PM)
  param = parametrization(PM);
  
  polys = unique(param.img_gens)
  polys = polys[findall(!is_zero, polys)]

  #equivalent_classes = Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}()
  equivalent_classes = Dict{MPolyRingElem, Vector{MPolyRingElem}}()
  for poly in polys
      eqv_class = sort([p for p in ps if param(p) == poly], rev = true)
      # equivalent_classes[Tuple(index(eqv_class[1]))] = eqv_class
      equivalent_classes[eqv_class[1]] = eqv_class
  end

  return equivalent_classes
end

@attr Tuple{
  MPolyRing, 
  Array
  } function reduced_model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
  eq_calasses = Oscar.equivalent_classes(PM)
  keys_eq_classes = sort(collect(keys(eq_calasses)), rev = true)

  return polynomial_ring(base_field(PM),
                        [Symbol(k) for k in keys_eq_classes];
                        cached=cached)
end

function reduced_parametrization(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel})
  _, p = model_ring(PM)
  S, _ = parameter_ring(PM)
  param = parametrization(PM)

  redR, _ = reduced_model_ring(PM)
  gens_redR = gens(redR)
  keys_eq_classes = index.(gens_redR)
  
  probabilities = [param(p[k...]) for k in keys_eq_classes]
  hom(redR, S, reduce(vcat, probabilities))
end

## Change of coordinates p ↔ q ##
# ----------------------------- #

function reduced_fourier_transform(PM::GroupBasedPhylogeneticModel)
  _, p = Oscar.model_ring(phylogenetic_model(PM))
  _, q = Oscar.model_ring(PM)

  Rp_red, _ = Oscar.reduced_model_ring(phylogenetic_model(PM))
  Rq_red, _ = Oscar.reduced_model_ring(PM)

  p_classes = Oscar.equivalent_classes(phylogenetic_model(PM))
  q_classes = Oscar.equivalent_classes(PM)

  np = length(p_classes)
  nq = length(q_classes)
  
  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  # keys_p_classes = collect(keys(p_classes))
  # sort!(keys_p_classes, rev = true)
  # keys_q_classes = collect(keys(q_classes))
  # sort!(keys_q_classes, rev = true)

  keys_p_classes = gens(Rp_red)
  keys_q_classes = gens(Rq_red)

  
  H = Rp_red.(hadamard(matrix_space(ZZ, n_states(PM), n_states(PM))))
  
  M = Rp_red.(Int.(zeros(nq, np)))
  for i in 1:nq
    q_key_index = index(keys_q_classes[i])
    current_fourier_class = q_classes[q[q_key_index...]]
    for j in 1:np
      p_key_index = index(keys_p_classes[j])
      current_prob_class = p_classes[p[p_key_index...]]
      current_entries_in_M = [prod([H[y,x] for (x,y) in zip(index(pp),index(qq))]) for pp in current_prob_class, qq in current_fourier_class]
      M[i,j] = Rp_red.(1//(length(current_prob_class)*length(current_fourier_class))*sum(current_entries_in_M))
    end
  end
  
  return M
end

function reduced_coordinate_change(PM::GroupBasedPhylogeneticModel)
  Rp_red, _ = Oscar.reduced_model_ring(phylogenetic_model(PM))
  Rq_red, _ = Oscar.reduced_model_ring(PM)
  
  M = reduced_fourier_transform(PM)
  hom(Rq_red, Rp_red, M*gens(Rp_red))

end

function inverse_reduced_fourier_transform(PM::GroupBasedPhylogeneticModel)
    _, p = Oscar.model_ring(phylogenetic_model(PM))
    _, q = Oscar.model_ring(PM)

    Rp_red, _ = Oscar.reduced_model_ring(phylogenetic_model(PM))
    Rq_red, _ = Oscar.reduced_model_ring(PM)

    p_classes = Oscar.equivalent_classes(phylogenetic_model(PM))
    q_classes = Oscar.equivalent_classes(PM)

    np = length(p_classes)
    nq = length(q_classes)
    
    ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
    # keys_p_classes = collect(keys(p_classes))
    # sort!(keys_p_classes, rev = true)
    # keys_q_classes = collect(keys(q_classes))
    # sort!(keys_q_classes, rev = true)

    keys_p_classes = gens(Rp_red)
    keys_q_classes = gens(Rq_red)

    
    H = Rq_red.(hadamard(matrix_space(ZZ, n_states(PM), n_states(PM))))
    Hinv = 1//n_states(PM) * H

    M = Rq_red.(Int.(zeros(np, nq)))
    for i in 1:np
      p_key_index = index(keys_p_classes[i])
      current_prob_class = p_classes[p[p_key_index...]]

      for j in 1:nq
        q_key_index = index(keys_q_classes[j])
        current_fourier_class = q_classes[q[q_key_index...]]
        
        current_entries_in_M = [prod([Hinv[x,y] for (x,y) in zip(index(pp),index(qq))]) for pp in current_prob_class, qq in current_fourier_class]

        M[i,j] = Rq_red.(sum(current_entries_in_M))
        
      end
    end
    M

    return M
end

function inverse_reduced_coordinate_change(PM::GroupBasedPhylogeneticModel)
  Rp_red, _ = Oscar.reduced_model_ring(phylogenetic_model(PM))
  Rq_red, _ = Oscar.reduced_model_ring(PM)
  
  M = inverse_reduced_fourier_transform(PM)
  hom(Rp_red, Rq_red, M*gens(Rq_red))

end


#########################
#### SPECIFIC MODELS ####
#########################

function cavender_farris_neyman_model(G::Graph{Directed})
  M = [:a :b;
       :b :a]
  x = [:x, :y] 

  group = collect(abelian_group(2))
  group = [group[1],group[2]]
  
  #PhylogeneticModel(G, M)
  GroupBasedPhylogeneticModel(G, M, x, group) 
end

function jukes_cantor_model(G::Graph{Directed})
  M = [:a :b :b :b;
       :b :a :b :b;
       :b :b :a :b;
       :b :b :b :a]
  # x = [Symbol("x[1]"), Symbol("x[2]"), Symbol("x[2]"), Symbol("x[2]")]
  # x = [:x1, :x2, :x2, :x2] #? Antony
  x = [:x, :y, :y, :y] 
  
  #PhylogeneticModel(G, M)
  GroupBasedPhylogeneticModel(G, M, x) 
end

function kimura2_model(G::Graph{Directed})
  M = [:a :b :c :b;
       :b :a :b :c;
       :c :b :a :b;
       :b :c :b :a]
  x = [:x, :y, :z, :z] 
  
  GroupBasedPhylogeneticModel(G, M, x)
end

function kimura3_model(G::Graph{Directed})
  M = [:a :b :c :d;
       :b :a :d :c;
       :c :d :a :b;
       :d :c :b :a]
  x = [:x, :y, :z, :t] 
  
  GroupBasedPhylogeneticModel(G, M, x)
end

function general_markov_model(G::Graph{Directed})

  M = [:m11 :m12 :m13 :m14;
       :m21 :m22 :m23 :m24;
       :m31 :m32 :m33 :m34;
       :m41 :m42 :m43 :m44]
  
  root_distr = [:π1, :π2, :π3, :π4]
  
  PhylogeneticModel(G, M)
end

function general_markov_model(G::Graph{Directed}, M::Matrix{<: VarName})
  PhylogeneticModel(G, M)
end
