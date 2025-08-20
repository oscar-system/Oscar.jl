################################################
#### PHYLOGENETIC DATA STRUCTURES & METHODS ####
################################################


###################################################################################
#
#       Phylogenetic Models
#
###################################################################################


@attributes mutable struct PhylogeneticModel{GT, M, L, T} <: GraphicalModel{GT, L}
  # do need to for T to be directed here? YES!
  base_field::Field
  graph::GT
  labelings::L
  trans_matrix_structure::Matrix{M}
  root_distribution::Vector{T}
  n_states::Int
  model_parameter_name::VarName

  
  function PhylogeneticModel(F::Field,
                             G::Graph{Directed},
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p")
    n_states = size(trans_matrix_structure)[1]
    if isnothing(root_distribution)
      root_distribution = F.(repeat([1//n_states], outer = n_states))
    end
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{Graph{Directed}, typeof(first(trans_matrix_structure)), 
               typeof(graph_maps), typeof(first(root_distribution))}(F, G,
                                                                    graph_maps,
                                                                    trans_matrix_structure,
                                                                    root_distribution,
                                                                    n_states,
                                                                    varname)
  end

  function PhylogeneticModel(G::Graph{Directed},
                             trans_matrix_structure::Matrix{M},
                             root_distribution::Vector{T},
                             varname::VarName="p") where {M <:  MPolyRing{<:MPolyRingElem}, T <: MPolyRing}
    trans_ring = parent(first(trans_matrix_structure))
    root_ring = base_ring(trans_ring)
    F = coefficient_ring(root_ring)
    
    if root_ring != parent(first(root_distribution))
      error("Rings for the root distribution parameters do not match. Expected ring `$(root_ring)` but got `$(parent(first(root_distribution)))`.")
    end

    PhylogeneticModel(F, G, trans_matrix_structure, root_distribution, varname)
  end

  function PhylogeneticModel(G::Graph{Directed},
                             trans_matrix_structure::Matrix,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p")
    return PhylogeneticModel(QQ, G, trans_matrix_structure, 
                             root_distribution, varname)
  end
end

n_states(PM::PhylogeneticModel) = PM.n_states
transition_matrix(PM::PhylogeneticModel) = PM.trans_matrix_structure
base_field(PM::PhylogeneticModel) = PM.base_field
varname(PM::PhylogeneticModel) = PM.model_parameter_name
root_distribution(PM::PhylogeneticModel) = PM.root_distribution

@attr Tuple{
  MPolyRing, 
  GraphTransDict,
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: VarName, L, T}; cached=false) where {L, T <: FieldElem} 
  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)
  
  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(graph(PM)))
      ), root_distribution(PM)
end

@attr Tuple{
  MPolyRing,
  GraphTransDict,
  Vector{<:MPolyRingElem}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: VarName}; cached=false)
  vars = unique(transition_matrix(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, r, x... = polynomial_ring(base_field(PM),
                               root_distribution(PM),
                               edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), 
      (j,e) in enumerate(sort_edges(graph(PM)))
      ), r
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: FieldElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)
  
  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM)))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, map)
  end

  R, dict_maps, root_distribution(PM)
end

@attr Tuple{
  MPolyRing, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: MPolyRingElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, rv, x... = polynomial_ring(base_field(PM), Symbol.(root_vars), edge_gens..., ; cached=cached)

  coef_map = hom(coefficient_ring(trans_ring), R, rv)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM)))
      map = [x[i][j] for i in 1:length(transition_vars)]
      dict_maps[e] = hom(trans_ring, R, coef_map, map)
  end

  R, dict_maps, rv

end

@attr Tuple{
  MPolyRing{T}, 
  Dict{Edge, Oscar.MPolyAnyMap}, 
  Vector{T}
} function parameter_ring(PM::PhylogeneticModel{Graph{Directed}, <: MPolyRingElem, L, T}; cached=false) where {L, T <: AbstractAlgebra.Generic.RationalFunctionFieldElem}
  trans_ring = parent(first(transition_matrix(PM)))
  transition_vars = gens(trans_ring)
  root_vars = gens(coefficient_ring(trans_ring))

  edge_gens = [x => 1:n_edges(graph(PM)) for x in Symbol.(transition_vars)]
  R, x... = polynomial_ring(coefficient_ring(trans_ring), edge_gens..., ; cached=cached)

  dict_maps = Dict{Edge, Oscar.MPolyAnyMap}()
  for (j,e) in enumerate(Oscar.sort_edges(graph(PM)))
    map = [x[i][j] for i in 1:length(transition_vars)]
    dict_maps[e] = hom(trans_ring, R, map)
  end
  R, dict_maps, root_vars
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
  print(io, "$(pm.trans_matrix_structure)")
end


###################################################################################
#
#       Group-based phylogenetic models
#
###################################################################################

@attributes mutable struct GroupBasedPhylogeneticModel{GT, L} <: GraphicalModel{GT, L}
  phylo_model::PhylogeneticModel{GT, <: VarName, L, <: RingElem} 
  fourier_param_structure::Vector{<: VarName}
  group::Vector{FinGenAbGroupElem}
  model_parameter_name::VarName
  
  function GroupBasedPhylogeneticModel(F::Field, 
                                       G::Graph{Directed},
                                       trans_matrix_structure::Matrix{<: VarName},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       root_distribution::Union{Nothing, Vector} = nothing,
                                       varname_phylo_model::VarName="p",
                                       varname_group_based::VarName="q")
    if isnothing(group)
      group = collect(abelian_group(2,2))
      group = [group[1],group[3],group[2],group[4]]
    end

    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{typeof(G), typeof(graph_maps)}(PhylogeneticModel(F, G, trans_matrix_structure, 
                                                                root_distribution,
                                                                varname_phylo_model),
                                              fourier_param_structure,
                                              group, varname_group_based)
  end

  function GroupBasedPhylogeneticModel(G::Graph{Directed},
                                       trans_matrix_structure::Matrix{<: VarName},
                                       fourier_param_structure::Vector{<: VarName},
                                       group::Union{Nothing, Vector{FinGenAbGroupElem}} = nothing,
                                       root_distribution::Union{Nothing, Vector} = nothing,
                                       varname_phylo_model::VarName="p",
                                       varname_group_based::VarName="q")
    return GroupBasedPhylogeneticModel(QQ, G, trans_matrix_structure, fourier_param_structure,
                                       group, root_distribution, 
                                       varname_phylo_model, varname_group_based)
  end

end

graph(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.graph # ? Antony

n_states(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.n_states
transition_matrix(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.trans_matrix_structure
base_field(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.base_field
root_distribution(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.root_distribution
varname_probabilities(PM::GroupBasedPhylogeneticModel) = PM.phylo_model.model_parameter_name # ? Antony

phylogenetic_model(PM::GroupBasedPhylogeneticModel) = PM.phylo_model
group(PM::GroupBasedPhylogeneticModel) = PM.group
fourier_parameters(PM::GroupBasedPhylogeneticModel) = PM.fourier_param_structure
varname_fourier(PM::GroupBasedPhylogeneticModel) = PM.model_parameter_name # ? Antony
varname(PM::GroupBasedPhylogeneticModel) = PM.model_parameter_name



@attr Tuple{
  MPolyRing,
  GraphTransDict
} function parameter_ring(PM::GroupBasedPhylogeneticModel; cached=false)
  vars = unique(fourier_parameters(PM))
  edge_gens = [x => 1:n_edges(graph(PM)) for x in vars]
  R, x... = polynomial_ring(base_field(PM), edge_gens...; cached=cached)

  R, Dict{Tuple{VarName, Edge}, MPolyRingElem}(
    (vars[i], e) => x[i][j] for i in 1:length(vars), (j,e) in enumerate(sort_edges(graph(PM)))
  )
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
  print(io, "$(pm.phylo_model.trans_matrix_structure) \n")

  print(io, " and fourier parameters of the form \n ")
  print(io, "$(pm.fourier_param_structure)")
end


###################################################################################
#
#       Parametrizations
#
###################################################################################

@attr Tuple{
  ModelRing{T, U}, 
  Dict{T, U}
} where {T, U} function full_model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
  l_indices = leaves_indices(PM)

  return model_ring(base_field(PM), varname(PM) => l_indices; cached=cached)
end

@attr MPolyAnyMap function full_parametrization(PM::PhylogeneticModel)
  gr = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = leaves_indices(PM::PhylogeneticModel)

  map = [leaves_probability(PM, Dict(i => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))

end

function full_parametrization(PM::GroupBasedPhylogeneticModel)
  gr = graph(PM)

  R, _ = full_model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = leaves_indices(PM)

  map = [leaves_fourier(PM, Dict(i => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))

end

#const EquivDict = Dict{T, Vector{T}} where T <: MPolyRingElem
@attr Dict{
      Tuple{Vararg{Int64}}, 
      Vector{MPolyRingElem}
} function equivalent_classes(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}) 
  _, ps = full_model_ring(PM)
  f = full_parametrization(PM);
  
  polys = unique(f.img_gens)
  polys = polys[findall(!is_zero, polys)]

  equivalent_classes = Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}()
  # equivalent_classes = Dict{MPolyRingElem, Vector{MPolyRingElem}}()
  for poly in polys
      eqv_class = sort([i for i in keys(ps) if f(ps[i]) == poly], rev = false)
      equivalent_classes[eqv_class[1]] = [ps[i...] for i in eqv_class]
  end
  return equivalent_classes
end

@attr Tuple{
  ModelRing{T, U}, 
  Dict{T, U}
} where {T, U} function model_ring(PM::Union{PhylogeneticModel, GroupBasedPhylogeneticModel}; cached=false)
  eq_classes = equivalent_classes(PM)
  ec_indices = sort(collect(keys(eq_classes)), rev = true)

  return model_ring(base_field(PM), varname(PM) => ec_indices; cached=cached)
end

function parametrization(PM::PhylogeneticModel)
  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)

  # need to double check the order of the keys here is consitent,
  # otherwise we need to add a sort
  lvs_indices = keys(gens(R))
  map = [leaves_probability(PM, Dict(i => k[i] for i in 1:n_leaves(graph(PM)))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))
end

function parametrization(PM::GroupBasedPhylogeneticModel)
  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)
  
  lvs_indices = keys(gens(R))

  map = [leaves_fourier(PM, Dict(i => k[i] for i in 1:n_leaves(graph(PM)))) for k in lvs_indices]
  hom(R, S, reduce(vcat, map))
end

function affine_parametrization(PM::PhylogeneticModel)
  R, p = model_ring(PM) 
  S, _ = parameter_ring(PM)
  f = parametrization(PM)
  ind = keys(gens(R))

  edgs = collect(edges(graph(PM)))
  vars = unique([entry_transition_matrix(PM, i, i, e) for i in 1:n_states(PM) for e in edgs])
  vals = unique([(1 - (sum([entry_transition_matrix(PM, i, j, e) for j in 1:n_states(PM)]) - entry_transition_matrix(PM, i, i, e)))  
          for i in 1:n_states(PM) for e in edgs])

  map = [evaluate(f(p[k...]), vars, vals) for k in ind]
  hom(R, S, reduce(vcat, map))

end

function affine_parametrization(PM::GroupBasedPhylogeneticModel)
  R, q = model_ring(PM)
  S, _ = parameter_ring(PM)
  f = parametrization(PM)
  ind = keys(gens(R))

  edgs = collect(edges(graph(PM)))
  vars = [entry_fourier_parameter(PM, 1, e) for e in edgs]
  vals = repeat([1], length(edgs))

  map = [evaluate(f(q[k...]), vars, vals) for k in ind]
  hom(R, S, reduce(vcat, map))
end


###################################################################################
#
#       Fourier - probabilities coordinate change
#
###################################################################################


function fourier_transform(PM::GroupBasedPhylogeneticModel)
  FRp, p = Oscar.full_model_ring(phylogenetic_model(PM))
  FRq, q = Oscar.full_model_ring(PM)

  Rp, rp = Oscar.model_ring(phylogenetic_model(PM))
  Rq, rq = Oscar.model_ring(PM)

  p_classes = Oscar.equivalent_classes(phylogenetic_model(PM))
  q_classes = Oscar.equivalent_classes(PM)

  np = length(p_classes)
  nq = length(q_classes)
  
  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  # keys_p_classes = collect(keys(p_classes))
  # sort!(keys_p_classes, rev = true)
  # keys_q_classes = collect(keys(q_classes))
  # sort!(keys_q_classes, rev = true)
  H = Rp.(hadamard(matrix_space(ZZ, n_states(PM), n_states(PM))))
  M = Rp.(Int.(zeros(nq, np)))
  
  for (i, q_index) in enumerate(sort(collect(keys(gens(Rq))); rev=true))
    current_fourier_class = q_classes[q_index]
    for (j, p_index) in enumerate(sort(collect(keys(gens(Rp))); rev=true))
      current_prob_class = p_classes[p_index]
      current_entries_in_M = [prod([H[y,x] for (x,y) in zip(FRp[pp], FRq[qq])]) for pp in current_prob_class, qq in current_fourier_class]
      M[i,j] = Rp.(1//(length(current_prob_class)*length(current_fourier_class))*sum(current_entries_in_M))
    end
  end
  
  
  return M
end

function coordinate_change(PM::GroupBasedPhylogeneticModel)
  Rp, _ = Oscar.model_ring(phylogenetic_model(PM))
  Rq, _ = Oscar.model_ring(PM)
  
  M = fourier_transform(PM)
  hom(Rq, Rp, M * gens(_ring(Rp)))
end

function inverse_fourier_transform(PM::GroupBasedPhylogeneticModel)
  FRp, p = Oscar.full_model_ring(phylogenetic_model(PM))
  FRq, q = Oscar.full_model_ring(PM)

  Rp, _ = Oscar.model_ring(phylogenetic_model(PM))
  Rq, _ = Oscar.model_ring(PM)

  p_classes = Oscar.equivalent_classes(phylogenetic_model(PM))
  q_classes = Oscar.equivalent_classes(PM)

  np = length(p_classes)
  nq = length(q_classes)
  
  ## We need to sort the equivalence classes: both inside each class as well as the collection of classes. 
  # keys_p_classes = collect(keys(p_classes))
  # sort!(keys_p_classes, rev = true)
  # keys_q_classes = collect(keys(q_classes))
  # sort!(keys_q_classes, rev = true)
  H = Rq.(hadamard(matrix_space(ZZ, n_states(PM), n_states(PM))))
  Hinv = 1//n_states(PM) * H

  M = Rq.(Int.(zeros(np, nq)))
  for (i, p_index) in enumerate(sort(collect(keys(gens(Rp))); rev=true))
    current_prob_class = p_classes[p_index]
    for (j, q_index) in enumerate(sort(collect(keys(gens(Rq))); rev=true))
      current_fourier_class = q_classes[q_index]
      current_entries_in_M = [
        prod([Hinv[x,y] for (x,y) in zip(FRp[pp],FRq[qq])
                ]) for pp in current_prob_class, qq in current_fourier_class]
      M[i,j] = Rq.(sum(current_entries_in_M))
    end
  end
  return M
end

function inverse_coordinate_change(PM::GroupBasedPhylogeneticModel)
  Rp, _ = Oscar.model_ring(phylogenetic_model(PM))
  Rq, _ = Oscar.model_ring(PM)
  
  M = inverse_fourier_transform(PM)
  hom(Rp, Rq, M*gens(_ring(Rq)))
end


###################################################################################
#
#       Construction of specific models
#
###################################################################################


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
  
  PhylogeneticModel(G, M, root_distr)
end


# TODO:
# function affine_phylogenetic_model!(PM::PhylogeneticModel)
# function affine_phylogenetic_model!(PM::GroupBasedPhylogeneticModel)
