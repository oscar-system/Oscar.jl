######################################
#### PHYLOGENETIC DATA STRUCTURES ####
######################################

@attributes mutable struct PhylogeneticModel{L} <: GraphicalModel{Directed, L}
  # do need to for T to be directed here?
  base_field::Field
  graph::Graph{Directed}
  labelings::L
  trans_mat_structure::Matrix{<: VarName}
  model_parameter_name::VarName
  n_states::Int
  root_distribution::Vector
  
  function PhylogeneticModel(F::Field,
                             G::Graph{Directed},
                             trans_mat_structure::Matrix{<: VarName},
                             n_states::Union{Nothing, Int} = nothing,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p")
    if isnothing(n_states)
      n_states = 4
    end
    if isnothing(root_distribution)
      root_distribution = repeat([1//n_states], outer = n_states)
    end
    graph_maps = NamedTuple(_graph_maps(G))
    graph_maps = isempty(graph_maps) ? nothing : graph_maps
    return new{typeof(graph_maps)}(F, G,
                                   graph_maps,
                                   trans_mat_structure,
                                   varname,
                                   n_states,
                                   root_distribution)
  end

  function PhylogeneticModel(G::Graph{Directed},
                             trans_mat_structure::Matrix{<: VarName},
                             n_states::Union{Nothing, Int} = nothing,
                             root_distribution::Union{Nothing, Vector} = nothing,
                             varname::VarName="p")
    return PhylogeneticModel(QQ, G, trans_mat_structure, n_states, root_distribution, varname)
  end
end

n_states(PM::PhylogeneticModel) = PM.n_states
transition_matrix(PM::PhylogeneticModel) = PM.trans_mat_structure
base_field(PM::PhylogeneticModel) = PM.base_field
varname(PM::PhylogeneticModel) = PM.model_parameter_name
root_distribution(PM::PhylogeneticModel) = PM.root_distribution

@attr Tuple{
  MPolyRing,
  GraphTransDict
} function parameter_ring(PM::PhylogeneticModel; cached=false)
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

function entry_transition_matrix(PM::PhylogeneticModel, e::Edge, i::Int, j::Int)
  tr_mat = transition_matrix(PM)
  parameter_ring(PM)[2][tr_mat[i,j], e]
end

@attr MPolyAnyMap function parametrization(PM::PhylogeneticModel)
  gr = graph(PM)

  R, _ = model_ring(PM)
  S, _ = parameter_ring(PM)

  lvs_indices = leaves_indices(PM::PhylogeneticModel)

  probabilities = [leaves_probability(PM, Dict(i => k[i] for i in 1:n_leaves(gr))) for k in lvs_indices]
  hom(R, S, reduce(vcat, probabilities))

end

@attr Dict{QQMPolyRingElem, Vector{QQMPolyRingElem}} function equivalent_classes(PM::PhylogeneticModel)
  _, ps = model_ring(PM)
  param = parametrization(PM);
  
  polys = unique(param.img_gens)
  polys = polys[findall(!is_zero, polys)]

  #equivalent_classes = Dict{Tuple{Vararg{Int64}}, Vector{QQMPolyRingElem}}()
  equivalent_classes = Dict{QQMPolyRingElem, Vector{QQMPolyRingElem}}()
  for poly in polys
      eqv_class = sort([p for p in ps if param(p) == poly], rev = true)
      # equivalent_classes[Tuple(index(eqv_class[1]))] = eqv_class
      equivalent_classes[eqv_class[1]] = eqv_class
  end

  return equivalent_classes
end

@attr Tuple{MPolyRing, Array} function reduced_model_ring(PM::PhylogeneticModel; cached=false)
    eq_calasses = Oscar.equivalent_classes(PM)
    keys_eq_classes = collect(keys(eq_calasses))

    return polynomial_ring(base_field(PM),
                         [Symbol(k) for k in keys_eq_classes];
                         cached=cached)
end

@attr MPolyAnyMap function reduced_parametrization(PM::PhylogeneticModel)

    _, p = model_ring(PM)
    S, _ = parameter_ring(PM)
    param = parametrization(PM)

    redR, _ = reduced_model_ring(PM)
    gens_redR = gens(redR)
    keys_eq_classes = index.(gens_redR)
    
    probabilities = [param(p[k...]) for k in keys_eq_classes]
    hom(redR, S, reduce(vcat, probabilities))
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



struct GroupBasedPhylogeneticModel
  phylo_model::PhylogeneticModel
  fourier_ring::MPolyRing{QQFieldElem}
  fourier_params::Dict{Edge, Vector{QQMPolyRingElem}}
  group::Vector{FinGenAbGroupElem}
end

# prob_ring::MPolyRing{QQFieldElem}
# root_distr::Vector{Any} #this need to become more precise
# trans_matrices::Dict{Edge, MatElem{QQMPolyRingElem}}

function Base.show(io::IO, pm::GroupBasedPhylogeneticModel)
  gr = graph(pm)
  nl = length(leaves(gr))
  ne = length(collect(edges(gr)))
  root_dist = join(Oscar.root_distribution(pm), ", ")
  c_edg = 2
  p_edg = inneighbors(gr, c_edg)[1]
  findall(x-> x==2, dst.(edges(gr)))
  M = transition_matrix(pm)[Edge(p_edg, c_edg)]
  idx = string(split(string(M[1,1]), "[")[2][1])

  print(io, "Group-based phylogenetic model on a tree with $(nl) leaves and $(ne) edges \n with distribution at the root [$(root_dist)]. \n")
  print(io, " The transition matrix associated to edge i is of the form \n ")
  print(io, replace(replace(string(M), "["*idx => "[i"), ";" => ";\n "))
  print(io, ", \n and the Fourier parameters are ")
  fp = transpose(fourier_parameters(pm)[Edge(p_edg, c_edg)])
  fp = replace(string(fp), "QQMPolyRingElem" => "")
  print(io, replace(replace(replace(string(fp), "["*idx => "[i"), ";" => ";\n "), "]]" => "]]."))
end


############################
#### GROUP-BASED MODELS ####
############################


function jukes_cantor_model(G::Graph{Directed})
  PhylogeneticModel(G, [:a :b :b :b;
                        :b :a :b :b;
                        :b :b :a :b;
                        :b :b :b :a]
                    )
  #G = collect(abelian_group(2,2))
  #group = [G[1],G[3],G[2],G[4]]
  #
  #pm = PhylogeneticModel(graph, ns, R, root_distr, matrices)
  #return GroupBasedPhylogeneticModel(pm, S, fourier_param, group)
end
