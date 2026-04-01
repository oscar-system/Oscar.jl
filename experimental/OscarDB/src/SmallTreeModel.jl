## OLD STRUCT
struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}
  model_type::String
end

## NEW STRUCTS
struct SmallGroupBasedModel
  _id::String # model encoding id, example GB-3-0-0-JC
  model::Union{GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}},
               GroupBasedPhylogeneticModel{PhylogeneticNetwork{QQFieldElem}}} ## Can this be just ::GroupBasedPhylogeneticModel?
  model_type::String # ex: jukes_cantor_model
  phylogenetic_model_id::String
  n_leaves::Int
  level::Int
  n_cycles::Int
  dimension::Int
  degree::Int
  n_coordinates::Int
  dim_sl::Union{Int, Nothing}
  deg_sl::Union{Int, Nothing}
  EDdeg::Union{Int, Nothing}
  parametrization::Oscar.MPolyAnyMap
  eq_classes::Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}
  vanishing_ideal::Union{MPolyIdeal{QQMPolyRingElem}, Nothing}
end


@doc raw"""
    small_group_based_model(name::String, model::Union{GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}},
                            GroupBasedPhylogeneticModel{PhylogeneticNetwork{QQFieldElem}}}, model_type::String, phylogenetic_model_id::String)

Creates a `SmallGroupBasedModel` which is the struct used to populate the collection "AlgebraicStatistics.SmallGroupBasedModels" of the `OscarDB`
"""
function small_group_based_model(name::String,
                          model::Union{GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}},
                                       GroupBasedPhylogeneticModel{PhylogeneticNetwork{QQFieldElem}}},
                          model_type::String,
                          phylogenetic_model_id::String)

  G = graph(model)
  f = parametrization(model)
  ec = equivalent_classes(model)
  I = vanishing_ideal(model)

  return SmallGroupBasedModel(
    name,
    model,
    model_type,
    phylogenetic_model_id,
    Oscar.n_leaves(G),
    Oscar.level_phylogenetic_network(graph_from_edges(Directed, edges(G))),
    sum([length(h)-1 for h in Oscar.hybrid_edges(graph_from_edges(Directed, edges(G)))]),
    dim(I),
    degree(I),
    length(ec),
    nothing,
    nothing,
    nothing,
    f,
    ec,
    I
  )
end

@doc raw"""
  group_based_phylogenetic_model(sgb::SmallGroupBasedModel)

Return the `GroupBasedPhylogeneticModel` of the small group-based model `sgb`.
"""
group_based_phylogenetic_model(sgb::SmallGroupBasedModel) = sgb.model

@doc raw"""
    phylogenetic_model(sgb::SmallGroupBasedModel)

Return the `PhylogeneticModel` of the small group-based model `sgb`.
"""
phylogenetic_model(sgb::SmallGroupBasedModel) = phylogenetic_model(group_based_phylogenetic_model(sgb))

@doc raw"""
    model_type(sgb::SmallGroupBasedModel)

Return the model type of the `sgb` as a small group-based model string.
"""
model_type(sgb::SmallGroupBasedModel) = sgb.model_type

#### ????
@doc raw"""
    phylogenetic_model_id(sgb::SmallGroupBasedModel)

Return the id of the corresponding `SmallPhylogeneticModel` from the collection "AlgebraicStatistics.SmallPhylogeneticModels" of the `OscarDB`
"""
phylogenetic_model_id(sgb::SmallGroupBasedModel) = sgb.phylogenetic_model_id

@doc raw"""
    graph(sgb::SmallGroupBasedModel)

Return the graph of the small group-based model `sgb`.
"""
graph(sgb::SmallGroupBasedModel) = graph(group_based_phylogenetic_model(sgb))

@doc raw"""
    n_leaves(sgb::SmallGroupBasedModel)

Return the number of leaves of the small group-based model `sgb`.
"""
n_leaves(sgb::SmallGroupBasedModel) = n_leaves(graph(sgb))

@doc raw"""
    dim(sgb::SmallGroupBasedModel)

Return the dimension of the small group-based model `sgb`.
"""
dim(sgb::SmallGroupBasedModel) = sgb.dimension

@doc raw"""
    degree(sgb::SmallGroupBasedModel)

Return the degree of the small group-based model `sgb`.
"""
degree(sgb::SmallGroupBasedModel) = sgb.degree

@doc raw"""
    n_coordinates(sgb::SmallGroupBasedModel)

Return the dimension of the smallest linear subspace containing the small group-based model `sgb` (in Fourier coordinates).
"""
n_coordinates(sgb::SmallGroupBasedModel) = sgb.n_coordinates

@doc raw"""
    dimension_singular_locus(sgb::SmallGroupBasedModel)

Return the dimension of the singular locus of the small group-based model `sgb`.
"""
dimension_singular_locus(sgb::SmallGroupBasedModel) = sgb.dim_sl

@doc raw"""
    degree_singular_locus(sgb::SmallGroupBasedModel)

Return the degree of the singular locus of the small group-based model `sgb`.
"""
degree_singular_locus(sgb::SmallGroupBasedModel) = sgb.deg_sl

@doc raw"""
    euclidean_distance_degree(sgb::SmallGroupBasedModel)

Return the Euclidean distance degree the small group-based model `sgb`.
"""
euclidean_distance_degree(sgb::SmallGroupBasedModel) = sgb.EDdeg

@doc raw"""
    parametrization(sgb::SmallGroupBasedModel)

Return the parametrization of the small group-based model `sgb` (in Fourier coordinates).
"""
parametrization(sgb::SmallGroupBasedModel) = sgb.parametrization

@doc raw"""
    equivalent_classes(sgb::SmallGroupBasedModel)

Return the equivalent classes of the small group-based model `sgb` (in Fourier coordinates).
"""
equivalent_classes(sgb::SmallGroupBasedModel) = sgb.eq_classes

@doc raw"""
    vanishing_ideal(sgb::SmallGroupBasedModel)

Return the vanishing ideal of the small group-based model `sgb` (in Fourier coordinates).
"""
vanishing_ideal(sgb::SmallGroupBasedModel) = sgb.vanishing_ideal

function Base.show(io::IO, sgb::SmallGroupBasedModel)
  print(io, "Small group-based phylogenetic model $(sgb._id)")
end




struct SmallPhylogeneticModel
  _id::String # model encoding id, example 3-0-0-JC
  model::Union{PhylogeneticModel{PhylogeneticTree{QQFieldElem}},
               PhylogeneticModel{PhylogeneticNetwork{QQFieldElem}}}
  model_type::String # ex: jukes_cantor_model
  extended_model_id::String ## Link to id of GBmodel (and ATRmodel in the future)
  n_leaves::Int
  level::Int
  n_cycles::Int
  dimension::Int
  degree::Int
  n_coordinates::Int
  dim_sl::Union{Int, Nothing}
  deg_sl::Union{Int, Nothing}
  MLdeg::Union{Int, Nothing}
  EDdeg::Union{Int, Nothing}
  parametrization::Oscar.MPolyAnyMap
  eq_classes::Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}
  vanishing_ideal::Union{MPolyIdeal{QQMPolyRingElem}, Nothing}
end
@doc raw"""
    small_phylogenetic_model(name::String, model::Union{PhylogeneticModel{PhylogeneticTree{QQFieldElem}},
                             PhylogeneticModel{PhylogeneticNetwork{QQFieldElem}}}, model_type::String, extended_model_id::Union{String, Nothing})

Creates a `SmallPhylogeneticModel` which is the struct used to populate the collection "AlgebraicStatistics.SmallPhylogeneticModels" of the `OscarDB`
"""
function small_phylogenetic_model(name::String,
                          model::Union{PhylogeneticModel{PhylogeneticTree{QQFieldElem}},
                                       PhylogeneticModel{PhylogeneticNetwork{QQFieldElem}}},
                          model_type::String,
                          extended_model_id::String="")

  G = graph(model)
  f = parametrization(model)
  ec = equivalent_classes(model)
  # I = vanishing_ideal(model)

  return SmallPhylogeneticModel(
    name,
    model,
    model_type,
    extended_model_id,
    Oscar.n_leaves(G),
    Oscar.level_phylogenetic_network(graph_from_edges(Directed, edges(G))),
    sum([length(h)-1 for h in Oscar.hybrid_edges(graph_from_edges(Directed, edges(G)))]),
    0,
    0,
    length(ec),
    nothing,
    nothing,
    nothing,
    nothing,
    f,
    ec,
    nothing
  )
end


@doc raw"""
    phylogenetic_model(spm::SmallPhylogeneticModel)

Return the `PhylogeneticModel` of the small phylogenetic model `spm`.
"""
phylogenetic_model(spm::SmallPhylogeneticModel) = spm.model

@doc raw"""
    model_type(spm::SmallPhylogeneticModel)

Return the model type of the `spm` as a small phylogenetic model string.
"""
model_type(spm::SmallPhylogeneticModel) = spm.model_type

#### ????
@doc raw"""
    extended_model_id(sgb::SmallPhylogeneticModel)

Return the id of the corresponding `GroupBasedModel` from the collection "AlgebraicStatistics.GroupBasedModels" of the `OscarDB`
"""
extended_model_id(sgb::SmallPhylogeneticModel) = sgb.extended_model_id

@doc raw"""
    graph(spm::SmallPhylogeneticModel)

Return the graph of the small phylogenetic model `spm`.
"""
graph(spm::SmallPhylogeneticModel) = graph(phylogenetic_model(spm))

@doc raw"""
    n_leaves(spm::SmallPhylogeneticModel)

Return the number of leaves of the small phylogenetic model `spm`.
"""
n_leaves(spm::SmallPhylogeneticModel) = n_leaves(graph(spm))

@doc raw"""
    dim(spm::SmallPhylogeneticModel)

Return the dimension of the small phylogenetic model `spm`.
"""
dim(spm::SmallPhylogeneticModel) = spm.dimension

@doc raw"""
    degree(spm::SmallPhylogeneticModel)

Return the degree of the small phylogenetic model `spm`.
"""
degree(spm::SmallPhylogeneticModel) = spm.degree

@doc raw"""
    n_coordinates(spm::SmallPhylogeneticModel)

Return the dimension of the smallest linear subspace containing the small phylogenetic model `spm`.
"""
n_coordinates(spm::SmallPhylogeneticModel) = spm.n_coordinates

@doc raw"""
    dimension_singular_locus(spm::SmallPhylogeneticModel)

Return the dimension of the singular locus of the small phylogenetic model `spm`.
"""
dimension_singular_locus(spm::SmallPhylogeneticModel) = spm.dim_sl

@doc raw"""
    degree_singular_locus(spm::SmallPhylogeneticModel)

Return the degree of the singular locus of the small phylogenetic model `spm`.
"""
degree_singular_locus(spm::SmallPhylogeneticModel) = spm.deg_sl

@doc raw"""
    euclidean_distance_degree(spm::SmallPhylogeneticModel)

Return the Euclidean distance degree the small phylogenetic model `spm`.
"""
euclidean_distance_degree(spm::SmallPhylogeneticModel) = spm.EDdeg

@doc raw"""
    maximum_likelihood_degree(spm::SmallPhylogeneticModel)

Return the maximum likelihood degree the small phylogenetic model `spm`.
"""
maximum_likelihood_degree(spm::SmallPhylogeneticModel) = spm.MLdeg

@doc raw"""
    parametrization(spm::SmallPhylogeneticModel)

Return the parametrization of the small phylogenetic model `spm`.
"""
parametrization(spm::SmallPhylogeneticModel) = spm.parametrization

@doc raw"""
    equivalent_classes(spm::SmallPhylogeneticModel)

Return the equivalent classes of the small phylogenetic model `spm`.
"""
equivalent_classes(spm::SmallPhylogeneticModel) = spm.eq_classes

@doc raw"""
    vanishing_ideal(spm::SmallPhylogeneticModel)

Return the vanishing ideal of the small phylogenetic model `spm`.
"""
vanishing_ideal(spm::SmallPhylogeneticModel) = spm.vanishing_ideal

function Base.show(io::IO, spm::SmallPhylogeneticModel)
  print(io, "Small phylogenetic model $(spm._id)")
end

