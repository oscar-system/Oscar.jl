## OLD STRUCT
struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}
  model_type::String
end

## NEW STRUCTS
struct SmallGroupBasedModel
  _id::String # model encoding id, example GB-3-0-0-JC
  model::GroupBasedPhylogeneticModel ## Before: Union{GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}},
                                     ## GroupBasedPhylogeneticModel{PhylogeneticNetwork}} --> not sure if necessary
  model_type::String # ex: jukes_cantor_model
  n_leaves::Int
  level::Int
  n_cycles::Int
  dimension::Union{Int, Nothing}
  degree::Union{Int, Nothing}
  n_coordinates::Int
  dim_sl::Union{Int, Nothing}
  deg_sl::Union{Int, Nothing}
  EDdeg::Union{Int, Nothing}
  parametrization::Oscar.MPolyAnyMap
  eq_classes::Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}
  vanishing_ideal::Union{MPolyIdeal{QQMPolyRingElem}, Nothing}
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

@doc raw"""
  graph(sgb::SmallGroupBasedModel)

Return the graph of the small group-based model `sgb`.
"""
graph(sgb::SmallGroupBasedModel) = graph(group_based_phylogenetic_model(sgb))

@doc raw"""
  n_leaves(sgb::SmallGroupBasedModel)

Return the number of leaves of the phylogenetic tree or network in the small group-based model `sgb`.
"""
n_leaves(sgb::SmallGroupBasedModel) = sgb.n_leaves

@doc raw"""
  level(sgb::SmallGroupBasedModel)

Return the level of the phylogenetic network `N` in the small group-based model `sgb`. If `N` is a tree, the ouptput is zero.
"""
level(sgb::SmallGroupBasedModel) = sgb.level

@doc raw"""
  n_cycles(sgb::SmallGroupBasedModel)

Return the number of cycles in the phylogenetic network `N` in the small group-based model `sgb`. If `N` is a tree, the ouptput is zero.
"""
n_cycles(sgb::SmallGroupBasedModel) = sgb.n_cycles

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
  model::PhylogeneticModel
  model_type::String # ex: jukes_cantor_model
  n_leaves::Int
  level::Int
  n_cycles::Int
  dimension::Union{Int, Nothing}
  degree::Union{Int, Nothing}
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
  phylogenetic_model(spm::SmallPhylogeneticModel)

Return the `PhylogeneticModel` of the small phylogenetic model `spm`.
"""
phylogenetic_model(spm::SmallPhylogeneticModel) = spm.model

@doc raw"""
  model_type(spm::SmallPhylogeneticModel)

Return the model type of the `spm` as a small phylogenetic model string.
"""
model_type(spm::SmallPhylogeneticModel) = spm.model_type

@doc raw"""
  graph(spm::SmallPhylogeneticModel)

Return the graph of the small phylogenetic model `spm`.
"""
graph(spm::SmallPhylogeneticModel) = graph(phylogenetic_model(spm))

@doc raw"""
  n_leaves(spm::SmallPhylogeneticModel)

Return the number of leaves of the phylogenetic tree or network in the small phylogenetic model `spm`.
"""
n_leaves(spm::SmallPhylogeneticModel) = spm.n_leaves

@doc raw"""
  level(spm::SmallPhylogeneticModel)

Return the level of the phylogenetic network `N` in the small phylogenetic model `spm`. If `N` is a tree, the ouptput is zero.
"""
level(spm::SmallPhylogeneticModel) = spm.level

@doc raw"""
  n_cycles(spm::SmallPhylogeneticModel)

Return the number of cycles in the phylogenetic network `N` in the small phylogenetic model `spm`. If `N` is a tree, the ouptput is zero.
"""
n_cycles(spm::SmallPhylogeneticModel) = spm.n_cycles

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


