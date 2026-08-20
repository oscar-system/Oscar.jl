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
    group_based_phylogenetic_model(sgb::SmallGroupBasedModel)

Return the `GroupBasedPhylogeneticModel` of the small group-based model `sgb`.
"""
group_based_phylogenetic_model(sgb::SmallGroupBasedModel) = sgb.model

@doc raw"""
    phylogenetic_model(sgb::SmallGroupBasedModel)
    phylogenetic_model(spm::SmallPhylogeneticModel)

Return the `PhylogeneticModel` of the small group-based model `sgb` or of the small phylogenetic model `spm`.
"""
phylogenetic_model(sgb::SmallGroupBasedModel) = phylogenetic_model(group_based_phylogenetic_model(sgb))
phylogenetic_model(spm::SmallPhylogeneticModel) = spm.model

@doc raw"""
    model_type(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the model type of the small model `sm` as a string.
"""
model_type(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.model_type

@doc raw"""
    graph(sgb::SmallGroupBasedModel)
    graph(spm::SmallPhylogeneticModel)

Return the graph of the small model `sm`.
"""
graph(sgb::SmallGroupBasedModel) = graph(group_based_phylogenetic_model(sgb))
graph(spm::SmallPhylogeneticModel) = graph(phylogenetic_model(spm))

@doc raw"""
    n_leaves(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the number of leaves of the phylogenetic tree or network in the small model `sm`.
"""
n_leaves(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.n_leaves

@doc raw"""
    level(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the level of the phylogenetic network in the small model `sm`. If the graph is a tree, the output is zero.
"""
level(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.level

@doc raw"""
    n_cycles(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the number of cycles in the phylogenetic network in the small model `sm`. If the graph is a tree, the output is zero.
"""
n_cycles(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.n_cycles

@doc raw"""
    dim(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the dimension of the small model `sm`.
"""
dim(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.dimension

@doc raw"""
    degree(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the degree of the small model `sm`.
"""
degree(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.degree

@doc raw"""
    n_coordinates(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the number of coordinates of the small model `sm`. 
For group-based models, this is the number of Fourier coordinates.
"""
n_coordinates(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.n_coordinates

@doc raw"""
    dimension_singular_locus(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the dimension of the singular locus of the small model `sm`.
"""
dimension_singular_locus(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.dim_sl

@doc raw"""
    degree_singular_locus(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the degree of the singular locus of the small model `sm`.
"""
degree_singular_locus(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.deg_sl

@doc raw"""
    euclidean_distance_degree(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the Euclidean distance degree of the small model `sm`.
"""
euclidean_distance_degree(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.EDdeg

@doc raw"""
    maximum_likelihood_degree(spm::SmallPhylogeneticModel)

Return the maximum likelihood degree of the small phylogenetic model `spm`.
"""
maximum_likelihood_degree(spm::SmallPhylogeneticModel) = spm.MLdeg

@doc raw"""
    parametrization(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the parametrization of the small model `sm`.
For group-based models, this is the parametrization in Fourier coordinates.
"""
parametrization(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.parametrization

@doc raw"""
    equivalent_classes(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the equivalent classes of the small model `sm`.
For group-based models, these are the classes in Fourier coordinates.
"""
equivalent_classes(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.eq_classes

@doc raw"""
    vanishing_ideal(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel})

Return the vanishing ideal of the small model `sm`.
For group-based models, this is the ideal in Fourier coordinates.
"""
vanishing_ideal(sm::Union{SmallGroupBasedModel, SmallPhylogeneticModel}) = sm.vanishing_ideal

function Base.show(io::IO, sgb::SmallGroupBasedModel)
  print(io, "Small group-based phylogenetic model $(sgb._id)")
end

function Base.show(io::IO, spm::SmallPhylogeneticModel)
  print(io, "Small phylogenetic model $(spm._id)")
end
