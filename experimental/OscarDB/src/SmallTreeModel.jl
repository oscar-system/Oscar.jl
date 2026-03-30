struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}} ## Includes the rings??
  model_type::String
  dimension::Int
  degree::Int
  np::Int
  nq::Int
  dim_sl::Int
  deg_sl::Int
  MLdeg::Int
  EDdeg::Int
  parametrization_p::Oscar.MPolyAnyMap
  parametrization_q::Oscar.MPolyAnyMap
  eq_classes_p::Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}
  eq_classes_q::Dict{Tuple{Vararg{Int64}}, Vector{MPolyRingElem}}
  coodinate_change_q_p::Oscar.MPolyAnyMap
  coodinate_change_p_q::Oscar.MPolyAnyMap
  vanishing_ideal::MPolyIdeal{QQMPolyRingElem}
end

@doc raw"""
    small_tree_model(name::String, model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}, model_type::String)

Creates a `SmallTreeModel` which is the struct used to populate the collection "AlgebraicStatistics.SmallTreeModels" of the `OscarDB`
"""
function small_tree_model(name::String,
                          model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}},
                          model_type::String)

    fq = Oscar.parametrization(model)
    fp = Oscar.parametrization(phylogenetic_model(model))
    ec_q = Oscar.equivalent_classes(model)
    ec_p = Oscar.equivalent_classes(phylogenetic_model(model))
    I = Oscar.vanishing_ideal(model)

  return SmallTreeModel(
    name,
    model,
    model_type,
    dim(I),
    Oscar.degree(I),
    length(ec_p),
    length(ec_q),
    0,
    0,
    0,
    0,
    fp, fq, 
    ec_p, ec_q, 
    Oscar.coordinate_change(model),
    Oscar.inverse_coordinate_change(model),
    I
  )
end

@doc raw"""
    group_based_phylogenetic_model(stm::SmallTreeModel)

Return the `GroupBasedPhylogeneticModel` of the small tree model `stm`.
"""
group_based_phylogenetic_model(stm::SmallTreeModel) = stm.model

@doc raw"""
    phylogenetic_model(stm::SmallTreeModel)

Return the `PhylogeneticModel` of the small tree model `stm`.
"""
phylogenetic_model(stm::SmallTreeModel) = phylogenetic_model(group_based_phylogenetic_model(stm))

@doc raw"""
    model_type(stm::SmallTreeModel)

Return the model type of the `stm` as a small tree model string.
"""
model_type(stm::SmallTreeModel) = stm.model_type

@doc raw"""
    graph(stm::SmallTreeModel)

Return the graph of the small tree model `stm`.
"""
graph(stm::SmallTreeModel) = graph(group_based_phylogenetic_model(stm))

@doc raw"""
    n_leaves(stm::SmallTreeModel)

Return the number of leaves of the small tree model `stm`.
"""
n_leaves(stm::SmallTreeModel) = n_leaves(graph(stm))

@doc raw"""
    dim(stm::SmallTreeModel)

Return the dimension of the small tree model `stm`.
"""
dim(stm::SmallTreeModel) = stm.dimension

@doc raw"""
    degree(stm::SmallTreeModel)

Return the degree of the small tree model `stm`.
"""
Oscar.degree(stm::SmallTreeModel) = stm.degree

@doc raw"""
    n_coordinates(stm::SmallTreeModel)

Return the dimension of the smallest linear subspace containing the small tree model `stm` (in Fourier coordinates).
"""
n_coordinates(stm::SmallTreeModel) = stm.nq

@doc raw"""
    n_coordinates_probabilities(stm::SmallTreeModel)

Return the dimension of the smallest linear subspace containing the small tree model `stm` (in probability coordinates).
"""
n_coordinates_probabilities(stm::SmallTreeModel) = stm.np

@doc raw"""
    parametrization(stm::SmallTreeModel)

Return the parametrization of the small tree model `stm` (in Fourier coordinates).
"""
Oscar.parametrization(stm::SmallTreeModel) = stm.parametrization_q

@doc raw"""
    parametrization_probabilities(stm::SmallTreeModel)

Return the parametrization of the small tree model `stm` (in probability coordinates).
"""
parametrization_probabilities(stm::SmallTreeModel) = stm.parametrization_p

@doc raw"""
    equivalent_classes(stm::SmallTreeModel)

Return the equivalent classes of the small tree model `stm` (in Fourier coordinates).
"""
Oscar.equivalent_classes(stm::SmallTreeModel) = stm.parametrization_q

@doc raw"""
    equivalent_classes_probabilities(stm::SmallTreeModel)

Return the equivalent classes of the small tree model `stm` (in probability coordinates).
"""
equivalent_classes_probabilities(stm::SmallTreeModel) = stm.parametrization_p

@doc raw"""
    coordinate_change(stm::SmallTreeModel)

Return the linear change of coordinates, from Fourier to probability coordinates of the small tree model `stm`.
"""
Oscar.coordinate_change(stm::SmallTreeModel) = stm.coodinate_change_q_p

@doc raw"""
    inverse_coordinate_change(stm::SmallTreeModel)

Return the linear change of coordinates, from probability to Fourier coordinates of the small tree model `stm`.
"""
Oscar.inverse_coordinate_change(stm::SmallTreeModel) = stm.coodinate_change_p_q

@doc raw"""
    vanishing_ideal(stm::SmallTreeModel)

Return the vanishing ideal of the small tree model `stm`.
"""
Oscar.vanishing_ideal(stm::SmallTreeModel) = stm.vanishing_ideal

function Base.show(io::IO, stm::SmallTreeModel)
  print(io, "Small tree phylogenetic model $(stm._id)")
end
