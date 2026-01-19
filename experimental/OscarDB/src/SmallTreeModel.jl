struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}
  model_type::String
end

@doc raw"""
    small_tree_model(name::String, model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}, model_type::String)

Creates a `SmallTreeModel` which is the struct used to populate the collection "AlgebraicStatistics.SmallTreeModels" of the `OscarDB`
"""
function small_tree_model(name::String,
                          model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}},
                          model_type::String)
  return SmallTreeModel(
    name,
    model,
    model_type
  )
end

@doc raw"""
    group_based_phylogenetic_model(stm::SmallTreeModel)

Returns the `GroupBasedPhylogeneticModel` of the `stm`.
"""
group_based_phylogenetic_model(stm::SmallTreeModel) = stm.model

@doc raw"""
    phylogenetic_model(stm::SmallTreeModel)

Returns the `PhylogeneticModel` of the `stm`.
"""
phylogenetic_model(stm::SmallTreeModel) = phylogenetic_model(group_based_phylogenetic_model(stm))

@doc raw"""
    model_type(stm::SmallTreeModel)

Returns the model type of the `stm` as a string.
"""
model_type(stm::SmallTreeModel) = stm.model_type

@doc raw"""
    graph(stm::SmallTreeModel)

Returns the graph of the `stm`.
"""
graph(stm::SmallTreeModel) = graph(group_based_phylogenetic_model(stm))

@doc raw"""
    n_leaves(stm::SmallTreeModel)

Returns the number of leaves of the `stm`.
"""
n_leaves(stm::SmallTreeModel) = n_leaves(graph(stm))

function Base.show(io::IO, stm::SmallTreeModel)
  print(io, "Small tree phylogenetic model $(stm._id)")
end
