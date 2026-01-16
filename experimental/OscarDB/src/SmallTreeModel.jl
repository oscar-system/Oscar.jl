struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}
  model_type::String
end

group_based_phylogenetic_model(stm::SmallTreeModel) = stm.model
phylogenetic_model(stm::SmallTreeModel) = phylogenetic_model(group_based_phylogenetic_model(stm))
model_type(stm::SmallTreeModel) = stm.model_type
graph(stm::SmallTreeModel) = graph(group_based_phylogenetic_model(stm))

n_leaves(stm::SmallTreeModel) = n_leaves(graph(stm))

function Base.show(io::IO, stm::SmallTreeModel)
  print(io, "Small tree phylogenetic model $(stm._id)")
end
