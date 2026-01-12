struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}
  model_type::String
end

group_based_phylogenetic_model(stm::SmallTreeModel) = stm.model
phylogenetic_model(stm::SmallTreeModel) = phylogenetic_model(group_based_phylogenetic_model(stm))
model_type(stm::SmallTreeModel) = stm.model_type
n_leaves(stm::SmallTreeModel) = stm.n_leaves
graph(stm::SmallTreeModel) = graph(group_based_phylogenetic_model(stm))

function Base.show(io::IO, stm::SmallTreeModel)
  println(io, "Small tree phylogenetic model")
  println(io, "Model type: $(model_type(stm))")
  println(io, "Number of leaves: $(n_leaves(stm))")
end
