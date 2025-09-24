@doc raw"""
    load_phylogenetic_model(tree_model_id::String)

Loads a phylogenetic model using the identifier given on the website. 
The identifier has the form tree_id-model_id, e.g. "3-0-0-JC" for the star3 graph
with Jukes-Cantor model.
"""
function load_phylogenetic_model(tree_model_id::String)
  file_path = joinpath(oscardir, "data", "AlgebraicStatistics",  "PhylogeneticModels", tree_model_id * ".mrdi")
  @req ispath(file_path) "Could not find model $tree_model_id, either the model has not been stored or the identifier was not entered correctly e.g. 3-0-0-JC"
  phylo_model = load(file_path)
  return phylo_model
end


