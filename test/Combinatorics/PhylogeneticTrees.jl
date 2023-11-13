@testset "Phylogenetic Trees" begin
  @testset "constructors" begin
    ptree1 = phylogenetic_tree(Float64, "((A:4,B:4):5,C:9);")
    ptree2 = phylogenetic_tree([0 4 9; 4 0 9; 9.0 9 0], ["a", "b", "c"])
  end
end
