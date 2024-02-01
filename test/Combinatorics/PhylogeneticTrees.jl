@testset "Phylogenetic Trees" begin
  ptree1 = phylogenetic_tree(Float64, "((A:4,B:4):5,C:9);")
  
  # waiting for dante and marcel for a graph constructor
  adjacency_tree(ptree1)

  @test equidistant(ptree1) == true

  # this will need a conversion
  M = matrix(QQ, [0 4 9; 4 0 9; 9 9//2 0])
  ptree2 = phylogenetic_tree(M, ["a", "b", "c"])

  cophenetic_matrix(ptree2)
  @test newick(ptree2) == "c:4.5,(a:2,b:2):2.5;"
  @test taxa(ptree2) == ["a", "b", "c"]
  
  ptree2_test = phylogenetic_tree([0 4 9; 4 0 9.0; 9 9 0], ["a", "b", "c"])
  tmc = tropical_median_consensus(ptree2_test, ptree2)
  tmc_v = tropical_median_consensus([ptree2_test, ptree2])
  @test newick(tmc) == "c:4.5,(a:2,b:2):2.5;"
  @test newick(tmc_v) == "c:4.5,(a:2,b:2):2.5;"
end
