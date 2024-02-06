@testset "Phylogenetic Trees" begin
  ptree1 = phylogenetic_tree(Float64, "((A:4,B:4):5,C:9);")
  
  # waiting for dante and marcel for a graph constructor
  adjacency_tree(ptree1)

  @test equidistant(ptree1) == true

  M = matrix(QQ, [0 4 9; 4 0 9; 9 9 0])
  ptree2 = phylogenetic_tree(M, ["a", "b", "c"])

  cophenetic_matrix(ptree2)
  @test newick(ptree2) == "c:9/2,(a:2,b:2):5/2;"
  @test taxa(ptree2) == ["a", "b", "c"]
  
  ptree1_test = phylogenetic_tree([0 4 9; 4 0 9.0; 9 9 0], ["A", "B", "C"])
  # tmc = tropical_median_consensus(ptree1_test, ptree1) -- Varargs method removed
  tmc_v = tropical_median_consensus([ptree1_test, ptree1])
  # @test newick(tmc) == "C:9,(A:6.5,B:6.5):2.5;"
  @test newick(tmc_v) == "C:9,(A:6.5,B:6.5):2.5;"
end
