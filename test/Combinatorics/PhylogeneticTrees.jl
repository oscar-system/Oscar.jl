@testset "Phylogenetic Trees" begin
  ptree1 = phylogenetic_tree(Float64, "((A:4,B:4):5,C:9);")
  test_tree = graph_from_edges(Directed, [[1, 2], [2, 3], [2, 4], [1, 5]])
  @test is_isomorphic(adjacency_tree(ptree1), test_tree)
  @test equidistant(ptree1) == true

  ptree1_test = phylogenetic_tree([0 4 9; 4 0 9.0; 9 9 0], ["A", "B", "C"])
  tmc = tropical_median_consensus(ptree1_test, ptree1)
  tmc_v = tropical_median_consensus([ptree1_test, ptree1])
  @test newick(tmc) == "C:9,(A:6.5,B:6.5):2.5;"
  @test newick(tmc_v) == "C:9,(A:6.5,B:6.5):2.5;"

  M = matrix(QQ, [0 4 9; 4 0 9; 9 9 0])
  ptree2 = phylogenetic_tree(M, ["a", "b", "c"])

  @test cophenetic_matrix(ptree2) == matrix(QQ, [0   4   9; 4   0   9; 9   9   0])
  @test newick(ptree2) == "c:9/2,(a:2,b:2):5/2;"
  @test taxa(ptree2) == ["a", "b", "c"]

  tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3], [5, 4], [5, 6]])
  pt = phylogenetic_tree(QQFieldElem, tree)
end
