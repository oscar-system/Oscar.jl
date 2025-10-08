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
  pt1 = phylogenetic_tree(QQFieldElem, tree)
  @test newick(pt1) == "(leaf 2:1,leaf 3:1,leaf 1:1):1,leaf 4:1;"
  @test collect(edges(tree)) == collect(edges(adjacency_tree(pt1)))

  label!(tree, nothing, Dict(1 => "a", 2 => "b", 3 => "c", 6 => "d"); name=:leaves)
  label!(tree, Dict((5, 6) => 1,
                    (5, 4) => 2,
                    (4, 1) => 3,
                    (4, 2) => 4,
                    (4, 3) => 5), nothing; name=:distance)

  pt2 = phylogenetic_tree(QQFieldElem, tree)
  
  @test newick(pt2) == "(b:4,c:5,a:3):2,d:1;"
  at = adjacency_tree(pt2)

  for v in 1:n_vertices(tree)
    @test at.leaves[v] == tree.leaves[v]
  end

  for e in edges(tree)
    @test at.distance[e] == tree.distance[e]
  end

  tree = graph_from_edges(Directed,[[4,1],[4,2],[4,3]])
  pt = phylogenetic_tree(QQFieldElem, tree)
  a_tree = adjacency_tree(pt)
  @test collect(edges(a_tree)) == collect(edges(tree))
  @test_throws ArgumentError a_tree.distance[1, 2]

  for l in leaves(pt)
    @test a_tree.leaves[l] == "leaf $l"
  end

  a_tree = adjacency_tree(pt; add_labels=false)
  @test !hasproperty(a_tree, :distance)
  @test !hasproperty(a_tree, :leaves)

  @test Oscar.interior_nodes(pt) == [4]
  @test collect(edges(pt)) == [Edge(4, 1), Edge(4, 2), Edge(4, 3)]
  @test n_edges(pt) == 3
  println(newick(pt))
  @test Oscar.is_isomorphic(pt, phylogenetic_tree(QQFieldElem, "a:1,b:2,c:3;"))
end
