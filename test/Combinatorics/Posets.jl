@testset "Posets" begin

  covrels = incidence_matrix(4, 4, [[2, 3], [4], [4], Int[]])
  a2_adj_cov = [
           0 1 1 0 0 0
           0 0 0 1 1 0
           0 0 0 1 1 0
           0 0 0 0 0 1
           0 0 0 0 0 1
           0 0 0 0 0 0
         ];
  sm = sparse_matrix(ZZ.(a2_adj_cov))
  bm = (!iszero).(a2_adj_cov)

  g = graph_from_edges(Directed, [[2,1],[2,4],[1,3],[4,3],[3,6],[2,5],[5,6],[7,1],[7,5],[5,8],[7,8]])
  ranks = Dict(2=>0, 1=>1, 4=>2, 3=>3, 5=>2, 6=>4, 7=>0, 8=>4)

  gcirc = graph_from_edges(Directed, [[1,2],[2,3],[3,1]])
  @test_throws ArgumentError partially_ordered_set(gcirc)

  vif = vertex_indices(facets(cube(3)))

  pos1 = partially_ordered_set(covrels)
  @test length(pos1) == 4
  pos2 = partially_ordered_set(a2_adj_cov)
  @test length(pos2) == 6
  pos2s = partially_ordered_set(sm)
  @test length(pos2s) == 6
  pos2b = partially_ordered_set(bm)
  @test length(pos2b) == 6

  posg = partially_ordered_set(g)
  @test length(posg) == 8

  posgr = partially_ordered_set(g, ranks)
  @test length(posgr) == 8

  poscubeinc = partially_ordered_set_from_inclusions(vif)
  @test length(poscubeinc) == 28

  poscube = face_poset(cube(3))
  @test length(poscube) == 28

  possq = face_poset(cube(2))
  @test length(possq) == 10

  posnfsq = face_poset(normal_fan(cube(2)))
  @test length(posnfsq) == 9

  pfm = lattice_of_flats(fano_matroid())
  @test length(pfm) == 16
  
  pcfm = lattice_of_cyclic_flats(fano_matroid())
  @test length(pcfm) == 9

  posmr = maximal_ranked_poset([2,4,3])
  @test length(posmr) == 11

  @testset "basics" begin
    @test rank(pos1) == 2
    @test rank(pos2) == 3
    @test rank(posg) == 3
    @test rank(posgr) == 4
    @test rank(poscubeinc) == 4
    @test rank(poscube) == 4
    @test rank(possq) == 3
    @test rank(posnfsq) == 2
    @test rank(pfm) == 3
    @test rank(pcfm) == 3
    @test rank(posmr) == 4

    @test n_atoms(posmr) == 2
    @test n_coatoms(posmr) == 3

    @test least_element(posmr) < first(atoms(posmr))
    @test least_element(posmr) < greatest_element(posmr)

    @test_throws ArgumentError least_element(posg)
    @test_throws ArgumentError greatest_element(posgr)

    @test issetequal(atoms(pos1),coatoms(pos1))

    @test rank(greatest_element(poscube)) == rank(poscube)

    @test rank(least_element(poscube)) == 0

    @test parent(least_element(pos1)) == pos1

    @test length(elements_of_rank(posmr, 2)) == 4
  end

  @testset "derived objects" begin
    @test comparability_graph(pos1) isa Graph{Undirected}
    cg = comparability_graph(possq)
    @test n_vertices(cg) == length(possq)
    @test n_edges(cg) == 25

    # the graph always includes the artificial top node
    @test is_isomorphic(graph(possq), graph(posnfsq))

    op = order_polytope(pos1)
    @test dim(op) == 2
    @test n_vertices(op) == 4

    cp = chain_polytope(pos1)
    @test dim(cp) == 2
    @test n_vertices(cp) == 4

    nr = Oscar.node_ranks(posmr)
    g = graph(posmr)
    posgnr = partially_ordered_set(g, nr)
    @test is_isomorphic(graph(posmr), graph(posgnr))
    @test Oscar.node_ranks(posmr) == Oscar.node_ranks(posgnr)
  end

  @testset "chains" begin
    mc = maximal_chains(pos1)
    @test length(mc) == 2
    @test maximal_chains(Int, pos1) isa Vector{Vector{Int}}
    @test maximal_chains(pos1) isa Vector{Vector{Oscar.PartiallyOrderedSetElement}}

    @test length(first(maximal_chains(possq))) == 4
    @test length(first(maximal_chains(posnfsq))) == 3

    mc = maximal_chains(possq)
    @test length(Set(reduce(vcat,mc))) == length(possq)
  end
end
