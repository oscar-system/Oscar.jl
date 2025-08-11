@testset "MixedGraphs" begin
  @testset "core functionality" begin
    g = mixed_graph(5)
    @test n_vertices(g) == 5
    @test n_edges(g) == 0
    @test vertices(g) == 1:5
    add_directed_edge!(g, 1, 2)
    @test n_edges(g) == 1
    @test has_directed_edge(g, 1, 2)
    rem_directed_edge!(g, 1, 2)
    @test n_edges(g) == 0
    @test !has_directed_edge(g, 1, 2)
    @test add_vertex!(g)
    @test n_vertices(g) == 6
    @test has_vertex(g, 6)
    rem_vertex!(g, 1)
    @test n_vertices(g) == 5
    @test has_vertex(g, 1)
    @test !has_vertex(g, 6)
    @test add_vertices!(g, 5) == 5
    @test n_vertices(g) == 10
    @test vertices(g) == 1:10
    @test rem_vertices!(g, [2, 4, 6, 11])
    @test n_vertices(g) == 7
    @test vertices(g) == 1:7
  end
end
