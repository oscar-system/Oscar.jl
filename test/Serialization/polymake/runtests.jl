@testset "polymake" begin
    
    @testset "loading polymake objects" begin
        s = load(joinpath(@__DIR__, "square.poly"))
        @test s isa Polyhedron{fmpq}
        @test nvertices(s) == 4
        @test nfacets(s) == 4

        nf = load(joinpath(@__DIR__, "nf_square.fan"))
        @test nf isa PolyhedralFan{fmpq}
        @test nrays(nf) == 4
        @test n_maximal_cones(nf) == 4

        g = load(joinpath(@__DIR__, "square_graph.graph"))
        @test g isa Graphs.Graph{Graphs.Undirected}
    end
end
