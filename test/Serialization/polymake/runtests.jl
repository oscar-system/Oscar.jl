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
        @test g isa Graph{Undirected}
        @test ne(g) == 4

        um = nothing
        @test_logs (:warn, r"No function for converting the deserialized Polymake type to Oscar") um = load(joinpath(@__DIR__, "um5.mat"))
        @test um isa Polymake.PropertyValueAllocated

        i = load(joinpath(@__DIR__, "ideal.mv"))
        @test i isa Ideal
        @test map(collect, map(coefficients, gens(i))) == [[1, 1], [1, -4]]
        @test map(collect, map(exponent_vectors, gens(i))) ==
            [[[2, 0], [0, 1]], [[3, 0], [0, 1]]]
    end
end
