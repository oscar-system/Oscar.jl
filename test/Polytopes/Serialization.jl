@testset "Serialization" begin
    
    mktempdir() do path
        # TODO: Enable the following testset, see #758
        # @testset "Graph" begin
        #     G = Graphs.complete_graph(4)
        #     Graphs.save_graph(G, joinpath(path, "c4.graph"))
        #     loaded = Graphs.load_graph(joinpath(path, "c4.graph"))
        #     @test loaded isa Graphs.Graph{Graphs.Undirected}
        #     @test Base.propertynames(G) == Base.propertynames(loaded)
        #     @test Graphs.nv(G) == Graphs.nv(loaded)
        #     @test Graphs.ne(G) == Graphs.ne(loaded)
        # end

        @testset "Cone" begin
            C = positive_hull([1 0; 0 1])
            save_cone(C, joinpath(path, "posorth.cone"))
            loaded = load_cone(joinpath(path, "posorth.cone"))
            @test loaded isa Cone
            @test Base.propertynames(C) == Base.propertynames(loaded)
            @test nrays(C) == nrays(loaded)
            @test dim(C) == dim(loaded)
            @test C == loaded
        end

        @testset "Polyhedron" begin
            square = cube(2)
            save_polyhedron(square, joinpath(path, "square.poly"))
            loaded = load_polyhedron(joinpath(path, "square.poly"))
            @test loaded isa Polyhedron
            @test Base.propertynames(square) == Base.propertynames(loaded)
            @test nvertices(square) == nvertices(loaded)
            @test dim(square) == dim(loaded)
            @test square == loaded
        end
        
        @testset "PolyhedralFan" begin
            nfsquare = normal_fan(cube(2))
            save_polyhedralfan(nfsquare, joinpath(path, "nfsquare.fan"))
            loaded = load_polyhedralfan(joinpath(path, "nfsquare.fan"))
            @test loaded isa PolyhedralFan
            @test Base.propertynames(nfsquare) == Base.propertynames(loaded)
            @test nrays(nfsquare) == nrays(loaded)
            @test nmaximal_cones(nfsquare) == nmaximal_cones(loaded)
            @test dim(nfsquare) == dim(loaded)
        end
        
        @testset "LinearProgram" begin
            P = cube(3)
            LP = LinearProgram(P,[3,-2,4];k=2,convention = :min)
            save_linearprogram(LP, joinpath(path, "lp.poly"))
            loaded = load_linearprogram(joinpath(path, "lp.poly"))
            @test loaded isa LinearProgram
            @test Base.propertynames(LP) == Base.propertynames(loaded)
            @test objective_function(LP) == objective_function(loaded)
            @test feasible_region(LP) == feasible_region(loaded)
        end
        
        @testset "SubdivisionOfPoints" begin
            moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
            moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])
            MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0)
            save_subdivisionofpoints(MOAE, joinpath(path, "moae.sop"))
            loaded = load_subdivisionofpoints(joinpath(path, "moae.sop"))
            @test loaded isa SubdivisionOfPoints
            @test Base.propertynames(MOAE) == Base.propertynames(loaded)
            @test nmaximal_cells(MOAE) == nmaximal_cells(loaded)
            @test npoints(MOAE) == npoints(loaded)
        end
    end

end
