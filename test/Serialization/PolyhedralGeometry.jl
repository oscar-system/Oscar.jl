@testset "PolyhedralGeometry" begin
    
    mktempdir() do path
        @testset "Graph" begin
            G = Graphs.complete_graph(4)
            save(G, joinpath(path, "c4.graph"))
            loaded = load(joinpath(path, "c4.graph"))
            @test loaded isa Graphs.Graph{Graphs.Undirected}
            @test Base.propertynames(G) == Base.propertynames(loaded)
            @test Graphs.nv(G) == Graphs.nv(loaded)
            @test Graphs.ne(G) == Graphs.ne(loaded)
        end

        @testset "Cone" begin
            C = positive_hull([1 0; 0 1])
            save(C, joinpath(path, "posorth.cone"))
            loaded = load(joinpath(path, "posorth.cone"))
            @test loaded isa Cone
            @test Base.propertynames(C) == Base.propertynames(loaded)
            @test nrays(C) == nrays(loaded)
            @test dim(C) == dim(loaded)
            @test C == loaded
        end

        @testset "Polyhedron" begin
            square = cube(2)
            save(square, joinpath(path, "square.poly"))
            loaded = load(joinpath(path, "square.poly"))
            @test loaded isa Polyhedron
            @test Base.propertynames(square) == Base.propertynames(loaded)
            @test nvertices(square) == nvertices(loaded)
            @test dim(square) == dim(loaded)
            @test square == loaded
        end

        @testset "PolyhedralComplex" begin
            IM = IncidenceMatrix([[1,2,3],[1,3,4]])
            vr = [0 0; 1 0; 1 1; 0 1]
            PC = PolyhedralComplex(IM, vr)
            save(PC, joinpath(path, "docu.pc"))
            loaded = load(joinpath(path, "docu.pc"))
            @test loaded isa PolyhedralComplex
            @test Base.propertynames(PC) == Base.propertynames(loaded)
            @test nrays(PC) == nrays(loaded)
            @test n_maximal_polyhedra(PC) == n_maximal_polyhedra(loaded)
            @test dim(PC) == dim(loaded)
        end

        @testset "PolyhedralFan" begin
            nfsquare = normal_fan(cube(2))
            save(nfsquare, joinpath(path, "nfsquare.fan"))
            loaded = load(joinpath(path, "nfsquare.fan"))
            @test loaded isa PolyhedralFan
            @test Base.propertynames(nfsquare) == Base.propertynames(loaded)
            @test nrays(nfsquare) == nrays(loaded)
            @test n_maximal_cones(nfsquare) == n_maximal_cones(loaded)
            @test dim(nfsquare) == dim(loaded)
        end

        @testset "LinearProgram" begin
            P = cube(3)
            LP = LinearProgram(P,[3,-2,4];k=2,convention = :min)
            save(LP, joinpath(path, "lp.poly"))
            loaded = load(joinpath(path, "lp.poly"))
            @test loaded isa LinearProgram
            @test Base.propertynames(LP) == Base.propertynames(loaded)
            @test objective_function(LP) == objective_function(loaded)
            @test feasible_region(LP) == feasible_region(loaded)
        end
        
        @testset "SubdivisionOfPoints" begin
            moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
            moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])
            MOAE = SubdivisionOfPoints(moaepts, moaeimnonreg0)
            save(MOAE, joinpath(path, "moae.sop"))
            loaded = load(joinpath(path, "moae.sop"))
            @test loaded isa SubdivisionOfPoints
            @test Base.propertynames(MOAE) == Base.propertynames(loaded)
            @test n_maximal_cells(MOAE) == n_maximal_cells(loaded)
            @test npoints(MOAE) == npoints(loaded)
        end


        @testset "SimplicialComplex" begin
            cpp = complex_projective_plane()
            save(cpp, joinpath(path, "cpp.top"))
            loaded = load(joinpath(path, "cpp.top"))
            @test loaded isa SimplicialComplex
            @test Base.propertynames(cpp) == Base.propertynames(loaded)
            @test euler_characteristic(cpp) == euler_characteristic(loaded)
            @test nvertices(cpp) == nvertices(loaded)
        end
    end

end
