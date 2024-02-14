using Oscar: _integer_variables
@testset "PolyhedralGeometry" begin
    
    mktempdir() do path
        @testset "Graph" begin
            G = complete_graph(4)
            test_save_load_roundtrip(path, G) do loaded
              @test loaded isa Graph{Undirected}
              @test Base.propertynames(G) == Base.propertynames(loaded)
              @test n_vertices(G) == n_vertices(loaded)
              @test n_edges(G) == n_edges(loaded)
            end
        end

        @testset "Cone" begin
            C = positive_hull([1 0; 0 1])
            test_save_load_roundtrip(path, C) do loaded
              @test loaded isa Cone
              @test Base.propertynames(C) == Base.propertynames(loaded)
              @test n_rays(C) == n_rays(loaded)
              @test dim(C) == dim(loaded)
              @test C == loaded
            end
        end

        @testset "Polyhedron" begin
            square = cube(2)
            test_save_load_roundtrip(path, square) do loaded
              @test loaded isa Polyhedron
              @test Base.propertynames(square) == Base.propertynames(loaded)
              @test n_vertices(square) == n_vertices(loaded)
              @test dim(square) == dim(loaded)
              @test square == loaded
            end
        end

        @testset "PolyhedralComplex" begin
            IM = IncidenceMatrix([[1,2,3],[1,3,4]])
            vr = [0 0; 1 0; 1 1; 0 1]
            PC = polyhedral_complex(IM, vr)
            test_save_load_roundtrip(path, PC) do loaded
              @test loaded isa PolyhedralComplex
              @test Base.propertynames(PC) == Base.propertynames(loaded)
              @test n_rays(PC) == n_rays(loaded)
              @test number_of_maximal_polyhedra(PC) == number_of_maximal_polyhedra(loaded)
              @test dim(PC) == dim(loaded)
            end
        end

        @testset "PolyhedralFan" begin
            nfsquare = normal_fan(cube(2))
            test_save_load_roundtrip(path, nfsquare) do loaded
              @test loaded isa PolyhedralFan
              @test Base.propertynames(nfsquare) == Base.propertynames(loaded)
              @test n_rays(nfsquare) == n_rays(loaded)
              @test number_of_maximal_cones(nfsquare) == number_of_maximal_cones(loaded)
              @test dim(nfsquare) == dim(loaded)
            end
        end

        @testset "LinearProgram" begin
            P = cube(3)
            LP = linear_program(P,[3,-2,4];k=2,convention = :min)
            test_save_load_roundtrip(path, LP) do loaded
              @test loaded isa LinearProgram
              @test Base.propertynames(LP) == Base.propertynames(loaded)
              @test objective_function(LP) == objective_function(loaded)
              @test feasible_region(LP) == feasible_region(loaded)
            end
        end
        
        @testset "MixedIntegerLinearProgram" begin
            P = cube(3)
            MILP = mixed_integer_linear_program(
                P,
                [3,-2,4];
                k=2,
                convention = :min,
                integer_variables=[1, 2]
            )
            test_save_load_roundtrip(path, MILP) do loaded
              @test loaded isa MixedIntegerLinearProgram
              @test Base.propertynames(MILP) == Base.propertynames(loaded)
              @test objective_function(MILP) == objective_function(loaded)
              @test feasible_region(MILP) == feasible_region(loaded)
              @test Oscar._integer_variables(MILP) == Oscar._integer_variables(loaded)
            end
        end

        @testset "SubdivisionOfPoints" begin
            moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
            moaeimnonreg0 = IncidenceMatrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])
            MOAE = subdivision_of_points(moaepts, moaeimnonreg0)
            test_save_load_roundtrip(path, MOAE) do loaded
              @test loaded isa SubdivisionOfPoints
              @test Base.propertynames(MOAE) == Base.propertynames(loaded)
              @test number_of_maximal_cells(MOAE) == number_of_maximal_cells(loaded)
              @test number_of_points(MOAE) == number_of_points(loaded)
            end
        end


        @testset "SimplicialComplex" begin
            cpp = complex_projective_plane()
            test_save_load_roundtrip(path, cpp) do loaded
              @test loaded isa SimplicialComplex
              @test Base.propertynames(cpp) == Base.propertynames(loaded)
              @test euler_characteristic(cpp) == euler_characteristic(loaded)
              @test n_vertices(cpp) == n_vertices(loaded)
            end
        end
    end

end
