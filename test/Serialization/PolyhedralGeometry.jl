using Oscar: _integer_variables
@testset "PolyhedralGeometry" begin
  Qx, x = QQ[:x]
  F, a = embedded_number_field(x^2 - 2, -1.0)

  mktempdir() do path
    @testset "Graph" begin
      G = complete_graph(4)
      test_save_load_roundtrip(path, G) do loaded
        @test n_vertices(G) == n_vertices(loaded)
        @test n_edges(G) == n_edges(loaded)
      end
    end

    @testset "Matroid" begin
      @testset "Fano" begin
        M = fano_matroid()
        test_save_load_roundtrip(path, M) do loaded
          @test sort(bases(M)) == sort(bases(loaded))
          @test length(M) == length(loaded)
          @test rank(M) == rank(loaded)
        end
      end
      @testset "uniform" begin
        M = uniform_matroid(2, 4)
        test_save_load_roundtrip(path, M) do loaded
          @test sort(bases(M)) == sort(bases(loaded))
          @test length(M) == length(loaded)
          @test rank(M) == rank(loaded)
        end
      end
    end

    @testset "Cone" begin
      C = positive_hull([1 0; 0 1])
      test_save_load_roundtrip(path, C) do loaded
        @test n_rays(C) == n_rays(loaded)
        @test dim(C) == dim(loaded)
        @test C == loaded
      end

      CF = positive_hull(F, [F(1) F(0); F(0) F(1)])
      test_save_load_roundtrip(path, CF) do loaded
        @test n_rays(CF) == n_rays(loaded)
        @test dim(CF) == dim(loaded)
        @test CF == loaded
      end
    end

    @testset "Polyhedron" begin
      square = cube(2)
      f_vector(square)
      test_save_load_roundtrip(path, square) do loaded
        @test n_vertices(square) == n_vertices(loaded)
        @test dim(square) == dim(loaded)
        @test square == loaded
        @test Polymake.exists(Oscar.pm_object(loaded), "HASSE_DIAGRAM.DECORATION")
      end

      n2 = (QQBarField()(5))^(QQ(4//5))
      c = cube(QQBarField(), 3, -1, n2)
      f_vector(c)
      lattice_points(c)
      test_save_load_roundtrip(path, c) do loaded
        @test n_vertices(c) == n_vertices(loaded)
        @test dim(c) == dim(loaded)
        @test c == loaded
        @test Polymake.exists(Oscar.pm_object(loaded), "HASSE_DIAGRAM.DECORATION")
      end

      d_hedron = dodecahedron()
      facets(d_hedron)
      vertices(d_hedron)
      f_vector(d_hedron)
      lattice_points(d_hedron)

      # this type needs to be any since the values have different type_params
      dict_ps = Dict{String, Any}(
        "unprecise" => polyhedron(
          Polymake.common.convert_to{Float64}(Oscar.pm_object(d_hedron))
        ),
        "precise" => d_hedron
      )

      test_save_load_roundtrip(path, dict_ps) do loaded
        @test dict_ps["precise"] == loaded["precise"]
      end
    end

    @testset "PolyhedralComplex" begin
      IM = incidence_matrix([[1,2,3],[1,3,4]])
      vr = [0 0; 1 0; 1 1; 0 1]
      PC = polyhedral_complex(IM, vr)
      test_save_load_roundtrip(path, PC) do loaded
        @test n_rays(PC) == n_rays(loaded)
        @test number_of_maximal_polyhedra(PC) == number_of_maximal_polyhedra(loaded)
        @test dim(PC) == dim(loaded)
      end

      vr_F = F.([0 0; 1 0; 1 1; 0 1])
      PC_F = polyhedral_complex(IM, vr_F)
      test_save_load_roundtrip(path, PC_F) do loaded
        @test n_rays(PC_F) == n_rays(loaded)
        @test number_of_maximal_polyhedra(PC_F) == number_of_maximal_polyhedra(loaded)
        @test dim(PC_F) == dim(loaded)
      end
    end

    @testset "PolyhedralFan" begin
      nfsquare = normal_fan(cube(2))
      test_save_load_roundtrip(path, nfsquare) do loaded
        @test n_rays(nfsquare) == n_rays(loaded)
        @test number_of_maximal_cones(nfsquare) == number_of_maximal_cones(loaded)
        @test dim(nfsquare) == dim(loaded)
      end

      nfdodecahedron = normal_fan(dodecahedron())
      # add some extra properties to check save / load
      Polymake.give(Oscar.pm_object(nfdodecahedron), :MAXIMAL_CONES_FACETS)
      test_save_load_roundtrip(path, nfdodecahedron) do loaded
        @test n_rays(nfdodecahedron) == n_rays(loaded)
        @test number_of_maximal_cones(nfdodecahedron) == number_of_maximal_cones(loaded)
        @test dim(nfdodecahedron) == dim(loaded)
      end
    end

    @testset "LinearProgram" begin
      P = cube(3)
      LP = linear_program(P,[3,-2,4];k=2,convention = :min)
      test_save_load_roundtrip(path, LP) do loaded
        @test objective_function(LP) == objective_function(loaded)
        @test feasible_region(LP) == feasible_region(loaded)
      end

      serializer=Oscar.Serialization.LPSerializer(joinpath(path, "original"))
      test_save_load_roundtrip(path, LP; serializer=serializer) do loaded
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
        @test objective_function(MILP) == objective_function(loaded)
        @test feasible_region(MILP) == feasible_region(loaded)
        @test Oscar._integer_variables(MILP) == Oscar._integer_variables(loaded)
      end
    end

    @testset "SubdivisionOfPoints" begin
      moaepts = [4 0 0; 0 4 0; 0 0 4; 2 1 1; 1 2 1; 1 1 2]
      moaeimnonreg0 = incidence_matrix([[4,5,6],[1,4,2],[2,4,5],[2,3,5],[3,5,6],[1,3,6],[1,4,6]])
      MOAE = subdivision_of_points(moaepts, moaeimnonreg0)
      test_save_load_roundtrip(path, MOAE) do loaded
        @test number_of_maximal_cells(MOAE) == number_of_maximal_cells(loaded)
        @test number_of_points(MOAE) == number_of_points(loaded)
      end

      MOAEF = subdivision_of_points(F, F.(moaepts), moaeimnonreg0)
      test_save_load_roundtrip(path, MOAEF) do loaded
        @test number_of_maximal_cells(MOAEF) == number_of_maximal_cells(loaded)
        @test number_of_points(MOAEF) == number_of_points(loaded)
      end
    end

    @testset "SimplicialComplex" begin
      cpp = complex_projective_plane()
      test_save_load_roundtrip(path, cpp) do loaded
        @test Base.propertynames(cpp) == Base.propertynames(loaded)
        @test euler_characteristic(cpp) == euler_characteristic(loaded)
        @test n_vertices(cpp) == n_vertices(loaded)
      end
    end
  end
end
