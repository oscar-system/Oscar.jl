@testset "LieTheory" begin
  mktempdir() do path
    @testset "RootSystem" begin
      for R in [
        root_system(zero_matrix(ZZ, 0, 0)), # rk 0
        root_system(:A, 1),
        root_system(:A, 6),
        root_system(:B, 2),
        root_system(:C, 5),
        root_system(:D, 4),
        root_system(:E, 6),
        root_system(:F, 4),
        root_system(:G, 2),
        begin # something mixed
          cm = cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2))
          for _ in 1:50
            i, j = rand(1:nrows(cm), 2)
            if i != j
              swap_rows!(cm, i, j)
              swap_cols!(cm, i, j)
            end
          end
          root_system(cm)
        end,
      ]
        test_save_load_roundtrip(path, R) do loaded
          # nothing, cause `R === loaded` anyway
        end

        if rank(R) >= 1
          r = positive_root(R, n_positive_roots(R))
          test_save_load_roundtrip(path, r) do loaded
            @test root_system(loaded) === R
            @test coefficients(loaded) == coefficients(r)
          end

          cr = positive_coroot(R, n_positive_roots(R))
          test_save_load_roundtrip(path, cr) do loaded
            @test root_system(loaded) === R
            @test coefficients(loaded) == coefficients(cr)
          end
        end

        test_save_load_roundtrip(path, positive_roots(R)) do loaded
          @test length(loaded) == n_positive_roots(R)
          @test all(
            coefficients(loaded[i]) == coefficients(positive_root(R, i)) for
            i in 1:n_positive_roots(R)
          )
        end

        test_save_load_roundtrip(path, negative_coroots(R)) do loaded
          @test length(loaded) == n_positive_roots(R)
          @test all(
            coefficients(loaded[i]) == coefficients(negative_coroot(R, i)) for
            i in 1:n_positive_roots(R)
          )
        end
      end
    end

    @testset "WeightLattice" begin
      @testset "simple saving and loading" begin
        for R in [
          root_system(zero_matrix(ZZ, 0, 0)), # rk 0
          root_system((:A, 2), (:B, 4)),
        ]
          P = weight_lattice(R)

          test_save_load_roundtrip(path, P) do loaded
            # nothing, cause `P === loaded` anyway
          end

          w = WeightLatticeElem(P, [rand(ZZ, -10:10) for _ in 1:rank(P)])
          test_save_load_roundtrip(path, w) do loaded
            @test parent(loaded) === P
            @test coefficients(loaded) == coefficients(w)
          end

          test_save_load_roundtrip(path, gens(P)) do loaded
            @test length(loaded) == rank(P)
            @test all(coefficients(loaded[i]) == coefficients(gen(P, i)) for i in 1:rank(R))
          end
        end
      end

      @testset "cyclic reference between R and P survives" begin
        Oscar.reset_global_serializer_state()

        R_filename = joinpath(path, "R.mrdi")
        P_filename = joinpath(path, "P.mrdi")

        R = root_system(:D, 5)
        P = weight_lattice(R)

        save(R_filename, R)
        save(P_filename, P)

        Oscar.reset_global_serializer_state()

        loaded_R = load(R_filename)
        loaded_P = load(P_filename)

        @test loaded_R === root_system(loaded_P)
        @test loaded_P === weight_lattice(loaded_R)

        loaded_R = loaded_P = nothing # unset all references

        Oscar.reset_global_serializer_state()

        loaded_P = load(P_filename)
        loaded_R = load(R_filename)

        @test loaded_R === root_system(loaded_P)
        @test loaded_P === weight_lattice(loaded_R)

        loaded_R = loaded_P = nothing # unset all references        
      end
    end

    @testset "WeylGroup" begin
      @testset "simple saving and loading" begin
        for W in [
          weyl_group(zero_matrix(ZZ, 0, 0)), # rk 0
          weyl_group((:A, 2), (:B, 4)),
        ]
          test_save_load_roundtrip(path, W) do loaded
            # nothing, cause `W === loaded` anyway
          end

          x = rand(W)
          test_save_load_roundtrip(path, x) do loaded
            @test parent(loaded) === W
            @test word(loaded) == word(x)
          end

          test_save_load_roundtrip(path, gens(W)) do loaded
            @test length(loaded) == ngens(W)
            @test all(
              word(loaded[i]) == word(gen(W, i)) for i in 1:ngens(W)
            )
          end
        end
      end

      @testset "cyclic reference between R and W survives" begin
        Oscar.reset_global_serializer_state()

        R_filename = joinpath(path, "R.mrdi")
        W_filename = joinpath(path, "W.mrdi")

        R = root_system(:D, 5)
        W = weyl_group(R)

        save(R_filename, R)
        save(W_filename, W)

        Oscar.reset_global_serializer_state()

        loaded_R = load(R_filename)
        loaded_W = load(W_filename)

        @test loaded_R === root_system(loaded_W)
        @test loaded_W === weyl_group(loaded_R)

        loaded_R = loaded_W = nothing # unset all references

        Oscar.reset_global_serializer_state()

        loaded_W = load(W_filename)
        loaded_R = load(R_filename)

        @test loaded_R === root_system(loaded_W)
        @test loaded_W === weyl_group(loaded_R)

        loaded_R = loaded_W = nothing # unset all references        
      end
    end
  end
end
