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

        if n_roots(R) >= 1
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
  end
end
