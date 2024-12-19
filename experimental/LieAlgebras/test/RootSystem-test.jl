@testset "LieAlgebras.RootSystem" begin
  @testset "property tests (only Serialization)" begin

    # merge this again with `root_system_property_tests` from test/LieTheory/RootSystem.jl once this is moved to src
    function root_system_property_tests(R::RootSystem, rk::Int, npositive_roots::Int)
      @testset "Serialization" begin
        mktempdir() do path
          test_save_load_roundtrip(path, R) do loaded
            # nothing, cause `R === loaded` anyway
          end

          if n_roots(R) >= 1
            r = positive_root(R, n_positive_roots(R))
            test_save_load_roundtrip(path, r) do loaded
              @test root_system(loaded) === R
              @test coefficients(loaded) == coefficients(r)
            end
          end

          test_save_load_roundtrip(path, positive_roots(R)) do loaded
            @test length(loaded) == n_positive_roots(R)
            @test all(
              coefficients(loaded[i]) == coefficients(root(R, i)) for
              i in 1:n_positive_roots(R)
            )
          end
        end
      end
    end

    @testset "rk 0" begin
      R = root_system(zero_matrix(ZZ, 0, 0))
      root_system_property_tests(R, 0, 0)
    end

    @testset "A_$n" for n in [1, 2, 6]
      R = root_system(:A, n)
      root_system_property_tests(R, n, binomial(n + 1, 2))
    end

    @testset "B_$n" for n in [2, 3, 6]
      R = root_system(:B, n)
      root_system_property_tests(R, n, n^2)
    end

    @testset "C_$n" for n in [2, 3, 6]
      R = root_system(:C, n)
      root_system_property_tests(R, n, n^2)
    end

    @testset "D_$n" for n in [4, 6]
      R = root_system(:D, n)
      root_system_property_tests(R, n, n^2 - n)
    end

    @testset "E_6" begin
      R = root_system(:E, 6)
      root_system_property_tests(R, 6, 36)
    end

    @testset "E_7" begin
      R = root_system(:E, 7)
      root_system_property_tests(R, 7, 63)
    end

    @testset "E_8" begin
      R = root_system(:E, 8)
      root_system_property_tests(R, 8, 120)
    end

    @testset "F_4" begin
      R = root_system(:F, 4)
      root_system_property_tests(R, 4, 24)
    end

    @testset "G_2" begin
      R = root_system(:G, 2)
      root_system_property_tests(R, 2, 6)
    end

    @testset "something mixed 1" begin
      cm = cartan_matrix((:A, 3), (:C, 3), (:E, 6), (:G, 2))
      for _ in 1:50
        i, j = rand(1:nrows(cm), 2)
        if i != j
          swap_rows!(cm, i, j)
          swap_cols!(cm, i, j)
        end
      end
      R = root_system(cm)
      root_system_property_tests(R, 3 + 3 + 6 + 2, binomial(3 + 1, 2) + 3^2 + 36 + 6)
    end

    @testset "something mixed 2" begin
      cm = cartan_matrix((:F, 4), (:B, 2), (:E, 7), (:G, 2))
      for _ in 1:50
        i, j = rand(1:nrows(cm), 2)
        if i != j
          swap_rows!(cm, i, j)
          swap_cols!(cm, i, j)
        end
      end
      R = root_system(cm)
      root_system_property_tests(R, 4 + 2 + 7 + 2, 24 + 2^2 + 63 + 6)
    end
  end
end
