@testset "Groups" begin
  mktempdir() do path
    @testset "Permutation groups" begin
      G = symmetric_group(5)

      # single element
      x = gen(G, 1)
      test_save_load_roundtrip(path, x) do loaded
        @test x == loaded
      end

      # full symmetric group
      test_save_load_roundtrip(path, G) do loaded
        @test G == loaded
      end

      # subgroup of a symmetric group
      P = sylow_subgroup(G, 2)[1]
      test_save_load_roundtrip(path, P) do loaded
        @test P == loaded
      end

      # element and group together
      v = [x, P, G]
      test_save_load_roundtrip(path, v) do loaded
        @test v == loaded
      end
    end

    @testset "Finitely presented group" begin
      F = free_group(2)

      # single element
      x = gen(F, 1)
      test_save_load_roundtrip(path, x) do loaded
        @test x == loaded
      end

      # full free group
      test_save_load_roundtrip(path, F) do loaded
        @test F == loaded
      end

      # subgroup of free group
      U = sub(F, [gen(F, 1)])[1]
      test_save_load_roundtrip(path, U) do loaded
        @test U == loaded
      end

      # full f.p. group with relators
      G = quo(F, [x^2])[1]
      test_save_load_roundtrip(path, G) do loaded
        @test G == loaded
      end
      
      # subgroup of f.p. group with relators
      S = sub(G, [gen(G, 1)])[1]
      test_save_load_roundtrip(path, S) do loaded
        @test gens(S) == gens(loaded)
      end

      # various objects together
      v = [x, F, U, G, gens(S)]
      test_save_load_roundtrip(path, v) do loaded
        @test v == loaded
      end

      # Load objects from another Oscar session,
      # and check consistency.
      # (The file contains the vector `v`.)
      x, F, U, G, Sgens = load(joinpath(@__DIR__, "Groups_input.json"))
      @test x == gen(F, 1)
      @test x == gen(U, 1)
      @test Sgens[1] == gen(G, 1)
    end
  end
end
