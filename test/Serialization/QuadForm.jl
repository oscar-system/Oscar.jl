@testset "QuadForm" begin
  mktempdir() do path
    @testset "ZZLat" begin
      L = integer_lattice(ZZ[1 2;],gram=ZZ[0 1;1 0])
      V = ambient_space(L)

      test_save_load_roundtrip(path, L) do loaded
        @test L == loaded
      end
      test_save_load_roundtrip(path, V) do loaded
        @test V == loaded
      end
    end

    @testset "QuadSpace" begin
      F,a = quadratic_field(2)
      V = quadratic_space(F, F[a;])

      test_save_load_roundtrip(path, V) do loaded
        @test gram_matrix(V) == gram_matrix(loaded)
      end
      R, x = polynomial_ring(F, :x)
      FF, b =extension_field(x^2 - a)
      VFF = quadratic_space(FF,FF[b;])
      test_save_load_roundtrip(path, VFF) do loaded
        @test gram_matrix(VFF) == gram_matrix(loaded)
      end
    end
  end
end
