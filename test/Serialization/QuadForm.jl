@testset "QuadForm" begin
    
  mktempdir() do path


    @testset "ZLat" begin
      L = Zlattice(ZZ[1 2;],gram=ZZ[0 1;1 0])
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
          F1 = base_ring(loaded)
          iso = hom(F1,F,gen(F))
          @test gram_matrix(V) == map_entries(iso,gram_matrix(loaded))
        end
      R, x = PolynomialRing(F, "x")
      FF, b =extension_field(x^2 - a)
      VFF = quadratic_space(FF,FF[b;])
      test_save_load_roundtrip(path, VFF) do loaded
          FF1 = base_ring(loaded)
          F1 = base_field(FF1)
          isoF = hom(F1,F,gen(F))
          isoFF = hom(FF1,FF,isoF,gen(FF))
          @test gram_matrix(VFF) == map_entries(isoFF,gram_matrix(loaded))
      end
    end
  end
end
