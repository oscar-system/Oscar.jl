@testset "Fields" begin
  mktempdir() do path
    Qx, x = QQ[:x]
    K1, _ = number_field(x^2 + 5)
    K2, _ = number_field([x^2 + 5, x^2 + 7])
    K3, _ = number_field(x^2 + 7)
    Kt, t = polynomial_ring(K1, cached = false)
    L1, = number_field(t^2 + 2)
    L2, = number_field([t^2 + 2, t^2 + 3])
    Kt, t = polynomial_ring(K2, cached = false)
    L3, = number_field(t^2 + 2)
    L4, = number_field([t^2 + 2, t^2 + 3])

    @testset "Field Embeddings" begin
      for K in (K1, K2, L1, L2, L3, L4)
        for E in complex_embeddings(K)
          test_save_load_roundtrip(path, E) do loaded
            @test loaded == E
          end

          L, = Oscar.Hecke.embedded_field(K, E)
          test_save_load_roundtrip(path, L) do loaded
            @test number_field(loaded) == K
            @test embedding(loaded) == E
          end
        end
      end
    end

    @testset "Vector of Fields" begin
      v = [K1, K3]
      test_save_load_roundtrip(path, v) do loaded
        @test loaded == v
      end
    end

    @testset "QQBar" begin
      test_save_load_roundtrip(path, QQBarField()) do loaded
        @test loaded isa QQBarField
      end
    end
  end
end
