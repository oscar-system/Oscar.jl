R, x = polynomial_ring(QQ, :x)
q = x^2 + 3//4
K, a = number_field(q)
Z7 = residue_ring(ZZ, 7)[1]
Z7t, t = polynomial_ring(GF(7), :t)
Fin, d = finite_field(t^2 + t + 3)
Frac = fraction_field(R)
A, (f,g) = free_associative_algebra(K, [:f, :g])

cases = [
  (QQ, [1 2; 3 4//5]),
  (Z7, [1 3; 4 2]),
  (K, [a a + 1; a - 1 0]),
  (Fin, [d d^3; d^2 0]),
  (Frac, [1 // x x^2; 3 0]),
  (A, [(f^3 + g^2) f^2; g (a*f + 1)]),
]

@testset "Matrices" begin
  mktempdir() do path
    @testset "Empty Matrix" begin
      m = zero_matrix(ZZ, 0, 2)
      test_save_load_roundtrip(path, m) do loaded
        @test m == loaded
      end
    end
    
    for case in cases
      @testset "Matrices over $(case[1])" begin
        m = matrix(case[1], case[2])
        test_save_load_roundtrip(path, m) do loaded
          @test loaded == m
        end

        test_save_load_roundtrip(path, m; params=parent(m)) do loaded
          @test loaded == m
        end
      end

      if case[1] isa Ring
        @testset "Sparse Matrices over $(case[1])" begin
          m = sparse_matrix(case[1], case[2])
          test_save_load_roundtrip(path, m) do loaded
            @test loaded == m
          end

          test_save_load_roundtrip(path, m; params=parent(m)) do loaded
            @test loaded == m
          end
        end
      end
    end
  end
end
