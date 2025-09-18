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

@testset "QuadFormAndIsom" begin
  mktempdir() do path
    @testset "ZZLatWithIsom" begin
      L = integer_lattice(gram = matrix(QQ, 0, 0, []))
      Lf = integer_lattice_with_isometry(L; neg=true)
      Vf = ambient_space(Lf)
      test_save_load_roundtrip(path, Vf) do loaded
        @test Vf == loaded
      end
      test_save_load_roundtrip(path, Lf) do loaded
        @test Lf == loaded
      end

      M = integer_lattice(; gram = QQ[1 2; 2 1])
      g = QQ[4 -1; 1 0]
      Mg = integer_lattice_with_isometry(M, g)
      test_save_load_roundtrip(path, Mg) do loaded
        @test Mg == loaded
      end

      B = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 1]);
      G = matrix(QQ, 8, 8, [-4 2 0 0 0 0 0 0; 2 -4 2 0 0 0 0 0; 0 2 -4 2 0 0 0 2; 0 0 2 -4 2 0 0 0; 0 0 0 2 -4 2 0 0; 0 0 0 0 2 -4 2 0; 0 0 0 0 0 2 -4 0; 0 0 2 0 0 0 0 -4]);
      N = integer_lattice(B; gram=G);
      h = matrix(QQ, 8, 8, [1 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; -2 -4 -6 -4 -3 -2 -1 -3; 2 4 6 5 4 3 2 3; -1 -2 -3 -3 -3 -2 -1 -1; 0 0 0 0 1 0 0 0; 1 2 3 3 2 1 0 2]);
      Nh = integer_lattice_with_isometry(N, h);
      C = coinvariant_lattice(Nh)
      test_save_load_roundtrip(path, C) do loaded
        @test C == loaded
      end
    end
  end
end
