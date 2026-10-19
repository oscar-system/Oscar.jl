cases = [
  (QQ, QQ(1//2), QQ(3//4), "Rational Coefficients"),
  (ZZ, ZZ(5), ZZ(7), "Integers Coefficients"),
  begin
    F, i = cyclotomic_field(4)
    (F, i, F(1//2), "Cyclotomic Field Coefficients")
  end,
]

@testset "Serialization.Algebras" begin
  mktempdir() do path
    for case in cases
      @testset "Free Associative Algebra over $(case[4])" begin
        A, g = free_associative_algebra(case[1], ["x","y"])
        f = case[2] * g[1] + case[3] * g[2] + case[3]
        test_save_load_roundtrip(path, f) do loaded
          @test loaded == f
        end

        # should add test for loading with parent

        I = ideal(A, [f, g[1] * g[2]])
        test_save_load_roundtrip(path, I) do loaded
          @test collect(gens(loaded)) == collect(gens(I))
        end
      end
    end

    # a quaternion algebra is saved as (a, b) and its elements as coordinates;
    # elements of one file share one parent
    @testset "Quaternion algebra over $(name)" for (name, K, a, b) in [
        ("QQ", QQ, QQ(-7), QQ(-97)),
        ("a number field", quadratic_field(5)..., QQ(-3)),
        ("a field of characteristic 2", GF(2, 2), gen(GF(2, 2)), 1)]
      A = quaternion_algebra(K, a, K(b))
      i, j, k = basis(A)[2:4]
      x = inv(K(3)) + 3 * i - a * j + 5 * k
      test_save_load_roundtrip(path, A) do loaded
        @test base_ring(loaded) === K
        @test Hecke.standard_form(loaded) == (a, K(b))
        li, lj = basis(loaded)[2:3]
        if characteristic(K) == 2
          @test li^2 + li == loaded(a) && lj^2 == loaded(b) && li * lj == lj * li + lj
        else
          @test li^2 == loaded(a) && lj^2 == loaded(b) && li * lj == -lj * li
        end
      end
      test_save_load_roundtrip(path, x) do loaded
        @test coefficients(loaded) == coefficients(x)
        @test Hecke.standard_form(parent(loaded)) == (a, K(b))
      end
      test_save_load_roundtrip(path, x; params=A) do loaded
        @test loaded == x
      end
      test_save_load_roundtrip(path, [x, x^2, one(A)]) do loaded
        @test parent(loaded[1]) === parent(loaded[2]) === parent(loaded[3])
        @test loaded[1]^2 == loaded[2]
      end
    end
  end
end

