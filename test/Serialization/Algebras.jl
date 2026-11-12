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
    # with uses_id, only check_func sees an algebra rebuilt from the file
    @testset "Quaternion algebra over $(name)" for (name, K, a, b) in [
        ("QQ", QQ, QQ(-7), QQ(-97)),
        ("a number field", quadratic_field(5)..., QQ(-3)),
        begin
          F, s = quadratic_field(5)
          Fy, y = polynomial_ring(F, :y)
          L, t = number_field(y^2 - s, :t)
          ("a relative number field", L, L(-1), t)
        end,
        begin
          Qz, z = polynomial_ring(QQ, :z)
          N, (r, s) = number_field([z^2 - 2, z^2 - 5])
          ("a non-simple number field", N, r, s)
        end,
        ("a field of characteristic 2", GF(2, 2), gen(GF(2, 2)), 1)]
      A = quaternion_algebra(K, a, K(b))
      i, j, k = basis(A)[2:4]
      x = inv(K(3)) + 3 * i - a * j + 5 * k
      function is_rebuilt(loaded)
        li, lj = basis(loaded)[2:3]
        relations = if characteristic(K) == 2
          li^2 + li == loaded(a) && lj^2 == loaded(b) && li * lj == lj * li + lj
        else
          li^2 == loaded(a) && lj^2 == loaded(b) && li * lj == -lj * li
        end
        return loaded !== A && base_ring(loaded) === K &&
          Hecke.standard_form(loaded) == (a, K(b)) && relations
      end
      test_save_load_roundtrip(path, A; check_func=is_rebuilt) do loaded
        # nothing to test here, since `loaded === A`
      end
      test_save_load_roundtrip(path, x) do loaded
        @test coefficients(loaded) == coefficients(x)
      end
      test_save_load_roundtrip(path, x; params=A) do loaded
        @test loaded == x
      end
    end
  end
end

