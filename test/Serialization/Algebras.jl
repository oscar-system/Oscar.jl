cases = [
  (QQ, QQ(1//2), QQ(3//4), "Rational Coefficients"),
  (ZZ, ZZ(5), ZZ(7), "Integers Coefficients"),
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
  end
end

