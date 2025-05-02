@testset "LieAlgebras.LieAlgebraIdeal" begin
  @testset "Constructor and basic properties" begin
    let
      L = general_linear_lie_algebra(QQ, 3)
      b = lie_algebra(
        sub(
          L, [basis(L, 1), basis(L, 2), basis(L, 3), basis(L, 5), basis(L, 6), basis(L, 9)]
        ),
      )[1]
      n = ideal(b, [basis(b, 2), basis(b, 3), basis(b, 5)])

      let L = b, I = n
        @test I == ideal(L, basis(I))
        @test base_lie_algebra(I) == L
        @test length(gens(I)) == ngens(I)
        @test length(basis(I)) == dim(I)
        @test all(in(I), gens(I))
        @test all(in(I), basis(I))

        LI, emb = lie_algebra(I)
        @test dim(LI) == dim(I)
        @test domain(emb) == LI
        @test codomain(emb) == L
        @test image(emb) == sub(L, I)
      end

      @test dim(b) == 6
      @test dim(n) == 3

      @test n == ideal(b, [basis(b, 2), basis(b, 3), basis(b, 5)])
      @test n == ideal(b, [basis(b, 2) + basis(b, 3), basis(b, 3), basis(b, 5)])
      @test n ==
        ideal(b, [basis(b, 2) + basis(b, 3), basis(b, 3), basis(b, 5)]; is_basis=true)
    end

    let # Example where ideal basis is only found after two steps
      sc = zeros(QQ, 4, 4, 4)
      sc[1, 2, 3] = 1
      sc[2, 1, 3] = -1
      sc[1, 3, 4] = 1
      sc[3, 1, 4] = -1
      L = lie_algebra(QQ, sc, ["a", "b", "c", "d"])
      a, b, c, d = basis(L)

      let I = ideal(L, b)
        @test I == ideal(L, [b])
        @test dim(I) == 3
        @test ngens(I) == 1

        @test I == ideal(L, basis(I))
        @test base_lie_algebra(I) == L
        @test length(gens(I)) == ngens(I)
        @test length(basis(I)) == dim(I)
        @test all(in(I), gens(I))
        @test all(in(I), basis(I))
      end

      let I = ideal(L, [a + b + c + d])
        @test dim(I) == 3
        @test ngens(I) == 1

        @test I == ideal(L, basis(I))
        @test base_lie_algebra(I) == L
        @test length(gens(I)) == ngens(I)
        @test length(basis(I)) == dim(I)
        @test all(in(I), gens(I))
        @test all(in(I), basis(I))
      end
    end
  end
end
