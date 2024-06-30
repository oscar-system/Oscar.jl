@testset "LieAlgebras.DirectSumLieAlgebra" begin
  @testset "conformance tests" begin
    @testset "sl_2(QQ)" begin
      S1 = special_linear_lie_algebra(QQ, 2)
      L = direct_sum(S1)

      lie_algebra_conformance_test(
        L, DirectSumLieAlgebra{QQFieldElem}, DirectSumLieAlgebraElem{QQFieldElem}
      )
      @test dim(L) == dim(S1)

      L = direct_sum([S1])
      lie_algebra_conformance_test(
        L, DirectSumLieAlgebra{QQFieldElem}, DirectSumLieAlgebraElem{QQFieldElem}
      )
      @test dim(L) == dim(S1)
    end

    @testset "sl_2(GF(7)) ⊕ gl_3(GF(7))" begin
      GF7 = GF(7)
      S1 = special_linear_lie_algebra(GF7, 2)
      S2 = general_linear_lie_algebra(GF7, 3)

      L = direct_sum(S1, S2)
      lie_algebra_conformance_test(
        L, DirectSumLieAlgebra{FqFieldElem}, DirectSumLieAlgebraElem{FqFieldElem}
      )
      @test dim(L) == dim(S1) + dim(S2)

      L = direct_sum([S1, S2])
      lie_algebra_conformance_test(
        L, DirectSumLieAlgebra{FqFieldElem}, DirectSumLieAlgebraElem{FqFieldElem}
      )
      @test dim(L) == dim(S1) + dim(S2)
    end

    @testset "sl_3(CF(4)) ⊕ B3(CF(4)) ⊕ G2(CF(4))" begin
      CF4 = cyclotomic_field(4)[1]
      S1 = special_linear_lie_algebra(CF4, 3)
      S2 = lie_algebra(CF4, :B, 3)
      S3 = lie_algebra(CF4, :G, 2)
      L = direct_sum(S1, S2, S3)

      lie_algebra_conformance_test(
        L, DirectSumLieAlgebra{AbsSimpleNumFieldElem},
        DirectSumLieAlgebraElem{AbsSimpleNumFieldElem},
      )
      @test dim(L) == dim(S1) + dim(S2) + dim(S3)

      L = direct_sum([S1, S2, S3])
      lie_algebra_conformance_test(
        L, DirectSumLieAlgebra{AbsSimpleNumFieldElem},
        DirectSumLieAlgebraElem{AbsSimpleNumFieldElem},
      )
      @test dim(L) == dim(S1) + dim(S2) + dim(S3)
    end

    @testset "empty sum over QQ" begin
      L = direct_sum(QQ, LieAlgebra{QQFieldElem}[])

      lie_algebra_conformance_test(
        L,
        DirectSumLieAlgebra{QQFieldElem},
        DirectSumLieAlgebraElem{QQFieldElem},
      )
      @test dim(L) == 0
      @test is_abelian(L)
    end
  end

  @testset "properties" begin
    S1 = special_linear_lie_algebra(QQ, 2)
    S2 = lie_algebra(QQ, :A, 1)
    L = direct_sum(S1, S2)

    @test length(canonical_injections(L)) == 2
    @test all(i -> canonical_injection(L, i) == canonical_injections(L)[i], 1:2)
    @test all(i -> codomain(canonical_injection(L, i)) == L, 1:2)
    @test domain(canonical_injection(L, 1)) == S1
    @test domain(canonical_injection(L, 2)) == S2

    @test length(canonical_projections(L)) == 2
    @test all(i -> canonical_projection(L, i) == canonical_projections(L)[i], 1:2)
    @test all(i -> domain(canonical_projection(L, i)) == L, 1:2)
    @test codomain(canonical_projection(L, 1)) == S1
    @test codomain(canonical_projection(L, 2)) == S2

    for _ in 1:5
      x = L(rand(-10:10, dim(L)))
      @test x ==
        sum([canonical_injection(L, i)(canonical_projection(L, i)(x)) for i in 1:2])

      x1 = S1(rand(-10:10, dim(S1)))
      x2 = S2(rand(-10:10, dim(S2)))
      @test x1 == canonical_projection(L, 1)(canonical_injection(L, 1)(x1))
      @test x2 == canonical_projection(L, 2)(canonical_injection(L, 2)(x2))
      @test iszero(canonical_projection(L, 2)(canonical_injection(L, 1)(x1)))
      @test iszero(canonical_projection(L, 1)(canonical_injection(L, 2)(x2)))
    end

    @test !is_abelian(L)
  end
end
