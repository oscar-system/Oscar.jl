@testset "LieAlgebras.AbstractLieAlgebra" begin
  function sl2_struct_consts(R::Field)
    sc = zeros(R, 3, 3, 3)
    sc[1, 2, 3] = R(1)
    sc[2, 1, 3] = R(-1)
    sc[3, 1, 1] = R(2)
    sc[1, 3, 1] = R(-2)
    sc[3, 2, 2] = R(-2)
    sc[2, 3, 2] = R(2)
    return sc
  end

  @testset "conformance tests" begin
    @testset "0-dim Lie algebra /QQ" begin
      L = lie_algebra(QQ, Matrix{SRow{QQFieldElem}}(undef, 0, 0), Symbol[])
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
      @test is_abelian(L)
    end

    @testset "sl_2(QQ) using structure constants" begin
      L = lie_algebra(QQ, sl2_struct_consts(QQ), ["e", "f", "h"])
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "sl_2(CF(4)) using structure constants" begin
      CF4 = cyclotomic_field(4)[1]
      L = lie_algebra(CF4, sl2_struct_consts(CF4), ["e", "f", "h"])
      lie_algebra_conformance_test(
        L,
        AbstractLieAlgebra{AbsSimpleNumFieldElem},
        AbstractLieAlgebraElem{AbsSimpleNumFieldElem},
      )
    end

    @testset "4-dim abelian Lie algebra /QQ" begin
      L = abelian_lie_algebra(AbstractLieAlgebra, QQ, 4)
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
      @test is_abelian(L)
    end

    @testset "A_4(QQ)" begin
      L = lie_algebra(QQ, :A, 4)
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "B_3(QQ)" begin
      L = lie_algebra(QQ, :B, 3)
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "A_4(CF(4))" begin
      L = lie_algebra(cyclotomic_field(4)[1], :A, 4)
      lie_algebra_conformance_test(
        L,
        AbstractLieAlgebra{AbsSimpleNumFieldElem},
        AbstractLieAlgebraElem{AbsSimpleNumFieldElem},
      )
    end

    @testset "B_3(CF(4))" begin
      L = lie_algebra(cyclotomic_field(4)[1], :B, 3)
      lie_algebra_conformance_test(
        L,
        AbstractLieAlgebra{AbsSimpleNumFieldElem},
        AbstractLieAlgebraElem{AbsSimpleNumFieldElem},
      )
    end
  end
end
