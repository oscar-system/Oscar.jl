@testset "LieAlgebras.AbstractLieAlgebra" begin
  function sl2_struct_consts(R::Oscar.Ring)
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
        L, AbstractLieAlgebra{nf_elem}, AbstractLieAlgebraElem{nf_elem}
      )
    end

    @testset "A_4(QQ)" begin
      L = lie_algebra(QQ, ('A', 4))
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "B_3(QQ)" begin
      L = lie_algebra(QQ, ('B', 3))
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{QQFieldElem}, AbstractLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "A_4(CF(4))" begin
      L = lie_algebra(cyclotomic_field(4)[1], ('A', 4))
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{nf_elem}, AbstractLieAlgebraElem{nf_elem}
      )
    end

    @testset "B_3(CF(4))" begin
      L = lie_algebra(cyclotomic_field(4)[1], ('B', 3))
      lie_algebra_conformance_test(
        L, AbstractLieAlgebra{nf_elem}, AbstractLieAlgebraElem{nf_elem}
      )
    end
  end
end
