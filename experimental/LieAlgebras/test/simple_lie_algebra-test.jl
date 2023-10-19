@testset "LieAlgebras.SimpleLieAlgebra" begin
  @testset "conformance tests" begin
    @testset "B2(QQ)" begin
      L = lie_algebra(QQ, :B, 2)
      lie_algebra_conformance_test(
        L, SimpleLieAlgebra{QQFieldElem}, SimpleLieAlgebraElem{QQFieldElem}
      )
    end

    @testset "A3(CF(4))" begin
      L = lie_algebra(cyclotomic_field(4)[1], :A, 3)
      lie_algebra_conformance_test(
        L, SimpleLieAlgebra{nf_elem}, SimpleLieAlgebraElem{nf_elem}
      )
    end
  end

  @testset "constructors and basic properties" begin
    L = lie_algebra(QQ, :A, 2)
    @test dim(L) == 8
    @test coefficient_ring(L) == QQ
    @test root_system(L) == RootSystem(:A, 2)
    @test root_system_type(L) == (:A, 2)
    @test characteristic(L) == 0
    @test chevalley_basis(L) == [
      [L([1, 0, 0, 0, 0, 0, 0, 0]), L([0, 1, 0, 0, 0, 0, 0, 0]), L([0, 0, 1, 0, 0, 0, 0, 0])],
      [L([0, 0, 0, 1, 0, 0, 0, 0]), L([0, 0, 0, 0, 1, 0, 0, 0]), L([0, 0, 0, 0, 0, 1, 0, 0])],
      [L([0, 0, 0, 0, 0, 0, 1, 0]), L([0, 0, 0, 0, 0, 0, 0, 1])],
      ]
    M = matrix_space(QQ, 8, 8)
    @test adjoint_matrix(L) == [M([0 0 0 0 0 0 -2 1; 0 0 0 0 0 0 0 0; 0 -1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0]),
                                M([0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 -2; 1 0 0 0 0 0 0 0; 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0]),
                                M([0 0 0 0 -1 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 -1 -1; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0; 0 0 0 0 0 1 0 0]),
                                M([0 0 0 0 0 0 0 0; 0 0 -1 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 2 -1; 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0; -1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0]),
                                M([0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 -1 2; 0 0 0 -1 0 0 0 0; 0 0 0 0 0 0 0 0; 0 -1 0 0 0 0 0 0]),
                                M([0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; -1 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 1; 0 0 -1 0 0 0 0 0; 0 0 -1 0 0 0 0 0]),
                                M([2 0 0 0 0 0 0 0; 0 -1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 -2 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0]),
                                M([-1 0 0 0 0 0 0 0; 0 2 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 -2 0 0 0; 0 0 0 0 0 -1 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0])]
  end
end
