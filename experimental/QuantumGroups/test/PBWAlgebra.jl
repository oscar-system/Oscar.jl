function make_Usl3()
  R, y = polynomial_ring(QQ, :y => 1:6)
  rels = [
    y[1] * y[2] + y[4], y[1] * y[3], y[1] * y[4], y[1] * y[5] + y[6], y[1] * y[6],
    y[2] * y[3] + y[5], y[2] * y[4], y[2] * y[5], y[2] * y[6],
    y[3] * y[4] - y[6], y[3] * y[5], y[3] * y[6],
    y[4] * y[5], y[4] * y[6],
    y[5] * y[6],
  ]

  return pbw_algebra(R, rels)
end

@testset "QuantumGroups" begin
  function AbstractAlgebra.ConformanceTests.generate_element(A::PBWAlgebra)
    R = coefficient_ring(A)
    len = rand(0:8)
    coeffs = filter!(
      !iszero, [AbstractAlgebra.ConformanceTests.generate_element(R) for _ in 1:len]
    )
    exps = [rand(0:4) for _ in 1:ngens(A) for _ in 1:length(coeffs)]
    return PBWAlgebraElem(A, MPolyRingElem(R, coeffs, exps, length(coeffs)))
  end

  @testset "Conformance Tests" begin
    U = make_Usl3()
    AbstractAlgebra.ConformanceTests.test_NCRing_interface(U)
  end

  @testset "multplication" begin
    # coefficient
    # q^2*(F[3]*F[1]) = (q^2*F[3])*F[1]

    # multiplication table growth
    # table computation order: (2,1), ... (2,4); (3,1) ... (6,3)
    _ = F[7]^2 * F[1]^4
    @test leading_monomial(F[7]^6 * F[1]^3) == F[1]^3 * F[3]^6

    U = make_Usl3()
    F = gens(U)
    # table computation order: (2, 1) ... (6, 1) ... (6, 3)
    @test leading_monomial(F[1]^6 * F[3]^3) == F[1]^3 * F[3]^6
  end
end
