function make_Usl3()
  R, y = polynomial_ring(QQ, :y => 1:6)
  rels = [
    y[1] * y[2] + y[4], y[1] * y[3], y[1] * y[4], y[1] * y[5] + y[6], y[1] * y[6],
    y[2] * y[3] + y[5], y[2] * y[4], y[2] * y[5], y[2] * y[6],
    y[3] * y[4] - y[6], y[3] * y[5], y[3] * y[6],
    y[4] * y[5], y[4] * y[6],
    y[5] * y[6],
  ]

  return Oscar.QuantumGroups.pbw_algebra(R, rels)
end

@testset "QuantumGroups.PBWAlgebra" begin
  @testset "Conformance Tests" begin
    U = make_Usl3()
    ConformanceTests.test_NCRing_interface(U)
  end

  @testset "multplication table growth" begin
    # this depends on test depends on the default size of the multiplication table
    # table computation order: (2,1), ... (2,4); (3,1) ... (6,3)
    U = make_Usl3()
    F = gens(U)
    _ = F[2]^2 * F[1]^4
    @test leading_monomial(F[2]^6 * F[1]^3) == F[1]^3 * F[2]^6

    # table computation order: (2, 1) ... (6, 1) ... (6, 3)
    U = make_Usl3()
    F = gens(U)
    @test leading_monomial(F[2]^6 * F[1]^3) == F[1]^3 * F[2]^6
  end
end
