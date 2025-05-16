function make_Usl3()
  rels = []

  return pbw_algebra(rels)
end

@testset "QuantumGroups" begin
  @testset "PBWAlgebra" begin
    function AbstractAlgebra.ConformanceTests.generate_element(A::PBWAlgebra)
    end

    AbstractAlgebra.ConformanceTests.test_NCRing_interface(PBWAlgebra)
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
