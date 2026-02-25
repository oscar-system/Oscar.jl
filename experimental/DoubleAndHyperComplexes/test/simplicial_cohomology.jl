@testset "simplicial_cohomology" begin
    K = torus()
    CZ = Oscar.SimplicialCoComplex(ZZ, K)

    @test typeof(CZ) <: Oscar.SimplicialCoComplex

    AKZ = Oscar.SimplicialCohomologyRing(CZ)

    @test AKZ <: Oscar.SimplicialCohomologyRing

    @test base_ring(AKZ) == ZZ

    CF_23 = Oscar.SimplicialCoComplex(GF(23), K)

    AK_23 = Oscar.SimplicialCohomologyRing(CZ)

    @test AK_23 <: Oscar.SimplicialCohomologyRing

    @test base_ring(AK_23) == GF(23)

    K = simplicial_complex([[1,2,3],[2,3,4]])
    Cmplx_p = Oscar.SimplicialCoComplex(padic_field(7, precision = 30),K)
    R_p = Oscar.SimplicialCohomologyRing(Cmplx_p)
    @test is_exact_type(typeof(one(R_p))) == false
    @test is_domain_type(typeof(one(R_p))) == false
    
    Cmplx = Oscar.SimplicialCoComplex(ZZ,K)
    R = Oscar.SimplicialCohomologyRing(Cmplx)
    @test is_exact_type(typeof(one(R))) == true
    @test is_domain_type(typeof(one(R))) == false

end
