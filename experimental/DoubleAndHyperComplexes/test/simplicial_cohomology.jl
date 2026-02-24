@testset "simplicial_cohomology" begin
    K = torus()
    CZ = Oscar.SimplicialCoComplex(ZZ, K)

    @test typeof(CZ) <: Oscar.SimplicialCoComplex

    AKZ = SimplicialCohomologyRing(CZ)

    @test AKZ <: Oscar.SimplicialCohomologyRing

    @test base_ring(AKZ) == ZZ

    CF_23 = Oscar.SimplicialCoComplex(GF(23), K)

    AK_23 = SimplicialCohomologyRing(CZ)

    @test AK_23 <: Oscar.SimplicialCohomologyRing

    @test base_ring(AK_23) == GF(23)
end
