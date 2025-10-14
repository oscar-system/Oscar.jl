@testset "hartshorne example" begin
    # this is Hartshorne's example from Section 20.5 in 24HLC (24 hours of local cohohomology)
    kQ = monoid_algebra([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
    a, b, c, d = gens(kQ)

    # M = k[Q] (as a k[Q]-module)
    I_M = ideal(kQ, [])
    M = Oscar.quotient_ring_as_module(I_M)
    inj_res = Oscar.injective_resolution(M, 3)

    I = ideal(kQ, [a, b])

    # cohomological degree 0
    H0 = Oscar.zeroth_local_cohomology(quotient_ring_as_module(I_M), I)
    @test is_zero(H0)

    # cohomological degree 1
    H1 = Oscar.local_cohomology(I_M, I, 1)
    H1_sectors = [h for h in H1.sectors if dim(h.H)>0] #sectors with non-zero local cohomomology
    @test !Oscar.is_zero(H1)
    @test length(H1_sectors) == 1 

    # cohomological degree 2 
    H2 = Oscar.local_cohomology(I_M, I, 2)
    H2_sectors = [h for h in H2.sectors if dim(h.H)>0] #sectors with non-zero local cohomology
    @test !Oscar.is_zero(H2)
    @test length(H2_sectors) == 1

    #cohomological degree 3
    H3 = Oscar.local_cohomology(I_M, I, 3)
    @test Oscar.is_zero(H3)
end
