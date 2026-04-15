@testset "first local cohomology example" begin
    # get MonoidAlgebra
    kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ);
    x, y, z = gens(kQ);

    #example of computing local cohomology modules 
    I_M = ideal(kQ, [x^2*z, x^4*y])

    @test base_ring(I_M) == kQ

    m = ideal(kQ, [x, y, z]) #maximal ideal
    @test is_subset(I_M,m)

    M = Oscar.quotient_ring_as_module(I_M)

    #compute local cohomology
    H0 = Oscar.zeroth_local_cohomology(quotient_ring_as_module(I_M),m)
    @test !is_zero(H0)

    H1 = Oscar.local_cohomology(I_M,m,1)
    @test !is_zero(H1)
    H1_sectors = [h for h in H1.sectors if dim(h.H) > 0]
    @test all([dim(h.H) == 1 for h in H1_sectors])
    @test all([ambient_dim(h.sector) == 2 for h in H1.sectors])

    H2 = Oscar.local_cohomology(I_M,m,2)
    @test Oscar.is_zero(H2)

    H3 = Oscar.local_cohomology(I_M,m,3)
    @test Oscar.is_zero(H3)

    H4 = Oscar.local_cohomology(I_M,m,4)
    @test Oscar.is_zero(H4)
end
