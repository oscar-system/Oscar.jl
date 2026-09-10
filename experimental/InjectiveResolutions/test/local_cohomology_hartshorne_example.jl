@testset "hartshorne example" begin
    # this is Hartshorne's example from Section 20.5 in 24HLC (24 hours of local cohohomology)
    kQ = monoid_algebra([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
    a, b, c, d = gens(kQ)

    # M = k[Q] (as a k[Q]-module)
    I_M = ideal(kQ, [])
    M = Oscar.quotient_ring_as_module(I_M)
    I = ideal(kQ, [a, b])

    # cohomological degree 0
    H0 = Oscar.zeroth_local_cohomology(M, I)
    @test is_zero(H0)

    # cohomological degrees 1 and 2 from one injective resolution
    H = local_cohomology_all(I_M, I, 2)
    H1_sectors = [h for h in sectors(H[1]) if !is_zero(h)] #sectors with non-zero local cohomology
    @test !is_zero(H[1])
    @test length(H1_sectors) == 1

    H2_sectors = [h for h in sectors(H[2]) if !is_zero(h)] #sectors with non-zero local cohomology
    @test !is_zero(H[2])
    @test length(H2_sectors) == 1
end
