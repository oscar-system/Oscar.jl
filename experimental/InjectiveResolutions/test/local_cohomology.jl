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

@testset "local cohomology of module" begin
    # get MonoidAlgebra
    kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)
    x,y,z = gens(kQ)

    # define ideals over monoid algebra
    I = ideal(kQ, [x^2*z, x^4*y])
    J = ideal(kQ, [x^5*y, z^3])

    M_I = quotient_ring_as_module(I)
    M_J = quotient_ring_as_module(J)

    _M = direct_sum(M_I, M_J; task=:none)
    M,_ = sub(_M, [y*_M[1]+y*_M[2], x^2*_M[2]])

    m = ideal(kQ,[x,z])
    H1 = Oscar.local_cohomology(M,m,1)
    @test !is_zero(H1)
end