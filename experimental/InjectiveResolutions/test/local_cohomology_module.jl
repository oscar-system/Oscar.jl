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
    @test length(sectors(H1)) > 0
    @test sprint(show, H1) == "Sector partition of the 1st local cohomology module supported on ideal (x_1, x_3)"
    @test sprint(show, H1; context = :supercompact => true) == "Sector partition of a local cohomology module"
    # all local cohomology modules at once, and the ideal-argument method
    H = local_cohomology_all(M, m, 2)
    @test length(H) == 2
    @test !is_zero(H[1])
    @test length(sectors(H[1])) == length(sectors(H1))
    @test [dim(s.H) for s in sectors(H[1])] == [dim(s.H) for s in sectors(H1)]
    H_ideal = local_cohomology_all(I, m, 1)
    @test length(H_ideal) == 1
    @test_throws ArgumentError local_cohomology(M, ideal(monoid_algebra([[1, 0], [0, 1]], QQ), []), 1)
    @test_throws ArgumentError local_cohomology(M, m, 0)
    @test_throws ArgumentError local_cohomology_all(M, m, 0)
end
