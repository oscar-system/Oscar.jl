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
