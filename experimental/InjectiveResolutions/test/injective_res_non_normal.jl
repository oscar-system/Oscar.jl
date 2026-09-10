@testset "injective resolution of one-dimensional non-normal affine semigroup" begin
   # one-dimensional affine semigroup ring
    kQ = monoid_algebra([[2],[3]],QQ)
    @test !is_normal(kQ)
    N = quotient_ring_as_module(ideal(kQ,[]))
    inj_res = injective_resolution(N,3)
    @test is_exact(Q_graded_part(inj_res))
end

@testset "injective resolution of two-dimensional non-normal affine semigroup" begin
    # 2-dimensional non-normal affine semigroup ring
    kQ = monoid_algebra([[2,0],[3,0],[0,2],[0,3]],QQ)
    @test !is_normal(kQ)
    KQ = quotient_ring_as_module(ideal(kQ,[]))
    inj_res = injective_resolution(KQ,3)
    @test is_exact(Q_graded_part(inj_res))
end

@testset "injective resolution of 3-dimensional non-normal affine semigroup" begin
    # Example 4.5 in Katthän
    kQ = monoid_algebra([[0,0,1],[1,0,1],[0,2,1],[0,3,1],[1,3,1],[1,2,1]],QQ)
    KQ = quotient_ring_as_module(ideal(kQ,[]))

    #shift K[Q] by hand
    _KQ = twist(KQ,-grading_group(kQ)([2,6,4]))
    M = _KQ
    irr_res = irreducible_resolution(M, 2)
    @test is_exact(irr_res)

    # Example 13.27 in Miller-Sturmfels
    kQ = monoid_algebra([[0,0,1],[1,0,1],[3,0,1],[0,1,1],[1,1,1]],QQ)
    KQ = quotient_ring_as_module(ideal(kQ,[]))
    M = twist(KQ,-grading_group(kQ)([3,3,3]))
    irr_res = irreducible_resolution(M, 2)
    @test is_exact(irr_res)
end


@testset "injective resolution of module over non-normal monoid algebra" begin
    # quotient of the numerical semigroup ring k[t^2, t^3] by a monomial ideal
    kQ = monoid_algebra([[2], [3]], QQ)
    x, y = gens(kQ)
    I = ideal(kQ, [x^3, y^3])
    M = quotient_ring_as_module(I)
    @test !is_normal(kQ)

    irr_res = irreducible_resolution(M)
    @test is_exact(irr_res)

    inj_res = injective_resolution(M, 2)
    @test inj_res.upto <= 2
    @test is_exact(Q_graded_part(inj_res))
end