@testset "injective resolution over k[x,y]" begin
    # get MonoidAlgebra
    kQ = monoid_algebra([[1, 0], [0, 1]],QQ)
    x,y = gens(kQ)

    # define ideal over monoid algebra
    I = ideal(kQ, [x^4, x^2*y^2, y^4])

    #irreducible decomposition 
    W = irreducible_decomposition(I)
    @test I == intersect(W)

    # irreducible resolution of M = kQ/I
    M = quotient_ring_as_module(I)
    irr_res = irreducible_resolution(M)
    @test is_exact(irr_res)

    # compute injective resolution up to cohomological degree 2
    inj_res = injective_resolution(I, 2)
    @test inj_res.upto <= 2

    ## irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
    inj_res_Q = q_graded_part(inj_res)
    @test is_exact(inj_res_Q)

    # getters, monomial matrices, minimality, injective hull
    @test length(injective_modules(inj_res)) == inj_res.upto + 1
    @test length(cochain_maps(inj_res)) == inj_res.upto + 1
    mm = monomial_matrix(0, inj_res)
    @test nrows(mm.matrix) == length(indecomposable_injectives(injective_modules(inj_res)[1]))
    @test ncols(mm.matrix) == length(indecomposable_injectives(injective_modules(inj_res)[2]))
    @test is_minimal(inj_res)
    @test degree_shift(inj_res) isa Vector{Int}
    @test monoid_algebra(injective_modules(inj_res)[1]) == kQ
    E, lambda = injective_hull(M)
    @test length(indecomposable_injectives(E)) == 2
    @test irreducible_sums(irr_res) == irr_res.irr_sums
    @test_throws ArgumentError injective_resolution(M, 1; shift=:foo)

    # irreducible hull and Bass numbers of M = k[x,y]/(x^4, x^2y^2, y^4):
    # the socle lives in degrees [1,3] and [3,1]
    W, lambda = irreducible_hull(M)
    @test length(indecomposable_injectives(W)) == 2
    @test sort(degrees_of_bass_numbers(M, 0)) == [[1, 3], [3, 1]]
    p = faces(kQ)[1]
    @test is_zero(p.prime) == false
    bass = graded_bass_numbers(M, p, 1)
    @test bass[1] == Dict([1, 3] => 1, [3, 1] => 1)
    @test bass[2] == Dict([-1, 3] => 1, [1, 1] => 1, [3, -1] => 1)

    # wrapping an ideal of the underlying algebra
    J = monoid_algebra_ideal(kQ, ideal(kQ.algebra, [kQ.algebra[1]^2]))
    @test J isa MonoidAlgebraIdeal
    @test ngens(J) == 1
end

@testset "injective resolutions over k[Q] = k[x,y,z]/(xz - y^2)" begin
    # get MonoidAlgebra
    kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)
    x,y,z = gens(kQ)


    ##first example
    # define ideal over monoid algebra
    I = ideal(kQ, [x^2*z, x^4*y])
    J = ideal(kQ, [x^5*y, z^3])

    M_I = quotient_ring_as_module(I)
    M_J = quotient_ring_as_module(J)

    _M = direct_sum(M_I, M_J; task=:none)
    M,_ = sub(_M, [y*_M[1]+y*_M[2], x^2*_M[2]])

    # compute irreducible resolution
    irr_res = irreducible_resolution(M)
    @test is_exact(irr_res)

    # minimal injective resolution of kQ/I up to cohomological degree 3
    inj_res = injective_resolution(M, 3)
    @test inj_res.upto <= 3

    # irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
    irr_res_Q = q_graded_part(inj_res)
    @test is_exact(irr_res_Q)



    ##second example
    I = ideal(kQ, [y^3, x^3 * z])

    #irreducible decomposition
    W = irreducible_decomposition(I)
    @test I == intersect(W...)

    # compute irreducible resolution of M = kQ/I
    M = quotient_ring_as_module(I)
    irr_res = irreducible_resolution(M)

    # compute injective resolution up to cohomological degree 3
    inj_res = injective_resolution(I, 3)
    @test inj_res.upto <= 3

    # irreducible resolution that is the Q-graded part of the minimal injective resolution above
    inj_res_Q = q_graded_part(inj_res)
    @test is_exact(inj_res_Q)



    ##third example
    I = ideal(kQ, [x^2*z])

    #irreducible decomposition
    W = irreducible_decomposition(I)
    @test I == intersect(W)

    # compute irreducible resolution of M = kQ/I
    M = quotient_ring_as_module(I)
    irr_res = irreducible_resolution(M)
    @test is_exact(irr_res)

    # compute injective resolution of kQ/I up to cohomological degree 3
    inj_res = injective_resolution(I, 3)
    @test inj_res.upto <= 3

    # get irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
    irr_res_Q = q_graded_part(inj_res)
    @test is_exact(irr_res_Q)



    ##fourth example 
    I = ideal(kQ, [z^2, z*y, x^4])

    #irreducible decomposition
    W = irreducible_decomposition(I)
    @test I == intersect(W...)

    # irreducible resolution of M = R/I
    M = quotient_ring_as_module(I)
    irr_res = irreducible_resolution(M)
    @test is_exact(irr_res)

    # injective resolution of M up to cohomological degree 3
    inj_res = injective_resolution(I, 3)
    @test inj_res.upto <= 3

    irr_res_Q = q_graded_part(inj_res)
    @test is_exact(irr_res_Q)
end

@testset "injective resolution over ZZ^3-graded MonoidAlgebra" begin
    # k[Q] for the cone over a square (Segre embedding of P^1 x P^1)
    kQ = monoid_algebra([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
    a, b, c, d = gens(kQ)

    # first example
    I = ideal(kQ, [a^2*b, c^2])
    W = irreducible_decomposition(I)
    @test I == intersect(W)

    inj_res = injective_resolution(I, 1)
    @test inj_res.upto <= 1
    @test is_exact(q_graded_part(inj_res))

    # second example
    I = ideal(kQ, [a^2*b, c^2, d*a^4])
    W = irreducible_decomposition(I)
    @test I == intersect(W)

    inj_res = injective_resolution(I, 1)
    @test inj_res.upto <= 1
    @test is_exact(q_graded_part(inj_res))
end
