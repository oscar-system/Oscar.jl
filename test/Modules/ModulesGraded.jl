using Random

RNG = Random.MersenneTwister(42)

"""
	randpoly(R::Ring,coeffs=0:9,max_exp=4,max_terms=8)

> Return a random Polynomial from the Polynomial Ring `R` with coefficients in `coeffs`
> with exponents between `0` and `max_exp` und between `0` and `max_terms` terms
"""
function randpoly(R::Oscar.Ring,coeffs=0:9,max_exp=4,max_terms=8)
	n = nvars(R)
	K = base_ring(R)
	E = [[Random.rand(RNG,0:max_exp) for i=1:n] for j=1:max_terms]
	C = [K(Random.rand(RNG,coeffs)) for i=1:max_terms]
	M = MPolyBuildCtx(R)
	for i=1:max_terms
		push_term!(M,C[i],E[i])
	end
	return finish(M)
end

@testset "Basic degree and homogeneity tests for free modules" begin
    Qx, x = PolynomialRing(QQ, 3)
    R, x = grade(Qx)
    F = FreeMod_dec(R, 3) #standard grading

    @test degree(F[1]) == decoration(R)[0]
    @test degree(x[2]*x[3]^2*F[1]+x[1]^3*F[3]) == 3*decoration(R)[1]
    @test !ishomogeneous((x[1]+x[2]^2)*F[1])
    @test !ishomogeneous(F[1]+x[3]*F[3])
    @test ishomogeneous(x[1]*x[2]*F[1]+x[3]^2*F[2])
    @test homogeneous_component(x[1]*F[1]+x[2]^2*F[2]+x[3]*F[3], decoration(R)[1]) == x[1]*F[1]+x[3]*F[3]
    D = homogeneous_components(x[1]*F[1]+x[2]^2*F[2]+x[3]*F[3])
    @test length(D) == 2
    @test D[2*decoration(R)[1]] == x[2]^2*F[2]

    g = abelian_group(2,0)
    Qx, x = PolynomialRing(QQ, 2)
    R, x = grade(Qx, [g[1], g[2]])
    F = FreeMod_dec(R, [g[2], g[1], g[1]-2*g[2]])

    @test degree(F[1]) == g[2]
    @test degree(x[1]*F[1]+x[2]^3*F[3]) == g[1]+g[2]
    @test !ishomogeneous((x[1]+x[2]^2)*F[1])
    @test !ishomogeneous(F[1]+x[2]*F[3])
    @test ishomogeneous(x[1]*x[2]*F[1]+x[2]^2*F[2])
    @test homogeneous_component(x[1]*F[1]+x[2]^2*F[2]+x[2]*F[3], g[1]+g[2]) == x[1]*F[1]
    D = homogeneous_components(x[1]^2*F[1]+x[2]^2*F[2]+x[1]*x[2]^3*F[3])
    @test length(D) == 2
    @test D[g[1]+2*g[2]] == x[2]^2*F[2]
    @test_throws ErrorException degree(x[1]*F[1]+x[2]^2*F[2]+x[2]*F[3])

    #TODO: tests for filtrated case
end