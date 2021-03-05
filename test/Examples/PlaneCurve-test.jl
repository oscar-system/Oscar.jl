#Oscar.example("PlaneCurve.jl")

@testset "AffinePlaneCurve constructors" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = y^3*x^6 - y^6*x^2
    C = Oscar.AffinePlaneCurve(F)

    @test Oscar.defining_equation(C) == F
    @test dim(C) == 1

    @test degree(C) == 9

    @test Oscar.curve_components(C) == Dict{Oscar.AffinePlaneCurve{fmpq}, Int64}(Oscar.AffinePlaneCurve(x) => 2, Oscar.AffinePlaneCurve(y) => 3, Oscar.AffinePlaneCurve(x^4 - y^3) => 1)

    @test C == Oscar.AffinePlaneCurve(2*F)
end

@testset "AffinePlaneCurve reducible functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = Oscar.AffinePlaneCurve((x^2+y^2))
    P = Oscar.Point([QQ(0), QQ(0)])

    @test Oscar.isirreducible(F)
    @test Oscar.isreduced(F)
    @test Oscar.reduction(F) == F

    G = Oscar.AffinePlaneCurve(y^2)
    @test !Oscar.isirreducible(G)
    @test !Oscar.isreduced(G)
    @test Oscar.reduction(G) == Oscar.AffinePlaneCurve(y)

    H = Oscar.AffinePlaneCurve(x*y)

    @test !Oscar.isirreducible(H)
    @test Oscar.isreduced(H)
    @test Oscar.reduction(H) == H

    @test Oscar.union(G, H) == Oscar.AffinePlaneCurve(x*y^3)
end

@testset "AffinePlaneCurve intersection functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = Oscar.AffinePlaneCurve(x*(x+y))
    G = Oscar.AffinePlaneCurve(x+y^2+1)
    H = Oscar.AffinePlaneCurve(x*(x+y)*y)
    M = Oscar.AffinePlaneCurve((x-y)*(x-2))

    P = Oscar.Point([QQ(0), QQ(0)])
    Q = Oscar.Point([QQ(2), QQ(-2)])

    @test Oscar.common_components(F, G) == []
    @test Oscar.common_components(F, H) == [Oscar.AffinePlaneCurve(x*(x+y))]

    @test Oscar.curve_intersect(F, G) == [[], []]
    @test Oscar.curve_intersect(F, H) == [[Oscar.AffinePlaneCurve(x*(x+y))], []]
    @test Oscar.curve_intersect(F, M) == [[], [P, Q]] || Oscar.curve_intersect(F, M) == [[], [Q, P]]

    @test !Oscar.aretransverse(F, H, P)
    @test !Oscar.aretransverse(F, G, P)
    @test Oscar.aretransverse(F, M, Q)
end

@testset "AffinePlaneCurve int_multiplicity functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = Oscar.AffinePlaneCurve((x^2+y^2)*(x^2 + y^2 + 2*y))
    G = Oscar.AffinePlaneCurve((x^2+y^2)*(y^3*x^6 - y^6*x^2))
    L = Oscar.curve_intersect(F, G)
    P = Oscar.Point([QQ(0), QQ(0)])
    Q = Oscar.Point([QQ(0), QQ(-2)])

    @test L == [[Oscar.AffinePlaneCurve(x^2+y^2)], [P, Q]] || L == [[Oscar.AffinePlaneCurve(x^2+y^2)], [Q, P]]
    @test Oscar.intersection_multiplicity(F, G, Q) == 2
    @test Oscar.intersection_multiplicity(F, G, P) == -1
end

@testset "AffinePlaneCurve singularity functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])

    F = Oscar.AffinePlaneCurve(x+y^2)

    @test Oscar.curve_singular_locus(F) == [[], []]
    @test Oscar.issmooth_curve(F)

    H = Oscar.AffinePlaneCurve(x*y*(x+y))
    @test !Oscar.issmooth_curve(H)

    G = Oscar.AffinePlaneCurve(x^2*(x+y)*(y^3-x^2))
    S = Oscar.curve_singular_locus(G)
    P1 = Oscar.Point([QQ(0), QQ(0)])
    P2 = Oscar.Point([QQ(-1), QQ(1)])
    P3 = Oscar.Point([QQ(2), QQ(-2)])
    P4 = Oscar.Point([QQ(1), QQ(2)])

    @test S == [[Oscar.AffinePlaneCurve(x)], [P1, P2]] || S == [[Oscar.AffinePlaneCurve(x)], [P2, P1]]
    @test !Oscar.issmooth(G, P1)
    @test Oscar.issmooth(G, P3)

    @test Oscar.tangent(G, P3) == Oscar.AffinePlaneCurve(x+y)
    @test Oscar.tangent_lines(G, P1) == Dict{Oscar.AffinePlaneCurve{fmpq}, Int64}(Oscar.AffinePlaneCurve(x) => 4, Oscar.AffinePlaneCurve(x+y) => 1)

    @test Oscar.multiplicity(G, P1) == 5
    @test Oscar.multiplicity(G, P3) == 1
    @test Oscar.multiplicity(G, P4) == 0
end

@testset "ProjPlaneCurve constructors" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
	F = T(y^3*x^6 - y^6*x^2*z)
    C = Oscar.ProjPlaneCurve(F)

    @test Oscar.defining_equation(C) == F.f
    @test dim(C) == 1

    @test degree(C) == 9

    @test Oscar.curve_components(C) == Dict{Oscar.ProjPlaneCurve{fmpq}, Int64}(Oscar.ProjPlaneCurve(T(x)) => 2, Oscar.ProjPlaneCurve(T(y)) => 3, Oscar.ProjPlaneCurve(T(x^4 - y^3*z)) => 1)

    @test C == Oscar.ProjPlaneCurve(2*F)
end

@testset "ProjPlaneCurve reducible functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
    F = Oscar.ProjPlaneCurve(T(x^2+y^2))
    P = Oscar.Point([QQ(0), QQ(0), QQ(1)])

    @test Oscar.isirreducible(F)
    @test Oscar.isreduced(F)
    @test Oscar.reduction(F) == F

    G = Oscar.ProjPlaneCurve(T(y^2))
    @test !Oscar.isirreducible(G)
    @test !Oscar.isreduced(G)
    @test Oscar.reduction(G) == Oscar.ProjPlaneCurve(T(y))

    H = Oscar.ProjPlaneCurve(T(x*y))

    @test !Oscar.isirreducible(H)
    @test Oscar.isreduced(H)
    @test Oscar.reduction(H) == H

    @test Oscar.union(G, H) == Oscar.ProjPlaneCurve(T(x*y^3))
end

@testset "ProjPlaneCurve intersection functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
    F = Oscar.ProjPlaneCurve(T(x*(x+y)))
    G = Oscar.ProjPlaneCurve(T(x*z+y^2+z^2))
    H = Oscar.ProjPlaneCurve(T(x*(x+y)*y))
    M = Oscar.ProjPlaneCurve((x-y)*(x-2*z))
	PP = projective_space(QQ, 2)

    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
	S = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
	Z = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(-1), QQ(0)])

    @test Oscar.common_components(F, G) == []
    @test Oscar.common_components(F, H) == [Oscar.ProjPlaneCurve(T(x*(x+y)))]

    @test Oscar.curve_intersect(PP[1], F, G) == [[], []]
    @test Oscar.curve_intersect(PP[1], F, H) == [[Oscar.ProjPlaneCurve(T(x*(x+y)))], []]
    @test Oscar.curve_intersect(PP[1], Oscar.ProjPlaneCurve(T(x+y+z)), Oscar.ProjPlaneCurve(T(z))) == [[], [Z]]

	L = Oscar.curve_intersect(PP[1], F, M)

	@test L[1] == []
	@test length(L[2]) == 3
	@test length(findall(x->x==P, L[2])) == 1
	@test length(findall(x->x==Q, L[2])) == 1
	@test length(findall(x->x==S, L[2])) == 1

    @test !Oscar.aretransverse(F, H, P)
    @test !Oscar.aretransverse(F, G, P)
    @test Oscar.aretransverse(F, M, Q)
end

@testset "ProjPlaneCurve int_multiplicity functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
	F = Oscar.ProjPlaneCurve(T((x^2+y^2)*(x^2 + y^2 + 2*y*z)))
    G = Oscar.ProjPlaneCurve(T((x^2+y^2)*(y^3*x^6 - y^6*x^2*z)))
	PP = projective_space(QQ, 2)

    L = Oscar.curve_intersect(PP[1], F, G)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-2), QQ(1)])

    @test L[1] == [Oscar.ProjPlaneCurve(T(x^2+y^2))]
	@test length(L[2]) == 2
	@test length(findall(x->x==P, L[2])) == 1
	@test length(findall(x->x==Q, L[2])) == 1

    @test Oscar.intersection_multiplicity(F, G, Q) == 2
    @test Oscar.intersection_multiplicity(F, G, P) == -1
end

@testset "ProjPlaneCurve singularity functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
    PP = projective_space(QQ, 2)
    F = Oscar.ProjPlaneCurve(T(x*z+y^2))
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(0), QQ(0)])

    @test Oscar.curve_singular_locus(F) == [[], []]
    @test Oscar.issmooth_curve(F)
    @test Oscar.multiplicity(F, P) == 1

    H = Oscar.ProjPlaneCurve(T(x*y*(x+y)))
    @test !Oscar.issmooth_curve(H)

    G = Oscar.ProjPlaneCurve(x^2*(x+y)*(y^3-x^2*z))
    S = Oscar.curve_singular_locus(G)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(1)])
    P3 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
    P4 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(2), QQ(1)])

    @test S[1] == [Oscar.ProjPlaneCurve(T(x))]
	@test length(S[2]) == 2
	@test length(findall(x->x==P1, S[2])) == 1
	@test length(findall(x->x==P2, S[2])) == 1
    @test !Oscar.issmooth(G, P1)
    @test Oscar.issmooth(G, P3)

    @test Oscar.tangent(G, P3) == Oscar.ProjPlaneCurve(T(x+y))
    @test Oscar.tangent_lines(G, P1) == Dict{Oscar.ProjPlaneCurve{fmpq}, Int64}(Oscar.ProjPlaneCurve(T(x)) => 4, Oscar.ProjPlaneCurve(T(x+y)) => 1)

    @test Oscar.multiplicity(G, P1) == 5
    @test Oscar.multiplicity(G, P3) == 1
    @test Oscar.multiplicity(G, P4) == 0
end

@testset "AffineCurveDivisor basic functions" begin
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	C = Oscar.AffinePlaneCurve(y^2 + y + x^2)
	P = Oscar.Point([QQ(0), QQ(0)])
	Q = Oscar.Point([QQ(0), QQ(-1)])
	D = Oscar.AffineCurveDivisor(C, Dict(P => 3, Q => -2))
	@test Oscar.AffineCurveDivisor(C, P, -3) + D == Oscar.AffineCurveDivisor(C, Q, -2)
	@test -2*D == Oscar.AffineCurveDivisor(C, Dict(P => -6, Q => 4))
	@test !Oscar.iseffective(D)
	@test Oscar.iseffective(Oscar.AffineCurveDivisor(C, P, 3))
	phi = y//x
	@test Oscar.multiplicity(C, phi, P) == 1
	@test Oscar.divisor(C, phi) == Oscar.AffineCurveDivisor(C, Dict(P => 1, Q => -1))
end

@testset "ProjCurveDivisor basic functions" begin
	S, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	T = grade(S)
	C = Oscar.ProjPlaneCurve(T(y^2 + y*z + x^2))
	PP = projective_space(QQ, 2)
	P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
	Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-1), QQ(1)])
	D = Oscar.ProjCurveDivisor(C, Dict(P => 3, Q => -2))
	@test Oscar.ProjCurveDivisor(C, P, -3) + D == Oscar.ProjCurveDivisor(C, Q, -2)
	@test -2*D == Oscar.ProjCurveDivisor(C, Dict(P => -6, Q => 4))
	@test !Oscar.iseffective(D)
	@test Oscar.iseffective(Oscar.ProjCurveDivisor(C, P, 3))
	F = T(x)
    phi = T(x)//T(y)
	@test Oscar.multiplicity(C, F, P) == 1
	@test Oscar.multiplicity(C, phi, P) == -1
	@test Oscar.divisor(PP[1], C, F) == Oscar.ProjCurveDivisor(C, Dict(P => 1, Q => 1))
	@test Oscar.divisor(PP[1], C, phi) == Oscar.ProjCurveDivisor(C, Dict(P => -1, Q => 1))
end

@testset "ProjCurveDivisor global sections" begin
	S, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
	T = grade(S)
	C = Oscar.ProjPlaneCurve(T(y^2*z - x*(x-z)*(x+3*z)))
	PP = projective_space(QQ, 2)
	P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
	R = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
	D = Oscar.ProjCurveDivisor(C, P, 4)
	L = Oscar.global_sections(D)
	@test length(L) == 4
	@test length(findall(a->a==S(1)//S(1), L)) == 1
	@test length(findall(a->a==S(y)//S(z), L)) == 1
	@test length(findall(a->a==S(x)//S(z), L)) == 1
	@test length(findall(a->a==S(x^2)//S(z^2), L)) == 1

	E = Oscar.ProjCurveDivisor(C, P)
	F = Oscar.ProjCurveDivisor(C, R)
	@test !Oscar.islinearly_equivalent(E, F)
	@test Oscar.islinearly_equivalent(2*E, 2*F)
	@test !Oscar.isprincipal(E)

	G = 2*E - 2*F
	@test Oscar.isprincipal(G)
	@test Oscar.islinearly_equivalent(G, Oscar.divisor(C, Oscar.principal_divisor(G)))
end
