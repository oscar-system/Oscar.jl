Oscar.example("PlaneCurve.jl")

@testset "AffinePlaneCurve constructors" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = y^3*x^6 - y^6*x^2
    C = AffinePlaneCurve(F)

    @test C.eq == F
    @test C.dimension == 1

    degree(C)
    @test C.degree == 9

    curve_components(C)
    @test C.components ==  Dict{AffinePlaneCurve{fmpq}, Int64}(AffinePlaneCurve(x) => 2, AffinePlaneCurve(y) => 3, AffinePlaneCurve(x^4 - y^3) => 1)

    @test C == AffinePlaneCurve(2*F)
end

@testset "AffinePlaneCurve reducible functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = AffinePlaneCurve((x^2+y^2))
    P = Point([QQ(0), QQ(0)])

    @test isirreducible(F)
    @test isreduced(F)
    @test reduction(F) == F

    G = AffinePlaneCurve(y^2)
    @test !isirreducible(G)
    @test !isreduced(G)
    @test reduction(G) == AffinePlaneCurve(y)

    H = AffinePlaneCurve(x*y)

    @test !isirreducible(H)
    @test isreduced(H)
    @test reduction(H) == H

    @test union(G, H) == AffinePlaneCurve(x*y^3)
end

@testset "AffinePlaneCurve intersection functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = AffinePlaneCurve(x*(x+y))
    G = AffinePlaneCurve(x+y^2+1)
    H = AffinePlaneCurve(x*(x+y)*y)
    M = AffinePlaneCurve((x-y)*(x-2))

    P = Point([QQ(0), QQ(0)])
    Q = Point([QQ(2), QQ(-2)])

    @test common_components(F, G) == []
    @test common_components(F, H) == [AffinePlaneCurve(x*(x+y))]

    @test curve_intersect(F, G) == [[], []]
    @test curve_intersect(F, H) == [[AffinePlaneCurve(x*(x+y))], []]
    @test curve_intersect(F, M) == [[], [P, Q]] || curve_intersect(F, M) == [[], [Q, P]]

    @test !aretransverse(F, H, P)
    @test !aretransverse(F, G, P)
    @test aretransverse(F, M, Q)
end

@testset "AffinePlaneCurve int_multiplicity functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])
    F = AffinePlaneCurve((x^2+y^2)*(x^2 + y^2 + 2*y))
    G = AffinePlaneCurve((x^2+y^2)*(y^3*x^6 - y^6*x^2))
    L = curve_intersect(F, G)
    P = Point([QQ(0), QQ(0)])
    Q = Point([QQ(0), QQ(-2)])

    @test L == [[AffinePlaneCurve(x^2+y^2)], [P, Q]] || L == [[AffinePlaneCurve(x^2+y^2)], [Q, P]]
    @test intersection_multiplicity(F, G, Q) == 2
    @test intersection_multiplicity(F, G, P) == -1
end

@testset "AffinePlaneCurve singularity functions" begin
    R, (x,y) = PolynomialRing(QQ, ["x", "y"])

    F = AffinePlaneCurve(x+y^2)

    @test curve_singular_locus(F) == [[], []]
    @test issmooth_curve(F)

    H = AffinePlaneCurve(x*y*(x+y))
    @test !issmooth_curve(H)

    G = AffinePlaneCurve(x^2*(x+y)*(y^3-x^2))
    S = curve_singular_locus(G)
    P1 = Point([QQ(0), QQ(0)])
    P2 = Point([QQ(-1), QQ(1)])
    P3 = Point([QQ(2), QQ(-2)])
    P4 = Point([QQ(1), QQ(2)])

    @test S == [[AffinePlaneCurve(x)], [P1, P2]] || S == [[AffinePlaneCurve(x)], [P2, P1]]
    @test !issmooth(G, P1)
    @test issmooth(G, P3)

    @test tangent(G, P3) == AffinePlaneCurve(x+y)
    @test tangent_lines(G, P1) == Dict{AffinePlaneCurve{fmpq}, Int64}(AffinePlaneCurve(x) => 4, AffinePlaneCurve(x+y) => 1)

    @test multiplicity(G, P1) == 5
    @test multiplicity(G, P3) == 1
    @test multiplicity(G, P4) == 0
end

@testset "ProjPlaneCurve constructors" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
	F = T(y^3*x^6 - y^6*x^2*z)
    C = ProjPlaneCurve(F)

    @test C.eq == F
    @test C.dimension == 1

    degree(C)
    @test C.degree == 9

    curve_components(C)
    @test C.components ==  Dict{ProjPlaneCurve{fmpq}, Int64}(ProjPlaneCurve(T(x)) => 2, ProjPlaneCurve(T(y)) => 3, ProjPlaneCurve(T(x^4 - y^3*z)) => 1)

    @test C == ProjPlaneCurve(2*F)
end

@testset "ProjPlaneCurve reducible functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
    F = ProjPlaneCurve(T(x^2+y^2))
    P = Point([QQ(0), QQ(0), QQ(1)])

    @test isirreducible(F)
    @test isreduced(F)
    @test reduction(F) == F

    G = ProjPlaneCurve(T(y^2))
    @test !isirreducible(G)
    @test !isreduced(G)
    @test reduction(G) == ProjPlaneCurve(T(y))

    H = ProjPlaneCurve(T(x*y))

    @test !isirreducible(H)
    @test isreduced(H)
    @test reduction(H) == H

    @test union(G, H) == ProjPlaneCurve(T(x*y^3))
end

@testset "ProjPlaneCurve intersection functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
    F = ProjPlaneCurve(T(x*(x+y)))
    G = ProjPlaneCurve(T(x*z+y^2+z^2))
    H = ProjPlaneCurve(T(x*(x+y)*y))
    M = ProjPlaneCurve((x-y)*(x-2*z))
	PP = projective_space(QQ, 2)

    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
	S = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
	Z = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(-1), QQ(0)])

    @test common_components(F, G) == []
    @test common_components(F, H) == [ProjPlaneCurve(T(x*(x+y)))]

    @test curve_intersect(PP[1], F, G) == [[], []]
    @test curve_intersect(PP[1], F, H) == [[ProjPlaneCurve(T(x*(x+y)))], []]
    @test curve_intersect(PP[1], ProjPlaneCurve(T(x+y+z)), ProjPlaneCurve(T(z))) == [[], [Z]]

	L = curve_intersect(PP[1], F, M)

	@test L[1] == []
	@test length(L[2]) == 3
	@test length(findall(x->x==P, L[2])) == 1
	@test length(findall(x->x==Q, L[2])) == 1
	@test length(findall(x->x==S, L[2])) == 1

    @test !aretransverse(F, H, P)
    @test !aretransverse(F, G, P)
    @test aretransverse(F, M, Q)
end

@testset "ProjPlaneCurve int_multiplicity functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
	F = ProjPlaneCurve(T((x^2+y^2)*(x^2 + y^2 + 2*y*z)))
    G = ProjPlaneCurve(T((x^2+y^2)*(y^3*x^6 - y^6*x^2*z)))
	PP = projective_space(QQ, 2)

    L = curve_intersect(PP[1], F, G)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-2), QQ(1)])

    @test L[1] == [ProjPlaneCurve(T(x^2+y^2))]
	@test length(L[2]) == 2
	@test length(findall(x->x==P, L[2])) == 1
	@test length(findall(x->x==Q, L[2])) == 1

    @test intersection_multiplicity(F, G, Q) == 2
    @test intersection_multiplicity(F, G, P) == -1
end

@testset "ProjPlaneCurve singularity functions" begin
    R, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    T = grade(R)
    PP = projective_space(QQ, 2)
    F = ProjPlaneCurve(T(x*z+y^2))
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(0), QQ(0)])

    @test curve_singular_locus(F) == [[], []]
    @test issmooth_curve(F)
    @test multiplicity(F, P) == 1

    H = ProjPlaneCurve(T(x*y*(x+y)))
    @test !issmooth_curve(H)

    G = ProjPlaneCurve(x^2*(x+y)*(y^3-x^2*z))
    S = curve_singular_locus(G)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(1)])
    P3 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
    P4 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(2), QQ(1)])

    @test S[1] == [ProjPlaneCurve(T(x))]
	@test length(S[2]) == 2
	@test length(findall(x->x==P1, S[2])) == 1
	@test length(findall(x->x==P2, S[2])) == 1
    @test !issmooth(G, P1)
    @test issmooth(G, P3)

    @test tangent(G, P3) == ProjPlaneCurve(T(x+y))
    @test tangent_lines(G, P1) == Dict{ProjPlaneCurve{fmpq}, Int64}(ProjPlaneCurve(T(x)) => 4, ProjPlaneCurve(T(x+y)) => 1)

    @test multiplicity(G, P1) == 5
    @test multiplicity(G, P3) == 1
    @test multiplicity(G, P4) == 0
end

@testset "AffineCurveDivisor basic functions" begin
	R, (x,y) = PolynomialRing(QQ, ["x", "y"])
	C = AffinePlaneCurve(y^2 + y + x^2)
	P = Point([QQ(0), QQ(0)])
	Q = Point([QQ(0), QQ(-1)])
	D = AffineCurveDivisor(C, Dict(P => 3, Q => -2))
	@test AffineCurveDivisor(C, P, -3) + D == AffineCurveDivisor(C, Q, -2)
	@test -2*D == AffineCurveDivisor(C, Dict(P => -6, Q => 4))
	@test !iseffective(D)
	@test iseffective(AffineCurveDivisor(C, P, 3))
	phi = y//x
	@test multiplicity(C, phi, P) == 1
	@test divisor(C, phi) == AffineCurveDivisor(C, Dict(P => 1, Q => -1))
end

@testset "ProjCurveDivisor basic functions" begin
	S, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
	T = grade(S)
	C = ProjPlaneCurve(T(y^2 + y*z + x^2))
	PP = projective_space(QQ, 2)
	P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
	Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-1), QQ(1)])
	D = ProjCurveDivisor(C, Dict(P => 3, Q => -2))
	@test ProjCurveDivisor(C, P, -3) + D == ProjCurveDivisor(C, Q, -2)
	@test -2*D == ProjCurveDivisor(C, Dict(P => -6, Q => 4))
	@test !iseffective(D)
	@test iseffective(ProjCurveDivisor(C, P, 3))
	F = T(x)
	@test multiplicity(C, F, P) == 1
	@test divisor(PP[1], C, F) == ProjCurveDivisor(C, Dict(P => 1, Q => 1))
end
