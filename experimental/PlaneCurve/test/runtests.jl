#Oscar.example("PlaneCurve.jl")

@testset "AffinePlaneCurve constructors" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = y^3 * x^6 - y^6 * x^2
    C = AffinePlaneCurve(F)

    @test Oscar.defining_equation(C) == F
    @test dim(C) == 1

    @test degree(C) == 9

    @test Oscar.curve_components(C) == Dict{AffinePlaneCurve{QQFieldElem},Int64}(
        AffinePlaneCurve(x) => 2,
        AffinePlaneCurve(y) => 3,
        AffinePlaneCurve(x^4 - y^3) => 1,
    )

    @test C == AffinePlaneCurve(2 * F)
end

@testset "AffinePlaneCurve reducible functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = AffinePlaneCurve((x^2 + y^2))
    P = Point([QQ(0), QQ(0)])

    @test Oscar.is_irreducible(F)
    @test Oscar.is_reduced(F)
    @test Oscar.reduction(F) == F

    G = AffinePlaneCurve(y^2)
    @test !Oscar.is_irreducible(G)
    @test !Oscar.is_reduced(G)
    @test Oscar.reduction(G) == AffinePlaneCurve(y)

    H = AffinePlaneCurve(x * y)

    @test !Oscar.is_irreducible(H)
    @test Oscar.is_reduced(H)
    @test Oscar.reduction(H) == H

    @test Oscar.union(G, H) == AffinePlaneCurve(x * y^3)
end

@testset "AffinePlaneCurve intersection functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = AffinePlaneCurve(x * (x + y))
    G = AffinePlaneCurve(x + y^2 + 1)
    H = AffinePlaneCurve(x * (x + y) * y)
    M = AffinePlaneCurve((x - y) * (x - 2))

    P = Point([QQ(0), QQ(0)])
    Q = Point([QQ(2), QQ(-2)])

    @test common_components(F, G) == []
    @test common_components(F, H) == [AffinePlaneCurve(x * (x + y))]

    @test curve_intersect(F, G) == [[], []]
    @test curve_intersect(F, H) == [[AffinePlaneCurve(x * (x + y))], []]
    @test curve_intersect(F, M) == [[], [P, Q]] ||
          curve_intersect(F, M) == [[], [Q, P]]

    @test !aretransverse(F, H, P)
    @test !aretransverse(F, G, P)
    @test aretransverse(F, M, Q)
end

@testset "AffinePlaneCurve int_multiplicity functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = AffinePlaneCurve((x^2 + y^2) * (x^2 + y^2 + 2 * y))
    G = AffinePlaneCurve((x^2 + y^2) * (y^3 * x^6 - y^6 * x^2))
    L = curve_intersect(F, G)
    P = Point([QQ(0), QQ(0)])
    Q = Point([QQ(0), QQ(-2)])

    @test L == [[AffinePlaneCurve(x^2 + y^2)], [P, Q]] ||
          L == [[AffinePlaneCurve(x^2 + y^2)], [Q, P]]
    @test intersection_multiplicity(F, G, Q) == 2
    @test intersection_multiplicity(F, G, P) == -1
end

@testset "AffinePlaneCurve singularity functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])

    F = AffinePlaneCurve(x + y^2)

    @test curve_singular_locus(F) == [[], []]
    @test is_smooth_curve(F)

    H = AffinePlaneCurve(x * y * (x + y))
    @test !is_smooth_curve(H)

    G = AffinePlaneCurve(x^2 * (x + y) * (y^3 - x^2))
    S = curve_singular_locus(G)
    P1 = Point([QQ(0), QQ(0)])
    P2 = Point([QQ(-1), QQ(1)])
    P3 = Point([QQ(2), QQ(-2)])
    P4 = Point([QQ(1), QQ(2)])

    @test S == [[AffinePlaneCurve(x)], [P1, P2]] ||
          S == [[AffinePlaneCurve(x)], [P2, P1]]
    @test !is_smooth(G, P1)
    @test is_smooth(G, P3)

    @test tangent(G, P3) == AffinePlaneCurve(x + y)
    @test tangent_lines(G, P1) == Dict{AffinePlaneCurve{QQFieldElem},Int64}(
        AffinePlaneCurve(x) => 4,
        AffinePlaneCurve(x + y) => 1,
    )

    @test multiplicity(G, P1) == 5
    @test multiplicity(G, P3) == 1
    @test multiplicity(G, P4) == 0
end

@testset "ProjPlaneCurve constructors" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = y^3 * x^6 - y^6 * x^2 * z
    C = ProjPlaneCurve(F)

    @test Oscar.defining_equation(C) == F.f
    @test dim(C) == 1

    @test degree(C) == 9

    @test Oscar.curve_components(C) == Dict{ProjPlaneCurve{QQFieldElem},Int64}(
        ProjPlaneCurve(x) => 2,
        ProjPlaneCurve(y) => 3,
        ProjPlaneCurve(x^4 - y^3 * z) => 1,
    )

    @test C == ProjPlaneCurve(2 * F)
end

@testset "ProjPlaneCurve reducible functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = ProjPlaneCurve(x^2 + y^2)
    P = Point([QQ(0), QQ(0), QQ(1)])

    @test Oscar.is_irreducible(F)
    @test Oscar.is_reduced(F)
    @test Oscar.reduction(F) == F

    G = ProjPlaneCurve(y^2)
    @test !Oscar.is_irreducible(G)
    @test !Oscar.is_reduced(G)
    @test Oscar.reduction(G) == ProjPlaneCurve(y)

    H = ProjPlaneCurve(x * y)

    @test !Oscar.is_irreducible(H)
    @test Oscar.is_reduced(H)
    @test Oscar.reduction(H) == H

    @test Oscar.union(G, H) == ProjPlaneCurve(x * y^3)
end

@testset "ProjPlaneCurve intersection functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = ProjPlaneCurve(x * (x + y))
    G = ProjPlaneCurve(x * z + y^2 + z^2)
    H = ProjPlaneCurve(x * (x + y) * y)
    M = ProjPlaneCurve((x - y) * (x - 2 * z))
    PP = proj_space(QQ, 2)

    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
    S = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
    Z = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(-1), QQ(0)])

    @test common_components(F, G) == []
    @test common_components(F, H) == [ProjPlaneCurve(x * (x + y))]

    @test curve_intersect(PP[1], F, G) == [[], []]
    @test curve_intersect(PP[1], F, H) == [[ProjPlaneCurve(x * (x + y))], []]
    @test curve_intersect(
        PP[1],
        ProjPlaneCurve(x + y + z),
        ProjPlaneCurve(z),
    ) == [[], [Z]]

    L = curve_intersect(PP[1], F, M)

    @test L[1] == []
    @test length(L[2]) == 3
    @test length(findall(x -> x == P, L[2])) == 1
    @test length(findall(x -> x == Q, L[2])) == 1
    @test length(findall(x -> x == S, L[2])) == 1

    @test !aretransverse(F, H, P)
    @test !aretransverse(F, G, P)
    @test aretransverse(F, M, Q)
end

@testset "ProjPlaneCurve int_multiplicity functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = ProjPlaneCurve((x^2 + y^2) * (x^2 + y^2 + 2 * y * z))
    G = ProjPlaneCurve((x^2 + y^2) * (y^3 * x^6 - y^6 * x^2 * z))
    PP = proj_space(QQ, 2)

    L = curve_intersect(PP[1], F, G)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-2), QQ(1)])

    @test L[1] == [ProjPlaneCurve(x^2 + y^2)]
    @test length(L[2]) == 2
    @test length(findall(x -> x == P, L[2])) == 1
    @test length(findall(x -> x == Q, L[2])) == 1

    @test intersection_multiplicity(F, G, Q) == 2
    @test intersection_multiplicity(F, G, P) == -1
end

@testset "ProjPlaneCurve singularity functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    PP = proj_space(QQ, 2)
    F = ProjPlaneCurve(x * z + y^2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(0), QQ(0)])

    @test curve_singular_locus(F) == [[], []]
    @test is_smooth_curve(F)
    @test multiplicity(F, P) == 1

    H = ProjPlaneCurve(x * y * (x + y))
    @test !is_smooth_curve(H)

    G = ProjPlaneCurve(x^2 * (x + y) * (y^3 - x^2 * z))
    S = curve_singular_locus(G)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(1)])
    P3 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
    P4 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(2), QQ(1)])

    @test S[1] == [ProjPlaneCurve(x)]
    @test length(S[2]) == 2
    @test length(findall(x -> x == P1, S[2])) == 1
    @test length(findall(x -> x == P2, S[2])) == 1
    @test !is_smooth(G, P1)
    @test is_smooth(G, P3)

    @test tangent(G, P3) == ProjPlaneCurve(x + y)
    @test tangent_lines(G, P1) == Dict{ProjPlaneCurve{QQFieldElem},Int64}(
        ProjPlaneCurve(x) => 4,
        ProjPlaneCurve(x + y) => 1,
    )

    @test multiplicity(G, P1) == 5
    @test multiplicity(G, P3) == 1
    @test multiplicity(G, P4) == 0
end

@testset "AffineCurveDivisor basic functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    C = AffinePlaneCurve(y^2 + y + x^2)
    P = Point([QQ(0), QQ(0)])
    Q = Point([QQ(0), QQ(-1)])
    D = Oscar.AffineCurveDivisor(C, Dict(P => 3, Q => -2))
    @test Oscar.AffineCurveDivisor(C, P, -3) + D == Oscar.AffineCurveDivisor(C, Q, -2)
    @test -2 * D == Oscar.AffineCurveDivisor(C, Dict(P => -6, Q => 4))
    @test !Oscar.is_effective(D)
    @test Oscar.is_effective(Oscar.AffineCurveDivisor(C, P, 3))
    phi = y // x
    @test multiplicity(C, phi, P) == 1
    @test Oscar.divisor(C, phi) == Oscar.AffineCurveDivisor(C, Dict(P => 1, Q => -1))
end

@testset "ProjCurveDivisor basic functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    C = ProjPlaneCurve(y^2 + y * z + x^2)
    PP = proj_space(QQ, 2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-1), QQ(1)])
    D = Oscar.ProjCurveDivisor(C, Dict(P => 3, Q => -2))
    @test Oscar.ProjCurveDivisor(C, P, -3) + D == Oscar.ProjCurveDivisor(C, Q, -2)
    @test -2 * D == Oscar.ProjCurveDivisor(C, Dict(P => -6, Q => 4))
    @test !Oscar.is_effective(D)
    @test Oscar.is_effective(Oscar.ProjCurveDivisor(C, P, 3))
    F = x
    phi = x // y
    @test multiplicity(C, F, P) == 1
    @test multiplicity(C, phi, P) == -1
    @test Oscar.divisor(PP[1], C, F) == Oscar.ProjCurveDivisor(C, Dict(P => 1, Q => 1))
    @test Oscar.divisor(PP[1], C, phi) == Oscar.ProjCurveDivisor(C, Dict(P => -1, Q => 1))
end

@testset "ProjCurveDivisor global sections" begin
    S, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    T, _ = grade(S)
    C = ProjPlaneCurve(y^2 * z - x * (x - z) * (x + 3 * z))
    PP = proj_space(QQ, 2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
    R = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    D = Oscar.ProjCurveDivisor(C, P, 4)
    L = Oscar.global_sections(D)
    @test length(L) == 4
    @test length(findall(a -> a == S(1) // S(1), L)) == 1
    @test length(findall(a -> a == S(y) // S(z), L)) == 1
    @test length(findall(a -> a == S(x) // S(z), L)) == 1
    @test length(findall(a -> a == S(x^2) // S(z^2), L)) == 1

    E = Oscar.ProjCurveDivisor(C, P)
    F = Oscar.ProjCurveDivisor(C, R)
    @test !Oscar.is_linearly_equivalent(E, F)
    @test Oscar.is_linearly_equivalent(2 * E, 2 * F)
    @test !Oscar.is_principal(E)

    G = 2 * E - 2 * F
    @test Oscar.is_principal(G)
    @test Oscar.is_linearly_equivalent(G, Oscar.divisor(C, Oscar.principal_divisor(G)))
end

@testset "Weierstrass form" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    PP = proj_space(QQ, 2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
    C = ProjPlaneCurve(y^2 * z - x^3 - x * z^2)
    D = ProjPlaneCurve(
            -x^3 - 3 * x^2 * y + 2 * x^2 * z - 3 * x * y^2 + 3 * x * y * z - 4 * x * z^2 - y^3 - y * z^2 + 6 * z^3,
    )
    @test Oscar.iselliptic(C)
    @test Oscar.toweierstrass(C, P) == y^2 * z - x^3 - x * z^2
    @test Oscar.toweierstrass(D, Q) ==
          y^2 * z + x * y * z + 3 * y * z^2 - x^3 - 2 * x^2 * z - 4 * x * z^2 - 6 * z^3
end

@testset "genus" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    C = ProjPlaneCurve(y^2 * z - x^3 - x * z^2)
    @test (Oscar.PlaneCurveModule.arithmetic_genus(C)) == ZZ(1)
    @test (geometric_genus(C)) == ZZ(1)
    R, (a, b) = polynomial_ring(GF(7), ["a", "b"])
    D = AffinePlaneCurve(b^9 - a^2 * (a - 1)^9)
    @test (Oscar.PlaneCurveModule.arithmetic_genus(D)) == ZZ(45)
    @test (geometric_genus(D)) == ZZ(0)
end

@testset "ProjEllipticCurve" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = y^2 * z - x^3 - x * z^2
    E = Oscar.ProjEllipticCurve(F)
    @test Oscar.discriminant(E) == -64
    @test Oscar.j_invariant(E) == 1728
end

@testset "Point addition" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = -x^3 - 3 * x^2 * y - 3 * x * y^2 - x * z^2 - y^3 + y^2 * z - y * z^2 - 4 * z^3
    PP = proj_space(QQ, 2)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
    E = Oscar.ProjEllipticCurve(F, P1)
    @test Oscar.weierstrass_form(E) == -x^3 - x * z^2 + y^2 * z - 4 * z^3
    Q1 = Oscar.Point_EllCurve(E, P1)
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-2), QQ(2), QQ(1)])
    Q2 = Oscar.Point_EllCurve(E, P2)
    @test Q1 + Q2 == Q2
    @test -Q2 == Oscar.Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)]))
end

@testset "Counting Points on Elliptic Curves" begin
    K = GF(5, 7)
    T, (x, y, z) = graded_polynomial_ring(K, ["x", "y", "z"])
    F = -x^3 - 3 * x^2 * y - 3 * x * y^2 - x * z^2 - y^3 + y^2 * z - y * z^2 - 4 * z^3
    PP = proj_space(K, 2)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [K(-1), K(1), K(0)])
    E = Oscar.ProjEllipticCurve(F, P1)
    @test order(E) == 78633
end

@testset "Torsion points on elliptic curves" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = -x^3 - 3 * x^2 * y - 3 * x * y^2 - x * z^2 - y^3 + y^2 * z - y * z^2 - 4 * z^3
    PP = proj_space(QQ, 2)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
    E = Oscar.ProjEllipticCurve(F, P1)
    Q1 = Oscar.Point_EllCurve(E, P1)
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-2), QQ(2), QQ(1)])
    Q2 = Oscar.Point_EllCurve(E, P2)
    @test Oscar.order(Q1) == 1
    @test Oscar.order(Q2) == 0
    @test is_torsion_point(Q1)
    @test !is_torsion_point(Q2)
    @test Oscar.torsion_points_lutz_nagell(E) == [Q1]
    @test Oscar.torsion_points_division_poly(E) == [Q1]
end

@testset "Elliptic Curve Method" begin
    n = Oscar.ECM(ZZ(4453))
    @test gcd(n, ZZ(4453)) == n

    A = residue_ring(ZZ, ZZ(4453))
    T, (x, y, z) = graded_polynomial_ring(A, ["x", "y", "z"])
    F = y^2 * z - x^3 - 10 * x * z^2 + 2 * z^3
    E = Oscar.ProjEllipticCurve(F)
    PP = proj_space(A, 2)
    P = Oscar.Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP[1], [A(1), A(3), A(1)]))
    Q = Oscar.Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP[1], [A(4332), A(3230), A(1)]))
    @test Oscar.sum_Point_EllCurveZnZ(P, P).Pt.v == Q.Pt.v
    @test Oscar.IntMult_Point_EllCurveZnZ(ZZ(2), P).Pt.v == Q.Pt.v
end

@testset "Primality Proving" begin
    n = ZZ(8051)
    m = Oscar.Pollard_rho(n)
    p = Oscar.Pollard_p_1(n)
    @test gcd(n, m) == m
    @test gcd(n, p) == p
end

@testset "ParaPlaneCurve" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    C1 = ProjPlaneCurve(1//2*x^5+x^2*y*z^2+x^3*y*z+1//2*x*y^2*z^2-2*x*y^3*z+y^5)
    I1 = Oscar.parametrization_plane_curve(C1)
    I2 = Oscar.adjoint_ideal(C1)
    R1 = parent(I1[1])
    s, t = gens(R1)
    @test I1 == [-4*s^4*t + 2*s^3*t^2, 2*s^2*t^3 - s*t^4, 4*s^5 - t^5]
    @test gens(I2) == [-x*y*z + y^3, -x^2*z + x*y^2, x^2*y + y^2*z, x^3 + x*y*z]
    C2 = ProjPlaneCurve(y^2 - x*z)
    P = Oscar.rational_point_conic(C2)
    I3 = Oscar.parametrization_conic(C2)
    R2 = parent(I3[1])
    s1, t1 = gens(R2)
    @test iszero(evaluate(C2.eq, P))
    @test I3 == [s1^2,  s1*t1, t1^2]
    C3 = ProjPlaneCurve(y^8-x^3*(z+x)^5)
    D = Oscar.map_to_rational_normal_curve(C3)
    I4 = Oscar.rat_normal_curve_anticanonical_map(D)
    Y = gens(parent(I4[1]))
    @test I4 == [Y[1], -Y[2], -Y[5], -Y[4], -Y[7]]
    C4 = Oscar.rat_normal_curve_It_Proj_Even(D)
    R3 = parent(C4[2].eq)
    U = gens(R3)
    @test C4[2] == ProjPlaneCurve(-U[1]*U[3] + U[2]^2)
    C5 = ProjPlaneCurve(-x^7-10*x^5*y^2-10*x^4*y^3-3*x^3*y^4+8*x^2*y^5+
    7*x*y^6+11*y^7+3*x^6*z+10*x^5*y*z+30*x^4*y^2*z+26*x^3*y^3*z-13*x^2*y^4*z-
    29*x*y^5*z-33*y^6*z-3*x^5*z^2-20*x^4*y*z^2-33*x^3*y^2*z^2-8*x^2*y^3*z^2+
    37*x*y^4*z^2+33*y^5*z^2+x^4*z^3+10*x^3*y*z^3+13*x^2*y^2*z^3-15*x*y^3*z^3-
    11*y^4*z^3)
    D2 = Oscar.map_to_rational_normal_curve(C5)
    I5 = Oscar.rat_normal_curve_It_Proj_Odd(D2)
    R4 = parent(I5[1])
    V = gens(R4)
    @test I5 == [121*V[3] + 77*V[4], -11*V[5] - 7*V[6]]
    C6 = ProjPlaneCurve(y^8 - x^3*(z+x)^5)
    I = Oscar.adjoint_ideal(C6)
    BM = Oscar.invert_birational_map(gens(I), C6)
    R5 = parent(BM["image"][1])
    W = gens(R5)
    @test BM["image"] == [-W[4]*W[7] + W[5]*W[6], -W[3]*W[7] + W[4]*W[6],
 -W[1]*W[7] + W[2]*W[6], -W[2]*W[7] + W[4]*W[5], -W[1]*W[7] + W[3]*W[5],
 W[1]*W[5] - W[7]^2, -W[1]*W[7] + W[4]^2, -W[1]*W[6] + W[3]*W[4],
 W[2]*W[4] - W[7]^2, W[1]*W[4] - W[6]*W[7], W[2]*W[3] - W[6]*W[7],
 W[1]*W[3] - W[6]^2, W[2]^2 - W[5]*W[7], W[1]*W[2] - W[4]*W[7],
 W[1]^2 - W[3]*W[7], W[1]*W[6]^2 - W[3]^2*W[7], -W[3]^3*W[7] + W[6]^4]
    @test BM["inverse"] == [-W[6]^2, -W[4]*W[7], -W[5]*W[7] + W[6]^2]
end

@testset "NonPlaneCurve" begin
    T, (x, y, z, t) = graded_polynomial_ring(QQ, ["x", "y", "z", "t"])
    I = ideal(T, [x^2, y^2*z, z^2])
    C = Oscar.ProjCurve(I)
    PP = proj_space(QQ, 3)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(2), QQ(0), QQ(5)])
    @test P in C
    @test Oscar.is_irreducible(C)
    J = Oscar.jacobi_ideal(C)
    L = gens(J)
    @test length(L) == 4
    @test length(findall(a -> a == 4*x*y*z, L)) == 1
    @test length(findall(a -> a == 2*x*y^2, L)) == 1
    @test length(findall(a -> a == 4*x*z, L)) == 1
    @test length(findall(a -> a == 4*y*z^2, L)) == 1
    C2 = Oscar.reduction(C)
    I2 = Oscar.defining_ideal(C2)
    @test I2 == ideal(T, [x, z])
end
