module PlaneCurveTest

using ..Oscar
using ..Test
import Oscar.PlaneCurveModule
const PC = Oscar.PlaneCurveModule


@testset "AffinePlaneCurve constructors" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = y^3 * x^6 - y^6 * x^2
    C = PC.AffinePlaneCurve(F)

    @test PC.defining_equation(C) == F
    @test dim(C) == 1

    @test PC.degree(C) == 9

    @test PC.curve_components(C) == Dict{PC.AffinePlaneCurve{QQFieldElem},Int64}(
        PC.AffinePlaneCurve(x) => 2,
        PC.AffinePlaneCurve(y) => 3,
        PC.AffinePlaneCurve(x^4 - y^3) => 1,
    )

    @test C == PC.AffinePlaneCurve(2 * F)
end

@testset "AffinePlaneCurve reducible functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = PC.AffinePlaneCurve((x^2 + y^2))
    P = PC.Point([QQ(0), QQ(0)])

    @test Oscar.is_irreducible(F)
    @test Oscar.is_reduced(F)
    @test PC.reduction(F) == F

    G = PC.AffinePlaneCurve(y^2)
    @test !Oscar.is_irreducible(G)
    @test !Oscar.is_reduced(G)
    @test PC.reduction(G) == PC.AffinePlaneCurve(y)

    H = PC.AffinePlaneCurve(x * y)

    @test !Oscar.is_irreducible(H)
    @test Oscar.is_reduced(H)
    @test PC.reduction(H) == H

    @test Oscar.union(G, H) == PC.AffinePlaneCurve(x * y^3)
end

@testset "AffinePlaneCurve intersection functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = PC.AffinePlaneCurve(x * (x + y))
    G = PC.AffinePlaneCurve(x + y^2 + 1)
    H = PC.AffinePlaneCurve(x * (x + y) * y)
    M = PC.AffinePlaneCurve((x - y) * (x - 2))

    P = PC.Point([QQ(0), QQ(0)])
    Q = PC.Point([QQ(2), QQ(-2)])

    @test PC.common_components(F, G) == []
    @test PC.common_components(F, H) == [PC.AffinePlaneCurve(x * (x + y))]

    @test PC.curve_intersect(F, G) == [[], []]
    @test PC.curve_intersect(F, H) == [[PC.AffinePlaneCurve(x * (x + y))], []]
    @test PC.curve_intersect(F, M) == [[], [P, Q]] ||
          PC.curve_intersect(F, M) == [[], [Q, P]]

    @test !PC.aretransverse(F, H, P)
    @test !PC.aretransverse(F, G, P)
    @test PC.aretransverse(F, M, Q)
end

@testset "AffinePlaneCurve int_multiplicity functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    F = PC.AffinePlaneCurve((x^2 + y^2) * (x^2 + y^2 + 2 * y))
    G = PC.AffinePlaneCurve((x^2 + y^2) * (y^3 * x^6 - y^6 * x^2))
    L = PC.curve_intersect(F, G)
    P = PC.Point([QQ(0), QQ(0)])
    Q = PC.Point([QQ(0), QQ(-2)])

    @test L == [[PC.AffinePlaneCurve(x^2 + y^2)], [P, Q]] ||
          L == [[PC.AffinePlaneCurve(x^2 + y^2)], [Q, P]]
    @test PC.intersection_multiplicity(F, G, Q) == 2
    @test PC.intersection_multiplicity(F, G, P) == -1
end

@testset "AffinePlaneCurve singularity functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])

    F = PC.AffinePlaneCurve(x + y^2)

    @test PC.curve_singular_locus(F) == [[], []]
    @test PC.is_smooth_curve(F)

    H = PC.AffinePlaneCurve(x * y * (x + y))
    @test !PC.is_smooth_curve(H)

    G = PC.AffinePlaneCurve(x^2 * (x + y) * (y^3 - x^2))
    S = PC.curve_singular_locus(G)
    P1 = PC.Point([QQ(0), QQ(0)])
    P2 = PC.Point([QQ(-1), QQ(1)])
    P3 = PC.Point([QQ(2), QQ(-2)])
    P4 = PC.Point([QQ(1), QQ(2)])

    @test S == [[PC.AffinePlaneCurve(x)], [P1, P2]] ||
          S == [[PC.AffinePlaneCurve(x)], [P2, P1]]
    @test !is_smooth(G, P1)
    @test is_smooth(G, P3)

    @test PC.tangent(G, P3) == PC.AffinePlaneCurve(x + y)
    @test PC.tangent_lines(G, P1) == Dict{PC.AffinePlaneCurve{QQFieldElem},Int64}(
        PC.AffinePlaneCurve(x) => 4,
        PC.AffinePlaneCurve(x + y) => 1,
    )

    @test multiplicity(G, P1) == 5
    @test multiplicity(G, P3) == 1
    @test multiplicity(G, P4) == 0
end

@testset "ProjPlaneCurve constructors" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = y^3 * x^6 - y^6 * x^2 * z
    C = PC.ProjPlaneCurve(F)

    @test PC.defining_equation(C) == F.f
    @test dim(C) == 1

    @test degree(C) == 9

    @test PC.curve_components(C) == Dict{PC.ProjPlaneCurve{QQFieldElem},Int64}(
        PC.ProjPlaneCurve(x) => 2,
        PC.ProjPlaneCurve(y) => 3,
        PC.ProjPlaneCurve(x^4 - y^3 * z) => 1,
    )

    @test C == PC.ProjPlaneCurve(2 * F)
end

@testset "ProjPlaneCurve reducible functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = PC.ProjPlaneCurve(x^2 + y^2)
    P = PC.Point([QQ(0), QQ(0), QQ(1)])

    @test Oscar.is_irreducible(F)
    @test Oscar.is_reduced(F)
    @test PC.reduction(F) == F

    G = PC.ProjPlaneCurve(y^2)
    @test !Oscar.is_irreducible(G)
    @test !Oscar.is_reduced(G)
    @test PC.reduction(G) == PC.ProjPlaneCurve(y)

    H = PC.ProjPlaneCurve(x * y)

    @test !Oscar.is_irreducible(H)
    @test Oscar.is_reduced(H)
    @test PC.reduction(H) == H

    @test Oscar.union(G, H) == PC.ProjPlaneCurve(x * y^3)
end

@testset "ProjPlaneCurve intersection functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = PC.ProjPlaneCurve(x * (x + y))
    G = PC.ProjPlaneCurve(x * z + y^2 + z^2)
    H = PC.ProjPlaneCurve(x * (x + y) * y)
    M = PC.ProjPlaneCurve((x - y) * (x - 2 * z))
    PP = proj_space(QQ, 2)

    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
    S = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
    Z = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(-1), QQ(0)])

    @test PC.common_components(F, G) == []
    @test PC.common_components(F, H) == [PC.ProjPlaneCurve(x * (x + y))]

    @test PC.curve_intersect(PP[1], F, G) == [[], []]
    @test PC.curve_intersect(PP[1], F, H) == [[PC.ProjPlaneCurve(x * (x + y))], []]
    @test PC.curve_intersect(
        PP[1],
        PC.ProjPlaneCurve(x + y + z),
        PC.ProjPlaneCurve(z),
    ) == [[], [Z]]

    L = PC.curve_intersect(PP[1], F, M)

    @test L[1] == []
    @test length(L[2]) == 3
    @test length(findall(x -> x == P, L[2])) == 1
    @test length(findall(x -> x == Q, L[2])) == 1
    @test length(findall(x -> x == S, L[2])) == 1

    @test !PC.aretransverse(F, H, P)
    @test !PC.aretransverse(F, G, P)
    @test PC.aretransverse(F, M, Q)
end

@testset "ProjPlaneCurve int_multiplicity functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = PC.ProjPlaneCurve((x^2 + y^2) * (x^2 + y^2 + 2 * y * z))
    G = PC.ProjPlaneCurve((x^2 + y^2) * (y^3 * x^6 - y^6 * x^2 * z))
    PP = proj_space(QQ, 2)

    L = PC.curve_intersect(PP[1], F, G)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-2), QQ(1)])

    @test L[1] == [PC.ProjPlaneCurve(x^2 + y^2)]
    @test length(L[2]) == 2
    @test length(findall(x -> x == P, L[2])) == 1
    @test length(findall(x -> x == Q, L[2])) == 1

    @test PC.intersection_multiplicity(F, G, Q) == 2
    @test PC.intersection_multiplicity(F, G, P) == -1
end

@testset "ProjPlaneCurve singularity functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    PP = proj_space(QQ, 2)
    F = PC.ProjPlaneCurve(x * z + y^2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(0), QQ(0)])

    @test PC.curve_singular_locus(F) == [[], []]
    @test PC.is_smooth_curve(F)
    @test multiplicity(F, P) == 1

    H = PC.ProjPlaneCurve(x * y * (x + y))
    @test !PC.is_smooth_curve(H)

    G = PC.ProjPlaneCurve(x^2 * (x + y) * (y^3 - x^2 * z))
    S = PC.curve_singular_locus(G)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(1)])
    P3 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)])
    P4 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(1), QQ(2), QQ(1)])

    @test S[1] == [PC.ProjPlaneCurve(x)]
    @test length(S[2]) == 2
    @test length(findall(x -> x == P1, S[2])) == 1
    @test length(findall(x -> x == P2, S[2])) == 1
    @test !is_smooth(G, P1)
    @test is_smooth(G, P3)

    @test PC.tangent(G, P3) == PC.ProjPlaneCurve(x + y)
    @test PC.tangent_lines(G, P1) == Dict{PC.ProjPlaneCurve{QQFieldElem},Int64}(
        PC.ProjPlaneCurve(x) => 4,
        PC.ProjPlaneCurve(x + y) => 1,
    )

    @test multiplicity(G, P1) == 5
    @test multiplicity(G, P3) == 1
    @test multiplicity(G, P4) == 0
end

@testset "AffineCurveDivisor basic functions" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    C = PC.AffinePlaneCurve(y^2 + y + x^2)
    P = PC.Point([QQ(0), QQ(0)])
    Q = PC.Point([QQ(0), QQ(-1)])
    D = PC.AffineCurveDivisor(C, Dict(P => 3, Q => -2))
    @test PC.AffineCurveDivisor(C, P, -3) + D == PC.AffineCurveDivisor(C, Q, -2)
    @test -2 * D == PC.AffineCurveDivisor(C, Dict(P => -6, Q => 4))
    @test !Oscar.is_effective(D)
    @test Oscar.is_effective(PC.AffineCurveDivisor(C, P, 3))
    phi = y // x
    @test multiplicity(C, phi, P) == 1
    @test PC.divisor(C, phi) == PC.AffineCurveDivisor(C, Dict(P => 1, Q => -1))
end

@testset "ProjCurveDivisor basic functions" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    C = PC.ProjPlaneCurve(y^2 + y * z + x^2)
    PP = proj_space(QQ, 2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(-1), QQ(1)])
    D = PC.ProjCurveDivisor(C, Dict(P => 3, Q => -2))
    @test PC.ProjCurveDivisor(C, P, -3) + D == PC.ProjCurveDivisor(C, Q, -2)
    @test -2 * D == PC.ProjCurveDivisor(C, Dict(P => -6, Q => 4))
    @test !PC.is_effective(D)
    @test PC.is_effective(PC.ProjCurveDivisor(C, P, 3))
    F = x
    phi = x // y
    @test multiplicity(C, F, P) == 1
    @test multiplicity(C, phi, P) == -1
    @test PC.divisor(PP[1], C, F) == PC.ProjCurveDivisor(C, Dict(P => 1, Q => 1))
    @test PC.divisor(PP[1], C, phi) == PC.ProjCurveDivisor(C, Dict(P => -1, Q => 1))
end

@testset "ProjCurveDivisor global sections" begin
    S, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    T, _ = grade(S)
    C = PC.ProjPlaneCurve(y^2 * z - x * (x - z) * (x + 3 * z))
    PP = proj_space(QQ, 2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
    R = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(0), QQ(1)])
    D = PC.ProjCurveDivisor(C, P, 4)
    L = PC.global_sections(D)
    @test length(L) == 4
    @test length(findall(a -> a == S(1) // S(1), L)) == 1
    @test length(findall(a -> a == S(y) // S(z), L)) == 1
    @test length(findall(a -> a == S(x) // S(z), L)) == 1
    @test length(findall(a -> a == S(x^2) // S(z^2), L)) == 1

    E = PC.ProjCurveDivisor(C, P)
    F = PC.ProjCurveDivisor(C, R)
    @test !PC.is_linearly_equivalent(E, F)
    @test PC.is_linearly_equivalent(2 * E, 2 * F)
    @test !PC.is_principal(E)

    G = 2 * E - 2 * F
    @test PC.is_principal(G)
    @test PC.is_linearly_equivalent(G, PC.divisor(C, PC.principal_divisor(G)))
end

@testset "Weierstrass form" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    PP = proj_space(QQ, 2)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(1), QQ(0)])
    Q = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
    C = PC.ProjPlaneCurve(y^2 * z - x^3 - x * z^2)
    D = PC.ProjPlaneCurve(
            -x^3 - 3 * x^2 * y + 2 * x^2 * z - 3 * x * y^2 + 3 * x * y * z - 4 * x * z^2 - y^3 - y * z^2 + 6 * z^3,
    )
    @test PC.iselliptic(C)
    @test PC.toweierstrass(C, P) == y^2 * z - x^3 - x * z^2
    @test PC.toweierstrass(D, Q) ==
          y^2 * z + x * y * z + 3 * y * z^2 - x^3 - 2 * x^2 * z - 4 * x * z^2 - 6 * z^3
end

@testset "genus" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    C = PC.ProjPlaneCurve(y^2 * z - x^3 - x * z^2)
    @test (Oscar.PlaneCurveModule.arithmetic_genus(C)) == ZZ(1)
    @test (PC.geometric_genus(C)) == ZZ(1)
    R, (a, b) = polynomial_ring(GF(7), ["a", "b"])
    D = PC.AffinePlaneCurve(b^9 - a^2 * (a - 1)^9)
    @test (Oscar.PlaneCurveModule.arithmetic_genus(D)) == ZZ(45)
    @test (PC.geometric_genus(D)) == ZZ(0)
end

@testset "ProjEllipticCurve" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = y^2 * z - x^3 - x * z^2
    E = PC.ProjEllipticCurve(F)
    @test PC.discriminant(E) == -64
    @test PC.j_invariant(E) == 1728
end

@testset "Point addition" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = -x^3 - 3 * x^2 * y - 3 * x * y^2 - x * z^2 - y^3 + y^2 * z - y * z^2 - 4 * z^3
    PP = proj_space(QQ, 2)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
    E = PC.ProjEllipticCurve(F, P1)
    @test PC.weierstrass_form(E) == -x^3 - x * z^2 + y^2 * z - 4 * z^3
    Q1 = PC.Point_EllCurve(E, P1)
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-2), QQ(2), QQ(1)])
    Q2 = PC.Point_EllCurve(E, P2)
    @test Q1 + Q2 == Q2
    @test -Q2 == PC.Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP[1], [QQ(2), QQ(-2), QQ(1)]))
end

@testset "Counting Points on Elliptic Curves" begin
    K = GF(5, 7)
    T, (x, y, z) = graded_polynomial_ring(K, ["x", "y", "z"])
    F = -x^3 - 3 * x^2 * y - 3 * x * y^2 - x * z^2 - y^3 + y^2 * z - y * z^2 - 4 * z^3
    PP = proj_space(K, 2)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [K(-1), K(1), K(0)])
    E = PC.ProjEllipticCurve(F, P1)
    @test order(E) == 78633
end

@testset "Torsion points on elliptic curves" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    F = -x^3 - 3 * x^2 * y - 3 * x * y^2 - x * z^2 - y^3 + y^2 * z - y * z^2 - 4 * z^3
    PP = proj_space(QQ, 2)
    P1 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-1), QQ(1), QQ(0)])
    E = PC.ProjEllipticCurve(F, P1)
    Q1 = PC.Point_EllCurve(E, P1)
    P2 = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(-2), QQ(2), QQ(1)])
    Q2 = PC.Point_EllCurve(E, P2)
    @test PC.order(Q1) == 1
    @test PC.order(Q2) == 0
    @test is_torsion_point(Q1)
    @test !is_torsion_point(Q2)
    @test PC.torsion_points_lutz_nagell(E) == [Q1]
    @test PC.torsion_points_division_poly(E) == [Q1]
end

@testset "Elliptic Curve Method" begin
    n = PC.ECM(ZZ(4453))
    @test gcd(n, ZZ(4453)) == n

    A = residue_ring(ZZ, ZZ(4453))
    T, (x, y, z) = graded_polynomial_ring(A, ["x", "y", "z"])
    F = y^2 * z - x^3 - 10 * x * z^2 + 2 * z^3
    E = PC.ProjEllipticCurve(F)
    PP = proj_space(A, 2)
    P = PC.Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP[1], [A(1), A(3), A(1)]))
    Q = PC.Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP[1], [A(4332), A(3230), A(1)]))
    @test PC.sum_Point_EllCurveZnZ(P, P).Pt.v == Q.Pt.v
    @test PC.IntMult_Point_EllCurveZnZ(ZZ(2), P).Pt.v == Q.Pt.v
end

@testset "Primality Proving" begin
    n = ZZ(8051)
    m = PC.Pollard_rho(n)
    p = PC.Pollard_p_1(n)
    @test gcd(n, m) == m
    @test gcd(n, p) == p
end

@testset "ParaPlaneCurve" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"])
    C1 = PC.ProjPlaneCurve(1//2*x^5+x^2*y*z^2+x^3*y*z+1//2*x*y^2*z^2-2*x*y^3*z+y^5)
    I1 = PC.parametrization_plane_curve(C1)
    I2 = PC.adjoint_ideal(C1)
    R1 = parent(I1[1])
    s, t = gens(R1)
    @test I1 == [-4*s^4*t + 2*s^3*t^2, 2*s^2*t^3 - s*t^4, 4*s^5 - t^5]
    @test gens(I2) == [-x*y*z + y^3, -x^2*z + x*y^2, x^2*y + y^2*z, x^3 + x*y*z]
    C2 = PC.ProjPlaneCurve(y^2 - x*z)
    P = PC.rational_point_conic(C2)
    I3 = PC.parametrization_conic(C2)
    R2 = parent(I3[1])
    s1, t1 = gens(R2)
    @test iszero(evaluate(C2.eq, P))
    @test I3 == [s1^2,  s1*t1, t1^2]
    C3 = PC.ProjPlaneCurve(y^8-x^3*(z+x)^5)
    D = PC.map_to_rational_normal_curve(C3)
    I4 = PC.rat_normal_curve_anticanonical_map(D)
    Y = gens(parent(I4[1]))
    @test I4 == [Y[1], -Y[2], -Y[5], -Y[4], -Y[7]]
    C4 = PC.rat_normal_curve_It_Proj_Even(D)
    R3 = parent(C4[2].eq)
    U = gens(R3)
    @test C4[2] == PC.ProjPlaneCurve(-U[1]*U[3] + U[2]^2)
    C5 = PC.ProjPlaneCurve(-x^7-10*x^5*y^2-10*x^4*y^3-3*x^3*y^4+8*x^2*y^5+
    7*x*y^6+11*y^7+3*x^6*z+10*x^5*y*z+30*x^4*y^2*z+26*x^3*y^3*z-13*x^2*y^4*z-
    29*x*y^5*z-33*y^6*z-3*x^5*z^2-20*x^4*y*z^2-33*x^3*y^2*z^2-8*x^2*y^3*z^2+
    37*x*y^4*z^2+33*y^5*z^2+x^4*z^3+10*x^3*y*z^3+13*x^2*y^2*z^3-15*x*y^3*z^3-
    11*y^4*z^3)
    D2 = PC.map_to_rational_normal_curve(C5)
    I5 = PC.rat_normal_curve_It_Proj_Odd(D2)
    R4 = parent(I5[1])
    V = gens(R4)
    @test I5 == [121*V[3] + 77*V[4], -11*V[5] - 7*V[6]]
    C6 = PC.ProjPlaneCurve(y^8 - x^3*(z+x)^5)
    I = PC.adjoint_ideal(C6)
    BM = PC.invert_birational_map(gens(I), C6)
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
    C = PC.ProjCurve(I)
    PP = proj_space(QQ, 3)
    P = Oscar.Geometry.ProjSpcElem(PP[1], [QQ(0), QQ(2), QQ(0), QQ(5)])
    @test P in C
    @test PC.is_irreducible(C)
    J = PC.jacobi_ideal(C)
    L = gens(J)
    @test length(L) == 4
    @test length(findall(a -> a == 4*x*y*z, L)) == 1
    @test length(findall(a -> a == 2*x*y^2, L)) == 1
    @test length(findall(a -> a == 4*x*z, L)) == 1
    @test length(findall(a -> a == 4*y*z^2, L)) == 1
    C2 = PC.reduction(C)
    I2 = PC.defining_ideal(C2)
    @test I2 == ideal(T, [x, z])
end

end
