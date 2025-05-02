@testset "ProjectivePlaneCurve" begin
  @testset "constructors" begin
      T, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
      F = y^3 * x^6 - y^6 * x^2 * z
      C = plane_curve(F)

      @test defining_equation(C) == y * x^5 - y^4 * x * z
      @test dim(C) == 1

      @test degree(C) == 6

      @test C == plane_curve(2 * F)
  end

  @testset "reducible functions" begin
      T, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
      F = plane_curve(x^2 + y^2)
      P = F([0, 0, 1])

      @test is_irreducible(F)
      @test is_reduced(F)

      G = plane_curve(y^2)
      @test is_irreducible(G)
      @test is_reduced(G)

      H = plane_curve(x * y)

      @test !is_irreducible(H)
      @test is_reduced(H)

      @test union(G, H) == plane_curve(x * y)
  end

  @testset "intersection functions" begin
      PP = projective_space(QQ,[:x,:y,:z])
      (x, y, z) = homogeneous_coordinates(PP)
      F = ProjectivePlaneCurve(x * (x + y))
      G = ProjectivePlaneCurve(x * z + y^2 + z^2)
      H = ProjectivePlaneCurve(x * (x + y) * y)
      M = ProjectivePlaneCurve((x - y) * (x - 2 * z))

      P = PP( [0, 0, 1])
      Q = PP( [2, -2, 1])
      S = PP( [0, 1, 0])
      Z = PP([1, -1, 0])

      @test common_components(F, G) == []
      @test common_components(F, H) == [ProjectivePlaneCurve(x * (x + y))]

      @test !is_transverse_intersection(F, H, P)
      @test !is_transverse_intersection(F, G, P)
      @test is_transverse_intersection(F, M, Q)
  end

  @testset "int_multiplicity functions" begin
      T, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
      F = ProjectivePlaneCurve((x^2 + y^2) * (x^2 + y^2 + 2 * y * z))
      G = ProjectivePlaneCurve((x^2 + y^2) * (y^3 * x^6 - y^6 * x^2 * z))
      PP = ambient_space(F)

      P = PP([0, 0, 1])
      Q = PP([0, -2, 1])

      @test intersection_multiplicity(F, G, Q) == 1
      @test_broken intersection_multiplicity(F, G, P) == inf
  end

  @testset "singularity functions" begin
      T, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
      F = ProjectivePlaneCurve(x * z + y^2)
      PP = ambient_space(F)
      P = PP([1, 0, 0])

      @test is_smooth(F)
      @test multiplicity(F, P) == 1

      H = ProjectivePlaneCurve(x * y * (x + y))
      @test !is_smooth(H)

      G = ProjectivePlaneCurve(x * (x + y) * (y^3 - x^2 * z))
      P1 = PP( [0, 0, 1])
      P2 = PP( [-1, 1, 1])
      P3 = PP( [2, -2, 1])
      P4 = PP( [1, 2, 1])

      @test is_smooth(G(P3))

      T = tangent_lines(G, P1)
      @test length(T) == 2
      @test T[ProjectivePlaneCurve(x)] == 3
      @test T[ProjectivePlaneCurve(x + y)] == 1

      @test multiplicity(G, P1) == 4
      @test multiplicity(G, P3) == 1
      @test multiplicity(G, P4) == 0
  end

  @testset "genus" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    C = plane_curve(y^2 * z - x^3 - x * z^2)
    @test arithmetic_genus(C) == 1
    @test geometric_genus(C) == 1
    R, (a, b) = polynomial_ring(GF(7), [:a, :b])
    D = plane_curve(b^9 - a^2 * (a - 1)^9)
    @test geometric_genus(D) == 0
    @test arithmetic_genus(projective_closure(D)) == 45
  end
end



@testset "ParametrizePlaneCurve" begin
    T, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    C1 = ProjectivePlaneCurve(1//2*x^5+x^2*y*z^2+x^3*y*z+1//2*x*y^2*z^2-2*x*y^3*z+y^5)
    I1 = parametrization(C1)
    I2 = Oscar.adjoint_ideal(C1)
    R1 = parent(I1[1])
    s, t = gens(R1)
    @test I1 == [-4*s^4*t + 2*s^3*t^2, 2*s^2*t^3 - s*t^4, 4*s^5 - t^5]
    @test gens(I2) == [-x*y*z + y^3, -x^2*z + x*y^2, x^2*y + y^2*z, x^3 + x*y*z]
    C2 = ProjectivePlaneCurve(y^2 - x*z)
    P = Oscar.rational_point_conic(C2)
    I3 = Oscar.parametrization_conic(C2)
    R2 = parent(I3[1])
    s1, t1 = gens(R2)
    @test iszero(evaluate(defining_equation(C2), P))
    @test I3 == [s1^2,  s1*t1, t1^2]
    C3 = ProjectivePlaneCurve(y^8-x^3*(z+x)^5)
    D = Oscar.map_to_rational_normal_curve(C3)
    I4 = Oscar.rat_normal_curve_anticanonical_map(D)
    Y = gens(parent(I4[1]))
    @test I4 == [Y[1], -Y[2], -Y[5], -Y[4], -Y[7]]
    C4 = Oscar.rat_normal_curve_It_Proj_Even(D)
    R3 = parent(defining_equation(C4[2]))
    U = gens(R3)
    @test C4[2] == ProjectivePlaneCurve(-U[1]*U[3] + U[2]^2)
    C5 = ProjectivePlaneCurve(-x^7-10*x^5*y^2-10*x^4*y^3-3*x^3*y^4+8*x^2*y^5+
    7*x*y^6+11*y^7+3*x^6*z+10*x^5*y*z+30*x^4*y^2*z+26*x^3*y^3*z-13*x^2*y^4*z-
    29*x*y^5*z-33*y^6*z-3*x^5*z^2-20*x^4*y*z^2-33*x^3*y^2*z^2-8*x^2*y^3*z^2+
    37*x*y^4*z^2+33*y^5*z^2+x^4*z^3+10*x^3*y*z^3+13*x^2*y^2*z^3-15*x*y^3*z^3-
    11*y^4*z^3)
    D2 = Oscar.map_to_rational_normal_curve(C5)
    I5 = Oscar.rat_normal_curve_It_Proj_Odd(D2)
    R4 = parent(I5[1])
    V = gens(R4)
    @test I5 == [121*V[3] + 77*V[4], -11*V[5] - 7*V[6]]
    C6 = ProjectivePlaneCurve(y^8 - x^3*(z+x)^5)
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

