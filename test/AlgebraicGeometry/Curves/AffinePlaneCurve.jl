@testset "AffinePlaneCurves" begin
  @testset "constructors" begin
      R, (x, y) = polynomial_ring(QQ, [:x, :y])
      F = y^3 * x^6 - y^6 * x^2
      C = plane_curve(F)

      @test defining_equation(C) == prod(i[1] for i in factor(F))
      @test dim(C) == 1

      @test C == plane_curve(2 * F)
  end

  @testset "reducible functions" begin
      R, (x, y) = polynomial_ring(QQ, [:x, :y])
      F = plane_curve((x^2 + y^2))
      P = F([0, 0])

      @test is_irreducible(F)
      @test is_reduced(F)

      G = plane_curve(y^2)
      @test G == plane_curve(y)

      H = plane_curve(x * y)

      @test !is_irreducible(H)
      @test is_reduced(H)

      @test union(G, H) == plane_curve(x * y)
  end

  @testset "intersection functions" begin
      R, (x, y) = polynomial_ring(QQ, [:x, :y])
      F = plane_curve(x * (x + y))
      G = plane_curve(x + y^2 + 1)
      H = plane_curve(x * (x + y) * y)
      M = plane_curve((x - y) * (x - 2))

      P = H([0, 0])
      Q = H([2, -2])

      @test common_components(F, G) == []
      @test common_components(F, H) == [plane_curve(x * (x + y))]

      @test !is_transverse_intersection(F, H, P)
      @test !is_transverse_intersection(F, G, P)
      @test is_transverse_intersection(F, M, Q)
  end

  @testset "int_multiplicity functions" begin
      R, (x, y) = polynomial_ring(QQ, [:x, :y])
      F = plane_curve((x^2 + y^2) * (x^2 + y^2 + 2 * y))
      G = plane_curve((x^2 + y^2) * (y^3 * x^6 - y^6 * x^2))
      P = F([0, 0])
      Q = F([0, -2])

      @test intersection_multiplicity(F, G, Q) == 1
      @test_throws AbstractAlgebra.InfiniteDimensionError intersection_multiplicity(F, G, P)
  end

  @testset "singularity functions" begin
      R, (x, y) = polynomial_ring(QQ, [:x, :y])

      F = plane_curve(x + y^2)

      @test is_smooth(F)

      H = plane_curve(x * y * (x + y))
      @test !is_smooth(H)

      G = plane_curve(x * (x + y) * (y^3 - x^2))
      PP = ambient_space(G)
      P1 = G([0, 0])
      P2 = G([-1, 1])
      P3 = G([2, -2])
      P4 = PP([1, 2])

      @test !is_smooth(P1)
      @test is_smooth(P3)

      @test tangent_lines(G, P1) == Dict{typeof(G),Int64}(
          plane_curve(x) => 3,
          plane_curve(x + y) => 1,
      )

      @test multiplicity(G, P1) == 4
      @test multiplicity(G, P3) == 1
      @test multiplicity(G, P4) == 0
  end
end

