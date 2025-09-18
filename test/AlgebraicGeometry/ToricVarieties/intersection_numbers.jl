@testset "Topological intersection numbers" begin
  antv = affine_normal_toric_variety(Oscar.positive_hull([1 1; -1 1]))

  antv2 = normal_toric_variety(
    incidence_matrix([[1, 2, 3, 4]]), [[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0]]
  )

  v = normal_toric_variety(incidence_matrix([[1], [2], [3]]), [[1, 0], [0, 1], [-1, -1]])

  dP1 = del_pezzo_surface(NormalToricVariety, 1)
  c0 = cohomology_class(dP1, gens(cohomology_ring(dP1))[1])

  dP3 = del_pezzo_surface(NormalToricVariety, 3)
  (x1, x2, x3, e1, e2, e3) = gens(cohomology_ring(dP3))
  c1 = cohomology_class(dP3, x1)
  c2 = cohomology_class(dP3, e1)

  product_space =
    hirzebruch_surface(NormalToricVariety, 5) * projective_space(NormalToricVariety, 2)

  @testset "Should fail" begin
    R, _ = polynomial_ring(QQ, 3)
    @test_throws ArgumentError cohomology_ring(antv)
    @test_throws ArgumentError chow_ring(antv2)
    @test_throws ArgumentError c1 - c0
    @test_throws ArgumentError c1 + c0
    @test_throws ArgumentError ideal_of_linear_relations(R, dP3)
  end

  @testset "Chow ring and volume form of direct product space" begin
    @test ngens(ideal_of_linear_relations(product_space)) == 4
    @test ngens(chow_ring(product_space).I) == 7
    @test is_trivial(volume_form(product_space)) == false
  end

  @testset "Chow ring for non-complete but simplicial varieties" begin
    @test is_simplicial(v) == true
    @test is_complete(v) == false
    @test length(gens(chow_ring(v).I)) == 5
  end

  @testset "Properties, attributes and arithmetic of cohomology classes" begin
    @test is_trivial(c1) == false
    @test length(exponents(c1)) == 3
    @test length(coefficients(c1)) == 3
    @test QQFieldElem(3) * c1 == ZZRingElem(3) * c1
    @test 2 * c1 != ZZRingElem(3) * c2
    @test (c0 == c1) == false
  end

  @testset "Intersection numbers on dP3" begin
    @test integrate(cohomology_class(dP3, e1*e1)) == -1
    @test integrate(cohomology_class(dP3, e2*e2)) == -1
    @test integrate(cohomology_class(dP3, e3*e3)) == -1
    @test integrate(cohomology_class(dP3, x1*x1)) == -1
    @test integrate(cohomology_class(dP3, x2*x2)) == -1
    @test integrate(cohomology_class(dP3, x3*x3)) == -1
    @test integrate(c1) == 0
    @test integrate(c1^2+c1-3//4*c1*c1) == -1//4
    @test length(intersection_form(dP3)) == 21
  end
end
