@testset "Normal toric varieties" begin
  ntv = normal_toric_variety(Oscar.normal_fan(Oscar.cube(2)))
  set_coordinate_names(ntv, ["x1", "x2", "y1", "y2"])
  ntv2 = normal_toric_variety(Oscar.cube(2))
  ntv3 = normal_toric_varieties_from_glsm(matrix(ZZ, [[1, 1, 1]]))
  ntv4 = normal_toric_varieties_from_star_triangulations(
    convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1])
  )

  @testset "Basic properties" begin
    @test is_complete(ntv) == true
    @test is_projective_space(ntv) == false
    @test torsion_free_rank(torusinvariant_cartier_divisor_group(ntv)) == 4
    @test torsion_free_rank(
      domain(
        map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(
          ntv
        ),
      ),
    ) == 4
    @test is_complete(ntv2) == true
    @test length(ntv3) == 1
    @test is_projective_space(ntv3[1]) == true
    @test length(ntv4) == 2
  end

  @testset "Speed test for Stanley-Reisner ideal (at most a few seconds)" begin
    success = false
    for i in 1:5
      ntv5 = normal_toric_variety(
        polarize(polyhedron(Polymake.polytope.rand_sphere(5, 60; seed=42)))
      )
      stats = @timed stanley_reisner_ideal(ntv5)
      duration = stats.time - stats.gctime
      if duration < 10
        success = true
        break
      else
        @warn "Stanley-Reisner ideal took $duration > 10 seconds (i=$i)"
      end
    end
    @test success == true
    ntv5 = normal_toric_variety(
      polarize(polyhedron(Polymake.polytope.rand_sphere(5, 60; seed=42)))
    )
    @test ngens(stanley_reisner_ideal(ntv5)) == 1648
  end

  p2 = projective_space(NormalToricVariety, 2)
  f2 = hirzebruch_surface(NormalToricVariety, 2)

  @testset "Equality of normal toric varieties" begin
    @test (p2 === f2) == false
    @test p2 === p2
    @test p2 != f2

    X = projective_space(NormalToricVariety, 2)
    X = domain(blow_up(X, [3, 4]))
    X = domain(blow_up(X, [-2, -3]))
    Y = weighted_projective_space(NormalToricVariety, [1, 2, 3])
    Y = domain(blow_up(Y, [-1, -1]))
    Y = domain(blow_up(Y, [3, 4]))
    @test X == Y

    Z = projective_space(NormalToricVariety, 2)
    X = domain(blow_up(Z, [1, 1]))
    Y = domain(blow_up(Z, [1, 2]))
    @test X != Y

    H = hirzebruch_surface(NormalToricVariety, 0)
    P1 = projective_space(NormalToricVariety, 1)
    ray_generators = [[1, 1], [1, 2]]
    max_cones = incidence_matrix([[1, 2]])
    X = normal_toric_variety(max_cones, ray_generators)
    @test length(Set([H, P1 * P1, X])) == 2

    @testset "Speed test hash (at most 0.5 seconds)" begin
      success = false
      ntv5 = normal_toric_variety(
        polarize(polyhedron(Polymake.polytope.rand_sphere(5, 60; seed=42)))
      )
      hash(ntv5)
      for i in 1:5
        stats = @timed hash(ntv5)
        duration = stats.time - stats.gctime
        if duration < 0.5
          success = true
          break
        else
          @warn "Hash took $duration > 0.5 seconds (i=$i)"
        end
      end
      @test success == true
    end
  end
end

@testset "Hodge and Betti numbers of toric varieties" begin
  X = domain(blow_up(projective_space(NormalToricVariety, 2), 1))
  @test betti_numbers(X) == ZZRingElem[1, 0, 2, 0, 1]
  @test hodge_numbers(X) == matrix(ZZ, [1 0 0; 0 2 0; 0 0 1])
  @test hodge_number(X, 1, 1) == 2
  @test hodge_number(X, 1, 0) == 0

  Y = weighted_projective_space(NormalToricVariety, [1, 1, 2])
  @test betti_numbers(Y) == ZZRingElem[1, 0, 1, 0, 1]
  @test hodge_numbers(Y) == matrix(ZZ, [1 0 0; 0 1 0; 0 0 1])
end

@testset "Hodge and Betti number exceptions" begin
  P2 = projective_space(NormalToricVariety, 2)
  for (p, q) in ((-1, 0), (3, 0), (0, -1), (0, 3))
    @test_throws ArgumentError hodge_number(P2, p, q)
  end

  noncomplete = affine_space(NormalToricVariety, 2)
  nonsimplicial = normal_toric_variety(normal_fan(cross_polytope(3)))
  for v in (noncomplete, nonsimplicial)
    @test_throws ArgumentError betti_numbers(v)
    @test_throws ArgumentError hodge_number(v, 0, 0)
    @test_throws ArgumentError hodge_numbers(v)
  end
end
