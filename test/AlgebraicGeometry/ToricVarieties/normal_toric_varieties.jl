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

# Test lazy Hodge-number caching and complete matrix construction.
@testset "Hodge numbers of projective space P3" begin
  P3 = projective_space(NormalToricVariety, 3)
  @test betti_number(P3, 2) == 1
  @test betti_numbers(P3) == ZZRingElem[1, 0, 1, 0, 1, 0, 1]
  @test !has_attribute(P3, :hodge_numbers)
  returned_hodge_number = hodge_number(P3, 1, 1)
  @test returned_hodge_number == 1
  cached_H = get_attribute(P3, :hodge_numbers)
  @test isassigned(cached_H, 2, 2)
  @test !isassigned(cached_H, 1, 1)
  @test_throws UndefRefError cached_H[1, 1]
  @test all(
    p == q || isassigned(cached_H, p, q) for p in axes(cached_H, 1), q in axes(cached_H, 2)
  )
  @test cached_H[1, 2] == cached_H[2, 1] == 0
  zero!(returned_hodge_number)
  @test hodge_number(P3, 1, 1) == 1
  @test hodge_number(P3, 1, 0) == 0
  @test isassigned(cached_H, 2, 1)
  @test hodge_numbers(P3) == matrix(ZZ, [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
  @test all(i -> isassigned(cached_H, i), eachindex(cached_H))
  returned_hodge_numbers = hodge_numbers(P3)
  @test returned_hodge_numbers isa ZZMatrix
  zero!(returned_hodge_numbers[2, 2])
  @test hodge_number(P3, 1, 1) == 1
  @test sprint(print_hodge_diamond, P3) ==
    "   1\n  0 0\n 0 1 0\n0 0 0 0\n 0 1 0\n  0 0\n   1\n"
  for p in -1:4, q in -1:4
    expected = p == q && 0 <= p <= dim(P3) ? betti_number(P3, 2 * p) : 0
    @test hodge_number(P3, p, q) == expected
  end
  @test sum((-1)^i * betti_numbers(P3)[i + 1] for i in 0:(2 * dim(P3))) ==
    euler_characteristic(P3)
end

@testset "Hodge numbers of a blowup of projective space P2" begin
  blown_up_P2 = domain(blow_up(projective_space(NormalToricVariety, 2), 1))
  @test betti_numbers(blown_up_P2) == ZZRingElem[1, 0, 2, 0, 1]
  @test !has_attribute(blown_up_P2, :hodge_numbers)
  @test hodge_numbers(blown_up_P2) == matrix(ZZ, [1 0 0; 0 2 0; 0 0 1])
  cached_H = get_attribute(blown_up_P2, :hodge_numbers)
  @test all(i -> isassigned(cached_H, i), eachindex(cached_H))
  @test sum(betti_numbers(blown_up_P2)) == euler_characteristic(blown_up_P2)
end

@testset "Hodge numbers of a weighted projective space P_112" begin
  singular_WP2 = weighted_projective_space(NormalToricVariety, [1, 1, 2])
  @test !is_smooth(singular_WP2)
  @test is_complete(singular_WP2) && is_simplicial(singular_WP2)
  @test betti_numbers(singular_WP2) == ZZRingElem[1, 0, 1, 0, 1]
  @test hodge_numbers(singular_WP2) == matrix(ZZ, [1 0 0; 0 1 0; 0 0 1])
  @test sum(betti_numbers(singular_WP2)) == euler_characteristic(singular_WP2)
end

@testset "Hodge numbers exceptions for toric varieties" begin
  affine_with_torus_factor = affine_normal_toric_variety(Oscar.positive_hull([1 0]))
  nonsimplicial = normal_toric_variety(normal_fan(cross_polytope(3)))
  @test has_torusfactor(affine_with_torus_factor)
  @test !is_simplicial(nonsimplicial)
  for v in (affine_with_torus_factor, nonsimplicial)
    @test_throws ArgumentError betti_numbers(v)
    @test_throws ArgumentError hodge_number(v, 0, 0)
    @test_throws ArgumentError hodge_numbers(v)
  end
end
