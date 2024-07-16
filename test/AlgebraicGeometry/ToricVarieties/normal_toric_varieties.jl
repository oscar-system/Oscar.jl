@testset "Normal toric varieties" begin
  
  ntv = normal_toric_variety(Oscar.normal_fan(Oscar.cube(2)))
  set_coordinate_names(ntv, ["x1", "x2", "y1", "y2"])
  ntv2 = normal_toric_variety(Oscar.cube(2))
  ntv3 = normal_toric_varieties_from_glsm(matrix(ZZ, [[1, 1, 1]]))
  ntv4 = normal_toric_varieties_from_star_triangulations(convex_hull([0 0 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]))
  
  @testset "Basic properties" begin
    @test is_complete(ntv) == true
    @test is_projective_space(ntv) == false
    @test torsion_free_rank(torusinvariant_cartier_divisor_group(ntv)) == 4
    @test torsion_free_rank(domain(map_from_torusinvariant_cartier_divisor_group_to_torusinvariant_weil_divisor_group(ntv))) == 4
    @test is_complete(ntv2) == true
    @test length(ntv3) == 1
    @test is_projective_space(ntv3[1]) == true
    @test length(ntv4) == 2
  end
  
  @testset "Speed test for Stanley-Reisner ideal (at most a few seconds)" begin
    success = false
    for i in 1:5
      ntv5 = normal_toric_variety(polarize(polyhedron(Polymake.polytope.rand_sphere(5, 60; seed=42))))
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
    ntv5 = normal_toric_variety(polarize(polyhedron(Polymake.polytope.rand_sphere(5, 60; seed=42))))
    @test ngens(stanley_reisner_ideal(ntv5)) == 1648
  end
end
