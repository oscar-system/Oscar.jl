@testset "ToricGeometry" begin
  mktempdir() do path
    @testset "NormalToricVariety" begin
      pp = projective_space(NormalToricVariety, 2)
      R = cox_ring(pp)
      check(x) = has_attribute(x, :cox_ring)

      test_save_load_roundtrip(path, pp; with_attrs=false, check_func=!check) do loaded
        @test rays(pp) == rays(loaded)
        @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
        @test coordinate_names(pp) == coordinate_names(loaded)
      end
      
      test_save_load_roundtrip(path, pp; with_attrs=true, check_func=check) do loaded
        @test rays(pp) == rays(loaded)
        @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
        @test coordinate_names(pp) == coordinate_names(loaded)
      end

      test_save_load_roundtrip(path, pp; check_func=check) do loaded
        @test rays(pp) == rays(loaded)
        @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
        @test coordinate_names(pp) == coordinate_names(loaded)
      end
    end

    @testset "ToricDivisor" begin
      pp = projective_space(NormalToricVariety, 2)
      td0 = toric_divisor(pp, [1,1,2])
      td1 = toric_divisor(pp, [1,1,3])
      vtd = [td0, td1]
      test_save_load_roundtrip(path, vtd) do loaded
        @test coefficients(td0) == coefficients(loaded[1])
        @test coefficients(td1) == coefficients(loaded[2])
        @test toric_variety(loaded[1]) === toric_variety(loaded[2])
      end
    end

    @testset "ToricDivisorClass" begin
      pp = projective_space(NormalToricVariety, 2)
      tdc0 = toric_divisor_class(toric_divisor(pp, [1,1,2]))
      tdc1 = toric_divisor_class(toric_divisor(pp, [1,1,3]))
      vtd = [tdc0, tdc1]
      test_save_load_roundtrip(path, vtd) do loaded
        @test divisor_class(tdc0) == divisor_class(loaded[1])
        @test divisor_class(tdc1) == divisor_class(loaded[2])
        @test toric_variety(loaded[1]) === toric_variety(loaded[2])
      end
    end

    @testset "ToricCohomologyClass" begin
      pp = projective_space(NormalToricVariety, 2)
      cc = volume_form(pp)
      cc_list = [cc]
      test_save_load_roundtrip(path, cc_list) do loaded
        @test cc == loaded[1]
        @test toric_variety(cc) === toric_variety(loaded[1])
      end
    end

    @testset "ToricChernClasses" begin
      pp = projective_space(NormalToricVariety, 2)
      cherns = chern_classes(pp)
      cc_list = [cherns]
      test_save_load_roundtrip(path, cc_list) do loaded
        @test cherns == loaded[1]
        @test integrate(loaded[1][2]) == 3
      end
    end

    @testset "Hodge and Betti numbers of toric varieties" begin
      expected_B = ZZRingElem[1, 0, 1, 0, 1]
      expected_H = matrix(ZZ, [1 0 0; 0 1 0; 0 0 1])

      before_computation = projective_space(NormalToricVariety, 2)
      function check_uncached(loaded)
        @test !has_attribute(loaded, :hodge_numbers)
        @test !has_attribute(loaded, :betti_numbers)
        return true
      end
      test_save_load_roundtrip(
        path, before_computation; check_func=check_uncached
      ) do loaded
        @test dim(loaded) == 2
      end

      after_computation = projective_space(NormalToricVariety, 2)
      original_H = hodge_numbers(after_computation)
      original_B = betti_numbers(after_computation)
      @test original_H == expected_H
      @test original_B == expected_B
      @test has_attribute(after_computation, :hodge_numbers)
      @test has_attribute(after_computation, :betti_numbers)
      @test !has_attribute(after_computation, :betti_number)
      function check_computed(loaded)
        @test !has_attribute(loaded, :hodge_numbers)
        @test has_attribute(loaded, :betti_numbers)
        @test get_attribute(loaded, :betti_numbers) == ZZRingElem[1, 1, 1]
        @test betti_numbers(loaded) == expected_B
        @test hodge_numbers(loaded) == expected_H
        @test has_attribute(loaded, :hodge_numbers)
        return true
      end
      test_save_load_roundtrip(path, after_computation; check_func=check_computed) do loaded
        @test has_attribute(loaded, :betti_numbers)
        @test get_attribute(loaded, :betti_numbers) == ZZRingElem[1, 1, 1]
        @test betti_numbers(loaded) == original_B
        @test hodge_numbers(loaded) == original_H
      end

      # Test serialization and defensive copies for a partial Betti-number cache.
      partial_computation = projective_space(NormalToricVariety, 2)
      returned_partial_betti = betti_number(partial_computation, 2)
      returned_partial_hodge = hodge_number(partial_computation, 1, 1)
      partial_H = get_attribute(partial_computation, :hodge_numbers)
      @test isassigned(partial_H, 2, 2)
      @test !isassigned(partial_H, 1, 1)
      function check_partial(loaded)
        @test has_attribute(loaded, :betti_numbers)
        @test get_attribute(loaded, :betti_numbers) == ZZRingElem[-1, 1, -1]
        @test !has_attribute(loaded, :hodge_numbers)
        @test hodge_numbers(loaded) == expected_H
        return true
      end
      test_save_load_roundtrip(
        path, partial_computation; check_func=check_partial
      ) do loaded
        @test has_attribute(loaded, :betti_numbers)
        @test get_attribute(loaded, :betti_numbers) == ZZRingElem[-1, 1, -1]
        @test betti_number(loaded, 2) == returned_partial_betti
        @test hodge_number(loaded, 1, 1) == returned_partial_hodge
      end
    end

  end
end
