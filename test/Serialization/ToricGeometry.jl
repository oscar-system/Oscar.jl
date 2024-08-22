@testset "ToricGeometry" begin
  mktempdir() do path
    @testset "NormalToricVariety" begin
      pp = projective_space(NormalToricVariety, 2)
      R = cox_ring(pp)
      check(x) = has_attribute(x, :cox_ring)

      test_save_load_roundtrip(path, pp; with_attrs=false, check_func=!check) do loaded
        @test rays(pp) == rays(loaded)
        @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
      end
      
      test_save_load_roundtrip(path, pp; with_attrs=true, check_func=check) do loaded
        @test rays(pp) == rays(loaded)
        @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
      end

      test_save_load_roundtrip(path, pp; check_func=check) do loaded
        @test rays(pp) == rays(loaded)
        @test ray_indices(maximal_cones(pp)) == ray_indices(maximal_cones(loaded))
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

  end
end
