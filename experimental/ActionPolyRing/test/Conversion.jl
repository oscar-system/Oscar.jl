@testset "all tests - Conversion.jl" verbose = true begin

  __algebraic_conversion_data = Oscar.__algebraic_conversion_data
  __algebraic_conversion_data_from_polys = Oscar.__algebraic_conversion_data_from_polys
  __idd = Oscar.__is_data_up_to_date

  @testset "Check conversion" begin
    dpr_diffnc, vars_diffnc = difference_polynomial_ring(QQ, [:u, :v], 2)
    dpr_diffntl, vars_diffntl = differential_polynomial_ring(QQ, [:u, :v], 2)

    for (dpr, (u, v)) in [(dpr_diffnc, vars_diffnc), (dpr_diffntl, vars_diffntl)]

      u_0 = u[0, 0]
      u_x = dpr[1, [1, 0]]
      u_y = dpr[1, [0, 1]]
      v_x = dpr[2, [1, 0]]
      v_y = dpr[2, [0, 1]]

      if dpr isa DifferencePolyRing
        other_dpr, w = difference_polynomial_ring(QQ, :w, 2)
      elseif dpr isa DifferentialPolyRing
        other_dpr, w = differential_polynomial_ring(QQ, :w, 2)
      end

      @testset "Construction & Initialization" begin
        @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x + 1]) # Not a jet variable
        @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x, u_y, u_x]) # Duplicate jet variable
        @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x, w[1, 0]]) # Distinct parents of jet variables

        conv = __algebraic_conversion_data(dpr, [u_x, u_y, v_x])
        @test conv.apr === dpr
        @test __idd(conv)
        @test conv.supported_vars == Set([u_x, u_y, v_x])
        @test ngens(conv.mpr) == 3
      end

      @testset "Empty conversion data" begin
        conv = __algebraic_conversion_data(dpr, elem_type(dpr)[])

        @test conv(conv(dpr(0))) == dpr(0)
        @test conv(conv(dpr(17))) == dpr(17)
        @test ngens(conv.mpr) == 0
        @test_throws ArgumentError conv(u)
        @test_throws ArgumentError conv(v)
      end

      @testset "Forward and Backward mappings" begin
        conv = __algebraic_conversion_data(dpr, [u_0, u_x, u_y, v_x])

        @test conv(conv(dpr(0))) == dpr(0)
        @test conv(conv(dpr(17))) == dpr(17)

        for jetvar in [u_0, u_x, u_y, v_x]
          @test conv(conv(jetvar)) == jetvar
        end

        p1 = u_x^2 * v_x - 3*u_y + u_0 + 1
        cp1 = conv(p1)
        @test parent(cp1) === conv.mpr

        ccp1 = conv(cp1)
        @test parent(ccp1) === dpr
        @test ccp1 == p1
      end

      @testset "Error handling" begin
        conv = __algebraic_conversion_data(dpr, [u_x, u_y])
        @test_throws ArgumentError conv(u_x^2 + v_y) # v_y not supported

        @test_throws ArgumentError conv(w[1, 0]) # wrong parent (fwd map)

        S, (x, y) = polynomial_ring(QQ, [:x, :y])
        @test_throws ArgumentError conv(x + y) # wrong parent (bwd map)

        _ = dpr[1, [10, 10]]
        @test_throws ArgumentError conv(u_x) # outdated data (fwd map)
        @test_throws ArgumentError conv(gen(conv.mpr, 1)) # outdated data (bwd map)
      end

      @testset "Variable ordering (revsort)" begin
        cs1 = __algebraic_conversion_data(dpr, [u_0, u_y, u_x])
        cs2 = __algebraic_conversion_data_from_polys(dpr, [u_0, u_y, u_x])
        for conv_sorted in [cs1, cs2]
          @test u_0 < u_y < u_x
          @test conv_sorted(u_x) == gen(conv_sorted.mpr, 1)
          @test conv_sorted(u_y) == gen(conv_sorted.mpr, 2)
          @test conv_sorted(u_0) == gen(conv_sorted.mpr, 3)
        end

        conv_unsorted = __algebraic_conversion_data(dpr, [u_0, u_y, u_x]; revsort=false)
        @test conv_unsorted(u_0) == gen(conv_unsorted.mpr, 1)
        @test conv_unsorted(u_y) == gen(conv_unsorted.mpr, 2)
        @test conv_unsorted(u_x) == gen(conv_unsorted.mpr, 3)
      end

      @testset "Extracting variables" begin
        p = v_y^2 + v_x - u_0
        q = v_x * u_x - u_y

        @test_throws ArgumentError __algebraic_conversion_data(dpr, [p, q])
        conv = __algebraic_conversion_data_from_polys(dpr, [p, q])
        @test conv.supported_vars == Set(union(vars(p), vars(q)))
        @test u_0 < v_y < u_y < v_x < u_x

        @test conv(u_0) == gen(conv.mpr, 5)
        @test conv(v_y) == gen(conv.mpr, 4)
        @test conv(u_y) == gen(conv.mpr, 3)
        @test conv(v_x) == gen(conv.mpr, 2)
        @test conv(u_x) == gen(conv.mpr, 1)
        @test_throws ArgumentError conv(v[0, 0])

        @test __algebraic_conversion_data_from_polys(dpr, [p]).supported_vars == Set(vars(p))
        @test __algebraic_conversion_data_from_polys(dpr, [p, p]).supported_vars == Set(vars(p))
        @test __algebraic_conversion_data_from_polys(dpr, [p, p + 1]).supported_vars == Set(vars(p))
        @test __algebraic_conversion_data_from_polys(dpr, [q]).supported_vars == Set(vars(q))
        @test __algebraic_conversion_data_from_polys(dpr, [zero(dpr)]).supported_vars == Set{elem_type(dpr)}()
        @test __algebraic_conversion_data_from_polys(dpr, elem_type(dpr)[]).supported_vars == Set{elem_type(dpr)}()
      end
    end # for
  end
end # all tests

