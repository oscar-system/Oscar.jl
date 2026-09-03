@testset "all tests - Conversion.jl" verbose = true begin
  
  __algebraic_conversion_data = Oscar.__algebraic_conversion_data
  __jtv = Oscar.__jtv

  @testset "Difference to algebraic" begin
    dpr, (u, v) = difference_polynomial_ring(QQ, [:u, :v], 2)
    
    u_0 = u[0, 0]
    u_x = dpr[1, [1, 0]]
    u_y = dpr[1, [0, 1]]
    v_x = dpr[2, [1, 0]]
    v_y = dpr[2, [0, 1]]

    @testset "Construction & Initialization" begin
      @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x + 1]) # Not a jet variable
      @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x, u_y, u_x]) # Duplicate jet variable
      
      dpr2, w = difference_polynomial_ring(QQ, :w, 2)
      @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x, w[1, 0]])

      conv = __algebraic_conversion_data(dpr, [u_x, u_y, v_x])
      @test conv.apr === dpr
      @test conv.original_ring_size == length(__jtv(dpr))
      @test conv.supported_vars == Set([u_x, u_y, v_x])
      @test ngens(conv.mpr) == 3
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
      
      @test_throws ArgumentError conv(u_x^2 + v_y)
      
      dpr2, u2 = difference_polynomial_ring(QQ, :u, 2)
      @test_throws ArgumentError conv(u2[1, 0])
      
      S, (x, y) = polynomial_ring(QQ, [:x, :y])
      @test_throws ArgumentError conv(x + y)
      
      _ = dpr[1, [10, 10]]
        
      @test_throws ArgumentError conv(u_x)
        
      mpr_gen = gen(conv.mpr, 1)
      @test_throws ArgumentError conv(gen(conv.mpr, 1))
    end
    @testset "Variable ordering (revsorted)" begin
      conv_sorted = __algebraic_conversion_data(dpr, [u_0, u_y, u_x])
      @test u_0 < u_y < u_x
      @test conv_sorted(u_x) == gen(conv_sorted.mpr, 1)
      @test conv_sorted(u_y) == gen(conv_sorted.mpr, 2)
      @test conv_sorted(u_0) == gen(conv_sorted.mpr, 3)
      
      conv_unsorted = __algebraic_conversion_data(dpr, [u_0, u_y, u_x]; revsorted=false)
      @test conv_unsorted(u_0) == gen(conv_unsorted.mpr, 1)
      @test conv_unsorted(u_y) == gen(conv_unsorted.mpr, 2)
      @test conv_unsorted(u_x) == gen(conv_unsorted.mpr, 3)
    end
  end
  
  @testset "Differential to algebraic" begin
    dpr, (u, v) = differential_polynomial_ring(QQ, [:u, :v], 2)
    
    u_0 = u[0, 0]
    u_x = dpr[1, [1, 0]]
    u_y = dpr[1, [0, 1]]
    v_x = dpr[2, [1, 0]]
    v_y = dpr[2, [0, 1]]

    @testset "Construction & Initialization" begin
      @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x + 1]) # Not a jet variable
      @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x, u_y, u_x]) # Duplicate jet variable
      
      dpr2, w = differential_polynomial_ring(QQ, :w, 2)
      @test_throws ArgumentError __algebraic_conversion_data(dpr, [u_x, w[1, 0]])

      conv = __algebraic_conversion_data(dpr, [u_x, u_y, v_x])
      @test conv.apr === dpr
      @test conv.original_ring_size == length(__jtv(dpr))
      @test conv.supported_vars == Set([u_x, u_y, v_x])
      @test ngens(conv.mpr) == 3
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
      
      @test_throws ArgumentError conv(u_x^2 + v_y)
      
      dpr2, u2 = differential_polynomial_ring(QQ, :u, 2)
      @test_throws ArgumentError conv(u2[1, 0])
      
      S, (x, y) = polynomial_ring(QQ, [:x, :y])
      @test_throws ArgumentError conv(x + y)
      
      _ = dpr[1, [10, 10]]
        
      @test_throws ArgumentError conv(u_x)
        
      mpr_gen = gen(conv.mpr, 1)
      @test_throws ArgumentError conv(gen(conv.mpr, 1))
    end
    @testset "Variable ordering (revsorted)" begin
      conv_sorted = __algebraic_conversion_data(dpr, [u_0, u_y, u_x])
      @test u_0 < u_y < u_x
      @test conv_sorted(u_x) == gen(conv_sorted.mpr, 1)
      @test conv_sorted(u_y) == gen(conv_sorted.mpr, 2)
      @test conv_sorted(u_0) == gen(conv_sorted.mpr, 3)
      
      conv_unsorted = __algebraic_conversion_data(dpr, [u_0, u_y, u_x]; revsorted=false)
      @test conv_unsorted(u_0) == gen(conv_unsorted.mpr, 1)
      @test conv_unsorted(u_y) == gen(conv_unsorted.mpr, 2)
      @test conv_unsorted(u_x) == gen(conv_unsorted.mpr, 3)
    end
  end
end
