@testset "all tests - ThomasDec" verbose = true begin
  @testset "Differential Reduction Methods" begin
    @testset "Single differential indeterminate" begin
      dpr, u = differential_polynomial_ring(QQ, :u, 2)
     
      u_x = u[1, 0]
      u_y = u[0, 1]
      u_xy = u[1, 1]
      u_xxy = u[2, 1]

      @testset "separant" begin
        q1 = u_x - u^2
        @test separant(q1) == dpr(1)
      
        q2 = u_xy^2 + u_x
        @test separant(q2) == 2*u_xy
        @test separant(dpr(5)) == dpr(5)
        @test separant(dpr(0)) == 0
      end

      @testset "pseudorem and pseudodivrem" begin
        p1 = u_x^2 + u
        q1 = u_x - u^2

        @test pseudorem(p1, q1) == u^4 + u
        quot1, rem1 = pseudodivrem(p1, q1)
        @test quot1 == u_x + u^2
        @test rem1 == u^4 + u

        p2 = u_x^2 + u_x
        q2 = u * u_x - 1

        @test pseudorem(p2, q2) == u + 1
        quot2, rem2 = pseudodivrem(p2, q2)
        @test quot2 == u * u_x + u + 1
        @test rem2 == u + 1

        @test pseudorem(p2, q2, u_x) == u + 1
        q3, r3 = pseudodivrem(p2, q2, u_x)
        @test q3 == quot2
        @test r3 == rem2
      end
    
      @testset "_leader_shift_for_partial_reduction" begin
        foo = Oscar._leader_shift_for_partial_reduction
        q = u_x - u^2
      
        p1 = u_xxy + u
        @test foo(p1, q) == [1, 1]
      
        p2 = u_xy + u
        @test foo(p2, q) == [0, 1]
      
        p3 = u_y + u
        @test foo(p3, q) === nothing
        @test foo(dpr(5), q) === nothing
      end
    
      @testset "partially_reduce" begin
        p = u_xxy + u
        q = u_x - u^2
      
        # Iteration 1: Q1 = apply_action(q, [1,1]) = u_xxy - 2*u_x*u_y - 2*u*u_xy
        # p1 = p - Q1 = 2*u*u_xy + 2*u_x*u_y + u
        # Iteration 2: Q2 = apply_action(q, [0,1]) = u_xy - 2*u*u_y
        # p2 = p1 - 2u*Q2 = 4*u^2*u_y + 2*u_x*u_y + u
        @test partially_reduce(p, q) == 4*u^2*u_y + 2*u_x*u_y + u
        @test partially_reduce(u_y + u, q) == u_y + u 
        
        @test partially_reduce(p, dpr(5)) == zero(p)
        @test_throws ArgumentError partially_reduce(p, dpr(0))
      end
    
      @testset "reduce" begin
        p = u_xxy + u_x^2 + u
        q = u_x - u^2
      
        # 1. Partial reduction eliminates u_xxy, leaving:
        # p_pred = 4*u^2*u_y + 2*u_x*u_y + u_x^2 + u
        # 2. Algebraic reduction pseudo-divides out u_x^2 and u_x using (u_x - u^2):
        # p_fullred = 6*u^2*u_y + u^4 + u
      
        @test reduce(p, q) == 6*u^2*u_y + u^4 + u
        @test reduce(u_y + u, q) == u_y + u 
        
        @test reduce(p, dpr(5)) == zero(p)
        @test_throws ArgumentError reduce(p, dpr(0))
      end
    end
    @testset "Two differential indeterminates" begin
      dpr, (u, v) = differential_polynomial_ring(QQ, [:u, :v], 2)

      u_x = dpr[1, [1, 0]]
      u_y = dpr[1, [0, 1]]
      u_xx = dpr[1, [2, 0]]
      u_xy = dpr[1, [1, 1]]
      u_yy = dpr[1, [0, 2]]
      u_xxy = dpr[1, [2, 1]]

      v_x = dpr[2, [1, 0]]
      v_y = dpr[2, [0, 1]]
      v_xx = dpr[2, [2, 0]]
      v_xy = dpr[2, [1, 1]]
      
      @testset "separant" begin
        q1 = u_xx + u_x^2
        q2 = u_x^3 + v_x*u_x
        q3 = -2*v_y*v_xx^5 - v_xy^2*v_xx^2 + 2*v_xx 
        q4 = dpr(1)
        @test separant(q1) == 1 
        @test separant(q2) == 3*u_x^2 + v_x
        @test separant(q3) == -10*v_y*v_xx^4 - 2*v_xy^2*v_xx + 2
        @test separant(q4) == dpr(1)
      end
      @testset "pseudorem and pseudodivrem" begin
        p1 = u_x^2 + v_x
        q1 = u_x - v_y
        
        @test pseudorem(p1, q1) == v_y^2 + v_x
        
        quot1, rem1 = pseudodivrem(p1, q1)
        @test quot1 == u_x + v_y
        @test rem1 == v_y^2 + v_x

        p2 = v_xx^2 + u_x
        q2 = u_y * v_xx - v_x
        
        @test pseudorem(p2, q2) == v_x^2 + u_y^2 * u_x
        
        quot2, rem2 = pseudodivrem(p2, q2)
        @test quot2 == u_y * v_xx + v_x
        @test rem2 == v_x^2 + u_y^2 * u_x

        @test pseudorem(p2, q2, v_xx) == v_x^2 + u_y^2 * u_x
        
        quot2_, rem2_ = pseudodivrem(p2, q2, v_xx)
        @test quot2_ == quot2
        @test rem2_ == rem2
        
        p3 = u_x * v_x + u_x
        q3 = v_x + 1
        
        @test pseudorem(p3, q3) == dpr(0)
        quot3, rem3 = pseudodivrem(p3, q3)
        @test quot3 == u_x
        @test rem3 == dpr(0)

        # Trivial Case 1: Degree of p is strictly less than degree of q
        p_triv1 = v_x
        q_triv1 = v_x^2 + u
        
        @test pseudorem(p_triv1, q_triv1) == v_x
        qt1, rt1 = pseudodivrem(p_triv1, q_triv1)
        @test qt1 == dpr(0)
        @test rt1 == v_x

        # Trivial Case 2: Dividend is exactly zero
        p_triv2 = dpr(0)
        q_triv2 = u_x + v
        
        @test pseudorem(p_triv2, q_triv2) == dpr(0)
        qt2, rt2 = pseudodivrem(p_triv2, q_triv2)
        @test qt2 == dpr(0)
        @test rt2 == dpr(0)

        # Trivial Case 3: Divisor is the zero polynomial 
        @test_throws DivideError pseudorem(u_x + v, dpr(0))
        @test_throws DivideError pseudodivrem(u_x + v, dpr(0))
        @test_throws DivideError pseudorem(u_x + v, dpr(0), u_x)
        @test_throws DivideError pseudodivrem(u_x + v, dpr(0), u_x)

        @testset "nonzero constants" begin
          R, x = differential_polynomial_ring(ZZ, :x, 1)
          @test pseudorem(x, R(2)) == zero(R)
          @test pseudorem(2*x, R(2)) == zero(R)
          @test pseudodivrem(x, R(2)) == (x, R(0))
          @test pseudodivrem(2*x, R(2)) == (x, R(0))
          @test pseudodivrem(R(6), R(2)) == (3, R(0))
          @test pseudodivrem(R(5), R(2)) == (R(5), R(0))
        end
      end

      @testset "partially_reduce" begin
        p1 = u_xxy + v_x
        q1 = u_x - v^2
        @test partially_reduce(p1, q1) == 2*v_x*v_y + 2*v*v_xy+v_x
        @test partially_reduce(p1, [q1]) == 2*v_x*v_y + 2*v*v_xy+v_x
        @test partially_reduce(p1, [q1, q1]) == 2*v_x*v_y + 2*v*v_xy+v_x
        @test partially_reduce(p1, [q1, dpr(1)]) == zero(p1)
        @test_throws ArgumentError partially_reduce(p1, [q1, dpr(1), dpr(0)])
        
        p5 = u_x^3 + v_x
        q5 = u_x^2 + v
        @test partially_reduce(p5, q5) == p5
      end

      @testset "reduce" begin
        p2 = u_yy + u_y*v_x + u
        q2 = u_y - v^2
        @test reduce(p2, q2) == v^2*v_x + 2*v*v_y + u

        p3 = u_y + v_x
        q3 = u_x - v_y
        @test reduce(p3, q3) == p3

        p4 = u_x + v_y
        q4 = u_x + v_y
        @test reduce(p4, q4) == dpr(0)

        p5 = u_x^3 + v_x
        q5 = u_x^2 + v
        @test reduce(p5, q5) == v_x - u_x * v
      end
    end # two diff indets
    @testset "reduce wrt a set" begin
      dpr, (u, v) = differential_polynomial_ring(QQ, [:u, :v], 1, index_ordering_name=:deglex)
      
      q1 = v[1] - u[0]
      q2 = u[1] - v[0]
      
      p = v[2]
      
      @test reduce(p, [q1, q2]) == v[0]
      @test partially_reduce(p, [q1, q2]) == u[1]

      dpr, (u, v) = difference_polynomial_ring(QQ, [:u, :v], 1, index_ordering_name=:deglex)

      q1 = v[1] - u[0]
      q2 = u[1] - v[0]

      p = v[2]

      @test reduce(p, [q1, q2]) == v[0]
      @test partially_reduce(p, [q1, q2]) == u[1]
    end
  end # Differential reduction methods
  @testset "Difference reduction methods" begin
    @testset "Single shift operator " begin
      dpr, u = difference_polynomial_ring(QQ, :u, 2)

      u_0 = u[0, 0]
      u_x = u[1, 0]
      u_y = u[0, 1]
      u_xx = u[2, 0]
      u_xy = u[1, 1]

      @testset "partially reduce" begin
        q = u_x^2 - u_0

        p1 = u_xx + u_0 
        @test partially_reduce(p1, q) == p1

        p2 = u_xx^2 + u_xy 
        @test partially_reduce(p2, q) == u_xy + u_x

        p3 = u_xx^3 + u_0 
        @test partially_reduce(p3, q) == u_xx * u_x + u_0
      end
      @testset "reduce" begin
        q = u_x^2 - u_0
        
        p1 = u_x^3 + u_x
        @test reduce(p1, q) == u_x * u_0 + u_x

        p2 = u_xx^2 + u_x^2
        @test reduce(p2, q) == u_x + u_0

        q3 = u_y * u_x^2 - dpr(1)
        p3 = u_x^2 + u_0
        @test reduce(p3, q3) == u_y * u_0 + dpr(1)

        p4 = u_xx + u_x + u_0
        @test reduce(p4, q) == p4
      end
    end # single shift operator

    @testset "Two shift operators" begin
      dpr, (u, v) = difference_polynomial_ring(QQ, [:u, :v], 2)

      u_0 = u[0, 0]
      u_x = u[1, 0]
      u_y = u[0, 1]
      u_xx = u[2, 0]

      v_0 = v[0, 0]
      v_x = v[1, 0]
      v_y = v[0, 1]
      v_xy = v[1, 1]

      @testset "partially reduce" begin
        q1 = u_x^2 - v_y
        p1 = u_xx^2 - u_0
        @test partially_reduce(p1, q1) == v_xy - u_0

        q2 = v_0 * u_x^2 - 1
        p2 = u_xx^2
        @test partially_reduce(p2, q2) == dpr(1)
      end
      @testset "reduce" begin
        q1 = u_x^2 - v_y
        p1 = v_x * u_x^3 + v_0
        @test reduce(p1, q1) == v_x * u_x * v_y + v_0

        q2 = u_x^2 - v_0
        p2 = u_xx^2 + u_x^2
        @test reduce(p2, q2) == v_x + v_0

        q3 = v_0 * u_x^2 - u_0
        p3 = u_x^3 + v_x
        @test reduce(p3, q3) == u_0 * u_x + v_0 * v_x

        q4 = u_x^2 - v_y
        p4 = u_xx + v_x * u_x + u_0
        @test reduce(p4, q4) == p4
      end
    end # two diff indets
  end # Difference reduction methods
  
  @testset "autoreduction" begin
    @testset "single differential indeterminate and single action map" begin
      dxr, x = difference_polynomial_ring(QQ, :x, 1)

      p1 = x[1] - x[0]^2
      p2 = x[2]^2 - x[0]
      @test autoreduce([p1, p2]) == [x[0]^8 - x[0],  x[1] - x[0]^2]
      @test autoreduce([p2, p1]) == [x[0]^8 - x[0],  x[1] - x[0]^2]
      @test autoreduce([p1, p1]) == [p1]
      @test autoreduce(DifferencePolyRingElem[]) == DifferencePolyRingElem[]
      @test autoreduce([dxr(0), dxr(0)]) == DifferencePolyRingElem[]
      @test autoreduce([dxr(0), dxr(-1), dxr(2)]) == [dxr(-1)]


      dpr, u = differential_polynomial_ring(QQ, :u, 1)
      q1 = u[1]^2 - u[0]
      q2 = u[2] - u[1]
      @test autoreduce([q1, q2]) == [u[0]]
      @test autoreduce([q2, q1]) == [u[0]]
      @test autoreduce([q1, q1]) == [q1]
      @test autoreduce(DifferentialPolyRingElem[]) == DifferentialPolyRingElem[]
      @test autoreduce([dpr(0), dpr(0)]) == DifferentialPolyRingElem[]
      @test autoreduce([dpr(0), dpr(-1), dpr(2)]) == [dpr(-1)]
    end
    @testset "two differential indeterminates and two action maps" begin
      dpr, (u, v) = differential_polynomial_ring(QQ, [:u, :v], 2; index_ordering_name=:degrevlex)

      p1 = u[1, 1] - u[0, 0]
      p2 = u[1, 0] - v[1, 0] 
      p3 = v[1, 1] - u[0, 0]

      @test autoreduce([p1, p2, p3]) == [p2, p3] 
      @test autoreduce([p3, p2, p1]) == [p2, p3] 
      @test autoreduce([p2, p3, p1]) == [p2, p3] 
      @test autoreduce([p1, p2]) == [p2, p3] 

      p4 = u[1, 1] - u[0, 1]
      p5 = p2
      p6 = v[1, 1] + v[0, 0]
      
      res = [-(u[0, 1] + v[0, 0]), u[1, 0] - v[1, 0], v[1, 1] + v[0, 0]] 
      
      @test autoreduce([p4, p5, p6]) == res 
      @test autoreduce([p6, p5, p4]) == res 
      @test autoreduce([p5, p6, p4]) == res 
    end
  end

end # all tests 
