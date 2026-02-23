using PuiseuxPolynomials
using Test
using Oscar

@testset "PuiseuxPolynomial.jl" begin
    @testset "TrivialTests" begin
        @test 1+1==2
        @test 1+2==3
        # @test 1+1==3 # this one fails
    end

    @testset "Construction" begin
        K, (t1,t2,t3) = Oscar.polynomial_ring(QQ, ["t1","t2","t3"])
        K_p,(tp1,tp2,tp3) = puiseux_polynomial_ring(QQ, ["t1","t2","t3"])
        @test K_p.underlyingPolynomialRing == K
        @test K_p == PuiseuxPolynomials.MPuiseuxPolyRing(K)

        h = 1+t1 + 2*t2+3*t1^4+t1*t2^4+t3^2
        g = PuiseuxPolynomials.MPuiseuxPolyRingElem(K_p,h)
        @test h != g
        @test g.scale == 1
        @test g.shift == [0,0,0]

        g = PuiseuxPolynomials.MPuiseuxPolyRingElem(K_p,h,[ZZ(1),ZZ(1),ZZ(1)],ZZ(3))
        @test g.scale==3
        @test g.shift==[ZZ(1),ZZ(1),ZZ(1)]

        h = t1^(2)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1
        @test g.shift == [2,0,0]
        @test normalize!(g) == false

        h = t1^2*(1 + t1)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1
        @test g.shift == [2,0,0]

        h = t1*t2^(2)*t3^(3)*(1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1
        @test g.shift == [1,2,3]

        h = t1*t2^(2)*t3^(3)*(1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h,skip_normalization=true)
        @test g.scale == 1
        @test g.shift == [0,0,0]
        @test normalize!(g) == true

        h = (1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h,[ZZ(1),ZZ(2),ZZ(3)],ZZ(3),skip_normalization=true)
        @test g.scale == 3
        @test g.shift == [1,2,3]
        
        K, _ = polynomial_ring(QQ, ["t1","t2","t3"])
        Kt, _ = puiseux_polynomial_ring(QQ,["t1","t2","t3"])
        @test Kt.underlyingPolynomialRing == K
    end
    
    @testset "Getters" begin
        K, (t1,t2,t3) = polynomial_ring(QQ, ["t1","t2","t3"])
        Kp, (tp1,tp2,tp3) = puiseux_polynomial_ring(QQ,["t1","t2","t3"])

	    @test K == PuiseuxPolynomial.underlying_polynomial_ring(Kp)
        @test QQ == PuiseuxPolynomials.base_ring(Kp)
        @test QQ == PuiseuxPolynomials.coefficient_ring(Kp)
        @test PuiseuxPolynomials.ngens(Kp) == 3
        @test PuiseuxPolynomials.gens(Kp) == [tp1,tp2,tp3]
        
        g = tp1^(1//2)+tp3^(1//3)
        @test PuiseuxPolynomials.elem_type(Kp) == typeof(g)
        @test PuiseuxPolynomials.parent_type(g) == typeof(Kp)
        @test PuiseuxPolynomials.base_ring_type(Kp) == typeof(QQ)
        @test PuiseuxPolynomials.parent(g) == Kp
        @test PuiseuxPolynomials.poly(g) == t1^3 + t3^2
        @test PuiseuxPolynomials.scale(g) == 6
        @test PuiseuxPolynomials.shift(g) == [0,0,0]

        g = tp1^(2//3)*tp1*tp2^(1//2)*tp3 + tp3^(3//7)*tp1*tp2^(1//2)*tp3 + tp2^(1//2)*tp1*tp2^(1//2)*tp3
        @test PuiseuxPolynomials.parent(g) == Kp
        @test PuiseuxPolynomials.poly(g) == t1^28 + t2^21 + t3^18
        @test PuiseuxPolynomials.shift(g) == [42,21,42]
        @test PuiseuxPolynomials.scale(g) == 2*3*7

        
        K, (t,) = puiseux_polynomial_ring(QQ,["t"])
        @test valuation(K(0)) == PosInf()
        @test valuation(t^(-1)) == -1
        @test valuation(t^(-2//3)+t^(-1//2)) == -2//3
        @test_throws AssertionError valuation(g)
    end

    @testset "Arithmetic" begin

        F, (up,vp,wp) = polynomial_ring(QQ,["u","v","w"])
        K, (u,v,w) = puiseux_polynomial_ring(QQ,["u","v","w"])
        g = v^(1//3)+u^(1//2)
        h = v^(1//3) + w^(1//3)

        @test monomials(g) == [u^(1//2),v^(1//3)]
        @test monomials(h) == [v^(1//3),w^(1//3)]
        @test collect(coefficients(g)) == [1,1]
        @test collect(coefficients(h)) == [1,1]
        @test collect(exponents(g)) == [[QQ(1//2),QQ(0),QQ(0)],[QQ(0),QQ(1//3),QQ(0)]]
        @test 0*g == 0
        @test 1*g == g
        @test g*1 == g
        @test 4*g == 4*u^(1//2) + 4*v^(1//3)
        @test g+0 == g
        @test 0+g == g
        @test h+g == u^(1//2) + 2*v^(1//3) + w^(1//3)
        @test h-g == w^(1//3)-u^(1//2)
        @test h*g == u^(1//2)*v^(1//3) + u^(1//2)*w^(1//3) + v^(2//3) + v^(1//3)*w^(1//3)
        @test (g)^3 == u^(3//2) + 3*u*v^(1//3) + 3*u^(1//2)*v^(2//3) + v
        @test (g)^1 == g
        @test (g)^0 == 1
        @test g//(2) == (1//2)*u^(1//2) + (1//2)*v^(1//3)

        g = u^(1//2)*v^(2//3) + w^(1//4)
        h = u^(2//3)
        @test g*h == u^(7//6)*v^(2//3)+w^(1//4)*u^(2//3)

        @test (u^(1//2)+u^(-1//2))^2 == u+2+u^(-1)

        g = puiseux_polynomial_ring_elem(K,up*vp*wp*((up^4)*(vp^4)*(wp^4)+1),[ZZ(5),ZZ(9),ZZ(13)],ZZ(4),skip_normalization=true)
        @test monomials(g) == [u^(5//2)*v^(7//2)*w^(9//2),u^(3//2)*v^(5//2)*w^(7//2)]
        @test collect(exponents(g)) == [[5//2,7//2,9//2],[3//2,5//2,7//2]]
        @test PuiseuxPolynomials.poly(g) == up*vp*wp*(up^4*vp^4*wp^4+1)
        @test PuiseuxPolynomials.scale(g) == 4
        @test PuiseuxPolynomials.shift(g) == [ZZ(5),ZZ(9),ZZ(13)]
        @test normalize!(g) == true
        @test PuiseuxPolynomials.poly(g) == 1+up^2*vp^2*wp^2
        @test PuiseuxPolynomials.scale(g) == 2
        g_c = PuiseuxPolynomials.rescale(g,ZZ(10))
        @test PuiseuxPolynomials.scale(g_c) == 10
        @test PuiseuxPolynomials.poly(g_c) == 1+up^10*vp^10*wp^10
        @test PuiseuxPolynomials.shift(g_c) == [ZZ(15),ZZ(25),ZZ(35)]
        @test collect(exponents(g_c)) == [[5//2,7//2,9//2],[3//2,5//2,7//2]]
        @test normalize!(K(0)) == false
    end

    @testset "Conversions" begin
        F, (up,vp,wp) = polynomial_ring(QQ,["u","v","w"])
        K, (u,v,w) = puiseux_polynomial_ring(QQ,["u","v","w"])
        @test typeof(3) == Int64
        @test K(3) == puiseux_polynomial_ring_elem(K,F(3))
        @test typeof(4//5) == Rational{Int64}
        @test K(4//5) == puiseux_polynomial_ring_elem(K,F(4//5))
        @test K(0) == zero(K)
        @test K(1) == one(K)
    end
end
