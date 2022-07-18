@testset "initial" begin
    Kx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(QQ,6)
    Cyclic5Homogenized = ideal([x1+x2+x3+x4+x5,
                                   x1*x2+x2*x3+x3*x4+x1*x5+x4*x5,
                                   x1*x2*x3+x2*x3*x4+x1*x2*x5+x1*x4*x5+x3*x4*x5,
                                   x1*x2*x3*x4+x1*x2*x3*x5+x1*x2*x4*x5+x1*x3*x4*x5+x2*x3*x4*x5,
                                   -x0^5+x1*x2*x3*x4*x5])
    Katsura5Homogenized = ideal([-x0+x1+2*x2+2*x3+2*x4+2*x5,
                                    -x0*x1+x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2,
                                    -x0*x2+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5,
                                    x2^2-x0*x3+2*x1*x3+2*x2*x4+2*x3*x5,
                                    2*x2*x3-x0*x4+2*x1*x4+2*x2*x5])
    w = [0,0,0,0,0,0]
    @testset "val_2" begin
        val_2 = TropicalSemiringMap(QQ,2)
        @testset "Cyclic5Homogenized" begin
            computed = gens(initial(Cyclic5Homogenized, val_2, w))
            # Apparently the Groebner basis lives in a different ring, so we have
            # to hack around this...
            computed_ring = parent(computed[1])
            (y1, y2, y3, y4, y5, y6) = gens(computed_ring)
            expected = [y2 + y3 + y4 + y5 + y6, y3^2 + y3*y5 + y4*y5 + y4*y6 + y6^2, y3*y4^2 + y3*y4*y5 + y3*y5*y6 + y3*y6^2 + y4^2*y5 + y4^2*y6 + y4*y5*y6 + y5^2*y6 + y5*y6^2 + y6^3, y3*y4*y5^2 + y3*y4*y5*y6 + y3*y4*y6^2 + y3*y5^2*y6 + y3*y5*y6^2 + y3*y6^3 + y4^2*y5*y6 + y4*y5^2*y6 + y4*y6^3 + y5^3*y6 + y5*y6^3 + y6^4, y3*y4*y5^2 + y3*y5*y6^2 + y4^3*y5 + y4^3*y6 + y4^2*y5^2 + y4*y5^2*y6 + y4*y6^3 + y5^2*y6^2, y3*y4*y5*y6^2 + y3*y5^3*y6 + y4^3*y6^2 + y4^2*y5^2*y6 + y4^2*y5*y6^2 + y4^2*y6^3 + y4*y5^3*y6 + y4*y5*y6^3 + y5^4*y6 + y5^3*y6^2, y3*y4*y6^3 + y3*y5^3*y6 + y3*y5^2*y6^2 + y3*y6^4 + y4^2*y5^3 + y4^2*y5*y6^2 + y4^2*y6^3 + y5^4*y6 + y5^3*y6^2 + y5^2*y6^3 + y5*y6^4 + y6^5, y1^5 + y3*y4*y6^3 + y3*y6^4 + y4^2*y5*y6^2 + y4*y5*y6^3 + y4*y6^4 + y5^2*y6^3 + y6^5, y3*y4*y6^4 + y3*y5^3*y6^2 + y3*y5^2*y6^3 + y3*y6^5 + y4^4*y6^2 + y4^3*y6^3 + y4^2*y5*y6^3 + y4*y6^5 + y5^2*y6^4 + y6^6, y3*y4*y6^4 + y3*y5^3*y6^2 + y3*y5^2*y6^3 + y3*y6^5 + y4^4*y6^2 + y4^2*y6^4 + y4*y5^4*y6 + y4*y5^3*y6^2 + y4*y6^5 + y5^4*y6^2 + y5^2*y6^4 + y6^6, y3*y4*y6^4 + y3*y5^4*y6 + y3*y5^3*y6^2 + y3*y5^2*y6^3 + y3*y5*y6^4 + y3*y6^5 + y4^3*y6^3 + y4^2*y6^4 + y4*y5^4*y6 + y4*y6^5 + y5^5*y6 + y5^3*y6^3 + y5*y6^5 + y6^6, y3*y4*y5*y6^4 + y3*y5^3*y6^3 + y4^2*y5*y6^4 + y4*y5^2*y6^4 + y4*y5*y6^5 + y5^5*y6^2 + y5^4*y6^3 + y5^3*y6^4]
            @test computed == expected
        end
        @testset "Katsura5Homogenized" begin
            computed = gens(initial(Katsura5Homogenized, val_2, w))
            # Apparently the Groebner basis lives in a different ring, so we have
            # to hack around this...
            computed_ring = parent(computed[1])
            (y1, y2, y3, y4, y5, y6) = gens(computed_ring)
            expected = [y1 + y2, y2*y6 + y4^2 + y5^2 + y6^2, y2*y5, y2*y4 + y2*y5 + y3^2, y2*y3, y2*y5^2, y2*y4*y5 + y2*y5*y6, y2^2*y5, y2^2*y4 + y2*y5^2, y2^3*y6 + y2^2*y6^2, y2*y5^2*y6, y2^2*y6^2 + y2*y4*y6^2 + y2*y5^2*y6 + y2*y5*y6^2 + y2*y6^3 + y3*y5^3 + y3*y5*y6^2 + y4*y5^2*y6 + y5^4 + y5^2*y6^2, y2*y5^3 + y2*y5^2*y6, y2^4*y6 + y2^3*y6^2 + y2*y5^2*y6^2, y2*y5^2*y6^2]
            @test computed == expected
        end
    end

    @testset "val_3" begin
        val_3 = TropicalSemiringMap(QQ,3)
        @testset "Cyclic5Homogenized" begin
            computed = gens(initial(Cyclic5Homogenized, val_3, w))
            # Apparently the Groebner basis lives in a different ring, so we have
            # to hack around this...
            computed_ring = parent(computed[1])
            (y1, y2, y3, y4, y5, y6) = gens(computed_ring)
            expected = [y2 + y3 + y4 + y5 + y6, y3^2 + y3*y5 + 2*y3*y6 + 2*y4*y5 + y4*y6 + y6^2, y3*y4^2 + 2*y3*y4*y5 + y3*y5*y6 + 2*y3*y6^2 + y4^2*y5 + 2*y4^2*y6 + y4*y5*y6 + y4*y6^2 + y5^2*y6 + y5*y6^2 + 2*y6^3, y3*y4*y5^2 + y3*y4*y5*y6 + 2*y3*y4*y6^2 + 2*y3*y5^2*y6 + 2*y3*y5*y6^2 + y3*y6^3 + y4^2*y5*y6 + y4*y5^2*y6 + y4*y6^3 + 2*y5^3*y6 + y5^2*y6^2 + 2*y5*y6^3 + y6^4, y3*y4*y5^2 + y3*y4*y6^2 + y3*y6^3 + y4^3*y5 + 2*y4^3*y6 + 2*y4^2*y5^2 + y4^2*y5*y6 + 2*y4^2*y6^2 + y4*y5^2*y6 + 2*y4*y5*y6^2 + y4*y6^3 + 2*y5*y6^3 + y6^4, 2*y3*y4*y5*y6^2 + y3*y4*y6^3 + 2*y3*y5^3*y6 + y3*y5^2*y6^2 + y3*y5*y6^3 + 2*y3*y6^4 + y4^3*y6^2 + y4^2*y5^2*y6 + 2*y4^2*y5*y6^2 + 2*y4^2*y6^3 + y4*y5^2*y6^2 + 2*y5^4*y6 + y5^3*y6^2 + 2*y6^5, y3*y4*y5*y6^2 + 2*y3*y5^3*y6 + 2*y3*y5^2*y6^2 + 2*y3*y5*y6^3 + 2*y3*y6^4 + y4^2*y5^3 + 2*y4^2*y5*y6^2 + y4^2*y6^3 + y4*y5^2*y6^2 + y4*y5*y6^3 + 2*y5^4*y6 + 2*y5^2*y6^3 + 2*y6^5, y1^5 + y3*y4*y5*y6^2 + y3*y4*y6^3 + 2*y3*y5*y6^3 + 2*y3*y6^4 + 2*y4^2*y5*y6^2 + y4*y5^2*y6^2 + y4*y5*y6^3 + 2*y4*y6^4 + y5^2*y6^3 + 2*y5*y6^4 + 2*y6^5, 2*y3*y4*y5*y6^3 + y3*y6^5 + y4^4*y6^2 + 2*y4^2*y6^4 + 2*y4*y5^3*y6^2 + 2*y4*y5^2*y6^3 + 2*y4*y5*y6^4 + y4*y6^5 + 2*y5^4*y6^2 + 2*y5^3*y6^3 + 2*y5^2*y6^4 + y5*y6^5 + y6^6, y3*y4*y5*y6^3 + 2*y3*y4*y6^4 + y3*y5^3*y6^2 + y3*y5^2*y6^3 + y3*y5*y6^4 + 2*y4^4*y6^2 + y4^3*y6^3 + 2*y4^2*y5*y6^3 + y4*y5^4*y6 + y4*y5^3*y6^2 + y4*y5^2*y6^3 + y4*y6^5 + y5^4*y6^2 + y5^3*y6^3 + 2*y5^2*y6^4 + 2*y5*y6^5, y3*y4*y5*y6^3 + y3*y5^4*y6 + y3*y5^3*y6^2 + 2*y3*y5*y6^4 + y3*y6^5 + 2*y4^3*y6^3 + 2*y4^2*y5*y6^3 + 2*y4*y5^3*y6^2 + 2*y4*y6^5 + y5^5*y6 + y5^4*y6^2 + 2*y5^3*y6^3 + 2*y5*y6^5 + y6^6, y3*y4*y5*y6^4 + 2*y3*y5^3*y6^3 + y3*y5^2*y6^4 + y3*y5*y6^5 + y3*y6^6 + 2*y4^2*y5*y6^4 + 2*y4^2*y6^5 + y4*y5^3*y6^3 + y4*y5*y6^5 + y5^5*y6^2 + y5^3*y6^4 + y6^7]
            @test computed == expected
        end
        @testset "Katsura5Homogenized" begin
            computed = gens(initial(Katsura5Homogenized, val_3, w))
            # Apparently the Groebner basis lives in a different ring, so we have
            # to hack around this...
            computed_ring = parent(computed[1])
            (y1, y2, y3, y4, y5, y6) = gens(computed_ring)
            expected = [y1 + 2*y2 + y3 + y4 + y5 + y6, y2*y6 + 2*y3*y5 + y4^2 + 2*y4*y5 + y5^2 + 2*y6^2, y2*y5 + 2*y3*y4 + y3*y5 + 2*y3*y6 + y4*y5 + y5^2 + y5*y6, y2*y4 + y2*y5 + y3^2 + 2*y3*y6 + y4^2 + 2*y4*y5 + y5^2 + y5*y6, y2*y3 + 2*y2*y4 + 2*y2*y5 + y3*y5 + 2*y3*y6 + 2*y4^2 + 2*y5^2 + y5*y6, y2^2*y6 + y2*y5^2 + y2*y5*y6 + y3*y5^2 + y3*y5*y6 + 2*y4*y5^2 + 2*y4*y6^2 + y5^3 + 2*y5^2*y6 + 2*y5*y6^2 + 2*y6^3, y2^2*y6 + y2*y4*y5 + 2*y2*y4*y6 + y2*y6^2 + y3*y5^2 + 2*y3*y6^2 + 2*y4*y5*y6 + 2*y4*y6^2 + 2*y5^3 + 2*y5^2*y6 + 2*y5*y6^2 + y6^3, y2^2*y5 + y2^2*y6 + 2*y2*y4*y6 + y2*y5^2 + 2*y2*y5*y6 + y2*y6^2 + 2*y3*y5^2 + y3*y5*y6 + 2*y3*y6^2 + y4*y5*y6 + 2*y4*y6^2 + 2*y5^3 + 2*y5^2*y6 + 2*y5*y6^2 + y6^3, y2^2*y4 + 2*y2^2*y6 + 2*y2*y4*y6 + 2*y2*y6^2 + 2*y3*y5^2 + y3*y5*y6 + 2*y3*y6^2 + y4*y5^2 + 2*y4*y5*y6 + y5^3 + y5^2*y6 + 2*y5*y6^2 + 2*y6^3, 2*y2^3*y6 + y2^2*y6^2 + 2*y2*y4*y6^2 + 2*y2*y5^2*y6 + 2*y2*y5*y6^2 + 2*y3*y5^2*y6 + y3*y6^3 + 2*y4*y6^3 + 2*y5^3*y6 + y5^2*y6^2 + 2*y5*y6^3, 2*y2^2*y6^2 + 2*y2*y4*y6^2 + y2*y5^2*y6 + y2*y6^3 + 2*y3*y5*y6^2 + 2*y3*y6^3 + 2*y4*y5^3 + y4*y5^2*y6 + 2*y4*y5*y6^2 + y4*y6^3 + 2*y5^2*y6^2 + 2*y5*y6^3, y2^2*y6^2 + 2*y2*y4*y6^2 + y2*y5^2*y6 + 2*y2*y5*y6^2 + y2*y6^3 + y3*y5^3 + y3*y5^2*y6 + 2*y3*y6^3 + 2*y4*y6^3 + 2*y5^4 + 2*y5^3*y6 + y5^2*y6^2 + y6^4, 2*y2^2*y6^2 + y2*y4*y6^2 + y2*y5^3 + y2*y5*y6^2 + 2*y2*y6^3 + y3*y5^3 + 2*y3*y5^2*y6 + 2*y3*y5*y6^2 + y3*y6^3 + 2*y4*y5^3 + y4*y5*y6^2 + y4*y6^3 + y5^4 + 2*y6^4, y2^4*y6 + y2^3*y6^2 + 2*y2*y4*y6^3 + y2*y5^2*y6^2 + 2*y2*y5*y6^3 + 2*y3*y5*y6^3 + 2*y4*y5^2*y6^2 + y4*y6^4 + 2*y5^4*y6 + y5^3*y6^2 + y5*y6^4 + y6^5, 2*y2*y4*y6^3 + y2*y6^4 + 2*y3*y5^2*y6^2 + 2*y3*y5*y6^3 + 2*y3*y6^4 + y4*y5*y6^3 + 2*y5^4*y6 + 2*y5^2*y6^3 + 2*y5*y6^4 + 2*y6^5]
            @test computed == expected
        end
    end
    @testset "val_t" begin
        Kt,t = RationalFunctionField(QQ,"t")
        Ktx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(Kt,6)
        Cyclic5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Cyclic5Homogenized)])
        Katsura5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Katsura5Homogenized)])
        val_t = TropicalSemiringMap(Kt,t)
        @testset "Cyclic5Homogenized_Kt" begin
            computed = gens(initial(Cyclic5Homogenized_Kt, val_t, w))
            # Apparently the Groebner basis lives in a different ring, so we have
            # to hack around this...
            computed_ring = parent(computed[1])
            (y1, y2, y3, y4, y5, y6) = gens(computed_ring)
            expected = [y2 + y3 + y4 + y5 + y6, y2*y6 - y3^2 - y3*y5 - y3*y6 + y4*y5 + y5*y6, -y2*y3*y6 + y2*y4*y6 - y2*y5*y6 + y3*y4^2 - y3*y4*y5 + y4^2*y5, -y2*y3*y4*y6 - y2*y4*y5*y6 + y2*y5^2*y6 + y3*y4*y5^2 - y3*y4*y5*y6 + y4*y5^2*y6, -y2*y3^2*y6 + y2*y3*y4*y6 - y2*y3*y5*y6 + y2*y4^2*y6 - y2*y4*y5*y6 - y3*y4^2*y6 + y3*y4*y5^2 + y3*y4*y5*y6 + y4^3*y5 - y4^2*y5^2 + y4^2*y5*y6 - y4*y5^2*y6, 2*y2*y3^2*y6^2 - 2*y2*y3*y4*y6^2 + 4*y2*y3*y5*y6^2 - y2*y4^2*y6^2 - y2*y4*y5*y6^2 + y2*y5^3*y6 + 2*y2*y5^2*y6^2 + y3*y4^2*y6^2 - 2*y3*y4*y5*y6^2 + y4^2*y5^2*y6 + y4^2*y5*y6^2 - 2*y4*y5^3*y6 + 3*y4*y5^2*y6^2 - 2*y5^3*y6^2, y2*y3^2*y4*y6 + 2*y2*y3*y4*y5*y6 - y2*y3*y5^2*y6 - y2*y5^3*y6 + y3^2*y4*y5*y6 + y3*y4*y5^2*y6 - y4^2*y5^3 - 2*y4*y5^3*y6, -y1^5 + y2*y3*y4*y5*y6, -3*y2*y3^2*y5*y6^2 + 2*y2*y3^2*y6^3 + 2*y2*y3*y4*y5*y6^2 + 2*y2*y3*y4*y6^3 - y2*y3*y5^2*y6^2 - 2*y2*y4^2*y5*y6^2 - 6*y2*y4^2*y6^3 + 4*y2*y4*y5^2*y6^2 + 3*y2*y4*y5*y6^3 + y2*y5^3*y6^2 + y3^2*y4*y5*y6^2 - 2*y3^2*y4*y6^3 + 2*y3^2*y5*y6^3 + 4*y3*y4^2*y5*y6^2 + 4*y3*y4^2*y6^3 - 5*y3*y4*y5^2*y6^2 - y3*y4*y5*y6^3 + 2*y3*y5^3*y6^2 + y4^4*y6^2 - 4*y4^3*y5*y6^2 + y4^3*y6^3 + 2*y4^2*y5^2*y6^2 - 3*y4^2*y5*y6^3 + y5^4*y6^2, y2*y3^2*y5*y6^2 - 2*y2*y3*y4^2*y6^2 - 3*y2*y3*y4*y5*y6^2 + 3*y2*y3*y5^2*y6^2 - y2*y4^3*y6^2 - 2*y2*y4^2*y5*y6^2 + y2*y4*y5^2*y6^2 + 4*y2*y5^3*y6^2 - 2*y3^2*y4*y5*y6^2 + y3*y4^3*y6^2 - y3*y4^2*y5*y6^2 - 5*y3*y4*y5^2*y6^2 + y4^3*y5*y6^2 + 2*y4^2*y5^2*y6^2 - y4*y5^4*y6 + 4*y4*y5^3*y6^2, -y2*y3^2*y4*y6^2 - 2*y2*y3^2*y5*y6^2 - 3*y2*y3*y5^2*y6^2 + y2*y4^2*y5*y6^2 + y2*y4*y5^2*y6^2 - y2*y5^3*y6^2 - y3^2*y4*y5*y6^2 - y3*y4^2*y5*y6^2 + y3*y4*y5^2*y6^2 + y3*y5^4*y6 - y4^2*y5^2*y6^2 + 3*y4*y5^4*y6 - y4*y5^3*y6^2 + y5^5*y6 + 3*y5^4*y6^2, -4*y2*y3^2*y5*y6^3 + 13*y2*y3*y4^2*y6^3 + 16*y2*y3*y4*y5*y6^3 - 11*y2*y3*y5^2*y6^3 + 6*y2*y4^3*y6^3 + 9*y2*y4^2*y5*y6^3 - 12*y2*y4*y5^2*y6^3 - 16*y2*y5^3*y6^3 + 11*y3^2*y4*y5*y6^3 + 2*y3^2*y5^2*y6^3 - 6*y3*y4^3*y6^3 + 10*y3*y4^2*y5*y6^3 + 23*y3*y4*y5^2*y6^3 + 3*y3*y5^3*y6^3 - 7*y4^3*y5*y6^3 - 12*y4^2*y5^2*y6^3 - 21*y4*y5^3*y6^3 - y5^5*y6^2 - 3*y5^4*y6^3]
            @test computed == expected
        end
    end
end
