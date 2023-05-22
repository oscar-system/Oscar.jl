@testset "groebner_fan" begin
    @testset "cyclic4" begin
        R,(x1,x2,x3,x4) = PolynomialRing(QQ,4)

        I = ideal([x1+x2+x3+x4,
                   x1*x2+x2*x3+x1*x4+x3*x4,
                   x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4,
                   x1*x2*x3*x4-1])

        Sigma = groebner_fan(I)
        @test isequal(f_vector(Sigma),[19,70,92,40])
        @test issimplicial(Sigma)
        @test !issmooth(Sigma)
    end

    @testset "nonregular Groebner fan" begin
        R,(a,b,c,d) = PolynomialRing(QQ,["a","b","c","d"])
        I = ideal([a*c*d + a^2*c - a*b,
                   a*d^2 - c,
                   a*d^4 + a*c])
        Sigma = groebner_fan(I)
        @test isequal(f_vector(Sigma),[63,206,225,81])
        @test !issimplicial(Sigma)
    end
end
