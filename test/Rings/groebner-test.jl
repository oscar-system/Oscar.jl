@testset "groebner" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"])
    I = ideal(R,[x*y^2-x, x^3-2*y^5])
    @test leading_ideal(I) == ideal(R,[x*y^2, x^4, y^5])
    @test leading_ideal(I, :lex) == ideal(R,[y^7, x*y^2, x^3])
end
