
@testset "homogenization" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"])
    S, (a, b, c) = PolynomialRing(QQ, ["a", "b","c"])
    F = x^2 + 1
    G = x*y + y^3 + x
    V = [F, G]
    I = ideal(R, V)
    @test homogenization(S, F, 3) == S(a^2 + c^2)
    W = homogenization(V, variable = "z")
    T = parent(W[1])
    E = gens(T)
    @test W == [E[2]^2 + E[1]^2, E[1]*E[2]*E[3] + E[3]^3 + E[1]^2*E[2]]
    J = homogenization(S, I, 1)
    JJ = ideal(S, [b^2 + a^2, b*c*a + b*a^2 + c^3])
    @test J == JJ
end
