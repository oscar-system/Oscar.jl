@testset "Cartier divisors and line bundles" begin
    IP = projective_space(QQ, ["x", "y", "z"])

    X = covered_scheme(IP)
    D = IdDict{AbsSpec, RingElem}()
    U = affine_charts(X)
    S = ambient_coordinate_ring(IP)
    for i in 1:3
        D[U[i]] = dehomogenize(IP, i-1)(S[2])^3
    end

    H = CartierDivisor(X, D)

    @test all(V->(parent(H(V)) === OO(V)), U)

    L = LineBundle(H)
    U21 = PrincipalOpenSubset(U[2], OO(U[2])[1])
    @test L(U[1], U21)(L(U[1])[1]) == inv(gens(OO(U21))[1]^3)*L(U21)[1]
end
