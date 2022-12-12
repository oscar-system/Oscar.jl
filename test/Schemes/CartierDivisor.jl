@testset "Cartier divisors" begin
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
end
