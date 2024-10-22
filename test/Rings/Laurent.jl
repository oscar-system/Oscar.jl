@testset "Laurent" begin
  for K in [QQ, GF(5)]
    Kx, x = laurent_polynomial_ring(K, 2, :x)
    I = ideal(Kx, [x[1]])
    @test gens(I) == [x[1]]
    @test one(Kx) in I
    @test x[1]^-1 in I
    @test base_ring(I) === Kx
    _Kx, _x = laurent_polynomial_ring(K, 3, :xx)
    @test !(_x[1] in I)

    f = hom(Kx, K, [K(2), K(3)])
    @test domain(f) === Kx
    @test codomain(f) === K
    @test f(x[1]^-1 + x[2]) == K(2)^-1 + K(3)

    @test_throws ArgumentError hom(Kx, ZZ, [ZZ(2), ZZ(3)])

    Ky, y = polynomial_ring(K, 2, :y)
    f = hom(Ky, Kx, gens(Kx))
    @test f(y[1]) == x[1]
    @test preimage(f, I) == ideal(Ky, [one(Ky)])

    f = hom(Ky, Kx, reverse(gens(Kx)))
    @test f(y[1]) == x[2]
    @test preimage(f, I) == ideal(Ky, [one(Ky)])

    Q, = quo(Ky, ideal(Ky, [y[1] * y[2] - 1]))
    f = hom(Kx, Q, gens(Q))
    for q in gens(Q)
      @test f(preimage(f, q)) == q
    end
  end
end
