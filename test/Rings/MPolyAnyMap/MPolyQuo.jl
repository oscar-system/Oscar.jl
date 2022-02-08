@testset "MPolyAnyMap/MPolyQuo" begin
  Qsqrt2, = quadratic_field(-1)
  Zx, _ = ZZ["x"]
  Zxy, _ = ZZ["x", "y"]
  
  for K in [GF(2), GF(fmpz(2)), GF(2, 2), GF(fmpz(2), 2), ZZ, QQ, Qsqrt2, Zx, Zxy]
    Kx, (x, y) = K["x", "y"]
    I = ideal(Kx, [x^2 - y])
    Q, = quo(Kx, I)

    h = @inferred hom(Q, Kx, [x, x^2])
    @test domain(h) === Q
    @test codomain(h) === Kx
    @test (@inferred h(Q(x))) == x
    @test (@inferred h(Q(x))) == x
    @test (@inferred h(Q(y))) == x^2
    @test (@inferred h(Q(y))) == x^2
    @test (@inferred h(x)) == x
    @test (@inferred h(x)) == x
    @test (@inferred h(y)) == x^2
    @test (@inferred h(y)) == x^2
    @test (@inferred h(K(1))) == K(1)
    @test (@inferred h(K(2))) == K(2)
    @test (@inferred h(1)) == K(1)
    @test (@inferred h(2)) == K(2)

    @test_throws ArgumentError hom(Q, Kx, [x, y])

    if K === QQ
      @test_throws MethodError hom(Q, Kx, [GF(2)(1), GF(2)(1)])
    end

    # noncommutative image
    S = MatrixAlgebra(K, 2)
    a = S([1 1; 0 1])
    b = S([0 1; 1 0])
    @test_throws ArgumentError hom(Q, S, [a, b])
    h = hom(Q, S, [a, b], check = false)

    a = S([2 0; 0 2])
    b = S([4 0; 0 4])
    h = hom(Q, S, [a, b])
    # inferred is broken because of AA bug
    @test h(x * y) == a * b
    @test h(x * y) == a * b
  end
end
