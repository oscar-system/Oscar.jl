@testset "MPolyAnyMap/MPolyRing" begin
  # Construction 
  Qsqrt2, = quadratic_field(-1)
  Zx, _ = ZZ["x"]
  Zxy, _ = ZZ["x", "y"]
  
  for K in [GF(2), GF(fmpz(2)), GF(2, 2), GF(fmpz(2), 2), ZZ, QQ, Qsqrt2, Zx, Zxy]
    Kx, (x, y) = K["x", "y"]
    h = @inferred hom(Kx, Kx, [y, x])
    @test sprint(show, h) isa String
    @test sprint(show, "text/plain", h) isa String
    @test domain(h) === Kx
    @test codomain(h) === Kx
    @test (@inferred h(x)) == y
    @test (@inferred h(x)) == y
    @test (@inferred h(y)) == x
    @test (@inferred h(y)) == x
    @test (@inferred h(K(1))) == K(1)
    @test (@inferred h(K(2))) == K(2)
    @test (@inferred h(1)) == K(1)
    @test (@inferred h(2)) == K(2)
    @test_throws ArgumentError hom(Kx, Kx, [x])
    @test_throws ArgumentError hom(Kx, Kx, [y, y, x])

    if K === QQ
      @test_throws MethodError hom(Kx, Kx, [GF(2)(1), GF(2)(1)])
    end

    h = @inferred hom(Kx, Kx, fmpz[1, 1])
    @test (@inferred h(x + y)) == Kx(2)
    @test (@inferred h(x + y)) == Kx(2)

    Kh, (z1, z2) = K["z1", "z2"]
    Kh, (zz1, zz2) = grade(Kh)
    h = @inferred hom(Kx, Kh, [z1 + z2, z2])
    @test (@inferred h(x)) == zz1 + zz2
    @test (@inferred h(x)) == zz1 + zz2
    @test (@inferred h(y)) == zz2
    @test (@inferred h(y)) == zz2
    @test_throws ArgumentError hom(Kx, Kh, [zz1 + zz1^2, zz2])
    h = @inferred hom(Kx, Kh, [zz1 + zz1^2, zz2], check = false)
  
    # julia-function  
    R, (x, y) = K["x", "y"]
    S, (u, v) = R["u", "v"]
    h = hom(S, S, let S = S; a -> S(a^2); end, gens(S))
    @test (@inferred h(u)) == u
    @test (@inferred h(u)) == u
    @test (@inferred h(x*u)) == x^2 * u
    @test (@inferred h(x*u)) == x^2 * u
    h = hom(S, S, let S = S; a -> S(a^2); end, gens(S))
    @test (h(u)) == u
    @test (h(u)) == u
    @test (h(x*u)) == x^2 * u
    @test (h(x*u)) == x^2 * u
    h = hom(S, S, a -> a^2, gens(S))
    @test (@inferred h(u)) == u
    @test (@inferred h(u)) == u
    @test (@inferred h(x*u)) == x^2 * u
    @test (@inferred h(x*u)) == x^2 * u
  
    A, (x,y) = K["x", "y"]
    f = hom(A, A, [2*x, 5*y])
    R, (u, v) = A["u", "v"]
    h = hom(R, R, f, [u+v, u*v])
    @test (@inferred h(x*u)) == 2*x*u + 2*x*v
    @test (@inferred h(x*u)) == 2*x*u + 2*x*v

    # Noncommutative image
    
    A, (x, y) = K["x", "y"]
    S = MatrixAlgebra(K, 2)
    a = S([1 1; 0 1])
    b = S([0 1; 1 0])
    @test_throws ArgumentError hom(A, S, [a, b])

    a = S([1 0; 0 2])
    b = S([2 0; 0 2])
    h = hom(A, S, [a, b])

    # bug in AA
    #@test (@inferred h(x * y)) == a * b
    #@test (@inferred h(x * y)) == a * b
    @test (h(x * y)) == a * b
    @test (h(x * y)) == a * b

    # another issue
    R, (x, y) = K["x", "y"]
    S, (u, v) = R["u", "v"]
    h = hom(S, S, a -> a^2, gens(S))
    @test (@inferred h(x)) == x^2
    @test (@inferred h(x)) == x^2

    @test_throws ErrorException isinjective(h)
    @test_throws ErrorException issurjective(h)
    @test_throws ErrorException isbijective(h)
    @test_throws ErrorException kernel(h)
  end
end
