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

    if K == Qsqrt2
      h = hom(Qsqrt2, Qsqrt2, -gen(Qsqrt2))
      a = S([2 0; 0 2])
      b = S([4 0; 0 4])
      h = hom(Q, S, h, [a, b])
      # inferred is broken because of AA bug
      @test h(x * y) == a * b
      @test h(x * y) == a * b
      @test h(gen(Qsqrt2) * x * y) == - gen(Qsqrt2) * a * b
    end

    # The following is a copy of the MPolyRing tests with trivial I
    # So we need the zero-test in the codomain, which means this is a nono
    # for non-fields.
    if !(K isa Oscar.Field)
      continue
    end
    Kx, (x, y) = K["x", "y"]
    I = ideal(Kx, elem_type(Kx)[])
    Kx, = quo(Kx, I)
    x = Kx(x)
    y = Kx(y)
    h = @inferred hom(Kx, Kx, [y, x])
    @test h isa Oscar.morphism_type(Kx, Kx)
    @test h isa Oscar.morphism_type(typeof(Kx), typeof(Kx))
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
    @test h isa Oscar.morphism_type(Kx, Kx)
    @test (@inferred h(x + y)) == Kx(2)
    @test (@inferred h(x + y)) == Kx(2)

    Kh, (z1, z2) = K["z1", "z2"]
    Kh, (zz1, zz2) = grade(Kh)
    h = @inferred hom(Kx, Kh, [z1 + z2, z2])
    @test h isa Oscar.morphism_type(Kx, Kh)
    @test h isa Oscar.morphism_type(typeof(Kx), typeof(Kh))
    @test (@inferred h(x)) == zz1 + zz2
    @test (@inferred h(x)) == zz1 + zz2
    @test (@inferred h(y)) == zz2
    @test (@inferred h(y)) == zz2
    @test_throws ArgumentError hom(Kx, Kh, [zz1 + zz1^2, zz2])
    h = @inferred hom(Kx, Kh, [zz1 + zz1^2, zz2], check = false)
    @test h isa Oscar.morphism_type(Kx, Kh)
    @test h isa Oscar.morphism_type(typeof(Kx), typeof(Kh))
  
    # julia-function  
    R, (x, y) = K["x", "y"]
    S, (u, v) = R["u", "v"]
    g = let S = S; a -> S(a^2); end
    h = hom(S, S, g, gens(S))
    @test h isa Oscar.morphism_type(S, S, g)
    @test h isa Oscar.morphism_type(typeof(S), typeof(S), typeof(g))
    @test (@inferred h(u)) == u
    @test (@inferred h(u)) == u
    @test (@inferred h(x*u)) == x^2 * u
    @test (@inferred h(x*u)) == x^2 * u
    g = let S = S; a -> S(a^2); end
    h = hom(S, S, g, gens(S))
    @test h isa Oscar.morphism_type(S, S, g)
    @test h isa Oscar.morphism_type(typeof(S), typeof(S), typeof(g))
    @test (h(u)) == u
    @test (h(u)) == u
    @test (h(x*u)) == x^2 * u
    @test (h(x*u)) == x^2 * u
    g = a -> a^2
    h = hom(S, S, g, gens(S))
    @test h isa Oscar.morphism_type(S, S, g)
    @test h isa Oscar.morphism_type(typeof(S), typeof(S), typeof(g))
    @test (@inferred h(u)) == u
    @test (@inferred h(u)) == u
    @test (@inferred h(x*u)) == x^2 * u
    @test (@inferred h(x*u)) == x^2 * u
  
    A, (x,y) = K["x", "y"]
    f = hom(A, A, [2*x, 5*y])
    @test f isa Oscar.morphism_type(A, A)
    @test f isa Oscar.morphism_type(typeof(A), typeof(A))
    R, (u, v) = A["u", "v"]
    h = hom(R, R, f, [u+v, u*v])
    @test h isa Oscar.morphism_type(R, R, f)
    @test h isa Oscar.morphism_type(typeof(R), typeof(R), typeof(f))
    @test (@inferred h(x*u)) == 2*x*u + 2*x*v
    @test (@inferred h(x*u)) == 2*x*u + 2*x*v
  end

  Qi, i = quadratic_field(-1)
  Qix, (x, y) = Qi["x", "y"]
  I = ideal(Qix, elem_type(Qix)[])
  Qix, = quo(Qix, I)
  x = Qix(x)
  y = Qix(y)
 
  # composition with coefficient map
  f = hom(Qix, Qix, hom(Qi, Qi, -i), [x^2, y^2])
  g = hom(Qix, Qix, hom(Qi, Qi, -i), [x + 1, y + 1])
  fg = @inferred f * g
  @test fg(i) == g(f(i))
  @test fg(x) == g(f(x))
  @test fg(y) == g(f(y))

  f = hom(Qix, Qix, z -> z + 1, [x^2, y^2])
  g = hom(Qix, Qix, z -> z^2, [x + 1, y + 1])
  fg = @inferred f * g
  @test fg(i) == g(f(i))
  @test fg(x) == g(f(x))
  @test fg(y) == g(f(y))

  # composition with arbitrary maps
  h = hom(Qi, Qi, -i)
  f = hom(Qix, Qi, [i, 0])
  fh = @inferred f*h
  @test fh(x) == h(f(x))
  f = hom(Qix, Qi, h, [i, 0])
  fh = @inferred f * h
  @test fh(x) == h(f(x)) 
  f = hom(Qix, Qi, x -> x, [i, 0])
  @test fh(x) == h(f(x)) 
  f = hom(Qix, Qi, x -> x, [i, 0])

  # Construct stacked domain
  R, (x, y) = PolynomialRing(QQ, ["x", "y"])
  S, (u, v) = PolynomialRing(QQ, ["u", "v"])
  I = ideal(S, [ u - v^2 ])
  Q, StoQ = quo(S, I)
  QtoR = hom(Q, R, [ x^2, x ])
  T, (a, b, c) = PolynomialRing(Q, [ "a", "b", "c" ])
  J = ideal(T, [ a*b - c^2 ])
  A, TtoA = quo(T, J)
  # The test is whether the following two lines work at all
  AtoR = @inferred hom(A, R, QtoR, [ x^2, y^2, x*y ])
  @test isone(AtoR(A(1)))
end
