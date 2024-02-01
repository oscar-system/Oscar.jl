@testset "MPolyAnyMap/MPolyRing" begin
  # Construction 
  Qsqrt2, = quadratic_field(-1)
  Zx, _ = ZZ["x"]
  Zxy, _ = ZZ["x", "y"]
  
  for K in [GF(2), GF(ZZRingElem(2)), GF(2, 2), GF(ZZRingElem(2), 2), ZZ, QQ, Qsqrt2, Zx, Zxy]
    Kx, (x, y) = K["x", "y"]
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

    h = @inferred hom(Kx, Kx, ZZRingElem[1, 1])
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

    # Noncommutative image
    
    A, (x, y) = K["x", "y"]
    S = matrix_ring(K, 2)
    a = S([1 1; 0 1])
    b = S([0 1; 1 0])
    @test_throws ArgumentError hom(A, S, [a, b])

    a = S([1 0; 0 2])
    b = S([2 0; 0 2])
    h = hom(A, S, [a, b])
    @test h isa Oscar.morphism_type(A, S)
    @test h isa Oscar.morphism_type(typeof(A), typeof(S))

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

    @test_throws ErrorException is_injective(h)
    @test_throws ErrorException is_surjective(h)
    @test_throws ErrorException is_bijective(h)
    # The following line works only in some cases. TODO: Adjust!
    #@test_throws ErrorException kernel(h)

    # composition
    Kx, (x, y) = K["x", "y"];
    Kxz, (z1, z2) = Kx["z1", "z2"];
    f = hom(Kxz, Kxz, hom(Kx, Kxz, [z1, z2]), [z2, z1])
    g = hom(Kxz, Kxz, hom(Kx, Kx, [y, x]), [z1 + 1, z2 + 1])
    fg = @inferred f * g
    @test fg(x) == g(f(x))
    @test fg(y) == g(f(y))
    @test fg(z1) == g(f(z1))
    @test fg(z2) == g(f(z2))

    f = hom(Kx, Kx, [y, x])
    g = hom(Kx, K, [K(1), K(1)])
    fg = @inferred f * g
    @test fg(x) == g(f(x))
    @test fg(y) == g(f(y))
  end

  # composition with coefficient map
  Qi, i = quadratic_field(-1)
  Qix, (x, y) = Qi["x", "y"]
  f = hom(Qix, Qix, hom(Qi, Qi, -i), [x^2, y^2])
  g = hom(Qix, Qix, hom(Qi, Qi, -i), [x + 1, y + 1])
  fg = @inferred f * g
  @test fg(i) == g(f(i))
  @test fg(x) == g(f(x))
  @test fg(y) == g(f(y))

  Qi, i = quadratic_field(-1)
  Qix, (x, y) = Qi["x", "y"]
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
end

@testset "coefficient maps" begin
  R, (x, y) = QQ[:x, :y]
  L, _ = QQ[:t]
  phi = hom(R, R, x->zero(L), [x, y])
  @test_throws ErrorException compose(phi, phi)
end
