@testset "pullbacks of Cartier divisors" begin 
  P = projective_space(QQ, 2)
  SP = homogeneous_coordinate_ring(P)
  (x, y, z) = gens(SP)

  A = [1 2 4; 1 3 9; 1 5 25]
  v = A*[x, y, z]
  psi = hom(SP, SP, v)
  g = ProjectiveSchemeMor(P, P, psi)
  g_cov = covered_scheme_morphism(g)

  X = covered_scheme(P)

  D = IdDict{AbsAffineScheme, RingElem}()
  H = x^2 + y^2 + z^2
  for U in affine_charts(X)
    D[U] = dehomogenization_map(P, U)(H)
  end

  C = Oscar.EffectiveCartierDivisor(X, D)
  D = pullback(g_cov)(C)
  for U in patches(trivializing_covering(D))
    @test length(D(U)) == 1
  end

  # Test pullbacks of coverings and restrictions of morphisms to them
  phi = restrict(g_cov, default_covering(X))
  @test phi === restrict(g_cov, default_covering(X)) # test caching
  ref = domain(phi)
  double_ref = restrict(g_cov, ref)
  @test double_ref === restrict(g_cov, ref)
  K = keys(gluings(domain(double_ref)))
  for k in K
      @test underlying_gluing(gluings(domain(double_ref))[k]) isa SimpleGluing
  end

  D2 = pullback(g_cov)(C) # The same thing
  @test trivializing_covering(D) === trivializing_covering(D2)
  for U in patches(trivializing_covering(D))
    @test D(U) == D2(U) # In particular, this will do a parent check
  end
end

@testset "Arithmetic of CartierDivisors" begin
  IP = projective_space(QQ, 3)
  (x,y,z,w) = gens(homogeneous_coordinate_ring(IP))
  f = x^4 + y^4 + z^4 + w^4
  IPX = subscheme(IP, f)
  S = homogeneous_coordinate_ring(IPX)
  (x, y, z, w) = gens(S)
  X = covered_scheme(IPX)
  h = (x+y+z+w)^3
  u = (x-y+4*w)
  C = cartier_divisor(IPX, h)
  D = cartier_divisor(IPX, u)
  @test 2*C == C+C
  @test iszero(D-D)
  @test 3*(C + D) == 3*C + 3*D
  Z = Oscar.cartier_divisor(IPX, S(1))
  @test iszero(Z)
  C2 = cartier_divisor(IPX, h^2)
  @test iszero(2*C - C2)
end

@testset "conversion of Cartier to Weil divisors" begin
  IP2 = projective_space(QQ, [:x, :y, :z])
  S = homogeneous_coordinate_ring(IP2)
  (x,y,z) = gens(S)
  I = ideal(S, x^2*y^3*(x+y+z))
  C = Oscar.effective_cartier_divisor(IP2, gen(I, 1))
  CW = Oscar.weil_divisor(C)
  CW = Oscar.irreducible_decomposition(CW)
  @test length(components(CW)) == 3
  @test 1 in collect(values(CW.C.coefficients))
  @test 2 in collect(values(CW.C.coefficients))
  @test 3 in collect(values(CW.C.coefficients))
  @test 3*weil_divisor(C) == weil_divisor(3*C)
end
