@testset "pullbacks of Cartier divisors" begin 
  P = projective_space(QQ, 2)
  SP = ambient_coordinate_ring(P)
  (x, y, z) = gens(SP)

  A = [1 2 4; 1 3 9; 1 5 25]
  v = A*[x, y, z]
  psi = hom(SP, SP, v)
  g = ProjectiveSchemeMor(P, P, psi)
  g_cov = covered_scheme_morphism(g)

  X = covered_scheme(P)

  D = IdDict{AbsSpec, RingElem}()
  H = x^2 + y^2 + z^2
  for U in affine_charts(X)
    D[U] = dehomogenize(P, U)(H)
  end

  C = oscar.EffectiveCartierDivisor(X, D)
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
  K = keys(glueings(domain(double_ref)))
  for k in K
      @test underlying_glueing(glueings(domain(double_ref))[k]) isa SimpleGlueing
  end

  D2 = pullback(g_cov)(C) # The same thing
  @test trivializing_covering(D) === trivializing_covering(D2)
  for U in patches(trivializing_covering(D))
    @test D(U) == D2(U) # In particular, this will do a parent check
  end
end

@testset "Arithmetic of CartierDivisors" begin
  IP = projective_space(QQ, 3)
  (x,y,z,w) = gens(ambient_coordinate_ring(IP))
  f = x^4 + y^4 + z^4 + w^4
  IPX = subscheme(IP, f)
  X = covered_scheme(IPX)
  h = (x+y+z+w)^3
  u = (x-y+4*w)
  C = oscar.cartier_divisor(IPX, h)
  D = oscar.cartier_divisor(IPX, u)
  @test 2*C == C+C
  @test iszero(D-D)
  @test 3*(C + D) == 3*C + 3*D
end
