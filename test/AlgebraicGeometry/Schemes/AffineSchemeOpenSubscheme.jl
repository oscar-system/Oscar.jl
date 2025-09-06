@testset "AffineSchemeOpenSubscheme1" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  f = x^2 + y^2 + z^2
  A = spec(R)
  X = spec(quo(R, ideal(R, [x*f, y*f]))[1])
  @test is_subscheme(X, A)

  U = AffineSchemeOpenSubscheme(A, [x, y])
  UX = intersect(X, U)

  Y = closure(UX)
  Z = closure(UX, A)
  @test Z==Y

  @test isempty(complement(Y,Y))
end

@testset "AffineSchemeOpenSubscheme2" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  X = spec(R)
  N = subscheme(X, [x,y])
  U = complement(X, N)
  V = AffineSchemeOpenSubscheme(X, [x^2, y-x])
  I = ideal(R, [x*y-z^2])
  Y = subscheme(X, I)
  Xx = hypersurface_complement(X, x)
  Yx = hypersurface_complement(Y, x)
  @test !is_closed_embedding(Y, Yx)
  U = AffineSchemeOpenSubscheme(Y, [x^5, y-7*x^7])
  f = maximal_extension(Y, x//z)
  intersections(U)
  set_name!(U, "U")
  @test name(U) == "U"

  S, (u, v) = polynomial_ring(QQ, [:u, :v])
  Z = spec(S)
  g = maximal_extension(Y, Z, [x//z, y])
  g = restrict(g, domain(g), AffineSchemeOpenSubscheme(Z, [v]))

  R, (x,y,z) = QQ[:x, :y, :z]
  X = hypersurface_complement(spec(R),x)
  F = fraction_field(R)
  f = x//y
  g = F(x,x)
  D = domain(maximal_extension(X, f))
  @test D == hypersurface_complement(X,y)
  @test X ==domain(maximal_extension(X, g))
end

@testset "AffineSchemeOpenSubscheme3" begin
  R, (x,y,u,v) = QQ[:x, :y, :u, :v]
  A = spec(R)
  O = subscheme(A, [x,y,u,v])
  U = complement(A, O)
  @test A==closure(U)
  Oscar.AffineSchemeOpenSubschemeRing(U)
  OU = OO(U)
  gens(OO(U))
  f = OO(U)(x)
  @test Oscar.scheme(f) === A
  @test affine_patches(f) == affine_patches(U)
  @test n_patches(f) == 4
  @test f == deepcopy(f)
  @test f^2 == f*f
  @test f^2 == f^ZZRingElem(2)
  @test isone(divexact(f,f))
  @test is_unit(one(OU))
  @test !is_unit(OU(x))
  g = OU(2)
  @test isone(g*inv(g))
  @test isone(OU(1))
  @test iszero(f-f)
  Axy = hypersurface_complement(A,x*y)
  V = domain(maximal_extension(Axy, 1//x))
  @test Axy == V

  M = matrix(R, [x u; y v])
  h = det(M)
  X = subscheme(A, h)
  XU = intersect(X, U)
  S, (a,b) = QQ[:a, :b]
  B = spec(S)
  maximal_extension(X, x//u)
  maximal_extension(A, x//u)
  maximal_extension(X, [x//u])
  f = maximal_extension(X, B, [x//y, x//u])
end

@testset "AffineSchemeOpenSubscheme4" begin
  R, (x,y) = QQ[:x, :y]
  S, (t) = polynomial_ring(QQ, [:t])
  t = t[1]

  A = spec(R)
  B = spec(S)
  f = x^2+x^3 -y^2
  C = subscheme(A, f)

  phi = maximal_extension(C, B, [x//y])
  psi = maximal_extension(B, C, [(1-t^2)//t^2, (1-t^2)//t^3])
  psi = restrict(psi, AffineSchemeOpenSubscheme(B, [t*(1-t^2)]), domain(phi))
  phi_res = restrict(phi, domain(phi), domain(psi))
  #psi_res = restrict(psi, domain(psi), domain(phi))
  p = compose(psi, phi)
  q = compose(phi_res, psi)
  phi1 = restrict(phi,domain(phi),codomain(phi))
  @test phi == phi1
  @test restrict(p, domain(p), domain(p)) == id_hom(domain(p))
  @test restrict(q, domain(q), domain(q)) == id_hom(domain(q))
end

@testset "restriction_maps" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  X = spec(R, ideal(R, x*y))
  U = AffineSchemeOpenSubscheme(X, [x^7, y^8])
  f = OO(U)([7*x//x^2, 5*y//y^2])
  h = 5*x - 7*y+ 98*x*y^4
  V = hypersurface_complement(X, h)
  pb = restriction_map(U, V, h)
  result = pb(f)
  LX = U[1]
  LY = U[2]
  @test f[LX] == OO(LX)(result) && f[U[2]] == OO(U[2])(result)

  X = spec(R, ideal(R, x*(x+1)-y*z))
  U = AffineSchemeOpenSubscheme(X, [y^9, (x+1)^3])
  f = OO(U)([x//y, z//(x+1)])
  h = 3*y + x+1
  V = hypersurface_complement(X, h)
  pb = restriction_map(U, V, h)
  result = pb(f)
  @test f[U[1]] == OO(U[1])(result) && f[U[2]] == OO(U[2])(result)

  V = AffineSchemeOpenSubscheme(X, [((x+1)-y)^7, ((x+1)+y)^5])
  resUV = restriction_map(U, V)
  resVU = restriction_map(V, U)
  resVU(resUV(f)) == f
  @test Oscar.is_identity_map(compose(resVU, resUV))
  @test Oscar.is_identity_map(compose(resUV, resVU))
end
@testset "pullbacks" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  X = spec(R, ideal(R, [x^2-y*z]))
  U = AffineSchemeOpenSubscheme(X, [x, y])
  V = AffineSchemeOpenSubscheme(X, [(x+y)^2, x^2 - y^2, (x-y)^2])
  f = AffineSchemeOpenSubschemeMor(U, V, [x, y, z])
  f[affine_patches(domain(f))[1]]
  pbf = pullback(f)
  @test pbf(OO(V)(x)) == OO(V)(x)
  @test pbf(OO(V)(y)) == OO(V)(y)
  @test pbf(OO(V)(z)) == OO(V)(z)
  g = AffineSchemeOpenSubschemeMor(V, U, [x, y, z])
  @test Oscar.is_identity_map(compose(pullback(f), pullback(g)))
  @test Oscar.is_identity_map(pullback(compose(f,g)))
  W = subscheme(X, x)
  preimage(f, W)
  @test U ==preimage(f, V)


  S, (u, v) = QQ[:u, :v]
  B = spec(S)
  p = morphism(X, B, [x, y])
  U = AffineSchemeOpenSubscheme(X, [x, y])
  V = AffineSchemeOpenSubscheme(B, [u, v])
  pres = restrict(p, U, V)
  @test U == preimage(p, V)

  g = pullback(pres)
  pullback(pres, u)
  @test g(OO(V)(u)) == OO(U)(x)
  @test g(OO(V)(v)) == OO(U)(y)
end

@testset "AffineSchemeOpenSubschemeRings" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  I = ideal(R, [x^2-y*z])
  X = spec(R, I)
  U = AffineSchemeOpenSubscheme(X, [x, y])
  subscheme(U, ideal(R,[x,y]))
  V = AffineSchemeOpenSubscheme(X, [x+y, x^2 - y^2, (x-y)^5])
  f = canonical_isomorphism(OO(U), OO(V))
  a = OO(U)([x//y, z//x])
  b = f(a)
  @test preimage(f, b) == a
end

@testset "restriction map" begin
  R, (x, y) = QQ[:x, :y, :z]
  X = spec(R)
  U = AffineSchemeOpenSubscheme(X,[x])
  restriction_map(U, hypersurface_complement(X,x))
end


@testset "is_subscheme" begin
  R, (x, y) = QQ[:x, :y, :z]
  # AffineScheme{..., MPolyRing}
  A2 = spec(R)
  # AffineScheme{..., MPolyQuoRing}
  Vxy = subscheme(A2, x*y)
  # AffineScheme{MPolyLocRing}
  A2_minus_Vxy = hypersurface_complement(A2, x*y)
  A2_minus_Vxy_p = PrincipalOpenSubset(A2, x*y)
  # AffineSchemeOpenSubscheme{AffineScheme{...,MPolyRing},...}
  A2_minus_origin = AffineSchemeOpenSubscheme(A2, [x,y])
  # AffineSchemeOpenSubscheme{AffineScheme{MPolyQuoRing},...}
  Vxy_minus_origin = AffineSchemeOpenSubscheme(Vxy, [x,y])
  # PrincipalOpenSubset{MPolyQuoLocalized....}
  Vxy_minus_origin_p = PrincipalOpenSubset(Vxy, OO(Vxy)(x+y))

  @test is_subscheme(A2, A2)
  @test is_subscheme(Vxy, A2)
  @test is_subscheme(A2_minus_origin, A2)
  @test is_subscheme(A2_minus_Vxy, A2)
  @test is_subscheme(A2_minus_Vxy_p, A2)
  @test is_subscheme(Vxy_minus_origin, A2)
  @test is_subscheme(Vxy_minus_origin_p, A2)

  @test is_subscheme(Vxy, Vxy)
  @test !is_subscheme(A2_minus_origin, Vxy)
  @test !is_subscheme(A2_minus_Vxy, Vxy)
  @test !is_subscheme(A2_minus_Vxy_p, Vxy)
  @test is_subscheme(Vxy_minus_origin, Vxy)
  @test is_subscheme(Vxy_minus_origin_p, Vxy)

  @test !is_subscheme(A2, A2_minus_origin)
  @test !is_subscheme(Vxy, A2_minus_origin)
  @test is_subscheme(A2_minus_origin, A2_minus_origin)
  @test is_subscheme(A2_minus_Vxy, A2_minus_origin)
  @test is_subscheme(A2_minus_Vxy_p, A2_minus_origin)
  @test is_subscheme(Vxy_minus_origin, A2_minus_origin)
  @test is_subscheme(Vxy_minus_origin_p, A2_minus_origin)

  @test !is_subscheme(A2, A2_minus_Vxy)
  @test !is_subscheme(Vxy, A2_minus_Vxy)
  @test !is_subscheme(A2_minus_origin, A2_minus_Vxy)
  @test is_subscheme(A2_minus_Vxy, A2_minus_Vxy)
  @test is_subscheme(A2_minus_Vxy_p, A2_minus_Vxy)
  @test !is_subscheme(Vxy_minus_origin, A2_minus_Vxy)
  @test !is_subscheme(Vxy_minus_origin_p, A2_minus_Vxy)

  @test !is_subscheme(A2, A2_minus_Vxy_p)
  @test !is_subscheme(Vxy, A2_minus_Vxy_p)
  @test !is_subscheme(A2_minus_origin, A2_minus_Vxy_p)
  @test is_subscheme(A2_minus_Vxy, A2_minus_Vxy_p)
  @test is_subscheme(A2_minus_Vxy_p, A2_minus_Vxy_p)
  @test !is_subscheme(Vxy_minus_origin, A2_minus_Vxy_p)
  @test !is_subscheme(Vxy_minus_origin_p, A2_minus_Vxy_p)

  @test !is_subscheme(A2, Vxy_minus_origin)
  @test !is_subscheme(Vxy, Vxy_minus_origin)
  @test_broken !is_subscheme(A2_minus_origin, Vxy_minus_origin)
  @test !is_subscheme(A2_minus_Vxy, Vxy_minus_origin)
  @test !is_subscheme(A2_minus_Vxy_p, Vxy_minus_origin)
  @test is_subscheme(Vxy_minus_origin, Vxy_minus_origin)
  @test is_subscheme(Vxy_minus_origin_p, Vxy_minus_origin)

  @test !is_subscheme(A2, Vxy_minus_origin_p)
  @test !is_subscheme(Vxy, Vxy_minus_origin_p)
  @test !is_subscheme(A2_minus_origin, Vxy_minus_origin_p)
  @test !is_subscheme(A2_minus_Vxy, Vxy_minus_origin_p)
  @test !is_subscheme(A2_minus_Vxy_p, Vxy_minus_origin_p)
  @test is_subscheme(Vxy_minus_origin, Vxy_minus_origin_p)
  @test is_subscheme(Vxy_minus_origin_p, Vxy_minus_origin_p)
end

@testset "is_open_embedding" begin
  R, (x, y) = QQ[:x, :y, :z]
  # AffineScheme{..., MPolyRing}
  A2 = spec(R)
  # AffineScheme{..., MPolyQuoRing}
  Vxy = subscheme(A2, x*y)
  # AffineScheme{MPolyLocRing}
  A2_minus_Vxy = hypersurface_complement(A2, x*y)
  A2_minus_Vxy_p = PrincipalOpenSubset(A2, x*y)
  # AffineSchemeOpenSubscheme{AffineScheme{...,MPolyRing},...}
  A2_minus_origin = AffineSchemeOpenSubscheme(A2, [x,y])
  # AffineSchemeOpenSubscheme{AffineScheme{MPolyQuoRing},...}
  Vxy_minus_origin = AffineSchemeOpenSubscheme(Vxy, [x,y])
  # PrincipalOpenSubset{MPolyQuoLocalized....}
  Vxy_minus_origin_p = PrincipalOpenSubset(Vxy, OO(Vxy)(x+y))

  @test is_open_embedding(A2, A2)
  @test !is_open_embedding(Vxy, A2)
  @test_broken is_open_embedding(A2_minus_origin, A2)
  @test is_open_embedding(A2_minus_Vxy, A2)
  @test is_open_embedding(A2_minus_Vxy_p, A2)
  @test_broken !is_open_embedding(Vxy_minus_origin, A2)
  @test !is_open_embedding(Vxy_minus_origin_p, A2)

  @test !is_open_embedding(A2, Vxy)
  @test is_open_embedding(Vxy, Vxy)
  @test_broken !is_open_embedding(A2_minus_origin, Vxy)
  @test !is_open_embedding(A2_minus_Vxy, Vxy)
  @test !is_open_embedding(A2_minus_Vxy_p, Vxy)
  @test_broken is_open_embedding(Vxy_minus_origin, Vxy)
  @test is_open_embedding(Vxy_minus_origin_p, Vxy)

  @test_broken !is_open_embedding(A2, A2_minus_origin)
  @test_broken !is_open_embedding(Vxy, A2_minus_origin)
  @test_broken is_open_embedding(A2_minus_origin, A2_minus_origin)
  @test_broken is_open_embedding(A2_minus_Vxy, A2_minus_origin)
  @test_broken is_open_embedding(A2_minus_Vxy_p, A2_minus_origin)
  @test_broken !is_open_embedding(Vxy_minus_origin, A2_minus_origin)
  @test_broken !is_open_embedding(Vxy_minus_origin_p, A2_minus_origin)

  @test !is_open_embedding(A2, A2_minus_Vxy)
  @test !is_open_embedding(Vxy, A2_minus_Vxy)
  #@test !is_open_embedding(A2_minus_origin, A2_minus_Vxy)
  @test is_open_embedding(A2_minus_Vxy, A2_minus_Vxy)
  @test is_open_embedding(A2_minus_Vxy_p, A2_minus_Vxy)
  @test_broken !is_open_embedding(Vxy_minus_origin, A2_minus_Vxy)
  @test !is_open_embedding(Vxy_minus_origin_p, A2_minus_Vxy)

  @test !is_open_embedding(A2, A2_minus_Vxy_p)
  @test !is_open_embedding(Vxy, A2_minus_Vxy_p)
  @test_broken !is_open_embedding(A2_minus_origin, A2_minus_Vxy_p)
  @test is_open_embedding(A2_minus_Vxy, A2_minus_Vxy_p)
  @test is_open_embedding(A2_minus_Vxy_p, A2_minus_Vxy_p)
  @test_broken !is_open_embedding(Vxy_minus_origin, A2_minus_Vxy_p)
  @test !is_open_embedding(Vxy_minus_origin_p, A2_minus_Vxy_p)

  @test_broken !is_open_embedding(A2, Vxy_minus_origin)
  @test_broken !is_open_embedding(Vxy, Vxy_minus_origin)
  @test_broken !is_open_embedding(A2_minus_origin, Vxy_minus_origin)
  @test_broken !is_open_embedding(A2_minus_Vxy, Vxy_minus_origin)
  @test_broken !is_open_embedding(A2_minus_Vxy_p, Vxy_minus_origin)
  @test_broken is_open_embedding(Vxy_minus_origin, Vxy_minus_origin)
  @test_broken is_open_embedding(Vxy_minus_origin_p, Vxy_minus_origin)

  @test !is_open_embedding(A2, Vxy_minus_origin_p)
  @test !is_open_embedding(Vxy, Vxy_minus_origin_p)
  @test_broken !is_open_embedding(A2_minus_origin, Vxy_minus_origin_p)
  @test !is_open_embedding(A2_minus_Vxy, Vxy_minus_origin_p)
  @test !is_open_embedding(A2_minus_Vxy_p, Vxy_minus_origin_p)
  @test_broken is_open_embedding(Vxy_minus_origin, Vxy_minus_origin_p)
  @test is_open_embedding(Vxy_minus_origin_p, Vxy_minus_origin_p)

end

@testset "is_closed_embedding" begin
  R, (x, y) = QQ[:x, :y, :z]
  # AffineScheme{..., MPolyRing}
  A2 = spec(R)
  # AffineScheme{..., MPolyQuoRing}
  Vxy = subscheme(A2, x*y)
  # AffineScheme{MPolyLocRing}
  A2_minus_Vxy = hypersurface_complement(A2, x*y)
  A2_minus_Vxy_p = PrincipalOpenSubset(A2, x*y)
  # AffineSchemeOpenSubscheme{AffineScheme{...,MPolyRing},...}
  A2_minus_origin = AffineSchemeOpenSubscheme(A2, [x,y])
  # AffineSchemeOpenSubscheme{AffineScheme{MPolyQuoRing},...}
  Vxy_minus_origin = AffineSchemeOpenSubscheme(Vxy, [x,y])
  # PrincipalOpenSubset{MPolyQuoLocalized....}
  Vxy_minus_origin_p = PrincipalOpenSubset(Vxy, OO(Vxy)(x+y))

  @test is_closed_embedding(A2, A2)
  @test is_closed_embedding(Vxy, A2)
  @test_broken !is_closed_embedding(A2_minus_origin, A2)
  @test !is_closed_embedding(A2_minus_Vxy, A2)
  @test !is_closed_embedding(A2_minus_Vxy_p, A2)
  @test_broken !is_closed_embedding(Vxy_minus_origin, A2)
  @test_broken !is_closed_embedding(Vxy_minus_origin_p, A2)

  @test !is_closed_embedding(A2, Vxy)
  @test is_closed_embedding(Vxy, Vxy)
  @test_broken !is_closed_embedding(A2_minus_origin, Vxy)
  @test_broken !is_closed_embedding(A2_minus_Vxy, Vxy)
  @test_broken !is_closed_embedding(A2_minus_Vxy_p, Vxy)
  @test_broken !is_closed_embedding(Vxy_minus_origin, Vxy)
  @test !is_closed_embedding(Vxy_minus_origin_p, Vxy)

  @test_broken !is_closed_embedding(A2, A2_minus_origin)
  @test_broken is_closed_embedding(Vxy, A2_minus_origin)
  @test_broken is_closed_embedding(A2_minus_origin, A2_minus_origin)
  @test_broken !is_closed_embedding(A2_minus_Vxy, A2_minus_origin)
  @test_broken !is_closed_embedding(A2_minus_Vxy_p, A2_minus_origin)
  @test_broken is_closed_embedding(Vxy_minus_origin, A2_minus_origin)
  @test_broken is_closed_embedding(Vxy_minus_origin_p, A2_minus_origin)

  @test_broken !is_closed_embedding(A2, A2_minus_Vxy)
  @test_broken !is_closed_embedding(Vxy, A2_minus_Vxy)
  @test_broken !is_closed_embedding(A2_minus_origin, A2_minus_Vxy)
  @test_broken is_closed_embedding(A2_minus_Vxy, A2_minus_Vxy)
  @test_broken is_closed_embedding(A2_minus_Vxy_p, A2_minus_Vxy)
  @test_broken !is_closed_embedding(Vxy_minus_origin, A2_minus_Vxy)
  @test_broken !is_closed_embedding(Vxy_minus_origin_p, A2_minus_Vxy)

  @test_broken !is_closed_embedding(A2, A2_minus_Vxy_p)
  @test_broken !is_closed_embedding(Vxy, A2_minus_Vxy_p)
  @test_broken !is_closed_embedding(A2_minus_origin, A2_minus_Vxy_p)
  @test_broken is_closed_embedding(A2_minus_Vxy, A2_minus_Vxy_p)
  @test_broken is_closed_embedding(A2_minus_Vxy_p, A2_minus_Vxy_p)
  @test_broken !is_closed_embedding(Vxy_minus_origin, A2_minus_Vxy_p)
  @test_broken !is_closed_embedding(Vxy_minus_origin_p, A2_minus_Vxy_p)

  @test_broken !is_closed_embedding(A2, Vxy_minus_origin)
  @test_broken !is_closed_embedding(Vxy, Vxy_minus_origin)
  @test_broken !is_closed_embedding(A2_minus_origin, Vxy_minus_origin)
  @test_broken !is_closed_embedding(A2_minus_Vxy, Vxy_minus_origin)
  @test_broken !is_closed_embedding(A2_minus_Vxy_p, Vxy_minus_origin)
  @test_broken is_closed_embedding(Vxy_minus_origin, Vxy_minus_origin)
  @test_broken is_closed_embedding(Vxy_minus_origin_p, Vxy_minus_origin)

  @test_broken !is_closed_embedding(A2, Vxy_minus_origin_p)
  @test_broken !is_closed_embedding(Vxy, Vxy_minus_origin_p)
  @test_broken !is_closed_embedding(A2_minus_origin, Vxy_minus_origin_p)
  @test_broken !is_closed_embedding(A2_minus_Vxy, Vxy_minus_origin_p)
  @test_broken !is_closed_embedding(A2_minus_Vxy_p, Vxy_minus_origin_p)
  @test_broken is_closed_embedding(Vxy_minus_origin, Vxy_minus_origin_p)
  @test is_closed_embedding(Vxy_minus_origin_p, Vxy_minus_origin_p)

end

@testset "Issue 2129" begin
  R, (x,y) = QQ[:x, :y]
  I = ideal(R, x)
  Q, _ = quo(R, I)
  X = spec(Q)
  U = AffineSchemeOpenSubscheme(X, [y])
  W = OO(U)

  P, (u,v) = W[:u, :v]
  PQ, _ = quo(P, ideal(P, [u]))

  @test parent(evaluate(v, gens(PQ))) === PQ

  PP, (uu, vv) = P[:uu, :vv]
  @test parent(one(W) * uu) === PP
end

