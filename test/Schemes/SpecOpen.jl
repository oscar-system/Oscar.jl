@testset "SpecOpen_1" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  f = x^2 + y^2 + z^2
  A = Spec(R)
  X = Spec(quo(R, ideal(R, [x*f, y*f]))[1])
  @test issubset(X, A)

  U = SpecOpen(A, [x, y])
  UX = intersect(X, U)

  Y = closure(UX)
  Z = closure(UX, A)
  @test Z==Y

  @test isempty(complement(Y,Y))
end

@testset "SpecOpen_2" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  X = Spec(R)
  N = subscheme(X, [x,y])
  U = complement(X, N)
  V = SpecOpen(X, [x^2, y-x])
  I = ideal(R, [x*y-z^2])
  Y = subscheme(X, I)
  Xx = hypersurface_complement(X, x)
  Yx = hypersurface_complement(Y, x)
  @test !is_closed_embedding(Y, Yx)
  U = SpecOpen(Y, [x^5, y-7*x^7])
  f = maximal_extension(Y, x//z)
  intersections(U)
  set_name!(U, "U")
  @test name(U) == "U"

  S, (u, v) = PolynomialRing(QQ, ["u", "v"])
  Z = Spec(S)
  g = maximal_extension(Y, Z, [x//z, y])
  g = restrict(g, domain(g), SpecOpen(Z, [v]))

  R, (x,y,z) = QQ["x", "y", "z"]
  X = hypersurface_complement(Spec(R),x)
  F = FractionField(R)
  f = x//y
  g = F(x,x)
  D = domain(maximal_extension(X, f))
  @test D == hypersurface_complement(X,y)
  @test X ==domain(maximal_extension(X, g))
end

@testset "SpecOpen_3" begin
  R, (x,y,u,v) = QQ["x", "y", "u","v"]
  A = Spec(R)
  O = subscheme(A, [x,y,u,v])
  U = complement(A, O)
  @test A==closure(U)
  Oscar.SpecOpenRing(U)
  OU = OO(U)
  gens(OO(U))
  f = OO(U)(x)
  @test Oscar.scheme(f) === A
  @test affine_patches(f) == affine_patches(U)
  @test npatches(f) == 4
  @test f == deepcopy(f)
  @test f^2 == f*f
  @test f^2 == f^fmpz(2)
  @test isone(divexact(f,f))
  @test isunit(one(OU))
  @test !isunit(OU(x))
  g = OU(2)
  @test isone(g*inv(g))
  @test isone(OU(1))
  @test iszero(f-f)
  Axy = hypersurface_complement(A,x*y)
  V = domain(maximal_extension(Axy, 1//x))
  @test Axy == V

  M = MatrixSpace(R,2,2)([x u; y v])
  h = det(M)
  X = subscheme(A, h)
  XU = intersect(X, U)
  S, (a,b) = QQ["a", "b"]
  B = Spec(S)
  maximal_extension(X, x//u)
  maximal_extension(A, x//u)
  maximal_extension(X, [x//u])
  f = maximal_extension(X, B, [x//y, x//u])
end

@testset "SpecOpen_4" begin
  R, (x,y) = QQ["x", "y"]
  S, (t) = PolynomialRing(QQ, ["t"])
  t = t[1]

  A = Spec(R)
  B = Spec(S)
  f = x^2+x^3 -y^2
  C = subscheme(A, f)

  phi = maximal_extension(C, B, [x//y])
  psi = maximal_extension(B, C, [(1-t^2)//t^2, (1-t^2)//t^3])
  psi = restrict(psi, SpecOpen(B, [t*(1-t^2)]), domain(phi))
  phi_res = restrict(phi, domain(phi), domain(psi))
  #psi_res = restrict(psi, domain(psi), domain(phi))
  p = compose(psi, phi)
  q = compose(phi_res, psi)
  phi1 = restrict(phi,domain(phi),codomain(phi))
  @test phi == phi1
  @test restrict(p, domain(p), domain(p)) == identity_map(domain(p))
  @test restrict(q, domain(q), domain(q)) == identity_map(domain(q))
end

@testset "restriction_maps" begin
  R, (x, y, z) = QQ["x", "y", "z"]
  X = Spec(R, ideal(R, x*y))
  U = SpecOpen(X, [x^7, y^8])
  f = OO(U)([7*x//x^2, 5*y//y^2])
  h = 5*x - 7*y+ 98*x*y^4
  V = hypersurface_complement(X, h)
  pb = restriction_map(U, V, h)
  result = pb(f)
  LX = U[1]
  LY = U[2]
  @test f[LX] == OO(LX)(result) && f[U[2]] == OO(U[2])(result)

  X = Spec(R, ideal(R, x*(x+1)-y*z))
  U = SpecOpen(X, [y^9, (x+1)^3])
  f = OO(U)([x//y, z//(x+1)])
  h = 3*y + x+1
  V = hypersurface_complement(X, h)
  pb = restriction_map(U, V, h)
  result = pb(f)
  @test f[U[1]] == OO(U[1])(result) && f[U[2]] == OO(U[2])(result)

  V = SpecOpen(X, [((x+1)-y)^7, ((x+1)+y)^5])
  resUV = restriction_map(U, V)
  resVU = restriction_map(V, U)
  resVU(resUV(f)) == f
  @test Oscar.is_identity_map(compose(resVU, resUV))
  @test Oscar.is_identity_map(compose(resUV, resVU))
end
@testset "pullbacks" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  X = Spec(R, ideal(R, [x^2-y*z]))
  U = SpecOpen(X, [x, y])
  V = SpecOpen(X, [(x+y)^2, x^2 - y^2, (x-y)^2])
  f = SpecOpenMor(U, V, [x, y, z])
  f[affine_patches(domain(f))[1]]
  pbf = pullback(f)
  @test pbf(OO(V)(x)) == OO(V)(x)
  @test pbf(OO(V)(y)) == OO(V)(y)
  @test pbf(OO(V)(z)) == OO(V)(z)
  g = SpecOpenMor(V, U, [x, y, z])
  @test Oscar.is_identity_map(compose(pullback(f), pullback(g)))
  @test Oscar.is_identity_map(pullback(compose(f,g)))
  W = subscheme(X, x)
  preimage(f, W)
  @test U ==preimage(f, V)


  S, (u, v) = QQ["u", "v"]
  B = Spec(S)
  p = SpecMor(X, B, [x, y])
  U = SpecOpen(X, [x, y])
  V = SpecOpen(B, [u, v])
  pres = restrict(p, U, V)
  @test U == preimage(p, V)

  g = pullback(pres)
  pullback(pres, u)
  @test g(OO(V)(u)) == OO(U)(x)
  @test g(OO(V)(v)) == OO(U)(y)
end

@testset "SpecOpenRings" begin
  R, (x, y, z) = QQ["x", "y", "z"]
  I = ideal(R, [x^2-y*z])
  X = Spec(R, I)
  U = SpecOpen(X, [x, y])
  subscheme(U, ideal(R,[x,y]))
  V = SpecOpen(X, [x+y, x^2 - y^2, (x-y)^5])
  f = canonical_isomorphism(OO(U), OO(V))
  a = OO(U)([x//y, z//x])
  b = f(a)
  @test preimage(f, b) == a
end

@testset "restriction map" begin
  R, (x, y) = QQ["x", "y", "z"]
  X = Spec(R)
  U = SpecOpen(X,[x])
  restriction_map(U, hypersurface_complement(X,x))
end


@testset "issubset" begin
  R, (x, y) = QQ["x", "y", "z"]
  # Spec{..., MPolyRing}
  A2 = Spec(R)
  # Spec{..., MPolyQuo}
  Vxy = subscheme(A2, x*y)
  # Spec{MPolyLocalizedRing}
  A2_minus_Vxy = hypersurface_complement(A2, x*y)
  A2_minus_Vxy_p = PrincipalOpenSubset(A2, x*y)
  # SpecOpen{Spec{...,MPolyRing},...}
  A2_minus_origin = SpecOpen(A2, [x,y])
  # SpecOpen{Spec{MPolyQuo},...}
  Vxy_minus_origin = SpecOpen(Vxy, [x,y])
  # PrincipalOpenSubset{MPolyQuoLocalized....}
  Vxy_minus_origin_p = PrincipalOpenSubset(Vxy, OO(Vxy)(x+y))

  @test issubset(A2, A2)
  @test issubset(Vxy, A2)
  @test issubset(A2_minus_origin, A2)
  @test issubset(A2_minus_Vxy, A2)
  @test issubset(A2_minus_Vxy_p, A2)
  @test issubset(Vxy_minus_origin, A2)
  @test issubset(Vxy_minus_origin_p, A2)

  @test issubset(Vxy, Vxy)
  @test !issubset(A2_minus_origin, Vxy)
  @test !issubset(A2_minus_Vxy, Vxy)
  @test !issubset(A2_minus_Vxy_p, Vxy)
  @test issubset(Vxy_minus_origin, Vxy)
  @test issubset(Vxy_minus_origin_p, Vxy)

  @test !issubset(A2, A2_minus_origin)
  @test !issubset(Vxy, A2_minus_origin)
  @test issubset(A2_minus_origin, A2_minus_origin)
  @test issubset(A2_minus_Vxy, A2_minus_origin)
  @test issubset(A2_minus_Vxy_p, A2_minus_origin)
  @test issubset(Vxy_minus_origin, A2_minus_origin)
  @test issubset(Vxy_minus_origin_p, A2_minus_origin)

  @test !issubset(A2, A2_minus_Vxy)
  @test !issubset(Vxy, A2_minus_Vxy)
  @test !issubset(A2_minus_origin, A2_minus_Vxy)
  @test issubset(A2_minus_Vxy, A2_minus_Vxy)
  @test issubset(A2_minus_Vxy_p, A2_minus_Vxy)
  @test !issubset(Vxy_minus_origin, A2_minus_Vxy)
  @test !issubset(Vxy_minus_origin_p, A2_minus_Vxy)

  @test !issubset(A2, A2_minus_Vxy_p)
  @test !issubset(Vxy, A2_minus_Vxy_p)
  @test !issubset(A2_minus_origin, A2_minus_Vxy_p)
  @test issubset(A2_minus_Vxy, A2_minus_Vxy_p)
  @test issubset(A2_minus_Vxy_p, A2_minus_Vxy_p)
  @test !issubset(Vxy_minus_origin, A2_minus_Vxy_p)
  @test !issubset(Vxy_minus_origin_p, A2_minus_Vxy_p)

  @test !issubset(A2, Vxy_minus_origin)
  @test !issubset(Vxy, Vxy_minus_origin)
  @test_broken !issubset(A2_minus_origin, Vxy_minus_origin)
  @test !issubset(A2_minus_Vxy, Vxy_minus_origin)
  @test !issubset(A2_minus_Vxy_p, Vxy_minus_origin)
  @test issubset(Vxy_minus_origin, Vxy_minus_origin)
  @test issubset(Vxy_minus_origin_p, Vxy_minus_origin)

  @test !issubset(A2, Vxy_minus_origin_p)
  @test !issubset(Vxy, Vxy_minus_origin_p)
  @test !issubset(A2_minus_origin, Vxy_minus_origin_p)
  @test !issubset(A2_minus_Vxy, Vxy_minus_origin_p)
  @test !issubset(A2_minus_Vxy_p, Vxy_minus_origin_p)
  @test issubset(Vxy_minus_origin, Vxy_minus_origin_p)
  @test issubset(Vxy_minus_origin_p, Vxy_minus_origin_p)
end

@testset "is_open_embedding" begin
  R, (x, y) = QQ["x", "y", "z"]
  # Spec{..., MPolyRing}
  A2 = Spec(R)
  # Spec{..., MPolyQuo}
  Vxy = subscheme(A2, x*y)
  # Spec{MPolyLocalizedRing}
  A2_minus_Vxy = hypersurface_complement(A2, x*y)
  A2_minus_Vxy_p = PrincipalOpenSubset(A2, x*y)
  # SpecOpen{Spec{...,MPolyRing},...}
  A2_minus_origin = SpecOpen(A2, [x,y])
  # SpecOpen{Spec{MPolyQuo},...}
  Vxy_minus_origin = SpecOpen(Vxy, [x,y])
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
  R, (x, y) = QQ["x", "y", "z"]
  # Spec{..., MPolyRing}
  A2 = Spec(R)
  # Spec{..., MPolyQuo}
  Vxy = subscheme(A2, x*y)
  # Spec{MPolyLocalizedRing}
  A2_minus_Vxy = hypersurface_complement(A2, x*y)
  A2_minus_Vxy_p = PrincipalOpenSubset(A2, x*y)
  # SpecOpen{Spec{...,MPolyRing},...}
  A2_minus_origin = SpecOpen(A2, [x,y])
  # SpecOpen{Spec{MPolyQuo},...}
  Vxy_minus_origin = SpecOpen(Vxy, [x,y])
  # PrincipalOpenSubset{MPolyQuoLocalized....}
  Vxy_minus_origin_p = PrincipalOpenSubset(Vxy, OO(Vxy)(x+y))

  @test is_closed_embedding(A2, A2)
  @test is_closed_embedding(Vxy, A2)s
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
