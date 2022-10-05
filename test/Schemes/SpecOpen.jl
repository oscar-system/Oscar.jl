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
  @test is_canonically_isomorphic(Z, Y)

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
end

@testset "SpecOpen_3" begin
  R, (x,y,u,v) = QQ["x", "y", "u","v"]
  A = Spec(R)
  O = subscheme(A, [x,y,u,v])
  U = complement(A, O)
  @test is_canonically_isomorphic(A, closure(U))
  Oscar.SpecOpenRing(U)
  OU = OO(U)
  @test_throws ErrorException gens(OO(U))
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
  @test f == OU(OO(affine_patches(f)[1])(f))
  Axy = hypersurface_complement(A,x*y)
  V = domain(maximal_extension(Axy, 1//x))
  @test is_canonically_isomorphic(Axy, V)

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
  @test is_canonically_isomorphic(U,preimage(f, V))


  S, (u, v) = QQ["u", "v"]
  B = Spec(S)
  p = SpecMor(X, B, [x, y])
  U = SpecOpen(X, [x, y])
  V = SpecOpen(B, [u, v])
  pres = restrict(p, U, V)
  @test is_canonically_isomorphic(U, preimage(p, V))

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
