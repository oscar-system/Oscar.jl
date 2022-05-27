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
end

@testset "SpecOpen_2" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  X = Spec(R)
  N = subscheme(X, [x,y])
  U = complement(X, N)
  V = SpecOpen(X, [x^2, y-x])
  I = ideal(R, [x*y-z^2])
  Y = subscheme(X, I)
  U = SpecOpen(Y, [x^5, y-7*x^7])
  f = maximal_extension(Y, x//z)

  S, (u, v) = PolynomialRing(QQ, ["u", "v"])
  Z = Spec(S)
  g = maximal_extension(Y, Z, [x//z, y])
  g = restriction(g, domain(g), SpecOpen(Z, [v]))
end

@testset "SpecOpen_3" begin
  R, (x,y,u,v) = QQ["x", "y", "u","v"]
  A = Spec(R)
  O = subscheme(A, [x,y,u,v])
  U = complement(A, O)
  M = MatrixSpace(R,2,2)([x u; y v])
  h = det(M)
  X = subscheme(A, h)
  XU = intersect(X, U)
  S, (a,b) = QQ["a", "b"]
  B = Spec(S)
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
  psi = restriction(psi, SpecOpen(B, [t*(1-t^2)]), domain(phi))
  phi_res = restriction(phi, domain(phi), domain(psi))
  #psi_res = restriction(psi, domain(psi), domain(phi))
  p = compose(psi, phi)
  q = compose(phi_res, psi)
  @test restriction(p, domain(p), domain(p)) == identity_map(domain(p))
  @test restriction(q, domain(q), domain(q)) == identity_map(domain(q))
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
  pbf = pullback(f)
  @test pbf(OO(V)(x)) == OO(V)(x)
  @test pbf(OO(V)(y)) == OO(V)(y)
  @test pbf(OO(V)(z)) == OO(V)(z)
  g = SpecOpenMor(V, U, [x, y, z])
  @test Oscar.is_identity_map(compose(pullback(f), pullback(g)))
  @test Oscar.is_identity_map(pullback(compose(f,g)))


  S, (u, v) = QQ["u", "v"]
  B = Spec(S)
  p = SpecMor(X, B, [x, y])
  U = SpecOpen(X, [x, y])
  V = SpecOpen(B, [u, v])
  pres = restrict(p, U, V)

  g = pullback(pres)
  @test g(OO(V)(u)) == OO(U)(x)
  @test g(OO(V)(v)) == OO(U)(y)
end


