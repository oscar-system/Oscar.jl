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
