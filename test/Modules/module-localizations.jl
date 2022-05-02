@testset "module localizations 1" begin
  kk = QQ
  R, (x,y) = QQ["x", "y"]
  U = MPolyComplementOfKPointIdeal(R, [0, 0])
  L = Localization(U)
  F = FreeMod(L, 3)
  A = L[x 0 1; 0 y y^2]
  B = L[x^2 0 1; 0 y^2 y^3]
  M = SubQuo(F, A, B)
  phi = hom(F, F, [L(1//(y-1))*F[2], 0*F[1], x*F[2]])
  f = hom(F, M, [L(x, y-2)*M[1], L(y, y^2-4)*M[2], M[1]])
  Kphi, iphi = kernel(phi)
  @test represents_element(F[2], Kphi)
  @test represents_element(x*(y-1)*F[1]-F[3], Kphi)
  kernel(f)

  T = MPolyPowersOfElement(R, [x, y])
  W = Localization(T)
  F = FreeMod(W, 2)
  A = W[x 0; 0 y^2]
  B = W[x^2//y y; 0 y^7]
  M = SubQuo(F, A, B)
end

@testset "module localizations 2" begin
  R, (x,y) = QQ["x", "y"]
  U = MPolyPowersOfElement(x+y)
  S = Localization(U)
  F = FreeMod(S, 2)
  Fb = base_ring_module(F)
  A = S[x//(x+y); y//(x+y)^2]
  B, D = Oscar.clear_denominators(A)

  b = MatrixSpace(S, 1, 1)((x+y)*x + 5*y//(x+y)^10)
  Oscar.has_solution(A, b)

  V = MPolyComplementOfPrimeIdeal(ideal(R, [x,y]))
  S = Localization(V)
  A = S[x//(x+y+1); y*(x-5)^3]
  b = MatrixSpace(S, 1, 1)((x+y)*x + 5*y//(x+y+2)^10)
  Oscar.has_solution(A, b)

  F = FreeMod(S, 2)
  G = FreeMod(S, 1)
  f = hom(F, G, S[x; y])
  K, inc = kernel(f)

  M, p = cokernel(f)
  v = (x^2 + y^3)*G[1]
  coordinates(v, M)

  g = MatrixSpace(S, 1, 1)((x^2+y^3)//(x+1))
  N, pr = quo(G, g)
  Kpr, incKpr = kernel(pr)

  L, incL = sub(F, S[x y])
  h = compose(incL, compose(f, pr))
  K3, incK3 = kernel(h)
end

@testset "module localizations 3" begin
  R, (x,y) = QQ["x", "y"]
  U = MPolyPowersOfElement(x^7)
  S = Localization(U)
  F = FreeMod(S, 1)
  A = S[x^4*y^2; x^2*y]
  B = S[y^8; y^9]
  M = SubQuo(F, A, B)
  Fb = base_ring_module(F)
  represents_element(y*Fb[1], Oscar.pre_saturated_module(M))
  represents_element(y*F[1], M)
  coordinates(y*F[1], M)
  represents_element(y*Fb[1], Oscar.pre_saturated_module(M))
  coordinates(y*F[1], M)
  c = coordinates(y*F[1], M)
end
