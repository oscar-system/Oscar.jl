@testset "localizations at points" begin
  kk = QQ
  R, (x,y) = QQ["x", "y"]
  U = MPolyComplementOfKPointIdeal(R, [1, -1])
  L = Localization(U)
  F = FreeMod(L, 3)
  @test 1//y * (x*F[1] + 2*y^2*F[2]) == x//y*F[1] + 2*y*F[2]

  f = hom(F, F, [x*F[2], 1//y*F[2], F[1]-1//(x*y)*F[2]])
  K, inc_K = kernel(f)
  @test inc_K(K[1]) == -F[1] + x*y*F[2]

  A = L[x 0 1//x; 0 (y-1)//x (y-1)^2//x]
  B = L[x^2 0 1; 0 (x+1)*(y-1) (x+1)*(y-1)^2]
  M = SubQuo(F, A, B)

  g = hom(F, M, [M[1], M[2], x*M[2]])
  Kg, inc_Kg = kernel(g)

  v = x*F[2]-F[3] 
  groebner_basis(M)

end
