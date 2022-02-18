@testset "spinor norms" begin
  g = ZZ[-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
  G =  diagonal_matrix(fmpq[2, 3//2, 4//3, 5//4])
  s = Oscar.spin(G, g)
  @test issquare(s//5)

  G =  QQ[4//5 0 0 0 0; 0 0 1//4 0 0; 0 1//4 0 0 0; 0 0 0 1//2 1//4; 0 0 0 1//4 1//2]
  Oscar.sigma_sharp(4, 1280, G, 5)
  Oscar.sigma_sharp(4, 1280, G, 2)


  diag = QQ[3//2;]
  g = QQ[22876792454960;]

  @test (-1, 3) == Oscar.det_spin(diag, g, 3, 25)
end
