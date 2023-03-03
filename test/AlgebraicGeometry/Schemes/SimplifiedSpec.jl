@testset "simplified spectra" begin
  IP = projective_space(QQ, 2)
  S = ambient_coordinate_ring(IP)
  (x,y,z) = gens(S)
  f = x^2 - y*z
  CP = subscheme(IP, f)
  C = covered_scheme(CP)
  simplify!(C)
  CC = coverings(C)[2]
  U = patches(CC)
  OC = OO(C)
  (u, v) = gens(OC(U[1]))
  W = PrincipalOpenSubset(U[1], u-1)
  @test !(W == U[1])
  OC(U[1], W)
end


