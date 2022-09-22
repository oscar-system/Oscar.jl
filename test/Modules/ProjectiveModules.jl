@testset "projective modules" begin
  R, (x,y,z) = QQ["x", "y", "z"]

  M = R[x y+1 z-2 x+y; y z-3 x-y+5 y-z; z x-y z-2 x+y+z]

  I = ideal(R, minors(M, 3))

  Q, _ = quo(R, I)

  A = map_entries(Q, M)

  X = Spec(Q)

  success, P = Oscar._is_projective(A, X)
  @test success
  @test P^2 == P
end
