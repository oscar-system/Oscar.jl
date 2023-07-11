@testset "Rational maps" begin
  IP3 = projective_space(QQ, [:x, :y, :z, :w])
  S = homogeneous_coordinate_ring(IP3)
  (x, y, z, w) = gens(S)
  I = ideal(S, x*w - y*z)
  X_proj = subscheme(IP3, I)
  X = covered_scheme(X_proj)

  P1 = projective_space(QQ, 1)
  IP1 = covered_scheme(P1)
  U = X[1][4]
  V = IP1[1][2]
  (x, y, z) = gens(ambient_coordinate_ring(U))
  Phi = oscar.RationalMap(X, IP1, U, V, [x//y])

  @test domain(Phi) === X 
  @test codomain(Phi) === IP1

  oscar.realize_on_patch(Phi, U)
  oscar.realize_on_patch(Phi, X[1][1])
  oscar.realize_on_patch(Phi, X[1][2])
  oscar.realize_on_patch(Phi, X[1][3])
  oscar.realize(Phi)
end
