@testset begin "MPolyAnyMap/MPolyRing"
  # Construction 
  
  Qx, (x, y) = QQ["x", "y"]
  @inferred hom(Qx, Qx, [y, x])
  @test_throws ArgumentError hom(Qx, Qx, [x])
  @test_throws ArgumentError hom(Qx, Qx, [y, y, x])

  R, (x, y) = QQ["x", "y"]
  S, (u, v) = R["u", "v"]
  h = hom(S, S, (a -> S(a^2)), gens(S))
  @test h(u) == u
  @test h(x*u) == x^2 * u
  h = hom(S, S, a -> a^2, gens(S))
  @test h(u) == u
  @test h(x*u) == x^2 * u

  #Qh, (z1, z2) = QQ["z1", "z2"]
  #Qhg = grade(Qh)
  
  A, (x,y) = ZZ["x", "y"]
  f = hom(A, A, [2*x, 5*y])
  R, (u, v) = A["u", "v"]
  h = hom(R, R, f, [u+v, u*v])
  @test (@inferred h(x*u)) == 2*x*u + 2*x*v
end
