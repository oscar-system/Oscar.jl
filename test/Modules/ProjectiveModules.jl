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

  W = MPolyQuoLocalizedRing(R, I, units_of(R))
  F3 = FreeMod(W, 3)
  F4 = FreeMod(W, 4)
  h = hom(F3, F4, map_entries(W, M))

  mod, _ = cokernel(h)
  
  success, pr, inc = is_projective(mod)
  @test success
  Palt = compose(pr, inc)
  Palt2 = compose(Palt, Palt)
  @test all(x->(Palt2(x) == Palt(x)), gens(domain(Palt)))

end
