@testset "homomorphisms of decorated rings" begin
  P, (x, y) = polynomial_ring(QQ, [:x, :y])
  Q, (t, ) = polynomial_ring(QQ, [:t])

  # Test the case where we do not need to specify a morphism of the decorations, 
  # because it can be inferred from the images of the generators.
  grP, (X, Y) = grade(P, [1, 1])
  phi0 = hom(grP, grP, [X^2, Y^2 + X*Y])
  @test_throws ErrorException hom(grP, grP, [X^2-1, Y^2 + X*Y])

  grP, (X, Y) = grade(P, [2, 3])
  grQ, (T, ) = grade(Q, [1])

  @test hom(grP, P, gens(P)) isa Oscar.MPolyAnyMap
  @test hom(P, grP, gens(grP)) isa Oscar.MPolyAnyMap

  H = grading_group(grQ)
  dec_map = hom(H, H, 2*gens(H))
  phi1 = hom(grQ, grQ, [T^2])
  phi2 = hom(grQ, grQ, [T^2], decoration_map = dec_map)
  @test phi1 == phi2

  G = grading_group(grP)
  dec_map = hom(G, H, gens(H))
  # We need to specify the grading map here:
  @test_throws ErrorException hom(grP, grQ, [T^2, T^3])
  hom(grP, grQ, [T^2, T^3], decoration_map = dec_map)

  _, u = polynomial_ring(QQ, :u)
  mipo = u^2 + 1
  QQi, imag = extension_field(mipo)

  Pi, (u, v) = polynomial_ring(QQi, [:u, :v])
  grPi, (U, V) = grade(Pi, [10, 15])

  G = grading_group(grP)
  H = grading_group(grPi)
  dec_map = hom(G, H, 5*gens(H))
  phi1 = hom(grP, grPi, gens(grPi), decoration_map = dec_map)
  phi2 = hom(grP, grPi, gens(grPi), decoration_map = dec_map, coefficient_map = QQi)
  @test phi1 == phi2
  coeff_map = MapFromFunc(QQ, QQi, QQi)
  phi3 = hom(grP, grPi, gens(grPi), decoration_map = dec_map, coefficient_map = coeff_map)
  @test phi2 == phi3

  psi1 = hom(grPi, grPi, gens(grPi), decoration_map = identity_map(grading_group(grPi)))
  psi2 = hom(grPi, grPi, gens(grPi), 
             decoration_map = identity_map(grading_group(grPi)), 
             coefficient_map = x->x)
  coeff_map = identity_map(QQi)
  psi3 = hom(grPi, grPi, gens(grPi), decoration_map = identity_map(grading_group(grPi)), coefficient_map = coeff_map)

  # test the various constellations for composition
  @test compose(phi1, psi1) == compose(phi1, psi2)
  @test compose(phi1, psi3) == compose(phi1, psi2)
  @test compose(phi2, psi1) == compose(phi2, psi2)
  @test compose(phi2, psi3) == compose(phi2, psi2)
  @test compose(phi3, psi1) == compose(phi3, psi2)
  @test compose(phi3, psi3) == compose(phi3, psi2)
end
