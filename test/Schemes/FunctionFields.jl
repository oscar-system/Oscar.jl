@testset "fraction fields of varieties" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  @test is_irreducible(Spec(R))
  @test is_irreducible(Spec(R, ideal(R, x)))
  @test !is_irreducible(Spec(R, ideal(R, x*y)))
  @test is_irreducible(Spec(Localization(R, units_of(R))[1]))
  @test !is_irreducible(Spec(R, ideal(R, x*y), units_of(R)))

  P = projective_space(QQ, 2)
  S = homogeneous_poly_ring(P)
  C = subscheme(P, ideal(S, S[1]*S[2]-S[3]^2))
  Ccov = as_covered_scheme(C)

  KK = VarietyFunctionField(Ccov)
  U2 = patches(Ccov)[2]
  a = 5*gens(OO(U2))[1]*gens(OO(U2))[2]

  b = KK(a, one(a))
  @test b^2 - b == KK(a^2-a, one(a))
end
