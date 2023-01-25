#@testset "kernels via flattenings" begin
  R, (x,y,z) = QQ["x", "y", "z"]

  S, _ = PolynomialRing(R, ["s", "t"])
  phi = oscar.flatten(S)
  @test vcat(phi.(gens(S)), oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))

  f = hom(S, S, gens(S))
  K = kernel(f)
  @test zero(S) in K

  f = hom(S, S, hom(R, R, [first(gens(R)), first(gens(R)), last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gens(R)[1] - gens(R)[2]) in K

  R, _ = quo(R, ideal(R, sum(gens(R))))
  S, _ = PolynomialRing(R, ["s", "t"])
  phi = oscar.flatten(S)
  @test vcat(phi.(gens(S)), oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))

  f = hom(S, S, gens(S))
  K = kernel(f)
  @test zero(S) in K

  f = hom(S, S, hom(R, R, [first(gens(R)) + gens(R)[2], zero(R), last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gens(R)[2]) in K

  R = base_ring(R)
  U = powers_of_element(sum(gens(R)))
  R, _ = localization(R, U)
  S, _ = PolynomialRing(R, ["s", "t"])
  phi = oscar.flatten(S)
  @test vcat(phi.(gens(S)), oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))

  f = hom(S, S, gens(S))
  @test oscar.flatten(f) isa Hecke.Map
  K = kernel(f) 
  @test zero(S) in K 
 
  f = hom(S, S, x->x, [S[1], S[1]])
  K = kernel(f)
  @test S[1] - S[2] in K
 
  I = ideal(R, first(gens(R))^2)
  R, _ = quo(R, I)
  S, _ = PolynomialRing(R, ["s", "t"])
  phi = oscar.flatten(S)
  @test vcat(phi.(gens(S)), oscar.map_from_coefficient_ring_to_flattening(phi).(gens(coefficient_ring(S)))) == gens(codomain(phi))
  @test gens(S) == inverse(phi).(phi.(gens(S)))
  @test oscar.map_from_coefficient_ring_to_flattening(phi).(gens(R)) == phi.(S.(gens(R)))
 
  f = hom(S, S, [S[1], S[1]])
  K = kernel(f)
  @test S[1] - S[2] in K
 
  f = hom(S, S, hom(R, R, [zero(R), first(gens(R)) + gens(R)[2], last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gens(R)[1]) in K
#end
