@testset "kernels via flattenings" begin
  R, (x,y,z) = QQ["x", "y", "z"]

  S, _ = PolynomialRing(R, ["s", "t"])
  S_flat = oscar._flatten(S)
  R_to_S_flat = oscar._base_ring_to_flat(S)
  S_to_S_flat = oscar._orig_to_flat(S)
  @test vcat(S_to_S_flat.(gens(S)), R_to_S_flat.(gens(R))) == gens(S_flat)
  S_flat_to_S = oscar._flat_to_orig(S)
  @test gens(S) == S_flat_to_S.(S_to_S_flat.(gens(S)))
  @test R_to_S_flat.(gens(R)) == S_to_S_flat.(S.(gens(R)))

  f = hom(S, S, gens(S))
  K = kernel(f)
  @test zero(S) in K

  f = hom(S, S, hom(R, R, [first(gens(R)), first(gens(R)), last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gens(R)[1] - gens(R)[2]) in K

  R, _ = quo(R, ideal(R, sum(gens(R))))
  S, _ = PolynomialRing(R, ["s", "t"])
  S_flat = oscar._flatten(S)
  R_to_S_flat = oscar._base_ring_to_flat(S)
  S_to_S_flat = oscar._orig_to_flat(S)
  @test vcat(S_to_S_flat.(gens(S)), R_to_S_flat.(gens(R))) == gens(S_flat)
  S_flat_to_S = oscar._flat_to_orig(S)
  @test gens(S) == S_flat_to_S.(S_to_S_flat.(gens(S)))
  @test R_to_S_flat.(gens(R)) == S_to_S_flat.(S.(gens(R)))

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
  S_flat = oscar._flatten(S)
  R_to_S_flat = oscar._base_ring_to_flat(S)
  S_to_S_flat = oscar._orig_to_flat(S)
  @test vcat(S_to_S_flat.(gens(S)), R_to_S_flat.(gens(R))) == gens(S_flat)
  S_flat_to_S = oscar._flat_to_orig(S)
  @test gens(S) == S_flat_to_S.(S_to_S_flat.(gens(S)))
  @test R_to_S_flat.(gens(R)) == S_to_S_flat.(S.(gens(R)))

  f = hom(S, S, gens(S))
  @test oscar._flatten(f) isa Hecke.Map
  #K = kernel(f) # currently not implemented
  #@test zero(S) in K 

  f = hom(S, S, x->x, gens(S))
  #K = kernel(f)
  #@test S(gens(R)[1] - gens(R)[2]) in K

  I = ideal(R, first(gens(R))^2)
  R, _ = quo(R, I)
  S, _ = PolynomialRing(R, ["s", "t"])
  S_flat = oscar._flatten(S)
  R_to_S_flat = oscar._base_ring_to_flat(S)
  S_to_S_flat = oscar._orig_to_flat(S)
  @test vcat(S_to_S_flat.(gens(S)), R_to_S_flat.(gens(R))) == gens(S_flat)
  S_flat_to_S = oscar._flat_to_orig(S)
  @test gens(S) == S_flat_to_S.(S_to_S_flat.(gens(S)))
  @test R_to_S_flat.(gens(R)) == S_to_S_flat.(S.(gens(R)))

  f = hom(S, S, gens(S))
  K = kernel(f)
  @test zero(S) in K

  f = hom(S, S, hom(R, R, [zero(R), first(gens(R)) + gens(R)[2], last(gens(R))]), gens(S))
  K = kernel(f)
  @test S(gens(R)[1]) in K
end
