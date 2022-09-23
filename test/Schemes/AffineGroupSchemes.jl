import Oscar._kk_star

function map_to_prod_from_comp(
    f::SpecMor{GroupType, GroupType, <:Any}, 
    g::SpecMor{GroupType, GroupType, <:Any}
  ) where {GroupType<:_kk_star}
  G = domain(f) 
  G == codomain(f) == domain(g) == codomain(g) || error("both maps must have equal domains and codomains")
  GxG = product_over_ground_field(G)
  return SpecMor(G, GxG, hom(OO(GxG), OO(G), vcat(pullback(f).(gens(OO(G))), pullback(g).(gens(OO(G))))))
end

@testset "kk_star_as_group_scheme" begin
  P, t = QQ["t"]
  h = t^7 - 1
  kk = splitting_field(h)
  G = _kk_star(kk)

  @test is_identity_map(compose(first_inclusion(G), multiplication_map(G)))
  GxG = product_over_ground_field(G)
  inv_map = SpecMor(G, G, pullback(inverse_map(G)))
  id_map = SpecMor(G, G, pullback(identity_map(G)))
  anti_diag = map_to_prod_from_comp(inv_map, id_map)
  h = compose(anti_diag, multiplication_map(G))
  g = map_to_prod_from_comp(identity_map(G), h)
  tmp = compose(g, multiplication_map(G))
  @test is_identity_map(tmp)
end

@testset "matrix groups" begin

  G = special_linear_group(2, QQ)
  rho = canonical_representation(G)
  rep = induced_representation_on_symmetric_power(G, 4)

  O = omega_process(G)
  I = nullcone_ideal(rep)
  rop = reynolds_operator_from_omega_process(rep)
  g = [rop(g) for g in gens(I)]

  check = check_invariance_function(rep)
  @test all(h->check(h), g)

  A, m = invariant_ring(rep)

  l = lift_to_invariant_polynomial_func(rep)
  [l(4*y) for y in g]
  [l(y^2) for y in g]
  p = g[1]^3-5*g[2]^2 - 7*g[3]*g[1]
  h = l(p)
  @test h == A[1]^3-5*A[2]^2 - 7*A[3]*A[1]
end
