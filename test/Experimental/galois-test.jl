@testset "Experimental.galois_group" begin

  Zx, x = ZZ["x"]
  k, a = number_field(x^5-2)
  G, C = galois_group(k)
  @test transitive_group_identification(G) == (5, 3)

  U = trivial_subgroup(G)[1]
  L = fixed_field(C, U)
  @test degree(L) == order(G)
  @test length(roots(k.pol, L)) == 5

  R, x = PolynomialRing(QQ, "x")
  pol = x^6 - 366*x^4 - 878*x^3 + 4329*x^2 + 14874*x + 10471
  g, C = galois_group(pol)
  @test order(g) == 18

  s = symmetric_group(4)
  g = sub(s, [s([2,1,4,3]), s([3,4,1,2])])[1]
  act = Oscar.GaloisGrp.action_on_blocks(g, [1, 2])
  @test order(image(act)[1]) == 2

  Qx, x = QQ["x"]
  G, C = galois_group((1//13)*x^2+2)
  @test order(G) == 2

  K, a = number_field(x^4-2)
  G, C = galois_group(K)
  Gc, Cc = galois_group(K, algo = :Complex)
  Gs, Cs = galois_group(K, algo = :Symbolic)
  @test is_isomorphic(G, Gc)
  @test is_isomorphic(G, Gs)
end

import Oscar.GaloisGrp: primitive_by_shape, an_sn_by_shape, cycle_structures
naive_is_giant(G::PermGroup) = is_natural_alternating_group(G) || is_natural_symmetric_group(G)

# compute the cycle types of a bunch of random elements of a permutation group
# the number of samples we take was made up on the spot and has no deeper scientific significance.
sample_cycle_structures(G::PermGroup) = Set(cycle_structure(rand_pseudo(G)) for i in 1:(5*degree(G)+10))

@testset "Cycle types" begin

  # fetch the groups only once, to benefit from caching the result of conjugacy_classes
  grps = [all_transitive_groups(degree => n) for n in 1:20]

  @testset "primitive_by_shape (exact) in degree $n" for n in 1:11
    @test all(G -> is_primitive(G) == primitive_by_shape(cycle_structures(G),n), grps[n] )
  end

  # additional tests but without computing conjugacy classes (slows down the tests too much);
  # note that we only check the implication, as there are cases in degree 15 and above
  # where primitive_by_shape fails to detect primitivity
  @testset "primitive_by_shape (randomized) in degree $n" for n in 12:length(grps)
    @test all(G -> is_primitive(G) || !primitive_by_shape(sample_cycle_structures(G),n), grps[n] )
  end

  @testset "an_sn_by_shape (exact) in degree $n" for n in 1:11
    @test all(G -> naive_is_giant(G) == an_sn_by_shape(cycle_structures(G),n), grps[n] )
  end

  @testset "an_sn_by_shape (randomized) in degree $n" for n in 12:length(grps)
    @test all(G -> naive_is_giant(G) == an_sn_by_shape(sample_cycle_structures(G),n), grps[n] )
  end

end
