@testset "galois_group" begin

  Zx, x = ZZ[:x]
  k, a = number_field(x^5-2)
  G, C = galois_group(k)
  @test transitive_group_identification(G) == (5, 3)

  U = trivial_subgroup(G)[1]
  L = fixed_field(C, U)
  @test degree(L) == order(G)
  @test length(roots(L, k.pol)) == 5

  R, x = polynomial_ring(QQ, :x; cached = false)
  pol = x^6 - 366*x^4 - 878*x^3 + 4329*x^2 + 14874*x + 10471
  g, C = galois_group(pol)
  @test order(g) == 18

  s = symmetric_group(4)
  g = sub(s, [s([2,1,4,3]), s([3,4,1,2])])[1]
  act = Oscar.GaloisGrp.action_on_blocks(g, [1, 2])
  @test order(image(act)[1]) == 2

  Qx, x = QQ[:x]
  G, C = galois_group((1//13)*x^2+2)
  @test order(G) == 2

  K, a = number_field(x^4-2; cached = false)
  G, C = galois_group(K)
  Gc, Cc = galois_group(K, algorithm = :Complex)
  Gs, Cs = galois_group(K, algorithm = :Symbolic)
  @test is_isomorphic(G, Gc)
  @test is_isomorphic(G, Gs)
  _G, _C = galois_group(K)
  # test caching
  @test C === _C 
  _G, _C = galois_group(K; redo = true)
  @test C !== _C

  # from the book
  K, a = number_field(x^9 - 3*x^8 + x^6 + 15*x^5 - 13*x^4 -
                      3*x^3 + 4*x - 1, "a")
  G, C = galois_group(K)
  @test order(G) == 216

  # Ehrhart
  G, C = galois_group((2*x+1)^2)
  @test order(G) == 1
  @test degree(G) == 2

  #from errors:
  G, C = galois_group((x^4 + 1)^3 * x^2 * (x^2 - 4*x + 1)^5)
  @test order(G) == 8
  @test degree(G) == 24
  k = fixed_field(C, sub(G, [one(G)])[1])
  @test degree(k) == 8

  G, C = galois_group((x^3-2)^2*(x^3-5)^2*(x^2-6))
  @test order(G) == 36
  @test degree(G) == 14

  #Fabian Gundlach...
  G, C = galois_group(x^8 - 2*x^7 - 48*x^6 + 58*x^5 + 846*x^4 - 4614*x^3 + 6609*x^2 + 48742*x + 493474)
  @test order(G) == 32
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
    # If primitive_by_shape returns true then is_primitive must return true. If
    # primitive_by_shape returns false then we can't say anything, as we just
    # might have sampled bad elements.
    # In other words, we test for the implication `primitive_by_shape => is_primitive`.
    @test all(G -> !primitive_by_shape(sample_cycle_structures(G),n) || is_primitive(G), grps[n] )
  end

  @testset "an_sn_by_shape (exact) in degree $n" for n in 1:11
    @test all(G -> naive_is_giant(G) == an_sn_by_shape(cycle_structures(G),n), grps[n] )
  end

  @testset "an_sn_by_shape (randomized) in degree $n" for n in 12:length(grps)
    # If an_sn_by_shape returns true then naive_is_giant must return true. If
    # an_sn_by_shape returns false then we can't say anything, as we just
    # might have sampled bad elements.
    # In other words, we test for the implication `an_sn_by_shape => naive_is_giant`.
    @test all(G -> !an_sn_by_shape(sample_cycle_structures(G),n) || naive_is_giant(G), grps[n] )
  end

  let # Ehrhart polynomial problems
    Qx, x = QQ[:x]
    f = 4//45*x^6 + 4//15*x^5 + 14//9*x^4 + 8//3*x^3 + 196//45*x^2 + 46//15*x + 1
    G, = galois_group(f)
    @test small_group_identification(G) == (4, 2)
  end
end

@testset "Galois group issue" begin
  # Contributed by "Lloyd" on slack
  R, s = QQ[:s]
  K, q = number_field(s^2 - 2, "q")
  Kw, w = polynomial_ring(K, :w)
  f = w^16 - 32*w^14 - 192*w^12 + 22720*w^10 + 23104*w^8 - 2580480*w^6 + 41287680*w^4 + 106168320*w^2 + 84934656
  g, s = galois_group(f)
  @test order(g) == 1536
  ss = map(representative, subgroup_classes(g))
  #should be of order 16, so field of degree 96
  H = ss[1000]
  f = fixed_field(s, H)
  @test degree(f) == order(g)//order(H)
end
