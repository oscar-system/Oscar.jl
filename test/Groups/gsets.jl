@testset "G-sets of permutation groups" begin

  # natural constructions (determined by the types of the seeds)
  G = symmetric_group(6)
  Omega = gset(G)
  @test isa(Omega, GSet)
  @test length(Omega) == 6
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! is_regular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, [Set([1, 2])])  # action on unordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 15
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! is_regular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, [[1, 2]])  # action on ordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! is_regular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, [(1, 2)])  # action on ordered pairs (repres. by tuples)
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! is_regular(Omega)
  @test ! issemiregular(Omega)

  # constructions by explicit action functions
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test isa(Omega, GSet)
  @test length(Omega) == 740
  @test length(orbits(Omega)) == 2
  @test ! istransitive(Omega)
  @test ! is_regular(Omega)
  @test ! issemiregular(Omega)

  R, x = PolynomialRing(QQ, ["x1", "x2", "x3"]);
  f = x[1]*x[2] + x[2]*x[3]
  G = symmetric_group(3)
  Omega = gset(G, on_indeterminates, [f])
  @test isa(Omega, GSet)
  @test length(Omega) == 3
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! is_regular(Omega)
  @test ! issemiregular(Omega)

  # seeds can be anything iterable
  G = symmetric_group(6)
  @test isa(gset(G, 1:6), GSet)
  @test isa(gset(G, collect(1:6)), GSet)
  @test isa(gset(G, Set(1:6)), GSet)

  # basic functionality
  G = symmetric_group(6)
  Omega = gset(G, [Set([1, 2])])
  @test representative(Omega) in Omega
  @test acting_domain(Omega) == G

  # wrapped elements of G-sets
  G = symmetric_group(4)
  omega = [1, 2]
  Omega = gset(G, Set([omega]))  # action on ordered pairs
  g = gens(G)[1]
  x = Omega(omega)
  @test x in Omega
  @test unwrap(x) == omega
  @test unwrap(omega) == omega
  @test x^g == Omega(omega^g)  # action via `^` is defined for both
  @test length(orbit(x)) == 12

  omega = [0,1,0,1]
  Omega = gset(G, permuted, Set([omega]))
  g = gens(G)[1]
  x = Omega(omega)
  @test x in Omega
  @test unwrap(x) == omega
  @test unwrap(omega) == omega
  @test x^g == Omega(permuted(omega, g))  # action via `^` is defined for `x`
  @test_throws ErrorException omega^g  # ... but not for the unwrapped object
  orb = orbit(x)
  @test length(orb) == 6
  @test isa(orb, GSet)

  # construction from a known set
  G = sylow_subgroup(symmetric_group(4), 3)[1]
  Omega = as_gset(G, Set(1:4))
  @test length(Omega) == 4
  @test length(orbits(Omega)) == 2
  @test isa(Omega, GSet)
  @test isa(orbit(Omega, 1), GSet)

  # orbit
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test length(orbit(Omega, [0,1,0,1,0,1])) == length(Oscar.orbit_via_Julia(Omega, [0,1,0,1,0,1]))

  # permutation
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  g = gens(G)[1]
  pi = permutation(Omega, g)
  @test order(pi) == order(g)
  @test degree(parent(pi)) == length(Omega)

  # action homomorphism
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  acthom = action_homomorphism(Omega)
  @test pi == g^acthom
  @test haspreimage(acthom, pi)[1]
  @test order(image(acthom)[1]) == 720

  # is_conjugate
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test is_conjugate(Omega, [0,1,0,1,0,1], [1,0,1,0,1,0])
  @test ! is_conjugate(Omega, [0,1,0,1,0,1], [1,2,3,4,5,6])

  # representative_action
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  rep = representative_action(Omega, [0,1,0,1,0,1], [1,0,1,0,1,0])
  @test rep[1]
  @test permuted([0,1,0,1,0,1], rep[2]) == [1,0,1,0,1,0]
  rep = representative_action(Omega, [0,1,0,1,0,1], [1,2,3,4,5,6])
  @test ! rep[1]

end

@testset "natural action of permutation groups" begin

  G8 = transitive_group(8, 3)
  S4 = symmetric_group(4)
  @test order(G8) == 8

  # all_blocks
  bl = all_blocks(G8)
  @test length(bl) == 14
  @test [1, 2] in bl
  @test length(all_blocks(S4)) == 0

  # blocks
  bl = blocks(G8)
  @test elements(bl) == [[1, 8], [2, 3], [4, 5], [6, 7]]
  @test elements(bl) == elements(blocks(G8, 1:degree(G8)))

  # is_primitive
  @test ! is_primitive(G8)
  @test ! is_primitive(G8, 1:degree(G8))
  @test is_primitive(S4)
  @test is_primitive(S4, 1:3)

  # is_regular
  @test is_regular(G8)
  @test ! is_regular(G8, 1:9)
  @test ! is_regular(S4)

  # issemiregular
  @test issemiregular(G8)
  @test ! issemiregular(G8, 1:9)
  @test ! issemiregular(S4)

  # istransitive
  @test istransitive(G8)
  @test ! istransitive(G8, 1:9)

  # maximal_blocks
  bl = maximal_blocks(G8)
  @test elements(bl) == [[1, 2, 3, 8], [4, 5, 6, 7]]
  @test elements(bl) == elements(maximal_blocks(G8, 1:degree(G8)))
  @test elements(maximal_blocks(S4)) == [[1, 2, 3, 4]]

  # minimal_block_reps
  bl = minimal_block_reps(G8)
  @test bl == [[1,i] for i in 2:8]
  @test bl == minimal_block_reps(G8, 1:degree(G8))
  @test minimal_block_reps(S4) == [[1, 2, 3, 4]]

  # transitivity
  @test transitivity(G8) == 1
  @test transitivity(S4) == 4
  @test transitivity(S4, 1:3) == 3
  @test transitivity(S4, 1:5) == 0

end
