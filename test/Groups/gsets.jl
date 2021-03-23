@testset "G-sets of permutation groups" begin

  # natural constructions (determined by the types of the seeds)
  G = symmetric_group(6)
  Omega = gset(G)
  @test isa(Omega, GSet)
  @test length(Omega) == 6
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, [Set([1, 2])])  # action on unordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 15
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, [[1, 2]])  # action on ordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, [(1, 2)])  # action on ordered pairs (repres. by tuples)
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  # constructions by explicit action functions
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test isa(Omega, GSet)
  @test length(Omega) == 740
  @test length(orbits(Omega)) == 2
  @test ! istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  R, x = PolynomialRing(QQ, ["x1", "x2", "x3"]);
  f = x[1]*x[2] + x[2]*x[3]
  G = symmetric_group(3)
  Omega = gset(G, on_indeterminates, [f])
  @test isa(Omega, GSet)
  @test length(Omega) == 3
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
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

  # isconjugate
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test isconjugate(Omega, [0,1,0,1,0,1], [1,0,1,0,1,0])
  @test ! isconjugate(Omega, [0,1,0,1,0,1], [1,2,3,4,5,6])

  # representative_action
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  rep = representative_action(Omega, [0,1,0,1,0,1], [1,0,1,0,1,0])
  @test rep[1]
  @test permuted([0,1,0,1,0,1], rep[2]) == [1,0,1,0,1,0]
  rep = representative_action(Omega, [0,1,0,1,0,1], [1,2,3,4,5,6])
  @test ! rep[1]

end
