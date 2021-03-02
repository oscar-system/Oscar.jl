@testset "G-sets of permutation groups" begin

  # natural constructions (determined by the types of the seeds)
  G = symmetric_group(6)
  Omega = gset(G)
  @test length(Omega) == 6
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, Set([Set([1, 2])]))  # action on unordered pairs
  @test length(Omega) == 15
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  Omega = gset(G, Set([[1, 2]]))  # action on ordered pairs
  @test length(Omega) == 30
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  # constructions by explicit action functions
  Omega = gset(G, permuted, Set([[0,1,0,1,0,1], [1,2,3,4,5,6]]))
  @test length(Omega) == 740
  @test length(orbits(Omega)) == 2
  @test ! istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  R, x = PolynomialRing(QQ, ["x1", "x2", "x3"]);
  f = x[1]*x[2] + x[2]*x[3]
  G = symmetric_group(3)
  Omega = gset(G, on_indeterminates, Set([f]))
  @test length(Omega) == 3
  @test length(orbits(Omega)) == 1
  @test istransitive(Omega)
  @test ! isregular(Omega)
  @test ! issemiregular(Omega)

  # wrapped elements of G-sets
  G = symmetric_group(4)
  omega = [1, 2]
  Omega = gset(G, Set([omega]))  # action on ordered pairs
  g = gens(G)[1]
  x = Omega(omega)
  @test x^g == Omega(omega^g)  # action via `^` is defined for both
  @test length(orbit(x)) == 12

  omega = [0,1,0,1]
  Omega = gset(G, permuted, Set([omega]))
  g = gens(G)[1]
  x = Omega(omega)
  @test x^g == Omega(permuted(omega, g))  # action via `^` is defined for `x`
  @test_throws ErrorException omega^g  # ... but not for the unwrapped object
  @test length(orbit(x)) == 6

  # construction from a known set
  G = sylow_subgroup(symmetric_group(4), 3)[1]
  Omega = as_gset(G, Set(1:4))
  @test length(Omega) == 4
  @test length(orbits(Omega)) == 2

end
