@testset "G-sets of permutation groups" begin

  # natural constructions (determined by the types of the seeds)
  G = symmetric_group(6)
  Omega = gset(G)
  @test isa(Omega, GSet)
  @test length(Omega) == 6
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test collect(Omega) == 1:6  # ordering is kept

  Omega = gset(G, [Set([1, 2])])  # action on unordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 15
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  Omega = gset(G, [[1, 2]])  # action on ordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  Omega = gset(G, [(1, 2)])  # action on ordered pairs (repres. by tuples)
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  # constructions by explicit action functions
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test isa(Omega, GSet)
  @test length(Omega) == 740
  @test length(orbits(Omega)) == 2
  @test ! is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  R, x = PolynomialRing(QQ, ["x1", "x2", "x3"]);
  f = x[1]*x[2] + x[2]*x[3]
  G = symmetric_group(3)
  Omega = gset(G, on_indeterminates, [f])
  @test isa(Omega, GSet)
  @test length(Omega) == 3
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

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
  given = 1:4
  Omega = as_gset(G, given)
  @test length(Omega) == 4
  @test length(orbits(Omega)) == 2
  @test isa(Omega, GSet)
  @test isa(orbit(Omega, 1), GSet)
  @test collect(Omega) == given  # ordering is kept

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
  rest = restrict_homomorphism(acthom, derived_subgroup(G)[1])
  @test ! is_bijective(rest)

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

  # block homomorphism
  blhom = action_homomorphism(bl)
  @test ! is_bijective(blhom)

  # is_primitive
  @test ! is_primitive(G8)
  @test ! is_primitive(G8, 1:degree(G8))
  @test is_primitive(S4)
  @test ! is_primitive(S4, 1:3)

  # is_regular
  @test is_regular(G8)
  @test ! is_regular(G8, 1:9)
  @test ! is_regular(S4)

  # is_semiregular
  @test is_semiregular(G8)
  @test ! is_semiregular(G8, 1:9)
  @test ! is_semiregular(S4)

  # is_transitive
  @test is_transitive(G8)
  @test ! is_transitive(G8, 1:9)

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
  @test_throws ArgumentError transitivity(S4, 1:3)
  @test transitivity(S4, 1:4) == 4
  @test transitivity(S4, 1:5) == 0

end

@testset "G-sets of matrix groups over finite fields" begin

  # natural constructions (determined by the types of the seeds)
  G = general_linear_group(2, 3)
  V = free_module(base_ring(G), degree(G))
  Omega = gset(G)
  @test isa(Omega, GSet)
  @test length(Omega) == 9
  @test length(orbits(Omega)) == 2
  @test ! is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test collect(Omega) == collect(V)  # ordering is kept

  Omega = orbit(G, gen(V, 1))
  @test isa(Omega, GSet)
  @test length(Omega) == 8
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  Omega = gset(G, [Set(gens(V))])  # action on unordered pairs of vectors
  @test isa(Omega, GSet)
  @test length(Omega) == 24
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  Omega = gset(G, [gens(V)])  # action on ordered pairs of vectors
  @test isa(Omega, GSet)
  @test length(Omega) == 48
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test is_regular(Omega)
  @test is_semiregular(Omega)

  # orbit
  Omega = gset(G)
  v = gen(V, 1)
  @test length(orbit(Omega, v)) == length(Oscar.orbit_via_Julia(Omega, v))

  # permutation
  Omega = gset(G)
  g = gens(G)[1]
  pi = permutation(Omega, g)
  @test order(pi) == order(g)
  @test degree(parent(pi)) == length(Omega)

  # action homomorphism
  Omega = gset(G)
  acthom = action_homomorphism(Omega)
  @test pi == g^acthom
  @test haspreimage(acthom, pi)[1]
  @test order(image(acthom)[1]) == 48

  # isconjugate
  Omega = gset(G)
  @test isconjugate(Omega, gen(V, 1), gen(V, 2))
  @test ! isconjugate(Omega, zero(V), gen(V, 1))

  # representative_action
  Omega = gset(G)
  rep = representative_action(Omega, gens(V)...)
  @test rep[1]
  @test gen(V, 1) * rep[2] == gen(V, 2)
  rep = representative_action(Omega, zero(V), gen(V, 1))
  @test ! rep[1]

end

@testset "orbits of matrix groups over finite fields" begin

  @testset for F in [ GF(2), GF(3), GF(2,2) ]
    for n in 2:4
      q = order(F)
      V = VectorSpace(F, n)
      GL = general_linear_group(n, F)
      S = sylow_subgroup(GL, 2)[1]
      for G in [GL, S]
#       for k in 0:n   # k = 0 is a problem in GAP 4.12.0
        for k in 1:n
          res = orbit_representatives_and_stabilizers(G, k)
          total = ZZ(0)
          for (U, stab) in res
            total = total + index(G, stab)
            @test length(orbit(stab, U)) == 1
          end
          num = ZZ(1)
          for i in 0:(k-1)
            num = num * (q^n - q^i)
          end
          for i in 0:(k-1)
            num = divexact(num, q^k - q^i)
          end
          @test total == num
        end
      end
    end
  end

end
