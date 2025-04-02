@testset "G-sets of permutation groups" begin

  # natural constructions (determined by the types of the seeds)
  G = symmetric_group(6)
  Omega = gset(G)
  @test AbstractAlgebra.PrettyPrinting.repr_terse(Omega) == "G-set"
  @test isa(Omega, GSet)
  @test (@inferred length(Omega)) == 6
  @test (@inferred length(@inferred orbits(Omega))) == 1
  @test is_transitive(Omega)
  @test is_primitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test collect(Omega) == 1:6  # ordering is kept
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
  @test order(stabilizer(Omega, 1)[1]) == 120
  @test order(stabilizer(Omega, Set([1, 2]))[1]) == 48
  @test order(stabilizer(Omega, [1, 2])[1]) == 24
  @test order(stabilizer(Omega, (1, 2))[1]) == 24

  Omega = gset(G, [Set([1, 2])])  # action on unordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 15
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test is_primitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
  @test order(stabilizer(Omega, Set([1, 3]))[1]) == 48
  @test order(stabilizer(Omega, Set([Set([1, 2]), Set([1, 3])]))[1]) == 12
  @test order(stabilizer(Omega, [Set([1, 2]), Set([1, 3])])[1]) == 6
  @test order(stabilizer(Omega, (Set([1, 2]), Set([1, 3])))[1]) == 6
  @test_throws MethodError stabilizer(Omega, [1, 2])

  Omega = gset(G, [[1, 2]])  # action on ordered pairs
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
  @test order(stabilizer(Omega, [1, 3])[1]) == 24
  @test order(stabilizer(Omega, Set([[1, 2], [1, 3]]))[1]) == 12
  @test order(stabilizer(Omega, [[1, 2], [1, 3]])[1]) == 6
  @test order(stabilizer(Omega, ([1, 2], [1, 3]))[1]) == 6
  @test_throws MethodError stabilizer(Omega, Set([1, 2]))
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_primitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  Omega = gset(G, [(1, 2)])  # action on ordered pairs (repres. by tuples)
  @test isa(Omega, GSet)
  @test length(Omega) == 30
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_primitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  # larger examples
  G = symmetric_group(100)
  Omega = gset(G)
  S1, _ = stabilizer(Omega, [1, 2, 3, 4, 5])
  @test order(S1) == factorial(big(95))
  S2, _ = stabilizer(Omega, (1, 2, 3, 4, 5))
  @test S2 == S1
  S3, _ = stabilizer(Omega, Set([1, 2, 3, 4, 5]))
  @test order(S3) == order(S1) * factorial(big(5))

  # constructions by explicit action functions
  G = symmetric_group(6)
  omega = [0,1,0,1,0,1]
  Omega = gset(G, permuted, [omega, [1,2,3,4,5,6]])
  @test isa(Omega, GSet)
  @test length(Omega) == 740
  @test order(stabilizer(Omega, omega)[1]) * length(orbit(Omega, omega)) == order(G)
  @test order(stabilizer(Omega, Set([omega, [1,0,0,1,0,1]]))[1]) == 8
  @test order(stabilizer(Omega, [omega, [1,0,0,1,0,1]])[1]) == 4
  @test order(stabilizer(Omega, (omega, [1,0,0,1,0,1]))[1]) == 4
  @test_throws MethodError stabilizer(Omega, Set(omega))
  @test length(orbits(Omega)) == 2
  @test ! is_transitive(Omega)
  @test ! is_primitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test transitivity(Omega) == 0
  @test_throws ArgumentError rank_action(Omega)
  @test_throws ArgumentError gset(G, permuted, omega)

  R, x = polynomial_ring(QQ, [:x1, :x2, :x3]);
  f = x[1]*x[2] + x[2]*x[3]
  G = symmetric_group(3)
  Omega = gset(G, on_indeterminates, [f])
  @test isa(Omega, GSet)
  @test length(Omega) == 3
  @test length(orbits(Omega)) == 1
  @test order(stabilizer(Omega)[1]) * length(orbit(Omega, f)) == order(G)
  @test is_transitive(Omega)
  @test is_primitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test transitivity(Omega) == 3
  @test rank_action(Omega) == 2

  # seeds can be anything iterable
  G = symmetric_group(6)
  @test isa(gset(G, 1:6), GSet)
  @test isa(gset(G, collect(1:6)), GSet)
  @test isa(gset(G, Set(1:6)), GSet)

  # basic functionality
  G = symmetric_group(6)
  Omega = gset(G, [Set([1, 2])])
  @test representative(Omega) in Omega
  @test acting_group(Omega) == G

  # wrapped elements of G-sets
  G = symmetric_group(4)
  omega = [1, 2]
  Omega = gset(G, Set([omega]))  # action on ordered pairs
  g = gen(G, 1)
  x = Omega(omega)
  @test x in Omega
  @test unwrap(x) == omega
  @test unwrap(omega) == omega
  @test x^g == Omega(omega^g)  # action via `^` is defined for both
  @test length(orbit(x)) == 12

  omega = [0,1,0,1]
  Omega = gset(G, permuted, Set([omega]))
  g = gen(G, 1)
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
  orb = orbit(Omega, [0,1,0,1,0,1])
  @test length(orb) == length(Oscar.orbit_via_Julia(Omega, [0,1,0,1,0,1]))
  @test orbits(orb) == [orb]

  # permutation
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  g = gen(G, 1)
  pi = permutation(Omega, g)
  @test order(pi) == order(g)
  @test degree(parent(pi)) == length(Omega)

  # action homomorphism
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  acthom = action_homomorphism(Omega)
  @test pi == g^acthom
  @test has_preimage_with_preimage(acthom, pi)[1]
  @test order(image(acthom)[1]) == 720
  rest = restrict_homomorphism(acthom, derived_subgroup(G)[1])
  @test ! is_bijective(rest)

  # is_conjugate
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test is_conjugate(Omega, [0,1,0,1,0,1], [1,0,1,0,1,0])
  @test ! is_conjugate(Omega, [0,1,0,1,0,1], [1,2,3,4,5,6])

  # is_conjugate_with_data
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  rep = is_conjugate_with_data(Omega, [0,1,0,1,0,1], [1,0,1,0,1,0])
  @test rep[1]
  @test permuted([0,1,0,1,0,1], rep[2]) == [1,0,1,0,1,0]
  rep = is_conjugate_with_data(Omega, [0,1,0,1,0,1], [1,2,3,4,5,6])
  @test ! rep[1]

  # stabilizer
  G = symmetric_group(6)
  Omega = gset(G, permuted, [[0,1,0,1,0,1], [1,2,3,4,5,6]])
  @test_throws ArgumentError stabilizer(Omega, [0,0,0,0,0,0])
  omega = representative(Omega)
  @test stabilizer(Omega) == stabilizer(Omega, omega)
  @test stabilizer(Omega) !== stabilizer(Omega, omega)
  @test stabilizer(Omega) === stabilizer(Omega)
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
  @test transitivity(gset(G8)) == 1
  @test transitivity(S4) == 4
  @test transitivity(gset(S4)) == 4
  @test_throws ArgumentError transitivity(S4, 1:3)
  @test transitivity(S4, 1:4) == 4
  @test transitivity(S4, 1:5) == 0

end

@testset "G-sets of matrix groups over finite fields" begin

  # natural constructions (determined by the types of the seeds)
  G = general_linear_group(2, 3)
  V = free_module(base_ring(G), degree(G))
  v = gen(V, 1)
  Omega = gset(G)
  @test isa(Omega, GSet)
  @test length(Omega) == 9
  @test order(stabilizer(Omega, v)[1]) * length(orbit(Omega, v)) == order(G)
  @test length(orbits(Omega)) == 2
  @test ! is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test collect(Omega) == collect(V)  # ordering is kept

  Omega = orbit(G, v)
  @test isa(Omega, GSet)
  @test length(Omega) == 8
  @test order(stabilizer(Omega, v)[1]) * length(Omega) == order(G)
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  Omega = gset(G, [Set(gens(V))])  # action on unordered pairs of vectors
  @test isa(Omega, GSet)
  @test length(Omega) == 24
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  Omega = gset(G, [gens(V)])  # action on ordered pairs of vectors
  @test isa(Omega, GSet)
  @test length(Omega) == 48
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
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
  g = gen(G, 1)
  pi = permutation(Omega, g)
  @test order(pi) == order(g)
  @test degree(parent(pi)) == length(Omega)

  # action homomorphism
  Omega = gset(G)
  acthom = action_homomorphism(Omega)
  @test pi == g^acthom
  @test has_preimage_with_preimage(acthom, pi)[1]
  @test order(image(acthom)[1]) == 48

  # is_conjugate
  Omega = gset(G)
  @test is_conjugate(Omega, gen(V, 1), gen(V, 2))
  @test ! is_conjugate(Omega, zero(V), gen(V, 1))

  # is_conjugate_with_data
  Omega = gset(G)
  rep = is_conjugate_with_data(Omega, gens(V)...)
  @test rep[1]
  @test gen(V, 1) * rep[2] == gen(V, 2)
  rep = is_conjugate_with_data(Omega, zero(V), gen(V, 1))
  @test ! rep[1]

end

@testset "orbits of matrix groups over finite fields" begin

  @testset for F in [ GF(2), GF(3), GF(2,2) ], n in 2:4
    q = order(F)
    V = vector_space(F, n)
    GL = general_linear_group(n, F)
    S = sylow_subgroup(GL, 2)[1]
    for G in [GL, S]
      for k in 0:n
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

@testset "G-sets of matrix groups in characteristic zero" begin

  # natural constructions (determined by the types of the seeds)
  G = matrix_group(permutation_matrix(QQ,[3,1,2]))
  R, (x,y,z) = polynomial_ring(QQ, [:x,:y,:z])
  f = x^2 + y
  orb = orbit(G, f)
  @test length(orb) == 3

  F = QQBarField()
  e = one(F)
  s, c = sincospi(2 * e / 3)
  mat_rot = matrix([c -s; s c])
  G = matrix_group(mat_rot)
  p = F.([1, 0])
  orb = orbit(G, *, p)
  @test length(orb) == 3


end

@testset "G-sets by right transversals" begin
  G = symmetric_group(5)
  H = sylow_subgroup(G, 2)[1]
  Omega = right_cosets(G, H)
  @test AbstractAlgebra.PrettyPrinting.repr_terse(Omega) == "Right cosets of groups"
  @test isa(Omega, GSet)
  @test acting_group(Omega) == G
  @test length(Omega) == index(G, H)
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
  @test Omega[end] == Omega[length(Omega)]
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)
  @test length(blocks(Omega)) == 5
  @test length(minimal_block_reps(Omega)) == 1
  @test length(all_blocks(Omega::GSet)[1]) == 3


  @test eltype(Omega) == typeof(representative(Omega))

  # iteration
  for i in 1:length(Omega)
    @test findfirst(is_equal(Omega[i]), Omega) == i
  end
  rep = representative(Omega)
  for omega in Omega
    @test one(G) in omega || omega != rep
  end

  # orbit
  g = gen(G, 1)
  pnt = right_coset(H, g)
  @test pnt in Omega
  @test length(orbit(Omega, pnt)) == length(Oscar.orbit_via_Julia(Omega, pnt))

  # permutation
  pi = permutation(Omega, g)
  @test order(pi) == order(g)
  @test degree(parent(pi)) == length(Omega)
  fun = Oscar.action_function(Omega)
  for i in 1:length(Omega)
    @test Omega[i^pi] == fun(Omega[i], g)
  end

  # action homomorphism
  acthom = action_homomorphism(Omega)
  @test pi == g^acthom
  flag, pre = has_preimage_with_preimage(acthom, pi)
  @test flag
  @test pre == g
  @test order(image(acthom)[1]) == order(G)
  rest = restrict_homomorphism(acthom, derived_subgroup(G)[1])
  @test ! is_bijective(rest)

  # is_conjugate
  x, y = [right_coset(H, g) for g in gens(G)]
  @test is_conjugate(Omega, x, y)

  # is_conjugate_with_data
  rep = is_conjugate_with_data(Omega, x, y)
  @test rep[1]
  @test x * rep[2] == y
  @test Oscar.action_function(Omega)(x, rep[2]) == y
end

@testset "G-sets by left transversals" begin
  G = symmetric_group(5)
  H = sylow_subgroup(G, 2)[1]
  Omega = left_cosets(G, H)
  @test AbstractAlgebra.PrettyPrinting.repr_terse(Omega) == "Left cosets of groups"
  @test isa(Omega, GSet)
  @test acting_group(Omega) == G
  @test length(Omega) == index(G, H)
  @test order(stabilizer(Omega)[1]) * length(Omega) == order(G)
  @test Omega[end] == Omega[length(Omega)]
  @test length(orbits(Omega)) == 1
  @test is_transitive(Omega)
  @test ! is_regular(Omega)
  @test ! is_semiregular(Omega)

  @test eltype(Omega) == typeof(representative(Omega))

  # iteration
  for i in 1:length(Omega)
    @test findfirst(is_equal(Omega[i]), Omega) == i
  end
  rep = representative(Omega)
  for omega in Omega
    @test one(G) in omega || omega != rep
  end

  # orbit
  g = gen(G, 1)
  pnt = left_coset(H, g)
  @test pnt in Omega
  @test length(orbit(Omega, pnt)) == length(Oscar.orbit_via_Julia(Omega, pnt))

  # permutation
  pi = permutation(Omega, g)
  @test order(pi) == order(g)
  @test degree(parent(pi)) == length(Omega)
  fun = Oscar.action_function(Omega)
  for i in 1:length(Omega)
    @test Omega[i^pi] == fun(Omega[i], g)
  end

  # action homomorphism
  acthom = action_homomorphism(Omega)
  @test pi == g^acthom
  flag, pre = has_preimage_with_preimage(acthom, pi)
  @test flag
  @test pre == g
  @test order(image(acthom)[1]) == order(G)
  rest = restrict_homomorphism(acthom, derived_subgroup(G)[1])
  @test ! is_bijective(rest)

  # is_conjugate
  x, y = [left_coset(H, g) for g in gens(G)]
  @test is_conjugate(Omega, x, y)

  # is_conjugate_with_data
  rep = is_conjugate_with_data(Omega, x, y)
  @test rep[1]
  @test inv(rep[2]) * x == y
  @test Oscar.action_function(Omega)(x, rep[2]) == y
end

@testset "G-sets of PcGroups" begin
  G = small_group(24, 12)
  Omega = orbit(G, gen(G, 1))
  S, mp = stabilizer(Omega)
  @test length(Omega) == index(G, S)
end

@testset "G-sets of FinGenAbGroups" begin
  # Define an action on class functions.
  function galois_conjugate(chi::Oscar.GAPGroupClassFunction,
                            sigma::QQAbAutomorphism)
    return Oscar.class_function(parent(chi), [x^sigma for x in values(chi)])
  end

  # Compute Galois orbits on irreducible characters.
  t = character_table("L2(8)")
  N = lcm(map(conductor, t))
  u, mu = unit_group(quo(ZZ, N)[1])
  @test u isa FinGenAbGroup
  f = function(chi, g)
    return galois_conjugate(chi, QQAbAutomorphism(Int(lift(mu(g)))))
  end

  orb = @inferred orbit(u, f, t[3])
  @test length(collect(orb)) == 3

  Omega = @inferred gset(u, f, t)
  orbs = @inferred orbits(Omega)
  @test (@inferred length(orbs)) == 5
  @test sort(map(length, orbs)) == [1, 1, 1, 3, 3]
  @test all(o -> conductor(sum(collect(o))) == 1, orbs)
  o = orbs[findfirst(o -> length(o) == 3, orbs)]
  @test [order(permutation(o, x)) for x in gens(u)] == [1, 3]
end
