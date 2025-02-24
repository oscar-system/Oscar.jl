@testset "fraction fields of varieties" begin
  P = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(P)
  C = subscheme(P, ideal(S, S[1]*S[2]-S[3]^2))
  Ccov = covered_scheme(C)
  KK = VarietyFunctionField(Ccov)

  ConformanceTests.test_Field_interface(KK)
  #ConformanceTests.test_Field_interface_recursive(KK)  # FIXME: lots of ambiguity errors
end

@testset "fraction fields of varieties" begin
  R, (x,y,z) = QQ[:x, :y, :z]
  @test is_irreducible(spec(R))
  @test is_irreducible(spec(R, ideal(R, x)))
  @test !is_irreducible(spec(R, ideal(R, x*y)))
  @test is_irreducible(spec(localization(R, units_of(R))[1]))
  @test !is_irreducible(spec(R, ideal(R, x*y), units_of(R)))

  P = projective_space(QQ, 2)
  S = homogeneous_coordinate_ring(P)
  C = subscheme(P, ideal(S, S[1]*S[2]-S[3]^2))
  Ccov = covered_scheme(C)

  KK = VarietyFunctionField(Ccov)
  U2 = patches(Ccov)[2]
  a = 5*gens(OO(U2))[1]*gens(OO(U2))[2]

  b = KK(a, one(a))
  @test b^2 - b == KK(a^2-a, one(a))
end

@testset "fraction fields of varieties II" begin
  # We construct by hand the projective bundle
  # P(O_{P^1}(4)+O_{P^1}(6)+O_{P^1}(1))
  kk = GF(29)

  # Set up the base ℙ¹ with coordinates s and t
  S, _ = graded_polynomial_ring(kk, [:s, :t])

  base_P1 = proj(S)

  # split this into the standard covering
  base_covering = standard_covering(base_P1)

  A1s = patches(base_covering)[1]
  A1t = patches(base_covering)[2]

  # Set up relative projective space of relative dimension 2
  # over both base patches
  P2_s = projective_space(OO(A1s), ["xs", "ys", "zs"])

  Cs = standard_covering(P2_s)

  P2_t = projective_space(OO(A1t), ["xt", "yt", "zt"])

  Ct = standard_covering(P2_t)

  # Join the resulting schemes in a disjoint union with two
  # components
  C = disjoint_union(Cs, Ct)

  # Manually glue two (dense) patches of the two components
  X = Cs[3]
  Y = Ct[3]
  x = gens(OO(X))
  y = gens(OO(Y))
  f = maximal_extension(X, Y, [x[1]//(x[3])^4, x[2]//(x[3])^6, 1//x[3]])
  g = maximal_extension(Y, X, [y[1]//(y[3])^4, y[2]//(y[3])^6, 1//y[3]])
  add_gluing!(C, Gluing(X, Y, restrict(f, domain(f), domain(g)), restrict(g, domain(g), domain(f))))

  # Extend the gluing to the whole covered scheme
  Oscar.fill_transitions!(C)

  X = CoveredScheme(C)

  KK = VarietyFunctionField(X)

  U = representative_patch(KK)
  V = C[5]
  R = ambient_coordinate_ring(V)
  x = gens(R)
  f = x[1]^2 - 2*x[2]^5*x[3]^3
  g = 4*x[3]^2 - 5*x[2]
  @test KK(f, g) + KK(f, g) == KK(2*f//g)
  h = KK(f, g)
  @test h[V] == f//g
  K = fraction_field(R)
  @test K(h) == f//g
  @test KK(f, g) == KK(f//g)

  @test KK(f+2*g-5, g) + KK(f+2*g-5, g) == KK(2*(f+2*g-5)//g)
  h = KK(f+2*g-5, g)
  @test deepcopy(h)==h
  @test KK(1) == one(KK)
  @test iszero(zero(KK))
  @test h[V] == (f+2*g-5)//g
  @test divexact(h,h) == one(KK)
  @test h*inv(h) == one(KK)
  @test isone(h*inv(h))
  @test h^2 == h*h
  @test h^ZZRingElem(2)==h^2
  @test h^ZZRingElem(2)==h^Int8(2)
  @test KK(h[V]+h[V]) == h + h
  @test coefficient_ring(KK) === kk
  @test -h == KK(-h[V])
  @test 2//h == KK(2)//h
  @test 2//h == ZZ(2)//h
  @test elem_type(typeof(KK))==typeof(h)
  @test parent_type(typeof(h)) == typeof(KK)
  @test base_ring(KK)==base_ring(h)
  @test is_domain_type(typeof(h))
  @test is_exact_type(typeof(h))
  @test iszero(KK())
  @test KK(h)==h
  @test !iszero(KK(numerator(h)))
  @test_throws ErrorException KK(f,0*f)

  K = fraction_field(R)
  @test K(h) == (f+2*g-5)//g

  P2 = projective_space(QQ,2)
  S = homogeneous_coordinate_ring(P2)
  s0 = gen(S, 1)
  X = subscheme(P2, ideal(S, s0))
  Xc = covered_scheme(X)
  KX = function_field(Xc)
  @test is_irreducible(Xc) #fails

end

@testset "pullbacks for function fields" begin
  P = projective_space(QQ, [:x, :y, :z])
  (x, y, z) = gens(homogeneous_coordinate_ring(P))
  Y = covered_scheme(P)
  II = ideal_sheaf(P, [x,y])
  p = blow_up(II)
  X = domain(p)
  KY = function_field(Y)
  KX = function_field(X)
  U1 = first(affine_charts(Y))
  R1 = OO(U1)
  a = KY(R1[1], R1[2])
  b = pullback(projection(p))(a)
  U2 = affine_charts(Y)[2]
  R2 = OO(U2)
  a2 = KY(R2[1], R2[2])
  b2 = pullback(projection(p))(a2)
  V = representative_patch(KX)
  c = 5*a-a2^2
  pbc = pullback(projection(p))(c)
  d = 5*b - b2^2
  @test OO(V)(numerator(pbc)*denominator(d)) == OO(V)(numerator(d)*denominator(pbc))
end

@testset "refinements" begin
  P = projective_space(QQ, [:x, :y, :z])
  S = homogeneous_coordinate_ring(P)
  (x, y, z) = gens(S)
  Y = covered_scheme(P)
  I = ideal(S,[x^2-y*z])
  Q = subscheme(P, I)
  X = covered_scheme(Q)
  C = Oscar.simplified_covering(X)
  KK = function_field(X)
  for U in patches(C)
    x = first(gens(OO(U)))
    a = lifted_numerator(x)
    b = KK(a, one(a))
    invb = KK(one(a), a)
    @test isone(b*invb)
  end
  P1 = projective_space(QQ, ["a", "b"])
  S1 = homogeneous_coordinate_ring(P1)
  (a, b) = gens(S1)
  Phi = ProjectiveSchemeMor(P1, Q, [a*b, a^2, b^2])
  phi = covered_scheme_morphism(Phi)
  phi_cov = covering_morphism(phi)
  rho = X[codomain(phi_cov), C]
  psi_cov = compose(phi_cov, rho)
  L = domain(phi)
  psi = CoveredSchemeMorphism(L, X, psi_cov)

  U = first(patches(C))
  a = first(gens(OO(U)))
  a = lifted_numerator(a)
  h = KK(a, one(a))
  pullback(psi)(h)
end

