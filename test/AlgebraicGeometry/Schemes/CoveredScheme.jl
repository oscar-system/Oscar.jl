@testset "Covered schemes 1" begin
  R, (x,y) = polynomial_ring(QQ, ["x", "y"])
  X = subscheme(Spec(R), [x^2+y^2])
  P = projective_space(X, 3)
  S = homogeneous_coordinate_ring(P)
  (u, v) = gens(S)[1], gens(S)[2]
  h = u^3 
  h = u^3 + u^2
  h = u^3 + (u^2)*v 
  h = u^3 + u^2*v - OO(X)(x)*v^3
  Z = subscheme(P, h)
  C = standard_covering(Z)
  f = dehomogenization_map(Z, 1)
  @test f(u) == gens(OO(C[2]))[1]
end

@testset "Covered schemes 2" begin
  P = projective_space(QQ, ["x", "y", "z", "w"])
  Pc = covered_scheme(P)
  S = homogeneous_coordinate_ring(P)
  (x,y,z,w) = gens(S)
  X = subscheme(P, [x*w - y*z])
  @test dim(Pc)==3
  @test dim(covered_scheme(X))==2
  Y = subscheme(P, [x*z,y*z])
  @test dim(covered_scheme(Y)) == 2
  C = standard_covering(X)
  D, i, j = simplify(C) 
  @test all( x->(ngens(ambient_coordinate_ring(x)) == 2), collect(D))
  @test transition_graph(Pc[1]) isa Graph
  @test transition_graph(C) isa Graph
end

@testset "standard_covering" begin
  R, t = polynomial_ring(QQ,["t"])
  T = Oscar.standard_spec(subscheme(Spec(R),t))
  Pt= projective_space(T, 2)
  X = covered_scheme(Pt)
  @test dim(X) == 2
end

@testset "covered schemes 3" begin
  IP2 = projective_space(QQ, 2, var_name="u")
  S = homogeneous_coordinate_ring(IP2)
  u0, u1, u2 = gens(S)
  C = subscheme(IP2, u0^2 - u1*u2)

  Ccov = covered_scheme(C)
  #IP2cov = covered_scheme(IP2)

  R = base_ring(OO(Ccov[1][2]))
  L, _ = localization(R, R[1])
  @test poly_type(Spec(R)) === poly_type(Spec(L)) === poly_type(Ccov[1][2])
  Lnew, f, g = simplify(L)
  @test !(L == Lnew)
  @test compose(f, g) == identity_map(L)
  @test compose(g, f) == identity_map(Lnew)

  C1 = default_covering(Ccov)
  C2, f, g = simplify(C1)
  tmp1 = compose(f, g)
  tmp2 = compose(g, f) 
  @test domain(tmp1) === codomain(tmp1) === C2
  @test domain(tmp2) === codomain(tmp2) === C1
  @test length(C2) == 3
  @test length(all_patches(C2)) == 3
  for (i, U) in zip(1:3, C2)
    @test U === C2[i]
  end

  U = C1[2]
  x = gens(ambient_coordinate_ring(U))[1]
  W = SpecOpen(U, [x, x-1])
  W1, W2 = affine_patches(W)
  #add_affine_refinement!(C1, W)

  @test sprint(show, C1) isa String

  @test base_ring_type(Ccov) == typeof(QQ)
  @test base_ring(Ccov) == QQ
  simplify!(Ccov)
  @test length(coverings(Ccov)) == 2
  C1, C3 = coverings(Ccov)
  @test codomain(Ccov[C1, C3]) === C3
  @test domain(Ccov[C1, C3]) === C1
  @test codomain(Ccov[C3, C1]) === C1
  @test domain(Ccov[C3, C1]) === C3
  @test glueings(Ccov) === glueings(C1)
  set_name!(Ccov, "C")
  @test name(Ccov) == "C"

  @test CoveredScheme(C1[1]) isa AbsCoveredScheme
  @test is_empty(empty_covered_scheme(QQ))

  E = subscheme(IP2, ideal(S, gens(S)))
  Ecov = covered_scheme(E)
  @test is_empty(Ecov)

  @attributes mutable struct DummyCoveredScheme{BRT}<:AbsCoveredScheme{BRT}
    C::CoveredScheme
    function DummyCoveredScheme(C::CoveredScheme)
      return new{base_ring_type(C)}(C)
    end
  end
  Oscar.underlying_scheme(D::DummyCoveredScheme) = D.C

  D = DummyCoveredScheme(Ccov)
  @test coverings(D) == coverings(Ccov)
  @test default_covering(D) == default_covering(Ccov)
  @test patches(D) == patches(Ccov)
  @test Ccov[1][1] in D
  @test D[Ccov[1]] == Ccov[Ccov[1]]
  @test D[Ccov[1], Ccov[2]] == Ccov[Ccov[1], Ccov[2]]
  @test refinements(D) == refinements(Ccov)
  @test glueings(D) == glueings(Ccov)
  @test base_ring(D) == base_ring(Ccov)
end

@testset "closed embeddings and singular loci" begin
  IP2 = projective_space(QQ, ["x", "y", "z"])
  S = homogeneous_coordinate_ring(IP2)
  (x, y, z) = gens(S)
  f = x^2*z + y^3 - y^2*z
  C = subscheme(IP2, ideal(S, f))
  Ccov = covered_scheme(C)
  Csing, inc = singular_locus_reduced(Ccov)
  for V in affine_charts(Csing)
    U = codomain(inc[V])
    @test V == subscheme(U, image_ideal(inc)(U))
  end
end


@testset "is_integral" begin
  P = projective_space(QQ, 2)

  S = homogeneous_coordinate_ring(P)
  u, v, w = gens(S)
  I = ideal(S, [u*(v^2 + w^2 + u^2)])

  X = subscheme(P, I)
  Xcov = covered_scheme(X)
  @test !is_irreducible(Xcov)
  @test is_reduced(Xcov)
  @test !is_integral(Xcov)
  U = Xcov[1][1]
  @test is_integral(U)
  @test is_irreducible(U)
  @test is_reduced(U)
  @test is_integral(hypersurface_complement(U, OO(U)[2]))

  Y = subscheme(P, ideal(S, v^2 + w^2 + u^2))
  Ycov = covered_scheme(Y)
  @test is_irreducible(Ycov)
  @test is_integral(Ycov)

  R = forget_grading(S)
  A = Spec(R)
  @test is_integral(A)
  @test is_integral(hypersurface_complement(A, R[1]))

  # The following produces a covered scheme consisting of three disjoint affine patches
  Ysep = CoveredScheme(Covering(affine_charts(Ycov)))
  @test !is_irreducible(Ysep)
  @test is_reduced(Ysep)
  @test !is_integral(Ysep)

  # Two points, one reduced, the other not
  J1 = ideal(S, [u^2, v])
  J2 = ideal(S, [u-w, v-w])

  Z = subscheme(P, J1*J2)
  Zcov = covered_scheme(Z)
  @test !is_integral(Zcov)
  @test !is_irreducible(Zcov)

  Z1 = subscheme(P, J1)
  Z1cov = covered_scheme(Z1)
  @test is_irreducible(Z1cov)
  @test !is_reduced(Z1cov)
  @test !is_integral(Z1cov)

  Z2 = subscheme(P, J2)
  Z2cov = covered_scheme(Z2)
  @test is_integral(Z2cov)
end

@testset "conversion of morphisms" begin
  P = projective_space(QQ, 2)
  SP = homogeneous_coordinate_ring(P)
  (x, y, z) = gens(SP)
  m = ideal(SP, gens(SP))
  m3 = m^3
  n = ngens(m3)
  Q = projective_space(QQ, n-1)
  SQ = homogeneous_coordinate_ring(Q)
  phi = hom(SQ, SP, gens(m3))
  f = ProjectiveSchemeMor(P, Q, phi)
  f_cov = covered_scheme_morphism(f)
  @test domain(f_cov) === covered_scheme(P)
  @test codomain(f_cov) === covered_scheme(Q)

  A = [1 2 4; 1 3 9; 1 5 25]
  v = A*[x, y, z]
  psi = hom(SP, SP, v)
  g = ProjectiveSchemeMor(P, P, psi)
  g_cov = covered_scheme_morphism(g)

  @test domain(g_cov) === codomain(g_cov) === covered_scheme(P)

  # Test the LazyGlueings:
  X = covered_scheme(P)

  gg = covering_morphism(g_cov)
  dom_cov = domain(gg)
  for k in keys(glueings(dom_cov))
      @test underlying_glueing(glueings(dom_cov)[k]) isa SimpleGlueing
  end
end
