@testset "CoveredSchemes" begin
  @testset "Covered schemes 1" begin
    R, (x,y) = polynomial_ring(QQ, [:x, :y])
    X = subscheme(spec(R), [x^2+y^2])
    P = projective_space(X, 3)
    S = homogeneous_coordinate_ring(P)
    (u, v) = gen(S, 1), gen(S, 2)
    h = u^3
    h = u^3 + u^2
    h = u^3 + (u^2)*v
    h = u^3 + u^2*v - OO(X)(x)*v^3
    Z = subscheme(P, h)
    C = standard_covering(Z)
    f = dehomogenization_map(Z, C[2])
    @test f(u) == gens(OO(C[2]))[1]
  end

  @testset "Covered schemes 2" begin
    P = projective_space(QQ, [:x, :y, :z, :w])
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
    @test Oscar.transition_graph(Pc[1]) isa Graph
    @test Oscar.transition_graph(C) isa Graph
  end

  @testset "standard_covering" begin
    R, t = polynomial_ring(QQ,[:t])
    T = Oscar.standard_spec(subscheme(spec(R),t))
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
    @test Oscar.poly_type(spec(R)) === Oscar.poly_type(spec(L)) === Oscar.poly_type(Ccov[1][2])
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
    @test length(Oscar.all_patches(C2)) == 3
    for (i, U) in zip(1:3, C2)
      @test U === C2[i]
    end

    U = C1[2]
    x = gens(ambient_coordinate_ring(U))[1]
    W = AffineSchemeOpenSubscheme(U, [x, x-1])
    W1, W2 = affine_patches(W)
    #add_affine_refinement!(C1, W)

    @test sprint(show, C1) isa String

    @test Oscar.base_ring_type(Ccov) == typeof(QQ)
    @test base_ring(Ccov) == QQ
    simplify!(Ccov)
    @test length(coverings(Ccov)) == 2
    C1, C3 = coverings(Ccov)
    @test codomain(Ccov[C1, C3]) === C3
    @test domain(Ccov[C1, C3]) === C1
    @test codomain(Ccov[C3, C1]) === C1
    @test domain(Ccov[C3, C1]) === C3
    @test gluings(Ccov) === gluings(C1)
    set_name!(Ccov, "C")
    @test name(Ccov) == "C"

    @test CoveredScheme(C1[1]) isa AbsCoveredScheme
    @test is_empty(Oscar.empty_covered_scheme(QQ))

    E = subscheme(IP2, ideal(S, gens(S)))
    Ecov = covered_scheme(E)
    @test is_empty(Ecov)

    @attributes mutable struct DummyCoveredScheme{BRT}<:AbsCoveredScheme{BRT}
      C::CoveredScheme
      function DummyCoveredScheme(C::CoveredScheme)
        return new{Oscar.base_ring_type(C)}(C)
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
    @test Oscar.refinements(D) == Oscar.refinements(Ccov)
    @test gluings(D) == gluings(Ccov)
    @test base_ring(D) == base_ring(Ccov)
  end

  @testset "closed embeddings and singular loci" begin
    IP2 = projective_space(QQ, [:x, :y, :z])
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
    A = spec(R)
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

  @testset "irreducible components" begin 
    P1 = projective_space(QQ,1)
    (s0,s1) = homogeneous_coordinates(P1)
    X = subscheme(P1,ideal(s0*s1))
    Xcov = covered_scheme(X)
    @test length(irreducible_components(Xcov))==2
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

    # Test the LazyGluings:
    X = covered_scheme(P)

    gg = covering_morphism(g_cov)
    dom_cov = domain(gg)
    for k in keys(gluings(dom_cov))
      @test underlying_gluing(gluings(dom_cov)[k]) isa SimpleGluing
    end
  end

  @testset "base change" begin
    kk, pr = quo(ZZ, 5)
    IP1 = covered_scheme(projective_space(ZZ, 1))
    IP1_red, red_map = base_change(pr, IP1)

    IP2 = projective_space(ZZ, 2)
    S = homogeneous_coordinate_ring(IP2)
    (x, y, z) = gens(S)
    I = ideal(S, x^2 + y^2 + z^2)
    IP2_cov = covered_scheme(IP2)
    II = IdealSheaf(IP2, I)

    inc_X = Oscar.CoveredClosedEmbedding(IP2_cov, II)
    (a, b, c) = base_change(pr, inc_X)
    @test compose(a, inc_X) == compose(b, c)
  end

  @testset "decomposition info" begin
    P3 = projective_space(ZZ, 3)
    X = covered_scheme(P3)
    kk = GF(29)
    X29, f = base_change(kk, X)
    @test Oscar.has_decomposition_info(default_covering(X29))
    orig_cov = default_covering(X)
    U = first(patches(orig_cov))
    x, y, z = gens(OO(U))
    V2 = PrincipalOpenSubset(U, x)
    V1 = PrincipalOpenSubset(U, x-1)
    new_cov = Covering(append!(AbsAffineScheme[V1, V2], patches(orig_cov)[2:end]))
    Oscar.inherit_gluings!(new_cov, orig_cov)
    Oscar.inherit_decomposition_info!(X, new_cov, orig_cov=orig_cov)
    @test Oscar.decomposition_info(new_cov)[V1] == [OO(V1)(x)]
  end

  @testset "fiber products of coverings" begin
    IP1 = projective_space(QQ, [:x, :y])
    S = homogeneous_coordinate_ring(IP1)
    (x, y) = gens(S)
    X = covered_scheme(IP1)
    cov = default_covering(X)
    f = identity_map(cov)
    cc, p1, p2 = fiber_product(f, f)
    Phi = hom(S, S, [x+y, x-y])
    phi = ProjectiveSchemeMor(IP1, IP1, Phi)
    g = covered_scheme_morphism(phi)
    g_cov = covering_morphism(g)
    cc, p1, p2 = fiber_product(g_cov, f)

    orig = default_covering(X)
    ref = domain(g_cov)

    inc_cc = Oscar.refinement_morphism(ref, orig)
    id_orig = identity_map(orig)

    fp_cc_orig, p1, p2 = fiber_product(inc_cc, id_orig)

    fp_orig_cc, p1, p2 = fiber_product(id_orig, inc_cc)

    fp_cc_cc = fiber_product(inc_cc, inc_cc)
  end

  @testset "composition and fiber products of morphisms of covered schemes" begin
    IP1 = projective_space(QQ, [:x, :y])
    S = homogeneous_coordinate_ring(IP1)
    (x, y) = gens(S)

    X = covered_scheme(IP1)
    id_X = identity_map(X)
    fiber_product(id_X, id_X)

    Phi = hom(S, S, [x+y, x-y])
    phi = ProjectiveSchemeMor(IP1, IP1, Phi)
    f = covered_scheme_morphism(phi)
    f2 = covered_scheme_morphism(ProjectiveSchemeMor(IP1, IP1, Phi)) # The same, but as a non-identical copy

    fiber_product(f, id_X)

    compose(f, id_X)
    compose(id_X, f)
    compose(f, f)
    f_cov = covering_morphism(f)
    f_cov2 = covering_morphism(f2)
    ref = Oscar.refinement_morphism(domain(f_cov), default_covering(X))
    ref2 = Oscar.refinement_morphism(domain(f_cov2), default_covering(X))
    _, _, f_cov_ref = fiber_product(f_cov, ref)
    _, _, f_cov_ref2 = fiber_product(f_cov, ref2)

    ff = CoveredSchemeMorphism(X, X, f_cov_ref)
    ff2 = CoveredSchemeMorphism(X, X, f_cov_ref2)
    compose(id_X, ff)
    compose(ff, id_X)
    compose(ff, ff)
    compose(ff, ff2)
    compose(ff2, ff)

    fiber_product(ff, f)
    fiber_product(f, ff)
    fiber_product(ff, ff)
    fiber_product(ff, ff2)
    fiber_product(ff2, ff)
  end

  @testset "normalization" begin
    # Example integral
    R, (x, y, z) = grade(QQ[:x, :y, :z][1])
    I = ideal(R, z*x^2 + y^3)
    X = covered_scheme(proj(R, I))
    @test !is_normal(X; check=false)
    N = normalization(X)
    # trigger the computation of some gluings
    Xnorm = N[1]
    @test is_normal(Xnorm; check=false)
    Cnorm = Xnorm[1] # a covering
    gluing_morphisms(Cnorm[1,2])
    gluing_morphisms(Cnorm[1,3])
    gluing_morphisms(Cnorm[2,3])
    gluing_morphisms(Cnorm[3,3])

    # Example non-integral, this also tests the function `disjoint_union`
    R, (x, y, z) = grade(QQ[:x, :y, :z][1])
    I = ideal(R, (z*x^2 + y^3)*(x))
    X = covered_scheme(proj(R, I))
    @test !is_normal(X; check=false)
    N = normalization(X)
    # trigger the computation of some gluings
    Xnorm = N[1]
    @test is_normal(Xnorm; check=false)
    Cnorm = Xnorm[1] # a covering
    gluing_morphisms(Cnorm[1,2])

    # A non-normal Enriques surface as constructed by Enriques himself
    S, (x0,x1,x2,x3) = graded_polynomial_ring(QQ,[:x0,:x1,:x2,:x3])
    J = ideal(S, [x1^2*x2^2*x3^2 + x0^2*x2^2*x3^2 + x0^2*x1^2*x3^2 + x0^2*x1^2*x2^2 + x0*x1*x2*x3*(x0^2+x1^2+2x0*x1+x2^2+x3^2)])
    X = proj(S, J)
    Xcov = covered_scheme(X)
    N = normalization(Xcov);
    Xnorm = N[1]
    Cnorm = Xnorm[1] # a covering
    gluing_morphisms(Cnorm[1,2])

    # Example non-reduced
    R, (x, y) = polynomial_ring(rational_field(), [:x, :y])
    I = ideal(R, x^2)
    X = covered_scheme(spec(R, I))
    @test !is_normal(X)
  end
end

