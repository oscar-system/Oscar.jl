@testset "algebra homomorphisms" begin
  r, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  s, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
  S = quo(s, ideal(s, [c-b^3]))[1]
  T, (X, Y, Z) = grade(r)
  U = [X, Y, Z]
  V = S.([2*a+b^6, 7*b-a^2, c^2])
  f = hom(r, S, V)
  @test f isa @inferred Oscar.affine_algebra_morphism_type(r, S)
  @test f isa @inferred Oscar.affine_algebra_morphism_type(typeof(r), typeof(S))
  g1 = hom(r, T, U)
  @test g1 isa @inferred Oscar.affine_algebra_morphism_type(r, T)
  @test g1 isa @inferred Oscar.affine_algebra_morphism_type(typeof(r), typeof(T))
  g2 = hom(T, T, U)
  @test g2 isa @inferred Oscar.affine_algebra_morphism_type(T)
  @test g2 isa @inferred Oscar.affine_algebra_morphism_type(typeof(T))
  K = kernel(f)
  R = quo(r, K)[1]
  phi = hom(R, S, V)
  psi = inverse(phi)
  W = psi.(gens(domain(psi)))
  prphi = preimage(phi, S(a))
  prpsi = preimage(psi, R(x))

  @test is_surjective(f) == true
  @test is_injective(f) == false
  @test is_bijective(f) == false
  @test is_surjective(phi) == true
  @test is_injective(phi) == true
  @test is_bijective(phi) == true
  @test prphi == R(1//2 * x-1//2*z)
  @test prpsi == S(2*a+c^2)
  @test W[1] == R(1//2 * x-1//2*z)
  @test W[2] == R(1//28*x^2 - 1//14*x*z + 1//7*y + 1//28*z^2)
  @test W[3] == R(1//21952*x^6 - 3//10976*x^5*z + 3//5488*x^4*y + 15//21952*x^4*z^2 - 3//1372*x^3*y*z - 5//5488*x^3*z^3 + 3//1372*x^2*y^2 + 9//2744*x^2*y*z^2 + 15//21952*x^2*z^4 - 3//686*x*y^2*z - 3//1372*x*y*z^3 - 3//10976*x*z^5 + 1//343*y^3 + 3//1372*y^2*z^2 + 3//5488*y*z^4 + 1//21952*z^6)

  # Something non-surjective
  f = hom(r, r, [x, x, x])
  @test !(@inferred is_surjective(f))

  r, (X, Y, Z) = polynomial_ring(QQ, [:X, :Y, :Z])
  s1, (a, b, c) = polynomial_ring(QQ, [:a, :b, :c])
  S1 = quo(s1, ideal(s1, [c-b^3]))[1]
  V = S1.([2*a+b^6, 7*b-a^2, c^2])
  f1 = hom(r, S1, V)

  S2 = quo(r, ideal(r, [X*Y]))[1]
  a, b, c = S2(X), S2(Y), S2(Z)
  f2 = hom(r, S2, [(a*b)^3+a^2+c, b^2-1, c^3])

  r2 , (x, y) = polynomial_ring(QQ, [:x, :y])
  s3 = quo(r, ideal(r, [X^2-Y^2*Z]))[1]
  f3 = hom(r2, s3, [gen(s3, 1), gen(s3, 2)])

  @test is_finite(f1) == true
  @test is_finite(f2) == true
  @test is_finite(f3) == false

  begin
    # #655
    R, vars = QQ[:x, :y]
    x = vars[1]
    y = vars[2]
    f = hom(R, R, vars)
    push!(vars, x)
    @test f(x) == x
  end

  let
    R, (xx, yy) = QQ[:xx, :yy]
    S, (x, y, z) = QQ[:x, :y, :z]
    f = hom(R, S, [x, y])
    @test preimage(f, ideal(S, [x, y, z])) == ideal(R, [xx, yy])
    @test preimage(f, ideal(S, [z])) == ideal(R, [zero(R)])
  end
end

@testset "cross-type kernels" begin
  R, (x, y) = QQ[:x, :y]
  S, _ = grade(R)
  A, _ = quo(R, ideal(R, zero(R)))
  Q, _ = quo(S, ideal(S, zero(S)))
  phi = hom(Q, A, gens(A))
  @test iszero(kernel(phi))
end

@testset "Laurent codomain" begin
  R, = polynomial_ring(QQ, 4)
  S, (x1, x2) = laurent_polynomial_ring(QQ, 2)
  f = hom(R, S, [x1, x2, x1^-2*x2^5, x1^-1*x2^3])
  K = kernel(f)
  @test all(is_zero, f.(gens(K)))
  Q, = quo(R, K)
  g = hom(Q, S, [x1, x2, x1^-2*x2^5, x1^-1*x2^3])
  @test is_zero(kernel(g))
end



@testset "Polynomial ring homomorphisms" begin
    ########################################################################
    #   Ungraded -> ungraded
    ########################################################################
    let
        R, (x, y) = polynomial_ring(QQ, [:x, :y])
        S, (u, v) = polynomial_ring(QQ, [:u, :v])

        F = hom(R, S, [u, v^2])

        @test domain(F) === R
        @test codomain(F) === S
        @test F(x + y) == u + v^2
        @test Oscar.hom_kind(F) == Oscar.HomUngraded
    end

    ########################################################################
    #   Graded -> graded, identity grading map, shift 0
    ########################################################################
    let
        Z = abelian_group(0)
        RG, (xg, yg) = graded_polynomial_ring(QQ, [:xg, :yg]; weights=[Z[1], Z[1]])
        SG, (ug, vg) = graded_polynomial_ring(QQ, [:ug, :vg]; weights=[Z[1], Z[1]])

        Gmap1 = hom(RG, SG, [ug, vg])

        @test domain(Gmap1) === RG
        @test codomain(Gmap1) === SG
        @test Gmap1(xg) == ug
        @test Gmap1(yg) == vg
        @test Oscar.hom_kind(Gmap1) == Oscar.HomGraded

        phi1 = get_attribute(Gmap1, :grading_group_hom, nothing)
        @test phi1 !== nothing
        @test phi1(degree(xg)) == degree(ug)
        @test phi1(degree(yg)) == degree(vg)

        dsh1 = get_attribute(Gmap1, :degree_shift, nothing)
        @test dsh1 !== nothing
        @test dsh1 == zero(grading_group(SG))
    end

    let
        RG, (xg, yg) = graded_polynomial_ring(QQ, [:xg, :yg])
        SG, (ug, vg) = graded_polynomial_ring(QQ, [:ug, :vg])

        Gmap1 = hom(RG, SG, [ug, vg])

        @test domain(Gmap1) === RG
        @test codomain(Gmap1) === SG
        @test Gmap1(xg) == ug
        @test Gmap1(yg) == vg
        @test Oscar.hom_kind(Gmap1) == Oscar.HomGraded

        phi1 = get_attribute(Gmap1, :grading_group_hom, nothing)
        @test phi1 !== nothing
        @test phi1(degree(xg)) == degree(ug)
        @test phi1(degree(yg)) == degree(vg)

        dsh1 = get_attribute(Gmap1, :degree_shift, nothing)
        @test dsh1 !== nothing
        @test dsh1 == zero(grading_group(SG))
    end

    ########################################################################
    #   Graded -> graded with nontrivial map, shift 0
    ########################################################################
    let
        R2, (x2,) = graded_polynomial_ring(QQ, [:x2])
        S2, (u2,) = graded_polynomial_ring(QQ, [:u2]; weights=[2])

        GR2 = grading_group(R2)
        GS2 = grading_group(S2)

        genS = gens(GS2)[1]
        phiG = hom(GR2, GS2, [2*genS])

        Fdeg = hom(R2, S2, [u2]; grading_group_hom = phiG)

        @test Oscar.hom_kind(Fdeg) == Oscar.HomGraded
        @test Fdeg(x2) == u2
        @test degree(Fdeg(x2)) == phiG(degree(x2))

        dsh = get_attribute(Fdeg, :degree_shift, nothing)
        @test dsh !== nothing
        @test dsh == zero(GS2)
    end

    ########################################################################
    #   Graded -> graded with non-zero degree shift
    ########################################################################
    let
        Z = abelian_group(0)
        Rsh, (x1, y1) = graded_polynomial_ring(QQ, [:x1, :y1]; weights=[Z[1], 2*Z[1]])
        Ssh, (u1,)    = graded_polynomial_ring(QQ, [:u1]; weights=[Z[1]])

        # deg(x1)=1 -> deg(u1^2)=2, deg(y1)=2 -> deg(u1^3)=3, shift=1
        Fsh = hom(Rsh, Ssh, [u1^2, u1^3])

        @test Oscar.hom_kind(Fsh) == Oscar.HomGraded
        @test Fsh(x1) == u1^2
        @test Fsh(y1) == u1^3

        dsh = get_attribute(Fsh, :degree_shift, nothing)
        @test dsh == degree(u1)

        # Incompatible images should error
        @test_throws ArgumentError hom(Rsh, Ssh, [u1^2, u1^4])
    end

    ########################################################################
    #   Graded -> ungraded
    ########################################################################
    let
        Rg, (xg2, yg2) = graded_polynomial_ring(QQ, [:xg2, :yg2])
        Su, (uu, vu)   = polynomial_ring(QQ, [:uu, :vu])

        F_GU = hom(Rg, Su, [uu, vu^2])

        @test Oscar.hom_kind(F_GU) == Oscar.HomGradedToUngraded
        @test F_GU(xg2) == uu
        @test F_GU(yg2) == vu^2
    end

    ########################################################################
    #   Ungraded -> graded
    ########################################################################
    let
        Ru, (xu, yu)   = polynomial_ring(QQ, [:xu, :yu])
        Sg, (ug2, vg2) = graded_polynomial_ring(QQ, [:ug2, :vg2])

        inhom = ug2 + vg2^2
        F_UG = hom(Ru, Sg, [ug2, inhom])

        @test Oscar.hom_kind(F_UG) == Oscar.HomUngradedToGraded
        @test F_UG(xu) == ug2
        @test F_UG(yu) == inhom
    end

    ########################################################################
    #   Graded -> graded with different grading groups
    ########################################################################
    let
        GZ  = abelian_group(0)
        GZ2 = abelian_group(ZZMatrix(2,2,[0,0,0,0]))

        RZ, (xz1,) = graded_polynomial_ring(QQ, [:xz1]; weights=[gen(GZ,1)])
        SZ2, (uz1,) = graded_polynomial_ring(QQ, [:uz1]; weights=[gen(GZ2,1)])

        @test_throws ErrorException hom(RZ, SZ2, [uz1])

        phiZ = hom(GZ, GZ2, [gen(GZ2, 1)])
        FGZ  = hom(RZ, SZ2, [uz1]; grading_group_hom = phiZ)

        @test Oscar.hom_kind(FGZ) == Oscar.HomGraded
        @test degree(FGZ(xz1)) == phiZ(degree(xz1))
        @test get_attribute(FGZ, :degree_shift, nothing) == zero(grading_group(SZ2))
    end

    ########################################################################
    #    Composition: graded->ungraded then ungraded->graded yields an ungraded map
    #    Kernel and image are illegal for such maps
    ########################################################################
    let
        Rg, (x,) = graded_polynomial_ring(QQ, [:x])
        Su, (s,) = polynomial_ring(QQ, [:s])
        Tg, (t,) = graded_polynomial_ring(QQ, [:t])

        f = hom(Rg, Su, [s])
        g = hom(Su, Tg, [t + t^2])

        @test Oscar.hom_kind(f) == Oscar.HomGradedToUngraded
        @test Oscar.hom_kind(g) == Oscar.HomUngradedToGraded

        h = compose(f, g)

        @test Oscar.hom_kind(h) == Oscar.HomUngraded
        @test domain(h) === Rg
        @test codomain(h) === Tg
        @test h(x) == t + t^2

        # direct graded->graded with inhomogeneous images is rejected
        @test_throws Exception hom(Rg, Tg, [t + t^2])
    end
end
