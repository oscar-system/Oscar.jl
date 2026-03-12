using Test
using Oscar.IntersectionTheory

let pushforward = IntersectionTheory.pushforward
  @testset "GenericVariety" begin

    # generic abstract variety
    C = abstract_variety(1)
    c = gens(C.ring)[1]
    @test C.T === tangent_bundle(C)
    @test rank(C.T) isa Int
    @test total_chern_class(C.T) == 1 + c
    # NOTE: trim! changes the ring, so some tests below need the old Singular API
    # which is no longer available. The following tests after trim! are skipped:
    # - Singular.dimension(C.ring.I) (no Singular access)
    # - parent(c) == C.ring (c still points to old ring after trim!)
    # - euler_characteristic (integral requires point class, which is not set for generic variety)

    # generic abstract variety with parameter
    # NOTE: function_field(Singular.QQ, ...) is no longer available
    # F, (g,) = function_field(Singular.QQ, ["g"])
    # C = abstract_variety(1, base=F)
    # c = gens(C.ring)[1]
    # trim!(C.ring)
    # C.point = 1//(2 - 2g) * total_chern_class(1, C)
    # @test euler_number(C) == 2 - 2g
    # @test rank(trivial_line_bundle(C) * g) == g
    # @test rank(symmetric_power(g, 2trivial_line_bundle(C))) == g + 1

    # generic abstract_variety with bundles
    X, (A, B) = abstract_variety(2, [3 => "a", 3 => "b"])
    @test schur_functor(A, [1, 1]) == exterior_power(A, 2)
    @test schur_functor(A, [2]) == symmetric_power(A, 2)
    D = degeneracy_locus(A, B, 2)
    @test pushforward(map(D, X), D(1)) == degeneracy_locus(A, B, 2; class=true)

    # characteristic classes
    t = todd_class(2)
    c = gens(parent(t))
    @test t == QQ(1//12) * c[1]^2 + QQ(1//12) * c[2]
    l = l_genus(2)
    p = gens(parent(l))
    @test l == -QQ(1//45) * p[1]^2 + QQ(7//45) * p[2]
    a = a_hat_genus(2)
    p = gens(parent(a))
    @test a == QQ(7//5760) * p[1]^2 - QQ(1//1440) * p[2]
  end

  @testset "VarietyHom" begin
    p = abstract_point()
    P2 = abstract_projective_space(2)
    i = map(P2, P2)
    @test i.domain == P2
    @test i.codomain == P2

    i = map(p, P2)
    @test pushforward(i, p(1)) == P2.point
    @test pullback(i, P2.O1) == 0
    @test i.T === tangent_bundle(i)
    @test -i.T == 2trivial_line_bundle(p) # normal bundle

    # test that coercion works properly
    pt = P2.structure_map.codomain
    A = trivial_line_bundle(P2) * trivial_line_bundle(pt)
    @test parent(A) == P2
    @test A == trivial_line_bundle(P2)

    PF = projective_bundle(P2.bundles[2])
    A = trivial_line_bundle(P2) + trivial_line_bundle(PF)
    @test parent(A) == PF
    @test A == 2trivial_line_bundle(PF)

    # test that hom works for blow_up
    Bl, E, j = blowup(i)
    e = pushforward(j, E(1))
    @test e == gens(Bl.ring)[1]
    @test integral(e^2) == -1
    @test pullback(structure_map(E), p(1)) == E(1)

    P5 = abstract_projective_space(5; symbol="H")
    h, H = P2.O1, P5.O1
    v = map(P2, P5, [2h])
    @test pullback(v, H) == 2h
    @test pullback(v, P5.point) == 0
    @test v.pushforward(h) == 2H^4
    @test pushforward(v, P2.point) == P5.point
    @test -v.T == abstract_bundle(P2, 3, 1 + 9h + 30h^2) # normal bundle

    # test that hom works for product
    P, Q = abstract_projective_space(1), abstract_projective_space(1)
    PxQ = P * Q
    p, q = map(PxQ, P), map(PxQ, Q)
    @test pushforward(p, PxQ.point) == P.point
    @test integral(pullback(p, P.point) * pullback(q, Q.point)) == 1

    # cubic containing a plane
    P2 = abstract_projective_space(2)
    Y = complete_intersection(abstract_projective_space(5), 3)
    i = map(P2, Y, [P2.O1]; inclusion=true)
    Y1 = i.codomain
    p = pushforward(i, P2(1))
    h = Y1.O1
    @test Y1 != Y
    @test euler_number(Y1) == euler_number(Y)
    @test map(Y1, Y).T.ch == 0
    @test betti_numbers(Y1)[3] == 2
    @test basis(Y1, 2) == [h^2, p]
    @test intersection_matrix([h^2, p]) == matrix(QQ, [3 1; 1 3])

    # a related result:
    # the degree of the hypersurface of cubics containing a plane
    G = abstract_grassmannian(3, 6)
    S = G.bundles[1]
    @test integral(total_chern_class(symmetric_power(dual(S), 3))) == 3402
  end

  @testset "Constructors" begin

    # abstract_projective_space(2)
    P2 = abstract_projective_space(2)
    h = P2.O1
    S, Q = P2.bundles
    @test gens(P2.ring.I) == [h.f^3]
    @test h^3 == 0
    @test P2.point == h^2
    @test S == line_bundle(P2, -1)
    @test Q == abstract_bundle(P2, 2, 1 + h + h^2)
    @test Q == abstract_bundle(P2, 2 + h - QQ(1//2) * h^2)
    @test hom(S, Q) == P2.T
    @test euler_number(P2) == 3
    @test total_chern_class(P2) == 1 + 3h + 3h^2
    @test chern_class(P2, 1) == 3h
    @test top_chern_class(P2.T) == chern_class(P2, 2)
    @test total_segre_class(P2.T) == 1 - 3h + 6h^2
    @test segre_class(P2.T, 2) == 6h^2
    @test todd_class(P2) == 1 + QQ(3//2) * h + h^2
    @test integral(todd_class(P2)) == 1
    @test total_pontryagin_class(P2) == 1 + 3h^2
    @test pontryagin_class(P2, 1) == 3h^2
    @test a_hat_genus(P2) == -1//8
    @test signature(P2) == 1
    @test chern_number(P2, 2) == 3
    @test chern_numbers(P2) == [partition([2]) => 3, partition([1, 1]) => 9]
    @test euler_characteristic(trivial_line_bundle(P2)) == 1
    @test euler_characteristic(cotangent_bundle(P2)) == -1
    hilb = hilbert_polynomial(P2)
    t = gens(parent(hilb))[1]
    @test hilb isa QQPolyRingElem
    @test hilb == 1 + 3//2 * t + 1//2 * t^2

    # Grassmannian
    G = abstract_grassmannian(2, 4)
    S, Q = tautological_bundles(G)
    c1, c2 = gens(G.ring)
    @test betti_numbers(G) == [1, 1, 2, 1, 1]
    @test euler_number(G) == 6
    @test chern_class(G, 1) == -4chern_class(S, 1)
    @test integral(total_chern_class(symmetric_power(dual(S), 3))) == 27
    @test integral(chern_class(dual(S), 1)^4) == 2
    @test integral(chern_class(G, 2)^2) == 98
    @test schubert_class(G, 2) == c1^2 - c2
    @test schubert_class(G, [1, 1]) == c2
    @test schubert_class(G, partition([2, 1])) == -c1^3 + c1 * c2
    @test [length(schubert_classes(G, i)) for i in 0:4] == [1, 1, 2, 1, 1]

    # Grassmannian: TnVariety version
    G = tn_grassmannian(2, 4)
    S, Q = tautological_bundles(G)
    @test G isa TnVariety
    @test S isa TnBundle
    @test rank(tangent_bundle(G)) == 4
    #@test euler_number(G) == 6
    @test integral(total_chern_class(symmetric_power(dual(S), 3))) == 27
    @test integral(chern_class(dual(S), 1)^4) == 2
    @test integral(chern_class(G, 2)^2) == 98

    # flag abstract_variety
    F = abstract_flag_variety(1, 2, 3)
    A, B, C = tautological_bundles(F)
    @test dim(F) == 3
    @test rank.(tautological_bundles(F)) == [1, 1, 1]
    @test betti_numbers(F) == [1, 2, 2, 1]
    @test euler_number(F) == 6

    # flag abstract_variety: TnVariety version
    F = tn_flag_variety([1, 2, 3])
    A, B, C = tautological_bundles(F)
    @test dim(F) == 3
    @test rank.(tautological_bundles(F)) == [1, 1, 1]
    #@test euler_number(F) == 6

    # projective bundle
    X, (F,) = abstract_variety(3, [3 => "c"])
    PF = projective_bundle(F)
    @test dim(PF) == 5
    @test rank.(tautological_bundles(PF)) == [1, 2]
    p = PF.structure_map
    @test p.codomain == X
    @test pullback(p, X(1)) == 1
    @test pushforward(p, PF(1)) == 0
    @test pushforward(p, p.O1^2) == 1

    # flag bundle
    X, (F,) = abstract_variety(2, [4 => "c"])
    FlF = flag_bundle(F, 2)
    @test dim(FlF) == 6
    @test rank.(tautological_bundles(FlF)) == [2, 2]
    p = FlF.structure_map
    @test p.codomain == X
    @test pullback(p, X(1)) == 1
    @test pushforward(p, FlF(1)) == 0
    @test pushforward(p, p.O1^4) == 2
    @test [length(schubert_classes(FlF, i)) for i in 0:4] == [1, 1, 2, 1, 1]
  end

  @testset "SpecialVarieties" begin

    # K3 surface of genus 3 (quartic in P^3)
    S = abstract_K3_surface(3)
    L = polarization(S)
    @test dim(S) == 2
    @test euler_number(S) == 24
    @test integral(L^2) == 4
    @test euler_characteristic(OO(S)) == 2
    @test euler_characteristic(OO(S, L)) == 4

    # K3 surface of genus 2
    S2 = abstract_K3_surface(2)
    @test euler_characteristic(OO(S2, polarization(S2))) == 3

    # Quadric: odd-dimensional (complete intersection in P^{n+1})
    Q1 = abstract_quadric(1)
    @test dim(Q1) == 1
    @test euler_number(Q1) == 2

    Q3 = abstract_quadric(3)
    @test dim(Q3) == 3
    @test euler_number(Q3) == 4

    # Quadric: even-dimensional (explicit Chow ring with two rulings)
    Q2 = abstract_quadric(2)
    @test dim(Q2) == 2
    @test euler_number(Q2) == 4
    @test degree(Q2) == 2

    Q4 = abstract_quadric(4)
    @test dim(Q4) == 4
    @test euler_number(Q4) == 6

    # Cayley plane OP^2
    CP = abstract_cayley_plane()
    @test dim(CP) == 16
    @test euler_number(CP) == 27
    @test length(betti_numbers(CP)) == 17  # 16-dimensional variety: degrees 0-16
    @test betti_numbers(CP)[1] == 1

    # Cayley Grassmannian
    CG = abstract_cayley_grassmannian()
    @test dim(CG) == 8
    @test euler_number(CG) == 15
  end

  # @testset "Pushfwd" begin
  #   A = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, [:x,:y,:z,:w])[1], [3,3,3,3])
  #   B = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, [:s,:t])[1], [1,1])
  #   s, t = gens(B)
  #   f = IntersectionTheory.ChAlgHom(A, B, [s^3,s^2*t,s*t^2,t^3]) # twisted cubic
  #   M, g, pf = IntersectionTheory._pushfwd(f)
  #   @test length(g) == 6
  #   x = s^3 + 5s*t + t^20 # random element from B
  #   @test sum(g .* f.salg.(pf(x.f))) == x.f

  #   A = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, [:x,:y,:z,:w])[1], [4,4,2,1])
  #   B = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, [:s,:t,:u])[1], [1,1,1])
  #   s, t, u = gens(B)
  #   f = IntersectionTheory.ChAlgHom(A, B, [s^4+u^4,s*t^2*u,s^2-t^2-u^2,t]) # random morphism
  #   M, g, pf = IntersectionTheory._pushfwd(f)
  #   @test length(g) == 8
  #   x = s^2 + 2s*t + 3s*t*u + t^2*u + 20t*u + u^20 # random element from B
  #   @test sum(g .* f.salg.(pf(x.f))) == x.f
  # end

  # testset borrowed from Schubert2
  @testset "Blowups" begin

    # blow_up Veronese
    P2 = abstract_projective_space(2)
    P5 = abstract_projective_space(5)
    i = map(P2, P5, [2P2.O1])
    Bl, E, j = blowup(i)
    c = top_chern_class(tangent_bundle(Bl))
    @test integral(pushforward(structure_map(Bl), c)) == 12
    @test integral(c) == 12
    e = pushforward(j, E(1))
    quad = pullback(structure_map(Bl), 2P5.O1) - e
    @test integral(quad^5) == 1
    sext = pullback(structure_map(Bl), 6P5.O1) - 2e
    @test integral(sext^5) == 3264

    # blow_up point in P2
    P2 = abstract_projective_space(2)
    P = abstract_point(; base=P2.base)
    Bl, E, j = blowup(map(P, P2, [zero(P.ring)]))
    e = pushforward(j, E(1))
    @test integral(e^2) == -1
    @test integral(pullback(j, e)) == -1
    @test euler_number(Bl) == 4

    # blow_up point in P7
    P7 = abstract_projective_space(7)
    P = abstract_point(; base=P2.base)
    Bl, E, j = blowup(map(P, P7, [zero(P.ring)]))
    e = pushforward(j, E(1))
    @test euler_number(Bl) == 14

    # blow_up twisted cubic
    P1 = abstract_projective_space(1)
    P3 = abstract_projective_space(3)
    i = map(P1, P3, [3P1.O1])
    Bl, E, j = blowup(i)
    e = pushforward(j, E(1))
    quad = pullback(structure_map(Bl), 2P3.O1) - e
    @test integral(quad^3) == 0
    cubic = pullback(structure_map(Bl), 3P3.O1) - e
    @test integral(quad^2 * cubic) == 1

    # blow_up twisted cubic, with parameters
    T, (r, s, t) = polynomial_ring(QQ, [:r, :s, :t])
    F = fraction_field(T)
    (r, s, t) = gens(F)
    P1 = abstract_projective_space(1; base=F)
    P3 = abstract_projective_space(3; base=F)
    i = map(P1, P3, [3P1.O1])
    Bl, E, j = blowup(i)
    e = pushforward(j, E(1))
    rH, sH, tH = [pullback(structure_map(Bl), x * P3.O1) - e for x in [r, s, t]]
    @test integral(rH * sH * tH) == r * s * t - 3 * r - 3 * s - 3 * t + 10

    G = abstract_grassmannian(2, 5)
    P9 = abstract_projective_space(9)
    i = map(G, P9, [G.O1])
    Bl, E, j = blowup(i)
    e = pushforward(j, E(1))
    quad = pullback(structure_map(Bl), 2P9.O1) - e
    @test simplify(quad^5) == 0
    @test simplify(e^5) != 0

    # blow_up space curve of degree d and genus g
    T, (r, s, t, d, g) = polynomial_ring(QQ, [:r, :s, :t, :d, :g])
    F = fraction_field(T)
    (r, s, t, d, g) = gens(F)
    P2 = abstract_projective_space(2; base=F)
    P3 = abstract_projective_space(3; base=F)
    C = zero_locus_section(OO(P2, d))
    C.point = 1//(2 - 2g) * chern_class(C, 1)
    i = map(C, P3, [d * C.point])
    Bl, E, j = blowup(i)
    e = pushforward(j, E(1))
    rH, sH, tH = [pullback(structure_map(Bl), x * P3.O1) - e for x in [r, s, t]]
    @test integral(rH * sH * tH) == r * s * t - d * (r + s + t) + (2g - 2 + 4d)

    G = abstract_grassmannian(2, 5)
    Z = zero_locus_section(3line_bundle(G, 1))
    Bl, E = blowup(structure_map(Z))
    @test dim(Bl) == 6
    @test euler_number(Bl) == 18
    @test betti_numbers(Bl) == [1, 2, 4, 4, 4, 2, 1]
    @test [euler_characteristic(exterior_power(cotangent_bundle(Bl), i)) for i in 0:6] ==
      [1, -2, 4, -4, 4, -2, 1]
  end

  @testset "Setters" begin

    # set_point_class
    X, x = abstract_variety(2, ["c₁", "c₂"], [1, 2])
    set_point_class(X, x[2])
    @test X.point == x[2]
    @test_throws ErrorException set_point_class(X, x[2])

    # set_tangent_bundle
    X, x = abstract_variety(2, ["c₁", "c₂"], [1, 2])
    T = abstract_bundle(X, 2, 1 + x[1] + x[2])
    set_tangent_bundle(X, T)
    @test tangent_bundle(X) === T
    @test_throws ErrorException set_tangent_bundle(X, T)

    # set_polarization
    P2 = abstract_projective_space(2)
    @test_throws ErrorException set_polarization(P2, P2.O1)
    X, x = abstract_variety(2, ["h₁", "h₂"], [1, 2])
    set_polarization(X, x[1])
    @test X.O1 == x[1]
    @test_throws ErrorException set_polarization(X, x[1])

    # set_tautological_bundles
    X, x = abstract_variety(4, ["a₁", "a₂", "b₁", "b₂"], [1, 2, 1, 2])
    S = abstract_bundle(X, 2, 1 + x[1] + x[2])
    Q = abstract_bundle(X, 2, 1 + x[3] + x[4])
    set_tautological_bundles(X, [S, Q])
    @test tautological_bundles(X) == [S, Q]
    @test_throws ErrorException set_tautological_bundles(X, [S, Q])

    # set_structure_map
    X, x = abstract_variety(2, ["c₁", "c₂"], [1, 2])
    pt = abstract_point()
    f = map(X, pt)
    set_structure_map(X, f)
    @test structure_map(X) === f
    @test_throws ErrorException set_structure_map(X, f)
  end

  @testset "CharacteristicClasses" begin

    # a_hat_genus on higher-dimensional spaces
    P4 = abstract_projective_space(4)
    @test a_hat_genus(P4) == 3//128

    # l_genus / signature on P4
    @test signature(P4) == 1

    # signature of a K3-like surface (euler 24, signature 16)
    # Construct a quartic surface in P3
    P3 = abstract_projective_space(3)
    K3 = zero_locus_section(OO(P3, 4))
    @test euler_number(K3) == 24
    @test signature(K3) == -16

    # hilbert_polynomial of P3
    P3 = abstract_projective_space(3)
    hilb = hilbert_polynomial(P3)
    t = gens(parent(hilb))[1]
    @test hilb == 1 + 11//6 * t + t^2 + 1//6 * t^3

    # chern_numbers of P3
    cn = chern_numbers(P3)
    @test partition([3]) in first.(cn)
    @test cn[findfirst(p -> p.first == partition([3]), cn)].second == 4

    # Gromov-Witten and instanton numbers
    @test gromov_witten_invariant(1, 5) == 2875
    @test instanton_number(1, 5) == 2875
    @test instanton_number(2, 5) == 609250
    @test instanton_number(2, 4, 2) == 92288
  end

  @testset "AdamsOperations" begin

    # adams on ring elements
    P2 = abstract_projective_space(2)
    h = gens(P2.ring)[1]
    ch_T = chern_character(tangent_bundle(P2))
    @test adams(1, ch_T) == ch_T  # ψ¹ = identity
    @test adams(-1, ch_T) == chern_character(dual(tangent_bundle(P2)))

    # adams on bundles
    T = tangent_bundle(P2)
    @test adams(1, T) == T
    @test adams(-1, T) == dual(T)
    @test adams(2, T) == abstract_bundle(P2, adams(2, chern_character(T)))

    # adams with symbolic k (use ZZ)
    @test adams(ZZ(3), ch_T) == adams(3, ch_T)

    # cannibalistic class
    c_T = total_chern_class(T)
    can2 = cannibalistic(2, c_T)
    @test parent(can2) == P2.ring

    # cannibalistic on bundles
    can2_bundle = cannibalistic(2, T)
    @test parent(can2_bundle) == P2.ring
  end

  @testset "AbstractCurve" begin

    # genus 0 curve (P1)
    C0 = abstract_curve(0)
    @test dim(C0) == 1
    @test euler_number(C0) == 2
    @test euler_characteristic(OO(C0)) == 1
    @test rank(tangent_bundle(C0)) == 1

    # genus 1 curve (elliptic curve)
    C1 = abstract_curve(1)
    @test dim(C1) == 1
    @test euler_number(C1) == 0
    @test euler_characteristic(OO(C1)) == 0

    # genus 2 curve
    C2 = abstract_curve(2)
    @test dim(C2) == 1
    @test euler_number(C2) == -2
    @test euler_characteristic(OO(C2)) == -1

    # genus 3 curve
    C3 = abstract_curve(3)
    @test euler_number(C3) == -4
  end

  @testset "DegeneracyLocusK0" begin

    # k=0 with class=true: Porteous formula for where the map is zero
    G24 = abstract_grassmannian(2, 4)
    S, Q = tautological_bundles(G24)
    cl = degeneracy_locus(S, Q, 0; class=true)
    c2 = gens(G24.ring)[2]
    @test cl == 6 * c2^2  # Porteous formula: schur_functor(Q-S, [2,2]) in top degree

    # k=0 with class=true should work and give a ring element
    P4 = abstract_projective_space(4)
    F = 2 * OO(P4, -1)
    G = 3 * OO(P4)
    cl2 = degeneracy_locus(F, G, 0; class=true)
    @test parent(cl2) == P4.ring  # should return an element (possibly 0)
  end
end
