using Oscar.IntersectionTheory

let pushforward = IntersectionTheory.pushforward
  @testset "GenericVariety" begin
    
    # # generic abstract_variety
    # C = abstract_variety(1)
    # c = gens(C.ring)[1]
    # @test C.T === tangent_bundle(C)
    # @test rank(C.T) isa Int
    # @test total_chern_class(C.T) == 1 + c
    # trim!(C.ring)
    # @test Singular.dimension(C.ring.I) == 0
    # @test parent(c) == C.ring
    # @test betti(C) == [1, 1]
    # @test basis(C) == [[C.ring(1)], [c]]
    # @test euler(C) == c
    # @test euler_characteristic(trivial_line_bundle(C)) == 1//2 * c

    # # generic abstract_variety with parameter
    # F, (g,) = function_field(Singular.QQ, ["g"])
    # C = abstract_variety(1, base=F)
    # c = gens(C.ring)[1]
    # trim!(C.ring)
    # C.point = 1//(2 - 2g) * total_chern_class(1, C)
    # @test euler(C) == 2 - 2g
    # @test rank(trivial_line_bundle(C) * g) == g
    # @test rank(symmetric_power(g, 2trivial_line_bundle(C))) == g + 1

    # generic abstract_variety with bundles
    X, (A, B) = abstract_variety(2, [3=>"a", 3=>"b"])
    @test schur_functor(A, [1,1]) == exterior_power(A, 2)
    @test schur_functor(A, [2]) == symmetric_power(A, 2)
    D = degeneracy_locus(2, A, B)
    @test pushforward(hom(D, X), D(1)) == degeneracy_locus(2, A, B, class=true)

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

    p = IntersectionTheory.point()
    P2 = IntersectionTheory.abstract_projective_space(2)
    i = hom(P2, P2)
    @test i.domain == P2
    @test i.codomain == P2

    i = hom(p, P2)
    @test pushforward(i, p(1)) == P2.point
    @test pullback(i, P2.O1) == 0
    @test i.T === tangent_bundle(i)
    @test -i.T == 2trivial_line_bundle(p) # normal bundle

    # test that coercion works properly
    pt = P2.struct_map.codomain
    A = trivial_line_bundle(P2) * trivial_line_bundle(pt)
    @test parent(A) == P2
    @test A == trivial_line_bundle(P2)

    PF = abstract_projective_bundle(P2.bundles[2])
    A = trivial_line_bundle(P2) + trivial_line_bundle(PF)
    @test parent(A) == PF
    @test A == 2trivial_line_bundle(PF)

    # # test that hom works for blowup
    # Bl, E = blowup(i)
    # e = pushforward(E → Bl, E(1))
    # @test e == gens(Bl.ring)[1]
    # @test integral(e^2) == -1
    # @test pullback(E → p, p(1)) == E(1)

    P5 = abstract_projective_space(5, symbol="H")
    h, H = P2.O1, P5.O1
    v = hom(P2, P5, [2h])
    @test pullback(v, H) == 2h
    @test pullback(v, P5.point) == 0
    @test v.pushforward(h) == 2H^4
    @test pushforward(v, P2.point) == P5.point
    @test -v.T == abstract_bundle(P2, 3, 1 + 9h + 30h^2) # normal bundle

    # test that hom works for product
    P, Q = abstract_projective_space(1), abstract_projective_space(1)
    PxQ = P * Q
    p, q = hom(PxQ, P), hom(PxQ, Q)
    @test pushforward(p, PxQ.point) == P.point
    @test integral(pullback(p, P.point) * pullback(q, Q.point)) == 1

    # # cubic containing a plane
    # P2 = abstract_projective_space(2)
    # Y = complete_intersection(abstract_projective_space(5), 3)
    # i = hom(P2, Y, [P2.O1], inclusion=true)
    # Y1 = i.codomain
    # p = pushforward(i, P2(1))
    # h = Y1.O1
    # @test Y1 != Y
    # @test euler(Y1) == euler(Y)
    # @test (Y1 → Y).T.ch == 0
    # @test betti(Y1)[3] == 2
    # @test basis(2, Y1) == [h^2, p]
    # @test intersection_matrix([h^2, p]) == Nemo.matrix(QQ, [3 1; 1 3])

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
    @test Q == abstract_bundle(P2, 2 + h - QQ(1//2)*h^2)
    @test hom(S, Q) == P2.T
    @test euler(P2) == 3
    @test total_chern_class(P2) == 1 + 3h + 3h^2
    @test chern_class(P2, 1) == 3h
    @test top_chern_class(P2.T) == chern_class(P2, 2)
    # @test segre_class(P2.T) == 1 - 3h + 6h^2
    # @test segre_class(P2.T, 2) == 6h^2
    @test todd_class(P2) == 1 + QQ(3//2)*h + h^2
    @test integral(todd_class(P2)) == 1
    @test total_pontryagin_class(P2) == 1 + 3h^2
    @test pontryagin_class(P2, 1) == 3h^2
    @test a_hat_genus(P2) == -1//8
    @test signature(P2) == 1
    @test chern_number(P2, 2) == 3
    @test chern_numbers(P2) == [partition([2]) => 3, partition([1,1]) => 9]
    @test euler_characteristic(trivial_line_bundle(P2)) == 1
    @test euler_characteristic(cotangent_bundle(P2)) == -1
    hilb = hilbert_polynomial(P2)
    t = gens(parent(hilb))[1]
    @test hilb isa QQPolyRingElem
    @test hilb == 1 + 3//2*t + 1//2*t^2

    # Grassmannian
    G = abstract_grassmannian(2, 4)
    S, Q = bundles(G)
    c1, c2 = gens(G.ring)
    @test betti(G) == [1,1,2,1,1]
    @test euler(G) == 6
    @test chern_class(G, 1) == -4chern_class(S, 1)
    @test integral(total_chern_class(symmetric_power(dual(S), 3))) == 27
    @test integral(chern_class(dual(S), 1)^4) == 2
    @test integral(chern_class(G, 2)^2) == 98
    @test schubert_class(G, 2) == c1^2-c2
    @test schubert_class(G, [1, 1]) == c2
    @test schubert_class(G, partition([2, 1])) == -c1^3 + c1 * c2
    @test [length(schubert_classes(G, i)) for i in 0:4] == [1,1,2,1,1]

    # Grassmannian: TnVariety version
    G = abstract_grassmannian(2, 4, bott=true)
    S, Q = bundles(G)
    @test G isa IntersectionTheory.TnVariety
    @test S isa IntersectionTheory.TnBundle
    @test rank(tangent_bundle(G)) == 4
    @test euler(G) == 6
    @test integral(total_chern_class(symmetric_power(dual(S), 3))) == 27
    @test integral(chern_class(dual(S), 1)^4) == 2
    @test integral(chern_class(G, 2)^2) == 98

    # flag abstract_variety
    F = abstract_flag_variety(1, 2, 3)
    A, B, C = bundles(F)
    @test dim(F) == 3
    @test rank.(bundles(F)) == [1, 1, 1]
    @test betti(F) == [1,2,2,1]
    @test euler(F) == 6

    # flag abstract_variety: TnVariety version
    F = abstract_flag_variety(1, 2, 3, bott=true)
    A, B, C = bundles(F)
    @test dim(F) == 3
    @test rank.(bundles(F)) == [1, 1, 1]
    @test euler(F) == 6

    # projective bundle
    X, (F,) = abstract_variety(3, [3=>"c"])
    PF = abstract_projective_bundle(F)
    @test dim(PF) == 5
    @test rank.(bundles(PF)) == [1, 2]
    p = PF.struct_map
    @test p.codomain == X
    @test pullback(p, X(1)) == 1
    @test pushforward(p, PF(1)) == 0
    @test pushforward(p, p.O1^2) == 1
    
    # flag bundle
    X, (F,) = abstract_variety(2, [4=>"c"])
    FlF = abstract_flag_variety(F, 2)
    @test dim(FlF) == 6
    @test rank.(bundles(FlF)) == [2, 2]
    p = FlF.struct_map
    @test p.codomain == X
    @test pullback(p, X(1)) == 1
    # @test pushforward(p, FlF(1)) == 0
    # @test pushforward(p, p.O1^4) == 2
    @test [length(schubert_classes(FlF, i)) for i in 0:4] == [1,1,2,1,1]

  end

  # @testset "Pushfwd" begin
  #   A = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, ["x","y","z","w"])[1], [3,3,3,3])
  #   B = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, ["s","t"])[1], [1,1])
  #   s, t = gens(B)
  #   f = IntersectionTheory.ChAlgHom(A, B, [s^3,s^2*t,s*t^2,t^3]) # twisted cubic
  #   M, g, pf = IntersectionTheory._pushfwd(f)
  #   @test length(g) == 6
  #   x = s^3 + 5s*t + t^20 # random element from B
  #   @test sum(g .* f.salg.(pf(x.f))) == x.f
     
  #   A = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, ["x","y","z","w"])[1], [4,4,2,1])
  #   B = IntersectionTheory.ChRing(polynomial_ring(Singular.QQ, ["s","t","u"])[1], [1,1,1])
  #   s, t, u = gens(B)
  #   f = IntersectionTheory.ChAlgHom(A, B, [s^4+u^4,s*t^2*u,s^2-t^2-u^2,t]) # random morphism
  #   M, g, pf = IntersectionTheory._pushfwd(f)
  #   @test length(g) == 8
  #   x = s^2 + 2s*t + 3s*t*u + t^2*u + 20t*u + u^20 # random element from B
  #   @test sum(g .* f.salg.(pf(x.f))) == x.f
  # end

  # # testset borrowed from Schubert2
  # @testset "Blowup" begin
    
  #   # blowup Veronese
  #   P2 = abstract_projective_space(2)
  #   P5 = abstract_projective_space(5)
  #   v = hom(P2, P5, [2P2.O1])
  #   Bl, E = blowup(v)
  #   c = top_chern_class(tangent_bundle(Bl))
  #   @test integral(pushforward(Bl → P5, c)) == 12
  #   @test integral(c) == 12
  #   e = pushforward(E → Bl, E(1))
  #   quad = pullback(Bl → P5, 2P5.O1) - e
  #   @test integral(quad^5) == 1
  #   sext = pullback(Bl → P5, 6P5.O1) - 2e
  #   @test integral(sext^5) == 3264
    
  #   # blowup point in P2
  #   P2 = abstract_projective_space(2)
  #   Bl, E = blowup(point() → P2)
  #   e = pushforward(E → Bl, E(1))
  #   @test integral(e^2) == -1
  #   @test integral(pullback(E → Bl, e)) == -1
  #   @test euler(Bl) == 4

  #   # blowup point in P7
  #   P7 = abstract_projective_space(7)
  #   Bl, E = blowup(point() → P7)
  #   e = pushforward(E → Bl, E(1))
  #   @test euler(Bl) == 14
    
  #   # blowup twisted cubic
  #   P1 = abstract_projective_space(1)
  #   P3 = abstract_projective_space(3)
  #   i = hom(P1, P3, [3P1.O1])
  #   Bl, E = blowup(i)
  #   e = pushforward(E → Bl, E(1))
  #   quad = pullback(Bl → P3, 2P3.O1) - e
  #   @test integral(quad^3) == 0
  #   cubic = pullback(Bl → P3, 3P3.O1) - e
  #   @test integral(quad^2 * cubic) == 1
    
  #   # blowup twisted cubic, with parameters
  #   F, (r, s, t) = function_field(Singular.QQ, ["r", "s", "t"])
  #   P1 = abstract_projective_space(1, base=F)
  #   P3 = abstract_projective_space(3, base=F)
  #   i = hom(P1, P3, [3P1.O1])
  #   Bl, E = blowup(i)
  #   e = pushforward(E → Bl, E(1))
  #   rH, sH, tH = [pullback(Bl → P3, x * P3.O1) - e for x in [r,s,t]]
  #   @test integral(rH * sH * tH) == r*s*t - 3*r - 3*s - 3*t + 10

  #   G = abstract_grassmannian(2, 5)
  #   P9 = abstract_projective_space(9)
  #   i = hom(G, P9, [G.O1])
  #   Bl, E = blowup(i)
  #   e = pushforward(E → Bl, E(1))
  #   quad = pullback(Bl → P9, 2P9.O1) - e
  #   @test simplify(quad^5) == 0
  #   @test simplify(e^5) != 0
    
  #   # blowup space curve of degree d and genus g
  #   F, (r,s,t,d,g) = function_field(Singular.QQ, ["r", "s", "t", "d", "g"])
  #   P3 = abstract_projective_space(3, base=F)
  #   C = abstract_variety(1, base=F)
  #   trim!(C.ring)
  #   C.point = 1//(2-2g) * chern_class(1, C)
  #   i = hom(C, P3, [d * C.point])
  #   Bl, E = blowup(i)
  #   e = pushforward(E → Bl, E(1))
  #   rH, sH, tH = [pullback(Bl → P3, x * P3.O1) - e for x in [r,s,t]]
  #   @test integral(rH * sH * tH) == r*s*t - d*(r+s+t) + (2g-2+4d)
    
  #   G = abstract_grassmannian(2, 5)
  #   Z = zero_locus_section(3line_bundle(G, 1))
  #   Bl, E = blowup(Z → G)
  #   @test dim(Bl) == 6
  #   @test euler(Bl) == 18
  #   @test betti(Bl) == [1,2,4,4,4,2,1]
  #   @test [euler_characteristic(exterior_power(i, cotangent_bundle(Bl))) for i in 0:6] == [1,-2,4,-4,4,-2,1]

  # end
end
