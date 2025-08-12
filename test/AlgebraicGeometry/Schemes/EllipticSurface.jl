
@testset "elliptic surfaces" begin
  @testset "trivial lattice" begin
    k = GF(29)
    # The generic fiber of the elliptic fibration
    # as an elliptic curve over k(t)
    kt, t = polynomial_ring(k, :t)
    kP1 = fraction_field(kt)
    t = gen(kP1)
    E = elliptic_curve(kP1, [3*t^8 + 10*t^7 + 6*t^6 + 17*t^5 + 25*t^4 + 4*t^3 + 23*t^2 + 9*t + 14, 5*t^12 + 25*t^11 + 2*t^10 + 28*t^9 + 28*t^8 + 19*t^7 + 3*t^6 + 17*t^5 + 19*t^4 + 12*t^3 + 25*t^2 + 12*t + 6])
    # A basis for the Mordell-Weil group of E
    mwl_basis = [E(collect(i)) for i in [
    (12*t^4 + 21*t^3 + 5*t^2 + 12*t + 18, 23*t^5 + 7*t^4 + 22*t^3 + 13*t^2),
    (12*t^4 + 20*t^3 + 22*t^2 + 27*t + 18, 15*t^4 + 12*t^3 + 12*t),
    (12*t^4 + 20*t^3 + 27*t^2 + 11*t + 18, -(16*t^4 + 24*t^3 + 3*t^2 + 24*t)),
    (4*t^4 + 5*t^3 + 5*t^2 + 6*t + 13,  9*t^6 + 21*t^5 + 17*t^4 + 12*t^2 + 3*t + 6)]]
    S = elliptic_surface(E, 2, mwl_basis)
    weierstrass_model(S)
    weierstrass_contraction(S)

    E = elliptic_curve(kP1, [0,0,0,1,t^10])
    X = elliptic_surface(E, 2)
    triv = trivial_lattice(X)
    @test det(triv[2])==-3
  end
    
  @testset "trivial lattice QQ" begin
    Qt, t = polynomial_ring(QQ, :t)
    Qtf = fraction_field(Qt)
    E = elliptic_curve(Qtf, [0,0,0,0,t^5*(t-1)^2])
    X3 = elliptic_surface(E, 2)
    weierstrass_contraction(X3)
    trivial_lattice(X3)
  end  
  
  # This test takes about 30 seconds
  @testset "mordel weil lattices" begin
    k = GF(29,2)
    # The generic fiber of the elliptic fibration
    # as an elliptic curve over k(t)
    kt, t = polynomial_ring(k, :t)
    kP1 = fraction_field(kt)
    t = gen(kP1)

    E = elliptic_curve(kP1, [0, 2*t^4 + 28*t^2, 0, 27*t^6 + 19, 10*t^2])
    P = E([0,sqrt(k(10))*t,1])
    X = elliptic_surface(E, 2, [P])
    triv = trivial_lattice(X)
    @test det(triv[2]) == 256
    @test length(Oscar.mordell_weil_torsion(X)) == 1
    alg = algebraic_lattice(X)
    @test det(alg[3]) == -192
    @test det(mordell_weil_sublattice(X)) == 3
    
    X1 = elliptic_surface(short_weierstrass_model(E)[1],2)
    Oscar.isomorphism_from_generic_fibers(X,X1)
  end




  @testset "elliptic parameter" begin
    k = GF(29)
    kt, _ = polynomial_ring(k, :t)
    kP1 = fraction_field(kt); t = gen(kP1)
    E = elliptic_curve(kP1, [(21*t^7+6*t^6+11*t^4),(21*t^10+15*t^9+17*t^7+18*t^6+t^5)])
    mwl_basis = E.([
    [t^3 + 24*t^2 + 22*t + 5, 10*t^5 + 6*t^4 + 6*t^3 + 8*t^2 + 14*t + 3],
    [7*t^3 + 24*t^2 + 9, 9*t^5 + 18*t^4 + 12*t^3 + 8*t^2 + 2],
    [13*t^3 + 9*t + 24, 2*t^5 + 7*t^4 + 6*t^3 + 16*t^2 + 16*t + 22],
    [(17*t^5 + 14*t^4 + 28*t^2 + 2*t + 1)//(t^2 + 11*t + 23), (28*t^8 + 10*t^7 + 27*t^6 + 28*t^5 + 3*t^4 + 27*t^3 + 3*t + 1)//(t^3 + 2*t^2 + 11*t + 25)],
    [(11*t^5 + 22*t^4 + 23*t^3 + 14*t^2 + 17*t + 16)//(t^2 + 24*t + 28), (22*t^8 + 6*t^7 + 21*t^6 + 24*t^5 + 4*t^4 + 23*t^3 + 25*t^2 + 15*t + 6)//(t^3 + 7*t^2 + 26*t + 17)]])

    X = elliptic_surface(E, 2, mwl_basis)
    ff = QQFieldElem[9, 4, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, -1]
    @test det(algebraic_lattice(X)[3])==-1183
    @test length(Oscar.mordell_weil_torsion(X)) == 0 # no torsion points
    # u = elliptic_parameter(X, ff)
    g, phi = two_neighbor_step(X, ff)


    mwl_gens_new = vcat([mwl_basis[1] + mwl_basis[2], mwl_basis[1] - mwl_basis[2]])
    set_mordell_weil_basis!(X, mwl_gens_new)
    @test det(algebraic_lattice(X)[3])==96
    Oscar.algebraic_lattice_primitive_closure!(X)
    @test det(algebraic_lattice(X)[3])==24

  end
end

@testset "normalize_quartic and transform_to_weierstrass" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  P, (u, v) = polynomial_ring(QQ, [:u, :v])
  f = (3*x^2 - 5)^2*(y - 5*x^3 + 30*x - 5)^2 - (7*x^3 + x^2 -5*x + 2)
  g, trans = Oscar._normalize_hyperelliptic_curve(f, parent=P)
  @test trans(f) == g

  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  f = y^2 - 4*x^4 + 5*x^3 - 3*x^2 + 7*x - 4
  f_trans, trans = Oscar.transform_to_weierstrass(f, x, y, QQ.([0, 2]))
  @test trans(f//1) == (16*x^6 + 1//8*x^5 + 14*x^4*y - 16383//4096*x^4 - 16*x^3*y^2 + 647//128*x^3*y)//y^4
end

@testset "normalize_quartic and transform_to_weierstrass over function fields" begin
  pt, t = QQ[:t]
  kt = fraction_field(pt)
  R, (x, y) = polynomial_ring(kt, [:x, :y])
  P, (u, v) = polynomial_ring(kt, [:u, :v])
  f = (3*x^2 - 5*t^2)^2*(y - 5*t*x^3 + 30*x - 5)^2 - (x^3 + t*x^2 -5*x + 2*t^2)
  g, trans = Oscar._normalize_hyperelliptic_curve(f, parent=P)
  @test trans(f) == g

  R, (x, y) = polynomial_ring(kt, [:x, :y])
  f = y^2 - 4*x^4 + 5*x^3 - 3*x^2 + 7*x - 4*t^8
  f_trans, trans = Oscar.transform_to_weierstrass(f, x, y, kt.([0, 2*t^4]))
  @test f_trans == x^3 + (-3//8*t^8 + 49//128)//t^16*x^2 + 7//8//t^8*x*y + (-1//4*t^24 + 9//256*t^16 - 147//2048*t^8 + 2401//65536)//t^32*x - y^2 + (5//16*t^16 - 21//128*t^8 + 343//2048)//t^24*y
  
  R, (x, y) = polynomial_ring(kt, [:x, :y])
  f = y^2 - 4*x^4 + 5*x^3 - 3*x^2 + 7*x - 4*t^8
  Y,trafo = elliptic_surface(f, kt.([0, 2*t^4]))
end


@testset "two neighbor steps" begin
  K = GF(7)
  Kt, t = polynomial_ring(K, :t)
  Ktf = fraction_field(Kt)
  E = elliptic_curve(Ktf, [0, -t^3, 0, t^3, 0])
  P = E([t^3, t^3])
  X2 = elliptic_surface(E, 2, [P]);
  KX2 = function_field(X2; check=false)
  U = weierstrass_chart(X2)
  (xx, yy, tt) = ambient_coordinates(U)
  u4 = KX2(yy, xx*tt)
  u5 = KX2(yy + tt^3, xx*tt - tt^4)
  u6 = KX2(yy + tt^3, xx - tt^3)
  fibers_in_X2 = [QQ.(vec(collect(i))) for i in [
                                                 [4   2   0   0   0   0   0   0   0   -4   -4   -8   -7   -6   -5   -4   -3   -2   -1   0],
                                                 [1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0],
                                                 [5   2   -2   -3   -4   -3   -2   -1   -2   -5   -4   -8   -7   -6   -5   -4   -3   -2   -1   0],
                                                 [4   2   -2   -4   -6   -9//2   -3   -3//2   -7//2   -3   -5//2   -5   -9//2   -4   -7//2   -3   -5//2   -2   -3//2   0],
                                                 [2   1   -1   -2   -3   -2   -1   0   -2   -1   -1   -2   -2   -2   -2   -2   -2   -1   0   1],
                                                 [2   1   0   0   0   0   0   0   0   -2   -2   -4   -4   -4   -4   -3   -2   -1   0   1]
                                                ]]
  NS = algebraic_lattice(X2)[3]
  if !(fibers_in_X2[4] in NS)
    # account for the order of the basis not being unique 
    f = fibers_in_X2[4]
    f[10:11] = reverse(f[10:11])
    fibers_in_X2[4] = f 
  end
  g4,_ = two_neighbor_step(X2, fibers_in_X2[4])  # this should be the 2-torsion case
  g5,_ = two_neighbor_step(X2, fibers_in_X2[5])  # the non-torsion case
  g6,_ = two_neighbor_step(X2, fibers_in_X2[6])  # the non-torsion case
end


#=
# The following tests take roughly 10 minutes which is too much for the CI testsuite.
# We keep them for long running tests to be checked every now and then.
@testset "translation on elliptic surfaces" begin
  P, t = GF(113)[:t]
  kt = fraction_field(P)

  R, (x, y) = kt[:x, :y]

  f = y^2 - (x^3 + (37*t^8 + 56*t^7 + 102*t^6 + 97*t^5 + 21*t^4 + 14*t^3 + 24*t^2 + 23*t + 59)*x + 18*t^12 + 111*t^11 + 17*t^10 + 25*t^9 + 108*t^8 + 91*t^7 + 90*t^6 + 21*t^5 + 68*t^4 + 8*t^3 + 92*t^2 + 66*t + 31)

  E = elliptic_curve(f, x, y)

  P = E([(10*t^6 + 23*t^5 + 94*t^4 + 32*t^3 + t^2 + 40*t + 52)//(t^2 + 77*t + 98), (22*t^9 + 40*t^8 + 63*t^7 + 42*t^6 + 34*t^5 + 48*t^4 + 72*t^3 + 92*t^2 + 85*t + 91)//(t^3 + 59*t^2 + 68*t + 44)])
  
  pts = E.([
  [111*t^4 + 94*t^3 + 57*t^2 + 55*t + 31, 106*t^6 + 2*t^5 + 61*t^4 + 92*t^3 + 105*t^2 + 42*t + 24],   
  [26*t^3 + 32*t^2 + 41*t + 95, 73*t^6 + 101*t^5 + 83*t^4 + 45*t^3 + 53*t^2 + 97], 
  [23*t^4 + 35*t^3 + 27*t^2 + 106*t + 40, 43*t^6 + 103*t^5 + 15*t^4 + 44*t^3 + 34*t^2 + 30*t + 88],
  [109*t^4 + 49*t^3 + 82*t^2 + 69*t + 80, 22*t^6 + 51*t^5 + 27*t^4 + 5*t^3 + 27*t^2 + 80*t + 62],
  [103*t^4 + 27*t^3 + 44*t^2 + 111*t + 1, 2*t^6 + 72*t^5 + 37*t^4 + 59*t^3 + 55*t^2 + 106*t + 59],
  [6*t^4 + 48*t^3 + 70*t^2 + 112*t + 61, 111*t^6 + 41*t^5 + 76*t^4 + 54*t^3 + 58*t^2 + 7*t + 54],
  [42*t^4 + 7*t^3 + 106*t^2 + 112*t + 69, 96*t^6 + 100*t^5 + 48*t^4 + 34*t^3 + 112*t^2 + 83*t + 74],
  [59*t^4 + 55*t^3 + 50*t^2 + 36*t + 8, 15*t^6 + 23*t^5 + 92*t^4 + 64*t^3 + 103*t^2 + 17*t + 87]
  ])

  X = elliptic_surface(E, 2, pts)

  #set_verbose_level(:EllipticSurface, 5)

  # The following should not take more than at most two minutes.
  # But it broke the tests at some point leading to timeout, 
  # so we put it here to indicate regression.
  #for (i, g) in enumerate(values(gluings(default_covering(X))))
  #  gluing_domains(g) # Trigger the computation once
  #end

  D_P = section(X, P)
 
  II = first(components(D_P));
 
  trans = Oscar.translation_morphism(X, P; divisor=D_P)
 
  JJ = Oscar._pushforward_section(trans, P; divisor=D_P)
 
  IIX = first(components(section(X, 2*P)));
  # We have little chance to get through with the computations on all charts. 
  # But it suffices to compare the result on the weierstrass charts.
  weier = weierstrass_chart_on_minimal_model(X)
  @test IIX(weier) == JJ(weier)

  
  # pushforward of the whole algebraic lattice with verification of the result.
  lat, _, A = algebraic_lattice(X)
  A = gram_matrix(ambient_space(A))
 
  ll = Oscar._pushforward_lattice_along_isomorphism(trans)
  res_mat = [Oscar.basis_representation(X, d) for d in ll]
  
  res_mat = matrix(QQ, res_mat)
  # Check that the base change matrix is indeed orthogonal for the given lattice
  @test res_mat*A*transpose(res_mat) == A
end
=#

@testset "moebius transformations" begin
  k = GF(113)
  kt,t = polynomial_ring(k,:t)
  Ft = fraction_field(kt)

  E = elliptic_curve(Ft, Ft.([0,0,0,(56*t^8 + 56),0]));
  basis_mwl = [[112*t^4 + 112*t^3 + 56*t^2 + 15*t + 1, 100*t^6 + 24*t^5 + 56*t^4 + 13*t^3 + 64*t^2 + 89*t + 31],
               [31*t^4 + 15*t^2 + 82, 44*t^5 + 14*t^3 + 69*t],
               [82*t^4 + 13, 37*t^4 + 10],
               [91*t^4 + 16*t^3 + 25*t^2 + 14*t + 22, 18*t^6 + 55*t^5 + 45*t^4 + 44*t^3 + 110*t^2 + 58*t + 69],
               [32*t^4 + 72*t^2 + 81, 78*t^6 + 79*t^4 + 58*t^2 + 40],
               [56*t^4 + 49*t^3 + 85*t^2 + 57*t + 57, 72*t^6 + 25*t^5 + 22*t^4 + 101*t^3 + 104*t^2 + 88*t + 50],
               [22*t^4 + 87*t^3 + 77*t^2 + 26*t + 22, 44*t^6 + 27*t^5 + 68*t^4 + 98*t^3 + 45*t^2 + 27*t + 69],
               [112*t^4 + 44*t^3 + 49*t^2 + 44*t + 112, 100*t^6 + 74*t^5 + 49*t^4 + 8*t^3 + 49*t^2 + 74*t + 100]]
  basis_mwl = [E(i) for i in basis_mwl];  

  X = elliptic_surface(E, 2, basis_mwl[1:0])  # speed up test by hiding the sections
  
  moeb = Oscar.admissible_moebius_transformations(X, X)

  @test length(moeb) == 16 # This must really be the number for this particular surface.

  for phi in moeb
    @test Oscar.check_isomorphism_on_generic_fibers(phi)
  end

  phi = moeb[5] # Just pick one which is not the identity

  W = weierstrass_chart_on_minimal_model(X)

  for (i, U) in enumerate(affine_charts(X))
    phi_loc = Oscar.cheap_realization(phi, W, U)
    @test all(iszero(pullback(phi_loc)(lifted_numerator(g))) for g in gens(modulus(OO(U))))
  end
  
  phi_star = Oscar.pushforward_on_algebraic_lattices(phi)
  # Check that the pushforward preserves the numerical lattice of X
  Triv = algebraic_lattice(X)[3]
  @test phi_star(Triv)==Triv
  
  # test translation
  P = Oscar.mordell_weil_torsion(X)[1]
  torsion_translation = Oscar.translation_morphism(X, P)
  tors_of_P = pushforward(torsion_translation,zero_section(X)) 
  Psect = section(X, P)
  @test tors_of_P == Psect
  
  # add a section so that there is something to compute (although rather trivial)
  set_mordell_weil_basis!(X, basis_mwl[1:1])
  P = Oscar.extract_mordell_weil_basis(torsion_translation)
  @test length(P) == 2
end

