@testset "elliptic surfaces" begin
  #=
  # The tests in this file are quite expensive
  # hence this one is disabled
  k = GF(29)
  # The generic fiber of the elliptic fibration
  # as an elliptic curve over k(t)
  kt, t = polynomial_ring(k, :t)
  kP1 = fraction_field(kt)
  t = gen(kP1)
  E = EllipticCurve(kP1, [3*t^8 + 10*t^7 + 6*t^6 + 17*t^5 + 25*t^4 + 4*t^3 + 23*t^2 + 9*t + 14, 5*t^12 + 25*t^11 + 2*t^10 + 28*t^9 + 28*t^8 + 19*t^7 + 3*t^6 + 17*t^5 + 19*t^4 + 12*t^3 + 25*t^2 + 12*t + 6])
  # A basis for the Mordell-Weil group of E
  mwl_basis = [E(collect(i)) for i in [
  (12*t^4 + 21*t^3 + 5*t^2 + 12*t + 18, 23*t^5 + 7*t^4 + 22*t^3 + 13*t^2),
  (12*t^4 + 20*t^3 + 22*t^2 + 27*t + 18, 15*t^4 + 12*t^3 + 12*t),
  (12*t^4 + 20*t^3 + 27*t^2 + 11*t + 18, -(16*t^4 + 24*t^3 + 3*t^2 + 24*t)),
  (4*t^4 + 5*t^3 + 5*t^2 + 6*t + 13,  9*t^6 + 21*t^5 + 17*t^4 + 12*t^2 + 3*t + 6)]]
  S = elliptic_surface(E, 2, mwl_basis)
  weierstrass_model(S)
  relatively_minimal_model(S)

  E = EllipticCurve(kP1, [0,0,0,1,t^10])
  X = elliptic_surface(E, 2)
  triv = trivial_lattice(X)
  =#

  #=
  # This test takes about 5 minutes
  @testset "mordel weil lattices" begin
    k = GF(29,2)
    # The generic fiber of the elliptic fibration
    # as an elliptic curve over k(t)
    kt, t = polynomial_ring(k, :t)
    kP1 = fraction_field(kt)
    t = gen(kP1)

    E = EllipticCurve(kP1, [0, 2*t^4 + 28*t^2, 0, 27*t^6 + 19, 10*t^2])
    P = E([0,sqrt(k(10))*t,1])
    X = elliptic_surface(E, 2, [P])
    triv = trivial_lattice(X)
    @test det(triv[2]) == 256
    @test length(Oscar.mordell_weil_torsion(X)) == 1
    alg = algebraic_lattice(X)
    @test det(alg[3]) == -192
    @test det(mordell_weil_lattice(X)) == 3
  end
  =#
  #=
  # this test is quite expensive
  # probably because it is over QQ
  Qt, t = polynomial_ring(QQ, :t)
  Qtf = fraction_field(Qt)
  E = EllipticCurve(Qtf, [0,0,0,0,t^5*(t-1)^2])
  X3 = elliptic_surface(E, 2)
  relatively_minimal_model(X3)
  trivial_lattice(X3)
  =#

  @testset "elliptic parameter" begin
    k = GF(29)
    kt, _ = polynomial_ring(k, :t)
    kP1 = fraction_field(kt); t = gen(kP1)
    E = EllipticCurve(kP1, [(21*t^7+6*t^6+11*t^4),(21*t^10+15*t^9+17*t^7+18*t^6+t^5)])
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
    elliptic_parameter(X, ff)
  end
end

@testset "normalize_quartic and transform_to_weierstrass" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  P, (u, v) = polynomial_ring(QQ, [:u, :v])
  f = (3*x^2 - 5)^2*(y - 5*x^3 + 30*x - 5)^2 - (x^2 -5*x + 2)
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
  f = (3*x^2 - 5*t^2)^2*(y - 5*t*x^3 + 30*x - 5)^2 - (t*x^2 -5*x + 2*t^2)
  g, trans = Oscar._normalize_hyperelliptic_curve(f, parent=P)
  @test trans(f) == g

  R, (x, y) = polynomial_ring(kt, [:x, :y])
  f = y^2 - 4*x^4 + 5*x^3 - 3*x^2 + 7*x - 4*t^8
  f_trans, trans = Oscar.transform_to_weierstrass(f, x, y, kt.([0, 2*t^4]))
  @test f_trans == x^3 + (-3//8*t^8 + 49//128)//t^16*x^2 + 7//8//t^8*x*y + (-1//4*t^24 + 9//256*t^16 - 147//2048*t^8 + 2401//65536)//t^32*x - y^2 + (5//16*t^16 - 21//128*t^8 + 343//2048)//t^24*y
end

