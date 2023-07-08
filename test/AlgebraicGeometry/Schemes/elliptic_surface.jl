@testset "elliptic surfaces" begin
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

  Qt, t = polynomial_ring(QQ, :t)
  Qtf = fraction_field(Qt)
  E = EllipticCurve(Qtf, [0,0,0,0,t^5*(t-1)^2])
  X3 = elliptic_surface(E, 2)
  relatively_minimal_model(X3)

  @test length(fiber_components(S, [-1,1]))
end
