@testset "elliptic surfaces: elliptic parameter" begin
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
