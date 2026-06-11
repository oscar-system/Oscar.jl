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
