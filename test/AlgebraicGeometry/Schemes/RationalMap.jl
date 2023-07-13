@testset "Rational maps" begin
  IP3 = projective_space(QQ, [:x, :y, :z, :w])
  S = homogeneous_coordinate_ring(IP3)
  (x, y, z, w) = gens(S)
  I = ideal(S, x*w - y*z)
  X_proj = subscheme(IP3, I)
  X = covered_scheme(X_proj)

  P1 = projective_space(QQ, 1)
  IP1 = covered_scheme(P1)
  U = X[1][4]
  V = IP1[1][2]
  (x, y, z) = gens(ambient_coordinate_ring(U))
  Phi = oscar.RationalMap(X, IP1, U, V, [x//y])

  @test domain(Phi) === X 
  @test codomain(Phi) === IP1

  oscar.realize_on_patch(Phi, U)
  oscar.realize_on_patch(Phi, X[1][1])
  oscar.realize_on_patch(Phi, X[1][2])
  oscar.realize_on_patch(Phi, X[1][3])
  oscar.realize(Phi)
end

@testset "The standard Cremona transformation" begin
  IP2_proj = projective_space(QQ, [:x, :y, :z])
  IP2 = covered_scheme(IP2_proj)
  S = homogeneous_coordinate_ring(IP2_proj)
  (x, y, z) = gens(S)
  I1 = ideal(S, [x, y])
  I2 = ideal(S, [x, z])
  I3 = ideal(S, [y, z])
  I = I1*I2*I3
  II = simplify!(IdealSheaf(IP2_proj, I))
  pr = blow_up(II)
  X = domain(projection(pr))
  set_attribute!(X, :is_irreducible, true)
  V = first(affine_charts(IP2))
  U = first(affine_charts(X))
  pr_cov = covering_morphism(projection(pr))
  U_simp = first(patches(oscar.simplified_covering(X)))
  pb_y, pb_z = pullback(pr_cov[U]).(gens(OO(V)))
  Phi = oscar.RationalMap(X, X, U, U, [inv(fraction(pb_y)), fraction(pb_z)//fraction(pb_y), fraction(pb_z)])
  oscar.realize_on_patch(Phi, U)
end

