#@testset "A K3 surface from an elliptic fibration" begin
  IP1 = projective_space(GF(29), ["s", "t"])

  O0 = twisting_sheaf(IP1, 0)
  O4 = twisting_sheaf(IP1, -4)
  O6 = twisting_sheaf(IP1, -6)

  E = direct_sum([O0, O4, O6])

  X_proj = projectivization(E, var_names=["z", "x", "y"])

  X = covered_scheme(X_proj)

  U = affine_charts(X)[1]
  (x, y, t) = gens(OO(U))
  ft = y^2 - (x^3 + 21*x + (28*t^7+18))
  I = IdealSheaf(X, U, [ft])

  # The fabulous K3-surface. Almost.
  S = subscheme(I)

  I_sing = oscar.ideal_sheaf_of_singular_locus(S)

  @test scheme(I_sing) === S

  Y1_proj = blow_up(I_sing)
  Y1 = covered_scheme(Y1_proj)
  simplify!(Y1)
  @test !is_smooth(Y1)
  @show "done"

  @show has_attribute(Y1, :simplified_covering)
  I_sing_Y1 = oscar.ideal_sheaf_of_singular_locus(Y1)
  @show "done"
  Y2_proj = blow_up(I_sing_Y1)
  @show "done"
  Y2 = covered_scheme(Y2_proj)
  @show "done"
  simplify!(Y2)
  @show "done"
  @test !is_smooth(Y2)
  
  @show "round 2"
  I_sing_Y2 = oscar.ideal_sheaf_of_singular_locus(Y2)
  @show "done"
  Y3_proj = blow_up(I_sing_Y2)
  @show "done"
  Y3 = covered_scheme(Y3_proj)
  @show "done"
  simplify!(Y3)
  @show "done"
  @test !is_smooth(Y3)
  
  @show "round 3"
  I_sing_Y3 = oscar.ideal_sheaf_of_singular_locus(Y3)
  @show "done"
  Y4_proj = blow_up(I_sing_Y3)
  @show "done"
  Y4 = covered_scheme(Y4_proj)
  @show "done"
  simplify!(Y4)
  @show "done"
  @test !is_smooth(Y4)

  @show "round 4"
  I_sing_Y4 = oscar.ideal_sheaf_of_singular_locus(Y4)
  @show "done"
  Y5_proj = blow_up(I_sing_Y4)
  @show "done"
  Y5 = covered_scheme(Y5_proj)
  @show "done"
  simplify!(Y5)
  @show "done"
  @test !is_smooth(Y5)

  @show "round 5"
  I_sing_Y5 = oscar.ideal_sheaf_of_singular_locus(Y5)
  @show "done"
  Y6_proj = blow_up(I_sing_Y5)
  @show "done"
  Y6 = covered_scheme(Y5_proj)
  @show "done"
  simplify!(Y6)
  @show "done"
  @test !is_smooth(Y6)

  @show "round 6"
  I_sing_Y6 = oscar.ideal_sheaf_of_singular_locus(Y6)
  @show "done"
  Y7_proj = blow_up(I_sing_Y6)
  @show "done"
  Y7 = covered_scheme(Y6_proj)
  @show "done"
  simplify!(Y7)
  @show "done"
  @test !is_smooth(Y7)
#end
