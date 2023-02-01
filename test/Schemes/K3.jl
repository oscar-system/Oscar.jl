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
  inc_S = oscar.CoveredClosedEmbedding(X, I)
  @test I === image_ideal(inc_S)
  S = domain(inc_S)

  I_sing = oscar.ideal_sheaf_of_singular_locus(S)
  I_sing_X = radical(pushforward(inc_S)(I_sing))

  @test scheme(I_sing) === S
  @test scheme(I_sing_X) === X

  Bl_X = blow_up(I_sing_X)
  E1 = oscar.exceptional_divisor(Bl_X)
  X1 = domain(Bl_X)
  inc_Y1 = strict_transform(Bl_X, inc_S)
  Y1 = domain(inc_Y1)
  #simplify!(Y1)
  #@test !is_smooth(Y1)
  #@show "done"

  #@show has_attribute(Y1, :simplified_covering)
  I_sing_Y1 = oscar.ideal_sheaf_of_singular_locus(Y1)
  I_sing_X1 = radical(pushforward(inc_Y1)(I_sing_Y1))
  @show "done"
  prX2 = blow_up(I_sing_X1, covering=oscar.simplified_covering(X1),
                var_name="t")
  @show "done2"
  E2 = exceptional_divisor(prX2)
  @show "done3"
  X2 = domain(prX2)
  @show "done4"
  inc_Y2 = strict_transform(prX2, inc_Y1)
  @show "done5"
  Y2 = domain(inc_Y2)
  @show "done6"
  simplify!(Y2)
  @show "done7"
  @show has_attribute(Y2, :simplified_covering)
  I_sing_Y2 = oscar.ideal_sheaf_of_singular_locus(Y2)
  I_sing_X2 = radical(pushforward(inc_Y2)(I_sing_Y2))
  @show "done"
  prX3 = blow_up(I_sing_X2, covering=oscar.simplified_covering(X2),
                var_name="u")
  @show "done3"
  E3 = exceptional_divisor(prX3)
  @show "done3"
  X3 = domain(prX3)
  @show "done4"
  inc_Y3 = strict_transform(prX3, inc_Y2)
  @show "done5"
  Y3 = domain(inc_Y3)
  @show "done6"
  simplify!(Y3)
  @show "done7"
  @show has_attribute(Y3, :simplified_covering)
  I_sing_Y3 = oscar.ideal_sheaf_of_singular_locus(Y3)
  I_sing_X3 = radical(pushforward(inc_Y3)(I_sing_Y3))
  @show "done"
  prX4 = blow_up(I_sing_X3, covering=oscar.simplified_covering(X3),
                var_name="v")
  @show "done4"
  E4 = exceptional_divisor(prX4)
  @show "done3"
  X4 = domain(prX4)
  @show "done4"
  inc_Y4 = strict_transform(prX4, inc_Y3)
  @show "done5"
  Y4 = domain(inc_Y4)
  @show "done6"
  I_sing_Y4 = oscar.ideal_sheaf_of_singular_locus(Y4)
  I_sing_X4 = radical(pushforward(inc_Y4)(I_sing_Y4))
  U = patches(oscar.simplified_covering(X4))
  @show gens.(I_sing_X4.(U))
#  @show "done"
#  prX5 = blow_up(I_sing_X4, covering=oscar.simplified_covering(X4),
#                var_name="w")
#  @show "done5"
#  E5 = exceptional_divisor(prX5)
#  @show "done3"
#  X5 = domain(prX5)
#  @show "done4"
#  inc_Y5 = strict_transform(prX5, inc_Y4)
#  @show "done5"
#  Y5 = domain(inc_Y5)
#  @show "done6"
#  simplify!(Y5)
#  @show "done7"
#  @show has_attribute(Y5, :simplified_covering)
#  I_sing_Y5 = oscar.ideal_sheaf_of_singular_locus(Y5)
#  @show "success"
#  I_sing_X5 = radical(pushforward(inc_Y5)(I_sing_Y5))
#  @show "done"
#  prX6 = blow_up(I_sing_X5, covering=oscar.simplified_covering(X5))
#  @show "done6"
#  E6 = exceptional_divisor(prX6)
#  @show "done3"
#  X6 = domain(prX6)
#  @show "done4"
#  inc_Y6 = strict_transform(prX6, inc_Y5)
#  @show "done5"
#  Y6 = domain(inc_Y6)
#  @show "done6"
#  simplify!(Y6)
#  @show "done7"
#  @show has_attribute(Y6, :simplified_covering)
#  I_sing_Y6 = oscar.ideal_sheaf_of_singular_locus(Y6)
#  I_sing_X6 = radical(pushforward(inc_Y6)(I_sing_Y6))
#  @show "done"
#  prX7 = blow_up(I_sing_X6, covering=oscar.simplified_covering(X6))
#  @show "done7"
#  E7 = exceptional_divisor(prX7)
#  @show "done3"
#  X7 = domain(prX7)
#  @show "done4"
#  inc_Y7 = strict_transform(prX7, inc_Y6)
#  @show "done5"
#  Y7 = domain(inc_Y7)
#  @show "done6"
#  simplify!(Y7)
#  @show "done7"
#  @test !is_smooth(Y2)
  
#  @show "round 2"
#  I_sing_Y2 = oscar.ideal_sheaf_of_singular_locus(Y2)
#  @show "done"
#  Y3_proj = blow_up(I_sing_Y2)
#  @show "done"
#  Y3 = covered_scheme(Y3_proj)
#  @show "done"
#  simplify!(Y3)
#  @show "done"
#  @test !is_smooth(Y3)
#  
#  @show "round 3"
#  I_sing_Y3 = oscar.ideal_sheaf_of_singular_locus(Y3)
#  @show "done"
#  Y4_proj = blow_up(I_sing_Y3)
#  @show "done"
#  Y4 = covered_scheme(Y4_proj)
#  @show "done"
#  simplify!(Y4)
#  @show "done"
#  @test !is_smooth(Y4)
# 
#  @show "round 4"
#  I_sing_Y4 = oscar.ideal_sheaf_of_singular_locus(Y4)
#  @show "done"
#  Y5_proj = blow_up(I_sing_Y4)
#  @show "done"
#  Y5 = covered_scheme(Y5_proj)
#  @show "done"
#  simplify!(Y5)
#  @show "done"
#  @test !is_smooth(Y5)
# 
#  @show "round 5"
#  I_sing_Y5 = oscar.ideal_sheaf_of_singular_locus(Y5)
#  @show "done"
#  Y6_proj = blow_up(I_sing_Y5)
#  @show "done"
#  Y6 = covered_scheme(Y5_proj)
#  @show "done"
#  simplify!(Y6)
#  @show "done"
#  @test !is_smooth(Y6)
# 
#  @show "round 6"
#  I_sing_Y6 = oscar.ideal_sheaf_of_singular_locus(Y6)
#  @show "done"
#  Y7_proj = blow_up(I_sing_Y6)
#  @show "done"
#  Y7 = covered_scheme(Y6_proj)
#  @show "done"
#  simplify!(Y7)
#  @show "done"
#  @test !is_smooth(Y7)
#end
