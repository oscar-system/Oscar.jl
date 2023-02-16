@testset "A K3 surface from an elliptic fibration" begin
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

  prX = blow_up(I_sing_X)
  E1 = oscar.exceptional_divisor(prX)
  X1 = domain(prX)
  Y1, inc_Y1, pr_Y1 = strict_transform(prX, inc_S)

  I_sing_Y1 = oscar.ideal_sheaf_of_singular_locus(Y1)
  I_sing_X1 = radical(pushforward(inc_Y1)(I_sing_Y1))
  prX2 = blow_up(I_sing_X1, covering=oscar.simplified_covering(X1),
                var_name="t")
  E2 = exceptional_divisor(prX2)
  X2 = domain(prX2)
  Y2, inc_Y2, pr_Y2 = strict_transform(prX2, inc_Y1)
  simplify!(Y2)
  I_sing_Y2 = oscar.ideal_sheaf_of_singular_locus(Y2)
  I_sing_X2 = radical(pushforward(inc_Y2)(I_sing_Y2))
  prX3 = blow_up(I_sing_X2, covering=oscar.simplified_covering(X2),
                var_name="u")
  E3 = exceptional_divisor(prX3)
  X3 = domain(prX3)
  Y3, inc_Y3, pr_Y3 = strict_transform(prX3, inc_Y2)
  simplify!(Y3)
  I_sing_Y3 = oscar.ideal_sheaf_of_singular_locus(Y3)
  I_sing_X3 = radical(pushforward(inc_Y3)(I_sing_Y3))
  prX4 = blow_up(I_sing_X3, covering=oscar.simplified_covering(X3),
                var_name="v")
  E4 = exceptional_divisor(prX4)
  X4 = domain(prX4)
  Y4, inc_Y4, pr_Y4 = strict_transform(prX4, inc_Y3)
  I_sing_Y4 = oscar.ideal_sheaf_of_singular_locus(Y4)
  I_sing_X4 = radical(pushforward(inc_Y4)(I_sing_Y4))
  U = patches(oscar.simplified_covering(X4))
  
  ref_patches = [x for x in U if !(x===V1) && !(x===V2)]

  id_dict = IdDict{AbsSpec, Ideal}()
  for x in ref_patches
    id_dict[x] = I_sing_X4(x)
  end

  V11 = PrincipalOpenSubset(V1, gens(OO(V1))[2])
  V12 = PrincipalOpenSubset(V1, gens(OO(V1))[2]-1)

  V21 = PrincipalOpenSubset(V2, gens(OO(V2))[1])
  V22 = PrincipalOpenSubset(V2, gens(OO(V2))[1]-1)

  id_dict[V11] = ideal(OO(V11), saturated_ideal(I_sing_X4(V11)))
  id_dict[V12] = ideal(OO(V12), saturated_ideal(I_sing_X4(V12)))
  id_dict[V21] = ideal(OO(V21), saturated_ideal(I_sing_X4(V21)))
  id_dict[V22] = ideal(OO(V22), saturated_ideal(I_sing_X4(V22)))

  ref_patches = vcat(ref_patches, [V11, V12, V21, V22])

  ref = Covering(ref_patches)
  oscar.inherit_glueings!(ref, oscar.simplified_covering(X4))
  push!(coverings(X4), ref)

  J = IdealSheaf(X4, id_dict, check=false)

  prX4 = blow_up(J, covering=ref)
  X5 = domain(prX4)
  E5 = exceptional_divisor(prX4)
  simplify!(X5)
  U = patches(oscar.simplified_covering(X5))

  Y5, inc_Y5, pr_Y5 = strict_transform(prX4, inc_Y4)
  Y5 = domain(inc_Y5)

  P5 = pr_Y5
  P4 = compose(P5, pr_Y4)
  P3 = compose(P4, pr_Y3)
  P2 = compose(P3, pr_Y2)
  P1 = compose(P2, pr_Y1)
  E11 = pullback(inc_Y1)(E1)
  E12 = pullback(pr_Y2)(E11)
  E13 = pullback(pr_Y3)(E12)
  E14 = pullback(pr_Y4)(E13)
  E15 = pullback(pr_Y5)(E14)

  E22 = pullback(inc_Y2)(E2)
  E23 = pullback(pr_Y3)(E22)
  E24 = pullback(pr_Y4)(E23)
  E25 = pullback(pr_Y5)(E24)
  
  E33 = pullback(inc_Y3)(E3)
  E34 = pullback(pr_Y4)(E33)
  E35 = pullback(pr_Y5)(E34)
  
  E44 = pullback(inc_Y4)(E4)
  E45 = pullback(pr_Y5)(E44)

  E55 = pullback(inc_Y5)(E5)
end
