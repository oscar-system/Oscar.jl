@testset "resolution of singularities of elliptic fibrations" begin
  IP1 = projective_space(GF(29), ["s", "t"])

  O0 = twisting_sheaf(IP1, 0)
  O4 = twisting_sheaf(IP1, -4)
  O6 = twisting_sheaf(IP1, -6)

  E = direct_sum([O0, O4, O6])

  X_proj = projectivization(E, var_names=[:z, :x, :y])

  X = covered_scheme(X_proj)

  U = affine_charts(X)[1]
  (x, y, t) = gens(OO(U))
  ft = y^2 - (x^3 + 21*x + (28*t^7+18))
  I = IdealSheaf(X, U, [ft])

  # The fabulous K3-surface. Almost.
  inc_S = Oscar.CoveredClosedEmbedding(X, I)
  @test I === image_ideal(inc_S)
  S = domain(inc_S)

  I_sing = Oscar.ideal_sheaf_of_singular_locus(S)
  I_sing_X = small_generating_set(pushforward(inc_S)(I_sing))
             # happens to be radical in this example -- no radical needed here
  @test scheme(I_sing) === S
  @test scheme(I_sing_X) === X

  prX = blow_up(I_sing_X)
  E1 = Oscar.exceptional_divisor(prX)
  X1 = domain(prX)
  Y1, inc_Y1, pr_Y1 = strict_transform(prX, inc_S)


  I_sing_Y1 = Oscar.ideal_sheaf_of_singular_locus(Y1)
  I_sing_X1 = small_generating_set(pushforward(inc_Y1)(I_sing_Y1))
             # happens to be radical in this example -- no radical needed here
  prX2 = blow_up(I_sing_X1, covering=Oscar.simplified_covering(X1),
                var_name="t")
  E2 = exceptional_divisor(prX2)
  X2 = domain(prX2)
  Y2, inc_Y2, pr_Y2 = strict_transform(prX2, inc_Y1)
  simplify!(Y2)
  # @show has_attribute(Y2, :simplified_covering)
  I_sing_Y2 = Oscar.ideal_sheaf_of_singular_locus(Y2)
  I_sing_X2 = small_generating_set(pushforward(inc_Y2)(I_sing_Y2))
             # happens to be radical in this example -- no radical needed here
  prX3 = blow_up(I_sing_X2, covering=Oscar.simplified_covering(X2),
                var_name="u")
  E3 = exceptional_divisor(prX3)
  X3 = domain(prX3)
  Y3, inc_Y3, pr_Y3 = strict_transform(prX3, inc_Y2)
  simplify!(Y3)
  # @show has_attribute(Y3, :simplified_covering)
  I_sing_Y3 = Oscar.ideal_sheaf_of_singular_locus(Y3)
  I_sing_X3 = pushforward(inc_Y3)(I_sing_Y3)

  # Now the singular locus consists of two points. 
  # We blow them up successively rather than at once. 
  #println("starting primary decomposition of singular locus")
  decomp = Oscar.maximal_associated_points(I_sing_X3)
  #println("finished primary decomposition of singular locus")
  l = simplify.([a for a in decomp])

  #@show "blowing up first center"
  #@show gens.(l[1].(patches(Oscar.simplified_covering(X3))))
  prX41 = blow_up(l[1], covering=Oscar.simplified_covering(X3),
                var_name="v")
  #@show "done blowing up"
  X41 = domain(prX41)
  E41 = exceptional_divisor(prX41)
  Y41, inc_Y41, pr_Y41 = strict_transform(prX41, inc_Y3)

  l2 = small_generating_set(strict_transform(prX41, l[2]))
          ## was radical before transformation, stays radical afterwards
  simplify!(X41)
  #@show scheme(l2) === X41
  #@show "blowing up second center"
  #@show gens.(l2.(patches(Oscar.simplified_covering(X41))))
  prX42 = blow_up(l2, var_name="vv")
 # prX42 = blow_up(l2, covering=Oscar.simplified_covering(X41),
 #               var_name="vv")
  E42 = exceptional_divisor(prX42)
  X42 = domain(prX42)
  Y42, inc_Y42, pr_Y42 = strict_transform(prX42, inc_Y41)
  I_sing_Y42 = Oscar.ideal_sheaf_of_singular_locus(Y42)
  I_sing_X42 = pushforward(inc_Y42)(I_sing_Y42)

  # Now the singular locus consists of three points. 
  # We blow them up successively rather than at once. 
  #println("starting primary decomposition of singular locus")
  centers = Oscar.maximal_associated_points(I_sing_X42)
  #println("finished primary decomposition of singular locus")

  prX51 = blow_up(small_generating_set(first(centers)))

  E51 = exceptional_divisor(prX51)

  #@show E51
  Y51, inc_Y51, prY51 = strict_transform(prX51, inc_Y42)
  centers = (x->strict_transform(prX51, x)).(centers[2:3])

  prX52 = blow_up(small_generating_set(first(centers)))
  E52 = exceptional_divisor(prX52)

  #@show E52
  Y52, inc_Y52, prY52 = strict_transform(prX52, inc_Y51)
  centers = (x->strict_transform(prX52, x)).(centers[2:2])

  prX53 = blow_up(small_generating_set(first(centers)))
  E53 = exceptional_divisor(prX53)

  #@show E53
  Y53, inc_Y53, prY53 = strict_transform(prX53, inc_Y52)

# error()
# U = patches(Oscar.simplified_covering(X4))
# println("Some context on how the architecture of the refinement:")
# @show gens.(I_sing_X4.(U))
#
# # We need to refine the covering so that we can separate the 
# # different points in the support of I_sing_X4.
#
# V1 = U[9]
# V2 = U[11]
# @show gens(I_sing_X4(V1))
# @show gens(I_sing_X4(V2))
# @show gens(OO(V1))
# @show gens(OO(V2))
# 
# ref_patches = [x for x in U if !(x===V1) && !(x===V2)]
#
# id_dict = IdDict{AbsAffineScheme, Ideal}()
# for x in ref_patches
#   id_dict[x] = I_sing_X4(x)
# end
#
# V11 = PrincipalOpenSubset(V1, gens(OO(V1))[2])
# V12 = PrincipalOpenSubset(V1, gens(OO(V1))[2]-1)
#
# V21 = PrincipalOpenSubset(V2, gens(OO(V2))[1])
# V22 = PrincipalOpenSubset(V2, gens(OO(V2))[1]-1)
#
# id_dict[V11] = ideal(OO(V11), saturated_ideal(I_sing_X4(V11)))
# id_dict[V12] = ideal(OO(V12), saturated_ideal(I_sing_X4(V12)))
# id_dict[V21] = ideal(OO(V21), saturated_ideal(I_sing_X4(V21)))
# id_dict[V22] = ideal(OO(V22), saturated_ideal(I_sing_X4(V22)))
#
# ref_patches = vcat(ref_patches, [V11, V12, V21, V22])
#
# ref = Covering(ref_patches)
# Oscar.inherit_gluings!(ref, Oscar.simplified_covering(X4))
# push!(coverings(X4), ref)
#
# J = IdealSheaf(X4, id_dict, check=false)
#
# prX5 = blow_up(J, covering=ref)
# X5 = domain(prX5)
# E5 = exceptional_divisor(prX4)
# simplify!(X5)
# U = patches(Oscar.simplified_covering(X5))

  Y53 = domain(inc_Y53)

  # @show is_smooth(Y53)

  # Pull all exceptional divisors up to the smooth model:
  E12 = strict_transform(prX2, E1)
  E13 = strict_transform(prX3, E12)
  E141 = strict_transform(prX41, E13)
  E142 = strict_transform(prX42, E141)
  E151 = strict_transform(prX51, E142)
  E152 = strict_transform(prX52, E151)
  E153 = strict_transform(prX53, E152)

  E23 = strict_transform(prX3, E2)
  E241 = strict_transform(prX41, E23)
  E242 = strict_transform(prX42, E241)
  E251 = strict_transform(prX51, E242)
  E252 = strict_transform(prX52, E251)
  E253 = strict_transform(prX53, E252)

  E341 = strict_transform(prX41, E3)
  E342 = strict_transform(prX42, E341)
  E351 = strict_transform(prX51, E342)
  E352 = strict_transform(prX52, E351)
  E353 = strict_transform(prX53, E352)

  E4142 = strict_transform(prX42, E41)
  E4151 = strict_transform(prX51, E4142)
  E4152 = strict_transform(prX52, E4151)
  E4153 = strict_transform(prX53, E4152)

  E4251 = strict_transform(prX51, E42)
  E4252 = strict_transform(prX52, E4251)
  E4253 = strict_transform(prX53, E4252)

  E5152 = strict_transform(prX52, E51)
  E5153 = strict_transform(prX53, E5152)

  E5253 = strict_transform(prX53, E52)

  E5353 = E53

  # Restrict them: 
  E1res = pullback(inc_Y53)(E153)
  E2res = pullback(inc_Y53)(E253)
  E3res = pullback(inc_Y53)(E353)
  E41res = pullback(inc_Y53)(E4153)
  E42res = pullback(inc_Y53)(E4253)
  E51res = pullback(inc_Y53)(E5153)
  E52res = pullback(inc_Y53)(E5253)
  E53res = pullback(inc_Y53)(E5353)

  # Compute the intersection matrix:
  # A = zero_matrix(ZZ, 8, 8)
  # E = [E1res, E2res, E3res, E41res, E42res, E51res, E52res, E53res]
  # F = weil_divisor.(E)   # too slow to be a test
  # for i in 1:8
  #   for j in 1:8
  #    if j == i 
  #      # avoid self intersection for now
  #      A[i, i] = -ZZ(2)
  #      continue
  #    end
  #    A[i, j] = integral(intersect(F[i], E[j]))
  #  end
  #end
  #@test A == ZZ[-2 0 0 0 2 0 0 0; 0 -2 0 0 0 0 0 2; 0 0 -2 0 2 2 0 0; 0 0 0 -2 0 2 2 2; 2 0 2 0 -2 0 0 0; 0 0 2 2 0 -2 0 0; 0 0 0 2 0 0 -2 0; 0 2 0 2 0 0 0 -2]
end
