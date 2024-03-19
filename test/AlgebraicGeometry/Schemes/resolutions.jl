@testset "non-snc-locus" begin
  R, (x, y, z) = QQ[:x, :y, :z]
  f = x^4 + y^2 + z^2
  I = ideal(R, f)
  IA3 = AffineScheme(R)
  inc_X = ClosedEmbedding(IA3, I)
  X = domain(inc_X)
  origin = ideal(R, [x, y, z])
  bl1 = blow_up(IA3, origin)
  X1, inc1, pr1 = strict_transform(bl1, inc_X)
  E1 = pullback(inc1, exceptional_divisor(bl1))
  I_sing = Oscar.ideal_sheaf_of_singular_locus(X1)
  new_centers = Oscar.maximal_associated_points(pushforward(inc1, I_sing))
  bl2 = blow_up(first(new_centers))
  X2, inc2, pr2 = strict_transform(bl2, inc1)

  E1 = pullback(inc2, strict_transform(bl2, exceptional_divisor(bl1)))
  E2 = pullback(inc2, exceptional_divisor(bl2))

  Oscar.non_snc_locus([E1, E2]) == Oscar.unit_ideal_sheaf(X2)


  E1 = pullback(inc2, pullback(bl2, exceptional_divisor(bl1)))
  E2 = pullback(inc2, exceptional_divisor(bl2))

  Oscar.non_snc_locus([E1, E2]) != Oscar.unit_ideal_sheaf(X2)
end
