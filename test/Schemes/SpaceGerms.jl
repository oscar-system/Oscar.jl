@testset "SpaceGerms" begin
  R, (x,y,z) = QQ["x", "y", "z"]

  X1 = Spec(R)
  P = ideal(R, [x, y])

  X1_loc, res = germ_at_point(X1, P)
  @test pullback(res)(x) == OO(X1_loc)(x)
  @test inv(pullback(res)(z)) == OO(X1_loc)(1//z)
  @test point(X1_loc) == P

  U = complement_of_ideal(R, [1, 2, 3])
  L, _ = Localization(R, U)

  X2 = Spec(L)
  X2_loc = SpaceGerm(X2)
  @test Oscar.underlying_scheme(X2_loc) == X2
  @test is_subset(representative(X2_loc), X1)
  @test point(X2_loc) == [1, 2, 3]

  X3 = subscheme(X1, [x^2 - y^2 + z^2])
  P = ideal(OO(X3), [z-1])
  X3_loc, res = germ_at_point(X3, P)
  @test z in inverted_set(OO(X3_loc))

  X4 = hypersurface_complement(X1, x)
  P = ideal(OO(X4), y)
  X4_loc, res = germ_at_point(X4, P)

  P2 = ideal(R, y)
  X4_loc2, res = germ_at_point(X4, P2)
  @test is_canonically_isomorphic(X4_loc, X4_loc2)
end
