@testset "singular locus" begin
  R, (x,y,z) = QQ["x", "y", "z"]

  I = ideal(R, [x^2 - y^2 + z^2])
  J = ideal(R, [x-1, y-2])
  X = Spec(R, I*J, units_of(R))
  Z = singular_locus(X)
  @test issubset(subscheme(X, [x,y,z]), Z)
  @test issubset(subscheme(X, [z^2-3, y-2, x-1]), Z)
  @test issmooth(X) == false
  Y = Spec(R, ideal(R, [x^2 - y^2 + z^2 - 1]))
  @test issmooth(Y)
end
