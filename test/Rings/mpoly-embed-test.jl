@testset "mpoly-embed.mpoly->mpoly" begin

  @test_throws Exception embed_names(QQ["x", "y"][1], QQ["x", "y", "x"][1])

  D, (x1, z1, y1) = PolynomialRing(QQ, ["x", "z", "y"])
  R, (y2, x2, t2, z2) = PolynomialRing(QQ, ["y", "x", "t", "z"])
  m = embed_names(D, R)
  @test domain(m) == D
  @test codomain(m) == R
  @test m(D(1)) == R(1)
  @test m(x1) == x2
  @test m(y1) == y2
  @test m(z1) == z2
  @test m(2 + x1 + y1^2//3 + z1^3) == 2 + x2 + y2^2//3 + z2^3
  @test m(ideal(D, [x1, y1^2, z1^3])) == ideal(R, [x2, y2^2, z2^3])
  @test_throws Exception m(x2)

  D, (x1, z1, y1) = PolynomialRing(QQ, ["x", "z", "y"])
  R, (y2, x2) = PolynomialRing(GF(5), ["y", "x"])
  m = embed_names(D, R)
  @test m(D(1)) == R(1)
  @test m(x1) == x2
  @test m(y1) == y2
  @test m(z1) == R(0)
  @test m(2 + x1 + y1^2//3) == 2 + x2 + y2^2//3
  @test m(2 + x1 + y1^2//3 + z1//2) == 2 + x2 + y2^2//3
  @test m(2 + x1 + y1^2//3 + z1*y1//2) == 2 + x2 + y2^2//3
  @test m(ideal(D, [x1, y1^2, z1^3])) == ideal(R, [x2, y2^2])

  D, (x1, y1) = PolynomialRing(QQ, ["x", "y"])
  R, (y2, x2) = PolynomialRing(FractionField(QQ["t"][1]), ["y", "x"])
  m = embed_names(D, R)
  @test m(D(1)) == R(1)
  @test m(x1) == x2
  @test m(y1) == y2
  @test m(x1 + y1^2//3 + x1^2*y1^3//2) == x2 + y2^2//3 + x2^2*y2^3//2
  @test m(ideal(D, [x1, y1^2])) == ideal(R, [x2, y2^2])

end
