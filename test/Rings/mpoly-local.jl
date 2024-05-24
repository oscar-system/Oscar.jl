@testset "mpoly-loc constructors" begin
  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  m = ideal(R, [y - 1, x - 2, z - 3])
  Q = localization(R, m)
  I = ideal(Q, [x - 2, (y - 1)^2*z])
  a = Q(x//(1 + x))

  @test gens(I) == Q.([x - 2, y^2*z - 2*y*z + z])
  @test fraction(a) == x // (1 + x)
end

@testset "mpoly-loc operations" begin
  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  m = ideal(R, [y - 1, x - 2, z - 3])
  Q = localization(R, m)
  I = ideal(Q, [x - 2, (y - 1)^2*z])
  J = ideal(Q, [x - 2, y - 1])
  a = Q(x//(1 + x))
  b = Q((x + y)//(3 + x - y))

  K1 = I+J
  K2 = I*J
  K3 = J^3
  c1 = a+b
  c2 = a*b

  @test iszero(a-a)
  @test a*deepcopy(a) == a^2
  @test det(matrix(Q, [1 x; y z])) == Q(z-x*y)

  @test gens(K1) == Q.([x - 2, y^2*z - 2*y*z + z, x - 2, y - 1])
  @test gens(K2) == Q.([x^2 - 4*x + 4, x*y - x - 2*y + 2, x*y^2*z - 2*x*y*z + x*z - 2*y^2*z + 4*y*z - 2*z, y^3*z - 3*y^2*z + 3*y*z - z])
  @test gens(K3) == Q.([x^3 - 6*x^2 + 12*x - 8, x^2*y - x^2 - 4*x*y + 4*x + 4*y - 4, x^2*y - x^2 - 4*x*y + 4*x + 4*y - 4, x*y^2 - 2*x*y + x - 2*y^2 + 4*y - 2, x^2*y - x^2 - 4*x*y + 4*x + 4*y - 4, x*y^2 - 2*x*y + x - 2*y^2 + 4*y - 2, x*y^2 - 2*x*y + x - 2*y^2 + 4*y - 2, y^3 - 3*y^2 + 3*y - 1])
  @test fraction(c1) == (x*(3 + x - y) + (1 + x)*(x + y))//((1 + x)*(3 + x - y))
  @test fraction(c2) == (x*(x + y))//((1 + x)*(3 + x - y))

  (x, y, z) = map(Q, (x, y, z))
  @test inv(x) + inv(y) == (x + y)*inv(x*y)
  @test_throws Exception x*inv(x - 2)
end

#=
# Disabled because of deprecation
@testset "mpoly-loc groebner" begin
  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  m = ideal(R, [y - 1, x - 2, z - 3])
  Q = localization(R, m)
  I = ideal(Q, [x - 2, (y - 1)^2*z])
  J = ideal(Q, [x - 1, y - 1])

  groebner_basis(I)
  GI  = collect(values(I.gb))[1]
  groebner_basis(J)
  GJ  = collect(values(J.gb))[1]
  @test GI.O == Q.([x - 2, (y - 1)^2])
  @test GJ.O == [one(Q)]
end
=#
