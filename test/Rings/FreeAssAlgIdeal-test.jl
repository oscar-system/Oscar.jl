@testset "FreeAssAlgIdeal.basic" begin
  Zt = PolynomialRing(ZZ, "t")[1]
  R, (x, y, z) = FreeAssociativeAlgebra(Zt, ["x", "y", "z", "w"])
  I = ideal(R, [x*y*x, y*z^2])
  @test base_ring(I) == R
  for p in gens(R)
    @test parent(p) == R
  end
end

@testset "FreeAssAlgIdeal.printing" begin
  R, (x, y, z) = FreeAssociativeAlgebra(GF(5), ["x", "y", "z", "w"])
  I = ideal(R, [x*y*x, y*z^2])
  @test length(string(I)) > 3
end

@testset "FreeAssAlgIdeal.membership" begin
  R, (x, y, z) = FreeAssociativeAlgebra(QQ, ["x", "y", "z"])
  I = ideal(R, [x*y - y*x, x*z - z*x])
  @test !in(x, I, 5)
  @test !in(x, I, 10)
  @test in(x*y*z - y*z*x, I, 9) # 9 should be enough
end
