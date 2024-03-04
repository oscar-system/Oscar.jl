@testset "FreeAssAlgIdeal.basic" begin
  Zt = polynomial_ring(ZZ, "t")[1]
  R, (x, y, z) = free_associative_algebra(Zt, ["x", "y", "z", "w"])
  I = ideal(R, [x*y*x, y*z^2])
  @test base_ring(I) == R
  for p in gens(R)
    @test parent(p) == R
  end
end

@testset "FreeAssAlgIdeal.printing" begin
  R, (x, y, z) = free_associative_algebra(GF(5), ["x", "y", "z", "w"])
  I = ideal(R, [x*y*x, y*z^2])
  @test length(string(I)) > 3
end

@testset "FreeAssAlgIdeal.membership" begin
  R, (x, y, z) = free_associative_algebra(QQ, ["x", "y", "z"])
  I = ideal(R, [x*y - y*x, x*z - z*x])
  @test !in(x, I, 5)
  @test !in(x, I, 10)
  @test in(x*y*z - y*z*x, I, 9) # 9 should be enough

  f1 = x*y + y*z
  I2 = ideal([f1])
  @test !ideal_membership(x*y, I2, 3)
  @test ideal_membership(f1, I2, 4) 
  @test ideal_membership(f1, I2) 
  @test !isdefined(I2, :gb)
  gb = groebner_basis(I2, 3; protocol=true)
  @test isdefined(I2, :gb)
end

@testset "FreeAssAlgIdeal.utils" begin
  R, (x, y, z) = free_associative_algebra(QQ, ["x", "y", "z"])
  I = Oscar.ideal(R, [x*y - y*x, x*z - z*x])
  @test isa(ngens(I),Int)
  @test isequal(ngens(I),2)
  @test isa(gen(I,ngens(I)),FreeAssAlgElem{QQFieldElem})
  @test isa(gens(I),Vector)

  lpring, _  = Oscar._to_lpring(R, 3)
  @test isa(lpring,NCRing)

  _, (x, y, z) = Singular.FreeAlgebra(QQ, ["x", "y","z"],6)
  free, _ = free_associative_algebra(QQ, ["x", "y", "z"])
  f1 = x*y + y*z

  F1 = free(f1)
  @test isa(F1,FreeAssAlgElem)
end 

@testset "FreeAssAlgIdeal.groebner_basis" begin
    free, (x,y,z) = free_associative_algebra(QQ, ["x", "y", "z"])
    f1 = x*y + y*z
    f2 = x^2 + y^2
    I = ideal([f1, f2])

    gb = groebner_basis(I, 3; protocol=true)
    @test maximum(total_degree.(gb))==3
    @test isdefined(I, :gb)
    gb2 = groebner_basis([f1, f2], 5; protocol=false)
    @test maximum(total_degree.(gb2))==5
end
