@testset "FreeAssociativeAlgebraIdeal.basic" begin
  Zt = polynomial_ring(ZZ, :t)[1]
  R, (x, y, z) = free_associative_algebra(Zt, [:x, :y, :z, :w])
  I = ideal(R, [x*y*x, y*z^2])
  @test base_ring(I) == R
  for p in gens(R)
    @test parent(p) == R
  end
end

@testset "FreeAssociativeAlgebraIdeal.printing" begin
  R, (x, y, z) = free_associative_algebra(GF(5), [:x, :y, :z, :w])
  I = ideal(R, [x*y*x, y*z^2])
  @test length(string(I)) > 3
end

@testset "FreeAssociativeAlgebraIdeal.membership" begin
  R, (x, y, z) = free_associative_algebra(QQ, [:x, :y, :z])
  I = ideal(R, [x*y - y*x, x*z - z*x])
  @test !ideal_membership(x, I, 5)
  @test !ideal_membership(x, I, 10)
  @test ideal_membership(x*y*z - y*z*x, I, 9) # 9 should be enough

  f1 = x*y + y*z
  I2 = ideal([f1])
  @test !ideal_membership(x*y, I2, 3)
  @test ideal_membership(f1, I2, 4) 
  @test ideal_membership(f1, I2.gens, 4)
  @test ideal_membership(f1, I2) 
  @test ideal_membership(f1, I2.gens)
  gb = groebner_basis(I2, 3; protocol=false)
  @test isdefined(I2, :gb)
  @test length(gens(I * I2)) == 2
end

@testset "FreeAssociativeAlgebraIdeal.utils" begin
  R, (x, y, z) = free_associative_algebra(QQ, [:x, :y, :z])
  I = ideal(R, [x*y - y*x, x*z - z*x])
  @test base_ring(I) == R
  @test isa(ngens(I),Int)
  @test isequal(ngens(I),2)
  @test isa(gen(I,ngens(I)),FreeAssociativeAlgebraElem{QQFieldElem})
  @test isa(gens(I),Vector)

  lpring, _  = Oscar._to_lpring(R, 3)
  @test isa(lpring,NCRing)

  _, (x, y, z) = Singular.FreeAlgebra(QQ, ["x", "y","z"],6)
  free, _ = free_associative_algebra(QQ, [:x, :y, :z])
  f1 = x*y + y*z

  F1 = free(f1)
  @test isa(F1,FreeAssociativeAlgebraElem)
end 

@testset "FreeAssociativeAlgebraIdeal.groebner_basis" begin
    free, (x,y,z) = free_associative_algebra(QQ, [:x, :y, :z])
    f1 = x*y + y*z
    f2 = x^2 + y^2
    I = ideal([f1, f2])

    gb = groebner_basis(I, 3; protocol=false)
    @test !is_groebner_basis(gb)
    @test maximum(total_degree.(gb))==3
    @test isdefined(I, :gb)
    gb2 = groebner_basis([f1, f2], 5; protocol=false)
    @test maximum(total_degree.(gb2))==5
    @test !is_groebner_basis(gb2)
end

@testset "FreeAssociativeAlgebraIdeal.groebner_basis.quantum_automorphism_group" begin
  M = uniform_matroid(3,4)
  qAut1 = gens(quantum_automorphism_group(M))
  gb1 = groebner_basis(qAut1; interreduce=true)
  gb2 = groebner_basis(qAut1, 4,ordering=:deglex,interreduce=true)
  @test sort(leading_monomial.(gb2)) == sort(leading_monomial.(gb1))
  @test is_groebner_basis(gb1)
  @test is_groebner_basis(gb2)
  I1 = ideal(gb1); I2 = ideal(gb2)
  @test I1 == I2

  qAut2 = quantum_symmetric_group(4)
  gena = gens(qAut2)

  gb3 = groebner_basis(gena; interreduce=false)
  I1 = ideal(copy(gb3))
  groebner_basis(I1)
  interreduce!(gb3)
  interreduce!(I1)
  @test length(I1.gb) == 78
  @test length(gb3) == 78


  gb4 = groebner_basis(gena, 4,ordering=:deglex,interreduce=true)
  @test is_groebner_basis(gb3)
  @test is_groebner_basis(gb4)
  @test sort(leading_monomial.(gb3)) == sort(leading_monomial.(gb4))

  I1 = ideal(gb3); I2 = ideal(gb4)
  @test gb4[end] in I1
  @test normal_form(gb4[end],I1) == 0
  @test I1 == I2

  x = base_ring(I1)[1]; y = base_ring(I1)[7]
  @test !(x*y - y*x in I1)

  gb5 = groebner_basis(gens(quantum_symmetric_group(3)); interreduce=true)
  A = parent(gb5[1])
  I3 = ideal(Oscar.IdealGens(A, gb5))
  @test is_groebner_basis(I3.gens)

  x = base_ring(ideal(gb5))[1]; y = base_ring(ideal(gb5))[7]
  @test x*y - y*x in ideal(gb5)
end

@testset "FreeAssociativeAlgebraIdeal.groebner_basis.f4" begin
    S1 = quantum_symmetric_group(4);
    S1_alt = quantum_symmetric_group(4);
    gb_f4 = groebner_basis(S1;algorithm=:f4)
    gb_bb = groebner_basis(S1_alt;algorithm=:buchberger);
    @test is_groebner_basis(gb_f4)
    @test is_groebner_basis(gb_bb)
    @test S1 == S1_alt
end
