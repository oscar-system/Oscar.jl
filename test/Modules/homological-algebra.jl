@testset "mpoly_affine_homological-algebra.fitting_ideal" begin
 R, (x, y) = polynomial_ring(QQ, [:x, :y]);
 F = free_module(R, 3)
 I = ideal(R, [zero(R)])
 @test fitting_ideal(F, 2) == I
 @test fitting_ideal(F, 3) == R
end

@testset "mpoly_affine_homological-algebra.is_flat" begin
 R, (x, y) = polynomial_ring(QQ, [:x, :y]);
 F = free_module(R, 3)
 @test is_flat(F) == true
end

@testset "mpoly_affine_homological-algebra.non_flat_locus" begin
 R, (x, y) = polynomial_ring(QQ, [:x, :y]);
 F = free_module(R, 3)
 @test non_flat_locus(F) == ideal(R, [one(R)])
end

@testset "mpoly_affine_homological-algebra.is_regular_sequence" begin
 R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);
 F = free_module(R, 1)
 U = matrix([x*y-y])
 M = quo(F, U)[1]
 V =  [x, x*z-z]
 @test is_regular_sequence(V, M) == true
 W = [x*z-z, x]
 @test is_regular_sequence(W, M) == false
end

@testset "mpoly_affine_homological-algebra.koszul_matrix" begin
 R, (x, y) = polynomial_ring(QQ, [:x, :y]);
 V = gens(R)
 KM = koszul_matrix(V, 1)
 KMS = Oscar._koszul_matrix_from_singular(V, 1)
 @test nrows(KM) == nrows(KMS) == 2
 @test ncols(KM) == ncols(KMS)
 #@test KM[1] == R[1] # Custom test for the deprecated Singular code
end

@testset "mpoly_affine_homological-algebra.koszul_complex" begin
 R, (x, y) = polynomial_ring(QQ, [:x, :y]);
 V = gens(R)
 K = koszul_complex(V)
 KS = Oscar._koszul_complex_from_singular(V)
 KM = matrix(map(K, 2))
 KMS = matrix(map(KS, 2))
 @test nrows(KM) == nrows(KMS) == 1
 @test ncols(KM) == ncols(KMS)
 #@test KM[1, 1] == -R[2]# Custom test for the deprecated Singular code
end

@testset "mpoly_affine_homological-algebra.koszul_homology" begin
 R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);
 F = free_module(R, 1)
 U = matrix([x*y-y])
 M = quo(F, U)[1]
 V =  [x, x*z-z]
 @test is_zero(koszul_homology(V, M, 0)) == false
 @test is_zero(Oscar._koszul_homology_from_singular(V, M, 0)) == false
 @test is_zero(koszul_homology(V, M, 1)) == true
 @test is_zero(Oscar._koszul_homology_from_singular(V, M, 1)) == true
end

@testset "mpoly_affine_homological-algebra.depth" begin
 R, (x, y) = polynomial_ring(QQ, [:x, :y]);
 F = free_module(R, 1)
 U = matrix([x*y])
 M = quo(F, U)[1]
 I = ideal(R, gens(R))
 @test depth(I, M) == Oscar._depth_from_singular(I, M) == 1

 R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);
 F = free_module(R, 1);
 I = ideal(R, [x*z-z, x*y-y, x])
 @test depth(I, F) == Oscar._depth_from_singular(I, F) == 3
end

@testset "singular vs oscar comparison" begin
  n = 3
  R, _ = polynomial_ring(QQ, :x=>1:n)
  FR = FreeMod(R, 1)
  IR = ideal(R, gens(R))
  IR = IR*IR
  f = gens(IR)
  K0 = koszul_complex(f)
  K1 = koszul_complex(f, FR)
  HK0 = [homology(K0, i) for i in 0:ngens(IR)]
  HK1 = [homology(K1, i) for i in 0:ngens(IR)]
  HK2 = [koszul_homology(f, FR, i) for i in 0:ngens(IR)]
  @test iszero.(HK1) == iszero.(HK2) == iszero.(HK0)
end

@testset "generic depth routine" begin
  P, _ = QQ[:u, :v]
  R, (x, y) = polynomial_ring(P, [:x, :y]);
  F = free_module(R, 1)
  U = matrix([x*y])
  M = quo(F, U)[1]
  I = ideal(R, gens(R))
  @test depth(I, M) == 1

  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);
  F = free_module(R, 1);
  I = ideal(R, [x*z-z, x*y-y, x])
  @test depth(I, F) == 3
end

@testset "twisting to degree zero" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  S3 = graded_free_module(S, [0, 0, 0])
  v = x*S3[1] + y*S3[2] + z*S3[3]
  C = koszul_complex(v)
  CC = Oscar._make_homogeneous(C)
  str = Oscar.strand(CC, 0)
  str[0:3]
end

