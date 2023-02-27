@testset "mpoly_affine_homological-algebra.fitting_ideal" begin
 R, (x, y) = polynomial_ring(QQ, ["x", "y"]);
 F = free_module(R, 3)
 I = ideal(R, [zero(R)])
 @test fitting_ideal(F, 2) == I
 @test fitting_ideal(F, 3) == R
end

@testset "mpoly_affine_homological-algebra.is_flat" begin
 R, (x, y) = polynomial_ring(QQ, ["x", "y"]);
 F = free_module(R, 3)
 @test is_flat(F) == true
end

@testset "mpoly_affine_homological-algebra.non_flat_locus" begin
 R, (x, y) = polynomial_ring(QQ, ["x", "y"]);
 F = free_module(R, 3)
 @test non_flat_locus(F) == ideal(R, [one(R)])
end

@testset "mpoly_affine_homological-algebra.is_regular_sequence" begin
 R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);
 F = free_module(R, 1)
 U = matrix([x*y-y])
 M = quo(F, U)[1]
 V =  [x, x*z-z]
 @test is_regular_sequence(V, M) == true
 W = [x*z-z, x]
 @test is_regular_sequence(W, M) == false
end

@testset "mpoly_affine_homological-algebra.koszul_matrix" begin
 R, (x, y) = polynomial_ring(QQ, ["x", "y"]);
 V = gens(R)
 KM = koszul_matrix(V, 1)
 @test nrows(KM) == 2
 @test KM[1] == R[1]
end

@testset "mpoly_affine_homological-algebra.koszul_complex" begin
 R, (x, y) = polynomial_ring(QQ, ["x", "y"]);
 V = gens(R)
 K = koszul_complex(V)
 KM = matrix(map(K, 2))
 @test ncols(KM) == 2
 @test KM[1, 1] == -R[2]
end

@testset "mpoly_affine_homological-algebra.koszul_homology" begin
 R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);
 F = free_module(R, 1)
 U = matrix([x*y-y])
 M = quo(F, U)[1]
 V =  [x, x*z-z]
 @test is_zero(koszul_homology(V, M, 0)) == false
 @test is_zero(koszul_homology(V, M, 1)) == true
end

@testset "mpoly_affine_homological-algebra.depth" begin
 R, (x, y) = polynomial_ring(QQ, ["x", "y"]);
 F = free_module(R, 1)
 U = matrix([x*y])
 M = quo(F, U)[1]
 I = ideal(R, gens(R))
 @test depth(I, M) == 1

 R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);
 F = free_module(R, 1);
 I = ideal(R, [x*z-z, x*y-y, x])
 @test depth(I, F) == 3
end
