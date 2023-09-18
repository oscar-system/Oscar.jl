using Oscar
using Test

@testset "Toric schemes (set_attributes = $set_attributes)" for set_attributes in [true, false]
  
  C = positive_hull([-1 1; 1 1])
  antv = affine_normal_toric_variety(C; set_attributes = set_attributes)
  
  @testset "A simplicial (and not smooth) affine toric scheme" begin
    @test is_smooth(underlying_scheme(antv)) == is_smooth(antv)
    @test dim(underlying_scheme(antv)) == dim(antv)
    @test polarize(cone(antv)) == weight_cone(antv)
    @test ngens(ambient_coordinate_ring(antv)) == nrows(hilbert_basis(antv)) == 3
  end
  
  X = hirzebruch_surface(NormalToricVariety, 3; set_attributes = set_attributes)
  
  @testset "Toric Scheme of Hirzebruch surface F3" begin
    @test is_smooth(underlying_scheme(X)) == is_smooth(X)
    @test dim(underlying_scheme(X)) == dim(X)
  end
  
  IP1 = projective_space(NormalToricVariety, 1; set_attributes = set_attributes)
  set_coordinate_names(IP1, ["x", "y"])
  Y = IP1*IP1
  
  @testset "Product of projective spaces" begin
    @test is_smooth(underlying_scheme(Y)) == is_smooth(Y)
    @test length(values(glueings(default_covering(Y)))) == 16
  end

  IP2 = projective_space(NormalToricVariety, 2; set_attributes = set_attributes)
  set_coordinate_names(IP2, ["x", "y", "z"])
  X, iso = Oscar.forget_toric_structure(IP2)

  @testset "Forget toric structure" begin
    @test X isa CoveredScheme
    @test codomain(iso) === IP2
    @test is_smooth(X) == is_smooth(IP2)
    @test length(keys(glueings(default_covering(X)))) == 9
    @test domain(iso) === X
    @test codomain(iso) === IP2
    @test domain(inverse(iso)) === IP2
    @test codomain(inverse(iso)) === X 
  end
  
  IP3 = projective_space(NormalToricVariety, 3; set_attributes = set_attributes)
  
  @testset "Toric ideal sheaves" begin
    S = cox_ring(IP3)
    I = ideal(S, gens(S)[2:3])
    II = IdealSheaf(IP3, I)
  end

  S = cox_ring(IP2)
  (x, y, z) = gens(S)
  I = IdealSheaf(IP2, ideal(S, [z*(x-y)]))
  J = IdealSheaf(IP2, ideal(S, [x, y]))
  bl = blow_up(IP2, [1, 1])
  pb_I = pullback(bl, I)
  pb_J = pullback(bl, J)

  @testset "toric blowdown morphism as morphism of covered schemes" begin
    @test scheme(I) === IP2
    @test length(Oscar.maximal_associated_points(I)) == 2
    @test length(Oscar.maximal_associated_points(pb_I)) == 3
    @test dim(J) == 0
    @test dim(pb_J) == 1
  end

  F2 = hirzebruch_surface(NormalToricVariety, 2; set_attributes = set_attributes)
  f = toric_morphism(F2, matrix(ZZ, [[1], [0]]), IP1; check = true)
  S = cox_ring(IP1)
  (x, y) = gens(S)
  K = IdealSheaf(IP1, ideal(S, [x]))
  pb_K = pullback(f, K)
  # The support of pb_K should be nothing but P1, it is the fiber
  # of the P1-fibration that defines the Hirebruch surface over the
  # point defined by K.
  # TODO: Add a test that verifies that indeed pb_K is P1.
  @testset "Hirzebruch surface as P1 fibration over P1" begin
    @test length(Oscar.maximal_associated_points(K)) == 1
    @test dim(K) == 0
    @test length(Oscar.maximal_associated_points(pb_K)) == 1
    @test dim(pb_K) == 1
    @test is_smooth(subscheme(pb_K)) == true
  end

end
