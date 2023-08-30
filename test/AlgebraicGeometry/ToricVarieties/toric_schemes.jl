using Oscar
using Test

@testset "Toric schemes (set_attributes = $set_attributes)" for set_attributes in [true, false]
  
  C = positive_hull([-1 1; 1 1])
  antv = affine_normal_toric_variety(C; set_attributes = set_attributes)
  
  @testset "A simplicial (and not smooth) affine toric scheme" begin
    @test is_smooth(underlying_scheme(antv)) == is_smooth(antv)
    @test dim(underlying_scheme(antv)) == dim(antv)
    @test polarize(cone(antv)) == dual_cone(antv)
    @test ngens(ambient_coordinate_ring(antv)) == nrows(hilbert_basis(antv)) == 3
  end
  
  X = hirzebruch_surface(NormalToricVariety, 3; set_attributes = set_attributes)
  
  @testset "Toric Scheme of Hirzebruch surface F3" begin
    @test is_smooth(underlying_scheme(X)) == is_smooth(X)
    @test dim(underlying_scheme(X)) == dim(X)
  end
  
  IP1 = projective_space(NormalToricVariety, 1; set_attributes = set_attributes)
  Y = IP1*IP1
  
  @testset "Product of projective spaces" begin
    @test is_smooth(underlying_scheme(Y)) == is_smooth(Y)
    @test length(values(glueings(default_covering(Y)))) == 16
  end

  IP2 = projective_space(NormalToricVariety, 2; set_attributes = set_attributes)
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
  
end
