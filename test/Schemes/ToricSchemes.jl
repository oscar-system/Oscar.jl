using Oscar
using Test

@testset "Toric schemes" begin
  
  C = positive_hull([-1 1; 1 1])
  antv = AffineNormalToricVariety(C)
  affine_toric_scheme = ToricSpec(antv)
  
  @testset "A simplicial (and not smooth) affine toric scheme" begin
    @test is_smooth(underlying_scheme(affine_toric_scheme)) == false
    @test is_simplicial(affine_normal_toric_variety(affine_toric_scheme)) == true
    @test is_simplicial(affine_toric_scheme) == true
    @test dim(cone(affine_toric_scheme)) == 2
    @test polarize(cone(affine_toric_scheme)) == dual_cone(affine_toric_scheme)
    @test nrows(hilbert_basis(affine_toric_scheme)) == 3
end
  
  S = hirzebruch_surface(3)
  X = ToricCoveredScheme(S)
  
  @testset "Toric Scheme of Hirzebruch surface F3" begin
    @test issmooth(X) == true
    @test issmooth(underlying_scheme(X)) == true
    @test issmooth(normal_toric_variety(X)) == true
    @test dim(fan(X)) == 2
  end
  
  IP1 = projective_space(NormalToricVariety, 1)
  IP1xIP1 = IP1*IP1
  Y = ToricCoveredScheme(IP1xIP1)
  
  @testset "Product of projective spaces" begin
    @test is_smooth(Y) == true
    @test is_simplicial(normal_toric_variety(Y)) == true
    @test length(values(glueings(default_covering(Y)))) == 16
  end
end
