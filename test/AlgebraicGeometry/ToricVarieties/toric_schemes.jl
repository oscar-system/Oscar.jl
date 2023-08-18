using Oscar
using Test

@testset "Toric schemes" begin
  
  C = positive_hull([-1 1; 1 1])
  antv = affine_normal_toric_variety(C)
  
  @testset "A simplicial (and not smooth) affine toric scheme" begin
    @test is_smooth(underlying_scheme(antv)) == is_smooth(antv)
    @test dim(underlying_scheme(antv)) == dim(antv)
    @test polarize(cone(antv)) == dual_cone(antv)
    @test ngens(ambient_coordinate_ring(antv)) == nrows(hilbert_basis(antv)) == 3
  end
  
  X = hirzebruch_surface(NormalToricVariety, 3)
  
  @testset "Toric Scheme of Hirzebruch surface F3" begin
    @test is_smooth(underlying_scheme(X)) == is_smooth(X)
    @test dim(underlying_scheme(X)) == dim(X)
  end
  
  IP1 = projective_space(NormalToricVariety, 1)
  Y = IP1*IP1
  
  @testset "Product of projective spaces" begin
    @test is_smooth(underlying_scheme(Y)) == is_smooth(Y)
    @test length(values(glueings(default_covering(Y)))) == 16
  end
end
