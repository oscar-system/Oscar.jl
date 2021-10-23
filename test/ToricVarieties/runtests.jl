using Oscar
using Test

C = Oscar.positive_hull([1 1; -1 1])
antv = AffineNormalToricVariety(C)

@testset "Affine toric varieties" begin
    @test issmooth( antv ) == false
    @test isorbifold( antv ) == true
    @test length( affine_open_covering( antv ) ) == 1
    # @test toric_ideal_binomial_generators( antv ) == [-1 -1 2]
    @test rank(torusinvariant_divisor_group( antv )) == 2
    @test rank(character_lattice( antv )) == 2
    map = map_from_character_to_principal_divisors( antv )
    @test rank(domain(map)) == 2
    @test rank(codomain(map)) == 2
    @test length(gens(codomain(map_from_weil_divisors_to_class_group(antv)))) == 2
    @test length(gens(class_group(antv))) == 2
    @test ngens(cox_ring(antv)) == 2
end

ntv = NormalToricVariety(C)

@testset "Affine toric varieties created as general normal toric varieties" begin
    @test isaffine( ntv ) == true
    associated_affine_variety = AffineNormalToricVariety( ntv )
    @test length( affine_open_covering( associated_affine_variety ) ) == 1
    @test length( list_of_variables_of_cox_ring( associated_affine_variety ) ) == 2
end

square = Oscar.cube(2)
nf = Oscar.normal_fan(square)
ntv = NormalToricVariety(nf)
ntv2 = NormalToricVariety(square)

@testset "Toric varieties from polyhedral fans" begin
    @test iscomplete( ntv ) == true
    @test iscomplete( ntv2 ) == true
end

H5 = hirzebruch_surface( 5 )

@testset "Hirzebruch surface" begin
    @test isnormal( H5 ) == true
    @test isaffine( H5 ) == false
    @test isprojective( H5 ) == true
    @test issmooth( H5 ) == true
    @test iscomplete( H5 ) == true
    @test hastorusfactor( H5 ) == false
    @test isorbifold( H5 ) == true
    @test issimplicial( H5 ) == true
    @test isgorenstein( H5 ) == true
    @test isq_gorenstein( H5 ) == true
    @test isfano( H5 ) == false
    nef_cone( H5 )
    mori_cone( H5 )
    @test dim( H5 ) == 2
    @test dim_of_torusfactor( H5 ) == 0
    @test euler_characteristic( H5 ) == 4
    @test ith_betti_number( H5, 0 ) == 1
    @test ith_betti_number( H5, 1 ) == 0
    @test ith_betti_number( H5, 2 ) == 2
    @test ith_betti_number( H5, 3 ) == 0
    @test ith_betti_number( H5, 4 ) == 1
    @test length( affine_open_covering( H5 ) ) == 4
    @test fan_of_variety( H5 ).pm_fan.FAN_DIM == 2
    @test rank(torusinvariant_divisor_group(H5)) == 4
    @test rank(character_lattice(H5)) == 2
    map = map_from_character_to_principal_divisors( H5 )
    @test rank(domain(map)) == 2
    @test rank(codomain(map)) == 4
    @test rank(class_group( H5 )) == 2
    @test rank(codomain(map_from_weil_divisors_to_class_group(H5))) == 2
    @test length( list_of_variables_of_cox_ring( H5 ) ) == 4
end

@testset "delPezzo surfaces" begin
    @test_throws ArgumentError del_pezzo( -1 )
    del_pezzo( 0 )
    del_pezzo( 1 )
    del_pezzo( 2 )
    del_pezzo( 3 )
    @test_throws ArgumentError del_pezzo( 4 )
end

P2 = toric_projective_space( 2 )

@testset "Projective space" begin
    @test isnormal( P2 ) == true
    @test isaffine( P2 ) == false
    @test isprojective( P2 ) == true
    @test issmooth( P2 ) == true
    @test iscomplete( P2 ) == true
    @test hastorusfactor( P2 ) == false
    @test isorbifold( P2 ) == true
    @test issimplicial( P2 ) == true
    @test ith_betti_number( P2, 0 ) == 1
    @test ith_betti_number( P2, 1 ) == 0
    @test ith_betti_number( P2, 2 ) == 1
    @test ith_betti_number( P2, 3 ) == 0
    @test ith_betti_number( P2, 4 ) == 1
    S = cox_ring( P2 )
    @test ngens( S ) == 3
end

D=ToricDivisor( [ 0,0,0,0 ], H5 )
        
@testset "Divisors" begin
    @test iscartier( D ) == true
    @test isprincipal( D ) == true
    @test isbasepoint_free( D ) == true
    @test isample( D ) == false
    @test isvery_ample( D ) == false
    @test isnef( D ) == true
    @test isintegral( D ) == true
    @test isq_cartier( D ) == true
end

@testset "Polytopes of divisors" begin
    p = polyhedron( D )
    @test dim( p ) == 0
    @test ambient_dim( p ) == 2
end
