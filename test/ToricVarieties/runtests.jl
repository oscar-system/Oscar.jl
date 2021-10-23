H5 = hirzebruch_surface( 5 )
P2 = toric_projective_space( 2 )

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
end

@testset "delPezzo surfaces" begin
    # Construct delPezzo surfaces
    @test_throws ArgumentError del_pezzo( -1 )
    del_pezzo( 0 )
    del_pezzo( 1 )
    del_pezzo( 2 )
    del_pezzo( 3 )
    @test_throws ArgumentError del_pezzo( 4 )
end

@testset "Projective space" begin
    # Perform tests for projective space
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
end

D=ToricDivisor( [ 0,0,0,0 ], H5 )
        
@testset "Divisors" begin
    # Compute properties of toric divisors on Hirzebruch surface
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

@testset "Affine toric varieties" begin
    C = Oscar.positive_hull([1 1; -1 1])
    antv = AffineNormalToricVariety(C)
    @test issmooth( antv ) == false
    @test isorbifold( antv ) == true
    @test length( affine_open_covering( antv ) ) == 1
    # @test toric_ideal_binomial_generators( antv ) == [-1 -1 2]
    ntv = NormalToricVariety(C)
    @test isaffine( ntv ) == true
    @test length( affine_open_covering( ntv ) ) == 1
end

@testset "Toric varieties from polyhedral fans" begin
    square = Oscar.cube(2)
    nf = Oscar.normal_fan(square)
    ntv = NormalToricVariety(nf)
    @test iscomplete( ntv ) == true
    ntv2 = NormalToricVariety(square)
    @test iscomplete( ntv2 ) == true
end
