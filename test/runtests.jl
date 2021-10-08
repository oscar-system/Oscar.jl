using Test
using JToric
import JToric: Polymake

H5 = hirzebruch_surface( 5 )
P2 = JToric.projective_space( 2 )

@testset "Hirzebruch surface" begin
    ntv_polymake2gap( H5.polymakeNTV )
    @test JToric.isnormal( H5 ) == true
    @test JToric.isaffine( H5 ) == false
    @test JToric.isprojective( H5 ) == true
    @test JToric.issmooth( H5 ) == true
    @test JToric.iscomplete( H5 ) == true
    @test JToric.has_torusfactor( H5 ) == false
    @test JToric.is_orbifold( H5 ) == true
    @test JToric.issimplicial( H5 ) == true
    @test JToric.is_gorenstein( H5 ) == true
    @test JToric.is_q_gorenstein( H5 ) == true
    @test JToric.is_fano( H5 ) == false
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
    @test JToric.isnormal( P2 ) == true
    @test JToric.isaffine( P2 ) == false
    @test JToric.isprojective( P2 ) == true
    @test JToric.issmooth( P2 ) == true
    @test JToric.iscomplete( P2 ) == true
    @test JToric.has_torusfactor( P2 ) == false
    @test JToric.is_orbifold( P2 ) == true
    @test JToric.issimplicial( P2 ) == true
    @test ith_betti_number( P2, 0 ) == 1
    @test ith_betti_number( P2, 1 ) == 0
    @test ith_betti_number( P2, 2 ) == 1
    @test ith_betti_number( P2, 3 ) == 0
    @test ith_betti_number( P2, 4 ) == 1
end

D=ToricDivisor( [ 0,0,0,0 ], H5 )
        
@testset "Divisors" begin
    # Compute properties of toric divisors on Hirzebruch surface
    @test JToric.iscartier( D ) == true
    @test JToric.isprincipal( D ) == true
    @test JToric.isbasepoint_free( D ) == true
    @test JToric.isample( D ) == false
    @test JToric.isvery_ample( D ) == false
    @test JToric.isnef( D ) == true
    @test JToric.isintegral( D ) == true
    @test JToric.isq_cartier( D ) == true
    @test JToric.iscartier( D2 ) == true
    @test JToric.isprincipal( D2 ) == true
    @test JToric.isbasepoint_free( D2 ) == true
    @test JToric.isample( D2 ) == false
    @test JToric.isvery_ample( D2 ) == false
    @test JToric.isnef( D2 ) == true
    @test JToric.isintegral( D2 ) == true
    @test JToric.isq_cartier( D2 ) == true
    @test JToric.iscartier( D3 ) == true
    @test JToric.isprincipal( D3 ) == false
    @test JToric.isbasepoint_free( D3 ) == true
    @test JToric.isample( D3 ) == true
    @test JToric.isvery_ample( D3 ) == true
    @test JToric.isnef( D3 ) == true
    @test JToric.isintegral( D3 ) == true
    @test JToric.isq_cartier( D3 ) == true
end

using Oscar

@testset "Polytopes of divisors" begin
    p = JToric.polyhedron_of_divisor( D )
    @test dim( p ) == 0
    @test ambient_dim( p ) == 2
    p2 = JToric.polyhedron_of_divisor( D2 )
    @test dim( p2 ) == 0
    @test ambient_dim( p2 ) == 2
    p3 = JToric.polyhedron_of_divisor( D3 )
    @test dim( p3 ) == 2
    @test ambient_dim( p3 ) == 2
end

@testset "Affine toric varieties" begin
    C = Oscar.positive_hull([1 1; -1 1])
    antv = AffineNormalToricVariety(C)
    JToric.show( antv )
    @test JToric.issmooth( antv ) == false
    @test JToric.is_orbifold( antv ) == true
    @test toric_ideal_binomial_generators( antv ) == [-1 -1 2]
    ntv = NormalToricVariety(C)
    @test JToric.isaffine( ntv ) == true
end

@testset "Toric varieties from polyhedral fans" begin
    square = Oscar.cube(2)
    nf = Oscar.normal_fan(square)
    ntv = NormalToricVariety(nf)
    @test JToric.iscomplete( ntv ) == true
    ntv2 = NormalToricVariety(square)
    @test JToric.iscomplete( ntv2 ) == true
end
