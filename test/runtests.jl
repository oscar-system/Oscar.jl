using JToric
using Test

@testset "JToric.jl" begin
    # Perform tests for Hirzebruch surface
    H5 = hirzebruch_surface( 5 )
    ntv_polymake2gap( H5.polymakeNTV )
    @test JToric.is_normal_variety( H5 ) == true
    @test JToric.is_affine( H5 ) == false
    @test JToric.is_projective( H5 ) == true
    @test JToric.is_smooth( H5 ) == true
    @test JToric.is_complete( H5 ) == true
    @test JToric.has_torusfactor( H5 ) == false
    @test JToric.is_orbifold( H5 ) == true
    @test JToric.is_simplicial( H5 ) == true
    @test JToric.is_isomorphic_to_projective_space( H5 ) == false
    @test JToric.is_direct_product_of_projective_spaces( H5 ) == false
    @test JToric.is_q_gorenstein( H5 ) == true
    @test size( affine_open_covering( H5 ) )[ 1 ] == 4
    cox_ring( H5 )
    cox_ring( H5, "u" )
    coordinate_ring_of_torus( H5, [ "u", "v", "w", "z" ] )
    @test size( list_of_variables_of_cox_ring( H5 ) )[ 1 ] == 4
    class_group( H5 )
    torus_invariant_divisor_group( H5 )
    map_from_character_to_principal_divisor( H5 )
    map_from_weil_divisors_to_class_group( H5 )
    @test dim( H5 ) == 2
    @test dim_of_torusfactor( H5 ) == 0
    coordinate_ring_of_torus( H5 )
    @test size( list_of_variables_of_coordinate_ring_of_torus( H5 ) )[ 1 ] == 4
    @test size( is_product_of( H5 ) )[ 1 ] == 1
    character_lattice( H5 )
    divisors = torus_invariant_prime_divisors( H5 )
    @test ( 0 in [ is_primedivisor( d ) for d in divisors ] ) == false
    irrelevant_ideal( H5 )
    stanley_reisner_ideal( H5 )
    morphism_from_cox_variety( H5 )
    cox_variety( H5 )
    fan_of_variety( H5 )
    fan( H5 )
    cartier_torus_invariant_divisor_group( H5 )
    picard_group( H5 )
    @test name_of_variety( H5 ) == "No name set for this variety"
    set_name_of_variety( H5, "Hirzebruch surface" )
    @test name_of_variety( H5 ) == "Hirzebruch surface"
    zariski_cotangent_sheaf( H5 )
    cotangent_sheaf( H5 )
    @test euler_characteristic( H5 ) == 4
    # UnderlyingSheaf( H5 ) <- Error in gap
    character_to_rational_function( [1,2,3,4], H5 )
    @test size( weil_divisors_of_variety( H5 ) )[ 1 ] == 4
    @test size( factors( H5 ) )[ 1 ] == 1
    zariski_cotangent_sheaf_via_euler_sequence( H5 )
    zariski_cotangent_sheaf_via_poincare_residue_map( H5 )
    @test ith_betti_number( H5, 0 ) == 1
    @test ith_betti_number( H5, 1 ) == 0
    @test ith_betti_number( H5, 2 ) == 2
    @test ith_betti_number( H5, 3 ) == 0
    @test ith_betti_number( H5, 4 ) == 1
    @test nr_of_q_rational_points( H5, 1 ) == 4
    
    # Construct delPezzo surfaces
    del_pezzo( -1 )
    del_pezzo( 0 )
    del_pezzo( 1 )
    del_pezzo( 2 )
    del_pezzo( 3 )
    del_pezzo( 4 )
    
    # Perform tests for projective space
    P2 = JToric.projective_space( 2 )
    @test JToric.is_normal_variety( P2 ) == true
    @test JToric.is_affine( P2 ) == false
    @test JToric.is_projective( P2 ) == true
    @test JToric.is_smooth( P2 ) == true
    @test JToric.is_complete( P2 ) == true
    @test JToric.has_torusfactor( P2 ) == false
    @test JToric.is_orbifold( P2 ) == true
    @test JToric.is_simplicial( P2 ) == true
    @test JToric.is_isomorphic_to_projective_space( P2 ) == true
    @test JToric.is_direct_product_of_projective_spaces( P2 ) == true
    @test ith_betti_number( P2, 0 ) == 1
    @test ith_betti_number( P2, 1 ) == 0
    @test ith_betti_number( P2, 2 ) == 1
    @test ith_betti_number( P2, 3 ) == 0
    @test ith_betti_number( P2, 4 ) == 1
    @test nr_of_q_rational_points( P2, 2 ) == 4

    # Perform tests for blowup of projective space
    blowup_variety = blowup_on_ith_minimal_torus_orbit( P2, 1 )
    @test JToric.is_normal_variety( blowup_variety ) == true
    @test JToric.is_affine( blowup_variety ) == false
    @test JToric.is_projective( blowup_variety ) == true
    @test JToric.is_smooth( blowup_variety ) == true
    @test JToric.is_complete( blowup_variety ) == true
    @test JToric.has_torusfactor( blowup_variety ) == false
    @test JToric.is_orbifold( blowup_variety ) == true
    @test JToric.is_simplicial( blowup_variety ) == true
    @test JToric.is_isomorphic_to_projective_space( blowup_variety ) == false
    @test JToric.is_direct_product_of_projective_spaces( blowup_variety ) == false
    @test ith_betti_number( blowup_variety, 0 ) == 1
    @test ith_betti_number( blowup_variety, 1 ) == 0
    @test ith_betti_number( blowup_variety, 2 ) == 2
    @test ith_betti_number( blowup_variety, 3 ) == 0
    @test ith_betti_number( blowup_variety, 4 ) == 1
    @test nr_of_q_rational_points( blowup_variety, 4 ) == 13
    
    # Perform tests for direct product of Hirzebruch surface and projective space
    v = H5 * P2;
    @test JToric.is_normal_variety( v ) == true
    @test JToric.is_affine( v ) == false
    @test JToric.is_projective( v ) == true
    @test JToric.is_smooth( v ) == true
    @test JToric.is_complete( v ) == true
    @test JToric.has_torusfactor( v ) == false
    @test JToric.is_orbifold( v ) == true
    @test JToric.is_simplicial( v ) == true
    @test JToric.is_isomorphic_to_projective_space( v ) == false
    @test JToric.is_direct_product_of_projective_spaces( v ) == false
    @test size( factors( v ) )[ 1 ] == 2
    @test ith_betti_number( v, 0 ) == 1
    @test ith_betti_number( v, 1 ) == 0
    @test ith_betti_number( v, 2 ) == 3
    @test ith_betti_number( v, 3 ) == 0
    @test ith_betti_number( v, 4 ) == 4
    @test ith_betti_number( v, 5 ) == 0
    @test ith_betti_number( v, 6 ) == 3
    @test ith_betti_number( v, 7 ) == 0
    @test ith_betti_number( v, 8 ) == 1
    @test nr_of_q_rational_points( v, 3 ) == 106
    
    # Compute properties of toric divisors on Hirzebruch surface
    D=create_divisor( [ 0,0,0,0 ], H5 )
    @test JToric.is_cartier( D ) == true
    @test JToric.is_principal( D ) == true
    @test JToric.is_primedivisor( D ) == false
    @test JToric.is_basepoint_free( D ) == true
    @test JToric.is_ample( D ) == false
    @test JToric.is_very_ample( D ) == "fail"
    @test JToric.is_nef( D ) == true
    D2 = divisor_of_character( [ 1,2 ], H5 )
    @test JToric.is_cartier( D2 ) == true
    @test JToric.is_principal( D2 ) == true
    @test JToric.is_primedivisor( D2 ) == false
    @test JToric.is_basepoint_free( D2 ) == true
    @test JToric.is_ample( D2 ) == false
    @test JToric.is_very_ample( D2 ) == "fail"
    @test JToric.is_nef( D2 ) == true
    D3 = divisor_of_class( H5, [ 1,2 ] )
    @test JToric.is_cartier( D3 ) == true
    @test JToric.is_principal( D3 ) == false
    @test JToric.is_primedivisor( D3 ) == false
    @test JToric.is_basepoint_free( D3 ) == true
    @test JToric.is_ample( D3 ) == true
    @test JToric.is_very_ample( D3 ) == true
    @test JToric.is_nef( D3 ) == true
end
