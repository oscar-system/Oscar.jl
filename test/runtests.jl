using JToric
using Test

@testset "JToric.jl" begin
    # load necessary gap-packages
    using CapAndHomalg
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    CapAndHomalg.LoadPackage( "JConvex" )
    CapAndHomalg.LoadPackage( "ToricV" )
    
    # Compute properties of toric varieties on the example of a Hirzebruch surface
    Rays = [[-1,5],[0,1],[1,0],[0,-1]]
    Cones = [[1,2],[2,3],[3,4],[4,1]]
    H5 = toric_variety( Rays, Cones )
    @test JToric.is_normal_variety( H5 ) == true
    @test JToric.is_affine( H5 ) == false
    @test JToric.is_projective( H5 ) == true
    @test JToric.is_smooth( H5 ) == true
    @test JToric.is_complete( H5 ) == true
    @test JToric.has_torusfactor( H5 ) == false
    @test JToric.has_no_torusfactor( H5 ) == true
    @test JToric.is_orbifold( H5 ) == true
    @test JToric.is_simplicial( H5 ) == true
    @test JToric.is_isomorphic_to_projective_space( H5 ) == false
    @test JToric.is_direct_product_of_projective_spaces( H5 ) == false
    
    # Compute properties of toric varieties on the example of a projective space
    P2 = JToric.projective_space( 2 )
    @test JToric.is_normal_variety( P2 ) == true
    @test JToric.is_affine( P2 ) == false
    @test JToric.is_projective( P2 ) == true
    @test JToric.is_smooth( P2 ) == true
    @test JToric.is_complete( P2 ) == true
    @test JToric.has_torusfactor( P2 ) == false
    @test JToric.has_no_torusfactor( P2 ) == true
    @test JToric.is_orbifold( P2 ) == true
    @test JToric.is_simplicial( P2 ) == true
    @test JToric.is_isomorphic_to_projective_space( P2 ) == true
    @test JToric.is_direct_product_of_projective_spaces( P2 ) == true

    # compute properties of toric divisors on the example of the trivial divisor
    D=create_divisor( [ 0,0,0,0 ], H5 )
    @test JToric.is_cartier( D ) == true
    @test JToric.is_principal( D ) == true
    @test JToric.is_primedivisor( D ) == false
    @test JToric.is_basepoint_free( D ) == true
    @test JToric.is_ample( D ) == false
    @test JToric.is_very_ample( D ) == "fail"
    @test JToric.is_numerically_effective( D ) == true

    # compute properties of toric divisors on the example of a non-trivial divisor
    D2 = divisor_of_character( [ 1,2 ], H5 )
    @test JToric.is_cartier( D2 ) == true
    @test JToric.is_principal( D2 ) == true
    @test JToric.is_primedivisor( D2 ) == false
    @test JToric.is_basepoint_free( D2 ) == true
    @test JToric.is_ample( D2 ) == false
    @test JToric.is_very_ample( D2 ) == "fail"
    @test JToric.is_numerically_effective( D2 ) == true
    
    # compute properties of toric divisors on the example of another non-trivial divisor
    D3 = divisor_of_given_class( H5, [ 1,2 ] )
    @test JToric.is_cartier( D3 ) == true
    @test JToric.is_principal( D3 ) == false
    @test JToric.is_primedivisor( D3 ) == false
    @test JToric.is_basepoint_free( D3 ) == true
    @test JToric.is_ample( D3 ) == true
    @test JToric.is_very_ample( D3 ) == true
    @test JToric.is_numerically_effective( D3 ) == true
    
    # compute attributes of toric varieties on the example of the Hirzebruch surface H5 defined above
    cover = AffineOpenCovering( H5 )
    @test size( cover )[ 1 ] == 4
    cox_ring = CoxRing( H5 )
    @test size( ListOfVariablesOfCoxRing( H5 ) )[ 1 ] == 4
    ClassGroup( H5 )
    TorusInvariantDivisorGroup( H5 )
    MapFromCharacterToPrincipalDivisor( H5 )
    MapFromWeilDivisorsToClassGroup( H5 )
    @test Dimension( H5 ) == 2
    @test DimensionOfTorusfactor( H5 ) == 0
    CoordinateRingOfTorus( H5 )
    @test size( ListOfVariablesOfCoordinateRingOfTorus( H5 ) )[ 1 ] == 4
    @test size( IsProductOf( H5 ) )[ 1 ] == 1
    CharacterLattice( H5 )
    divisors = TorusInvariantPrimeDivisors( H5 )
    @test ( 0 in [ is_primedivisor( divisors[ i ] ) for i in 1 : size(divisors)[1] ] ) == false
    IrrelevantIdeal( H5 )
    SRIdeal( H5 )
    MorphismFromCoxVariety( H5 )
    CoxVariety( H5 )
    FanOfVariety( H5 )
    Fan( H5 )
    CartierTorusInvariantDivisorGroup( H5 )
    PicardGroup( H5 )
    @test NameOfVariety( H5 ) == "No name set for this variety"
    SetNameOfVariety( H5, "Hirzebruch surface" )
    @test NameOfVariety( H5 ) == "Hirzebruch surface"
    ZariskiCotangentSheaf( H5 )
    CotangentSheaf( H5 )
    @test EulerCharacteristic( H5 ) == 0
    # UnderlyingSheaf( H5 ) <- Error in gap
    
    # apply methods to toric varieties on the example of the Hirzebruch surface H5 and projective space P2 defined above
    CoordinateRingOfTorus( H5, [ "u", "v", "w", "z" ] )
    CoxRing( H5, "u" )
    v = H5 * P2;
    @test JToric.is_normal_variety( v ) == true
    @test JToric.is_affine( v ) == false
    @test JToric.is_projective( v ) == true
    @test JToric.is_smooth( v ) == true
    @test JToric.is_complete( v ) == true
    @test JToric.has_torusfactor( v ) == false
    @test JToric.has_no_torusfactor( v ) == true
    @test JToric.is_orbifold( v ) == true
    @test JToric.is_simplicial( v ) == true
    @test JToric.is_isomorphic_to_projective_space( v ) == false
    @test JToric.is_direct_product_of_projective_spaces( v ) == false
end
