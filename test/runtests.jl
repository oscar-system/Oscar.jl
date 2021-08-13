using JToric
using Test

@testset "JToric.jl" begin
    # Write your tests here.
    
    # load necessary gap-packages
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    CapAndHomalg.LoadPackage( "JConvex" )
    CapAndHomalg.LoadPackage( "ToricV" )
    
    # Compute properties of toric varieties on the example of a Hirzebruch surface
    Rays = [[-1,5],[0,1],[1,0],[0,-1]]
    Cones = [[1,2],[2,3],[3,4],[4,1]]
    H5 = JToricVariety( Rays, Cones )
    @test JToric.IsNormalVariety( H5 ) == true
    @test JToric.IsAffine( H5 ) == false
    @test JToric.IsProjective( H5 ) == true
    @test JToric.IsSmooth( H5 ) == true
    @test JToric.IsComplete( H5 ) == true
    @test JToric.HasTorusfactor( H5 ) == false
    @test JToric.HasNoTorusfactor( H5 ) == true
    @test JToric.IsOrbifold( H5 ) == true
    @test JToric.IsSimplicial( H5 ) == true
    @test JToric.IsIsomorphicToProjectiveSpace( H5 ) == false
    @test JToric.IsDirectProductOfPNs( H5 ) == false
    
    # Compute properties of toric varieties on the example of a projective space
    P2 = JToric.ProjectiveSpace( 2 )
    @test JToric.IsNormalVariety( P2 ) == true
    @test JToric.IsAffine( P2 ) == false
    @test JToric.IsProjective( P2 ) == true
    @test JToric.IsSmooth( P2 ) == true
    @test JToric.IsComplete( P2 ) == true
    @test JToric.HasTorusfactor( P2 ) == false
    @test JToric.HasNoTorusfactor( P2 ) == true
    @test JToric.IsOrbifold( P2 ) == true
    @test JToric.IsSimplicial( P2 ) == true
    @test JToric.IsIsomorphicToProjectiveSpace( P2 ) == true
    @test JToric.IsDirectProductOfPNs( P2 ) == true

    # compute properties of toric divisors on the example of the trivial divisor
    D=CreateDivisor( [ 0,0,0,0 ], H5 )
    @test JToric.IsCartier( D ) == true
    @test JToric.IsPrincipal( D ) == true
    @test JToric.IsPrimedivisor( D ) == false
    @test JToric.IsBasepointFree( D ) == true
    @test JToric.IsAmple( D ) == false
    @test JToric.IsVeryAmple( D ) == "fail"
    @test JToric.IsNumericallyEffective( D ) == true

    # compute properties of toric divisors on the example of a non-trivial divisor
    D2 = DivisorOfCharacter( [ 1,2 ], H5 )
    @test JToric.IsCartier( D2 ) == true
    @test JToric.IsPrincipal( D2 ) == true
    @test JToric.IsPrimedivisor( D2 ) == false
    @test JToric.IsBasepointFree( D2 ) == true
    @test JToric.IsAmple( D2 ) == false
    @test JToric.IsVeryAmple( D2 ) == "fail"
    @test JToric.IsNumericallyEffective( D2 ) == true
    
    # compute properties of toric divisors on the example of another non-trivial divisor
    D3 = DivisorOfGivenClass( H5, [ 1,2 ] )
    @test JToric.IsCartier( D3 ) == true
    @test JToric.IsPrincipal( D3 ) == false
    @test JToric.IsPrimedivisor( D3 ) == false
    @test JToric.IsBasepointFree( D3 ) == true
    @test JToric.IsAmple( D3 ) == true
    @test JToric.IsVeryAmple( D3 ) == true
    @test JToric.IsNumericallyEffective( D3 ) == true
    
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
    #@assert IsPrimedivisor.(divisors)
    
end
