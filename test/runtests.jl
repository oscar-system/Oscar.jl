using JToric
using Test

@testset "JToric.jl" begin
    # Write your tests here.
    
    # load GAP-Julia interface
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    
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
    
end
