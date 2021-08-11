using JToric
using Test
using GAP

# load and export CapAndHomalg
import CapAndHomalg
export CapAndHomalg

@testset "JToric.jl" begin
    # Write your tests here.

    # load necessary gap packages
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    CapAndHomalg.LoadPackage( "JConvex" )
    CapAndHomalg.LoadPackage( "ToricV" )
    
    # Is P2 smooth?
    P2 = JToric.ProjectiveSpace( 2 )
    @test GAP.Globals.GAPToJulia( JToric.IsSmooth( P2 ) ) == true
    
    # Construct Hirzebruch surface
    Rays = [[-1,5],[0,1],[1,0],[0,-1]]
    Cones = [[1,2],[2,3],[3,4],[4,1]]
    H5 = JToric.Fan( Rays, Cones )
    
end
