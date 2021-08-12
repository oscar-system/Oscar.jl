using GAP
using CapAndHomalg


######################
# 1: The Julia type for ToricDivisors
######################

struct JToricDivisor
           bar
           GapToricDivisor
end


######################
# 2: Generic constructors
######################

#function JToricDivisor( rays::Array{Array{Int64,1},1}, cones::Array{Array{Int64,1},1} )

    # construct the toric variety in GAP
    #gap_rays = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( rays )
    #gap_cones = CapAndHomalg.GAP.Globals.ConvertJuliaToGAP( cones )
    #fan = CapAndHomalg.GAP.Globals.Fan( gap_rays, gap_cones )
    #variety = CapAndHomalg.GAP.Globals.ToricVariety( fan )
    
    # wrap it into a struct and return
    #return JToricVariety( 1, variety )

#end
#export JToricVariety

#CreateDivisor( ConvertJuliaToGAP( [ 0,0,0,0 ] ),H5 )
#DivisorOfCharacter( ConvertJuliaToGAP([ 1,2 ]),H5 )
#DivisorOfGivenClass( IsToricVariety, IsList )
