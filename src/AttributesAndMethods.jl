using GAP
using CapAndHomalg

include("ToricDivisors.jl")

######################
# 1: Attributes of ToricVarieties
######################

function AffineOpenCovering( v::JToricVariety )
    
    gap_cover = CapAndHomalg.GAP.Globals.AffineOpenCovering( v.GapToricVariety )
    len = GAP.Globals.GAPToJulia( GAP.Globals.Length( gap_cover ) )
    j_cover = [ JToricVariety( 1, gap_cover[ i ] ) for i in 1:len ]
    return j_cover
    
end
export AffineOpenCovering


struct JCoxRing
           bar
           GapCoxRing
end
export JCoxRing

function CoxRing( v::JToricVariety )
    
    gap_ring = CapAndHomalg.GAP.Globals.CoxRing( v.GapToricVariety )
    return JCoxRing( 1, gap_ring )
    
end
export CoxRing


function ListOfVariablesOfCoxRing( v::JToricVariety,  )
    
    gap_variables = CapAndHomalg.GAP.Globals.ListOfVariablesOfCoxRing( v.GapToricVariety )
    len = GAP.Globals.GAPToJulia( GAP.Globals.Length( gap_variables ) )
    julia_variables = [  GAP.Globals.GAPToJulia( gap_variables[ i ] ) for i in 1:len ]
    return julia_variables
    
end
export ListOfVariablesOfCoxRing


struct JClassGroup
           bar
           GapClassGroup
end
export JClassGroup

function ClassGroup( v::JToricVariety )
    
    gap_class_group = CapAndHomalg.GAP.Globals.ClassGroup( v.GapToricVariety )
    return JClassGroup( 1, gap_class_group )
    
end
export ClassGroup


struct JTorusInvariantDivisorGroup
           bar
           GapTorusInvariantDivisorGroup
end
export JTorusInvariantDivisorGroup

function TorusInvariantDivisorGroup( v::JToricVariety )
    
    gap_TorusInvariantDivisorGroup = CapAndHomalg.GAP.Globals.TorusInvariantDivisorGroup( v.GapToricVariety )
    return JTorusInvariantDivisorGroup( 1, gap_TorusInvariantDivisorGroup )
    
end
export TorusInvariantDivisorGroup


struct JMapFromCharacterToPrincipalDivisor
           bar
           GapMapFromCharacterToPrincipalDivisor
end
export JMapFromCharacterToPrincipalDivisor

function MapFromCharacterToPrincipalDivisor( v::JToricVariety )
    
    gap_MapFromCharacterToPrincipalDivisor = CapAndHomalg.GAP.Globals.MapFromCharacterToPrincipalDivisor( v.GapToricVariety )
    return JMapFromCharacterToPrincipalDivisor( 1, gap_MapFromCharacterToPrincipalDivisor )
    
end
export MapFromCharacterToPrincipalDivisor


struct JMapFromWeilDivisorsToClassGroup
           bar
           GapMapFromWeilDivisorsToClassGroup
end
export JMapFromWeilDivisorsToClassGroup

function MapFromWeilDivisorsToClassGroup( v::JToricVariety )
    
    gap_MapFromWeilDivisorsToClassGroup = CapAndHomalg.GAP.Globals.MapFromWeilDivisorsToClassGroup( v.GapToricVariety )
    return JMapFromWeilDivisorsToClassGroup( 1, gap_MapFromWeilDivisorsToClassGroup )
    
end
export MapFromWeilDivisorsToClassGroup


function Dimension( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.Dimension( v.GapToricVariety ) )
    
end
export Dimension


function DimensionOfTorusfactor( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.DimensionOfTorusfactor( v.GapToricVariety ) )
    
end
export DimensionOfTorusfactor


struct JCoordinateRingOfTorus
           bar
           GapCoordinateRingOfTorus
end
export JCoordinateRingOfTorus

function CoordinateRingOfTorus( v::JToricVariety )
    
    gap_CoordinateRingOfTorus = CapAndHomalg.GAP.Globals.CoordinateRingOfTorus( v.GapToricVariety )
    return JCoordinateRingOfTorus( 1, gap_CoordinateRingOfTorus )
    
end
export CoordinateRingOfTorus


function ListOfVariablesOfCoordinateRingOfTorus( v::JToricVariety )
    
    gap_variables = CapAndHomalg.GAP.Globals.ListOfVariablesOfCoordinateRingOfTorus( v.GapToricVariety )
    len = GAP.Globals.GAPToJulia( GAP.Globals.Length( gap_variables ) )
    julia_variables = [  GAP.Globals.GAPToJulia( gap_variables[ i ] ) for i in 1:len ]
    return julia_variables
    
end
export ListOfVariablesOfCoordinateRingOfTorus


function IsProductOf( v::JToricVariety )
    
    factors = CapAndHomalg.GAP.Globals.IsProductOf( v.GapToricVariety )
    len = GAP.Globals.GAPToJulia( GAP.Globals.Length( factors ) )
    j_factors = [ JToricVariety( 1, factors[ i ] ) for i in 1:len ]
    return j_factors
    
end
export IsProductOf


struct JCharacterLattice
           bar
           GapCharacterLattice
end
export JCharacterLattice

function CharacterLattice( v::JToricVariety )
    
    gap_CharacterLattice = CapAndHomalg.GAP.Globals.CharacterLattice( v.GapToricVariety )
    return JCharacterLattice( 1, gap_CharacterLattice )
    
end
export CharacterLattice


function TorusInvariantPrimeDivisors( v::JToricVariety )
    
    divisors = CapAndHomalg.GAP.Globals.TorusInvariantPrimeDivisors( v.GapToricVariety )
    len = GAP.Globals.GAPToJulia( GAP.Globals.Length( divisors ) )
    j_divisors = [ JToricDivisor( 1, divisors[ i ] ) for i in 1:len ]
    return j_divisors
    
end
export TorusInvariantPrimeDivisors


struct JIrrelevantIdeal
           bar
           GapIrrelevantIdeal
end
export JIrrelevantIdeal

function IrrelevantIdeal( v::JToricVariety )
    
    gap_IrrelevantIdeal = CapAndHomalg.GAP.Globals.IrrelevantIdeal( v.GapToricVariety )
    return JIrrelevantIdeal( 1, gap_IrrelevantIdeal )
    
end
export IrrelevantIdeal


struct JSRIdeal
           bar
           GapSRIdeal
end
export JSRIdeal

function SRIdeal( v::JToricVariety )
    
    gap_SRIdeal = CapAndHomalg.GAP.Globals.SRIdeal( v.GapToricVariety )
    return JSRIdeal( 1, gap_SRIdeal )
    
end
export SRIdeal


struct JMorphismFromCoxVariety
           bar
           GapMorphismFromCoxVariety
end
export JMorphismFromCoxVariety

function MorphismFromCoxVariety( v::JToricVariety )
    
    gap_MorphismFromCoxVariety = CapAndHomalg.GAP.Globals.MorphismFromCoxVariety( v.GapToricVariety )
    return JMorphismFromCoxVariety( 1, gap_MorphismFromCoxVariety )
    
end
export MorphismFromCoxVariety


function CoxVariety( v::JToricVariety )
    
    gap_CoxVariety = CapAndHomalg.GAP.Globals.CoxVariety( v.GapToricVariety )
    return JToricVariety( 1, gap_CoxVariety )
    
end
export CoxVariety


struct JuliaFan
           bar
           gap_fan
           rays::Vector{Vector{Int64}}
           cones::Vector{Vector{Int64}}
end
export JuliaFan

function FanOfVariety( v::JToricVariety )

    # collect data
    gap_fan = CapAndHomalg.GAP.Globals.FanOfVariety( v.GapToricVariety )
    rays = GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.RayGenerators( gap_fan ) )
    cones = GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.RaysInMaximalCones( gap_fan ) )
    cones = [ findall( x -> x == 1, cones[i] ) for i in 1 : size( cones )[ 1 ] ]

    # return the fan
    return JuliaFan( 1, gap_fan, rays, cones )
    
end
export FanOfVariety


function Fan( v::JToricVariety )
    
    return FanOfVariety( v )
    
end
export Fan


struct JCartierTorusInvariantDivisorGroup
           bar
           GapCartierTorusInvariantDivisorGroup
end
export JCartierTorusInvariantDivisorGroup

function CartierTorusInvariantDivisorGroup( v::JToricVariety )
    
    gap_CartierTorusInvariantDivisorGroup = CapAndHomalg.GAP.Globals.CartierTorusInvariantDivisorGroup( v.GapToricVariety )
    return JCartierTorusInvariantDivisorGroup( 1, gap_CartierTorusInvariantDivisorGroup )
    
end
export CartierTorusInvariantDivisorGroup


struct JPicardGroup
           bar
           GapPicardGroup
end
export JPicardGroup

function PicardGroup( v::JToricVariety )
    
    gap_PicardGroup = CapAndHomalg.GAP.Globals.PicardGroup( v.GapToricVariety )
    return JPicardGroup( 1, gap_PicardGroup )
    
end
export PicardGroup


function NameOfVariety( v::JToricVariety )
    
    if ! ( GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.HasNameOfVariety( v.GapToricVariety ) ) )
            return "No name set for this variety"
    end
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.NameOfVariety( v.GapToricVariety ) )
    
end
export NameOfVariety


function SetNameOfVariety( v::JToricVariety, s::String )
    
    CapAndHomalg.GAP.Globals.SetNameOfVariety( v.GapToricVariety, GAP.Globals.JuliaToGAP( GAP.Globals.IsString, s ) )
    return true
    
end
export SetNameOfVariety


struct JZariskiCotangentSheaf
           bar
           GapZariskiCotangentSheaf
end
export JZariskiCotangentSheaf

function ZariskiCotangentSheaf( v::JToricVariety )
    
    gap_ZariskiCotangentSheaf = CapAndHomalg.GAP.Globals.ZariskiCotangentSheaf( v.GapToricVariety )
    return JZariskiCotangentSheaf( 1, gap_ZariskiCotangentSheaf )
    
end
export ZariskiCotangentSheaf


struct JCotangentSheaf
           bar
           GapCotangentSheaf
end
export JCotangentSheaf

function CotangentSheaf( v::JToricVariety )
    
    gap_CotangentSheaf = CapAndHomalg.GAP.Globals.CotangentSheaf( v.GapToricVariety )
    return JCotangentSheaf( 1, gap_CotangentSheaf )
    
end
export CotangentSheaf


function EulerCharacteristic( v::JToricVariety )
    
    return GAP.Globals.GAPToJulia( CapAndHomalg.GAP.Globals.EulerCharacteristic( v.GapToricVariety ) )
    
end
export EulerCharacteristic


#struct JUnderlyingSheaf
#          bar
#          GapUnderlyingSheaf
#end
#export JUnderlyingSheaf

#function UnderlyingSheaf( v::JToricVariety )
    
#   gap_Underlying = CapAndHomalg.GAP.Globals.UnderlyingSheaf( v.GapToricVariety )
#    return JUnderlyingSheaf( 1, gap_UnderlyingSheaf )
    
#end
#export UnderlyingSheaf


######################
# 2: Methods of ToricVarieties
######################

function CoordinateRingOfTorus( v::JToricVariety, names::Vector{String} )
    
    gap_names = [ GAP.Globals.JuliaToGAP( GAP.Globals.IsString, names[ i ] ) for i in 1 : size(names)[1] ]
    gap_names = GAP.Globals.JuliaToGAP( GAP.Globals.IsList, gap_names )
    gap_CoordinateRingOfTorus = CapAndHomalg.GAP.Globals.CoordinateRingOfTorus( v.GapToricVariety, gap_names )
    return JCoordinateRingOfTorus( 1, gap_CoordinateRingOfTorus )
    
end
export CoordinateRingOfTorus


function CoxRing( v::JToricVariety, names::String )
    
    gap_names = GAP.Globals.JuliaToGAP( GAP.Globals.IsString, names )
    gap_CoxRing = CapAndHomalg.GAP.Globals.CoxRing( v.GapToricVariety, gap_names )
    return JCoxRing( 1, gap_CoxRing )
    
end
export CoordinateRingOfTorus


function Base.:*( v::JToricVariety, w::JToricVariety )
    
    gap_ToricVariety = v.GapToricVariety * w.GapToricVariety
    return JToricVariety( 1, gap_ToricVariety )
    
end
export *

#= 

#! @Description
#!  Computes the rational function corresponding to the character grid element <A>elem</A> or to the list of 
#!  integers <A>elem</A>. This computation needs to know the coordinate ring of the torus of the variety <A>vari</A>. By
#!  default this ring is introduced with variables <A>x1</A> to <A>xn</A> where <A>n</A> is the dimension of the variety. If
#!  different variables should be used, then <A>CoordinateRingOfTorus</A> has to be set accordingly before calling this method.
#! @Returns a homalg element
#! @Arguments elem, vari
DeclareOperation( "CharacterToRationalFunction",
                  [ IsList, IsToricVariety ] );

#! @Description
#!  Returns a list of the currently defined Divisors of the toric variety.
#! @Returns a list
#! @Arguments vari
DeclareOperation( "WeilDivisorsOfVariety",
                  [ IsToricVariety ] );

#! @Description
#!  
#! @Arguments vari
DeclareOperation( "Factors",
                  [ IsToricVariety ] );

#! @Description
#!  
#! @Arguments vari, p
DeclareOperation( "BlowUpOnIthMinimalTorusOrbit",
                  [ IsToricVariety, IsInt ] );

#! @Description
#!
DeclareGlobalFunction( "ZariskiCotangentSheafViaEulerSequence" );

#! @Description
#!
DeclareGlobalFunction( "ZariskiCotangentSheafViaPoincareResidueMap" );

#! @Description
#!  
#! @Arguments vari, p
DeclareOperation( "ithBettiNumber",
                  [ IsToricVariety, IsInt ] );

#! @Description
#!  
#! @Arguments vari, p
DeclareOperation( "NrOfqRationalPoints",
                  [ IsToricVariety, IsInt ] );

=#
