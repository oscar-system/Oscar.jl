using GAP
using CapAndHomalg

include("ToricDivisors.jl")

######################
# 1: Attributes of ToricVarieties
######################

function AffineOpenCovering( v::JToricVariety )
    gap_cover = GAP.Globals.AffineOpenCovering( v.GapToricVariety )
    return [ JToricVariety( v ) for v in gap_cover ]
end
export AffineOpenCovering


struct JCoxRing
           GapCoxRing
end
export JCoxRing

function CoxRing( v::JToricVariety )
    
    gap_ring = GAP.Globals.CoxRing( v.GapToricVariety )
    return JCoxRing( gap_ring )
    
end
export CoxRing


function ListOfVariablesOfCoxRing(v::JToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoxRing(v.GapToricVariety)
    return Vector{String}(vars)
end
export ListOfVariablesOfCoxRing


struct JClassGroup
           GapClassGroup
end
export JClassGroup

function ClassGroup( v::JToricVariety )
    
    gap_class_group = GAP.Globals.ClassGroup( v.GapToricVariety )
    return JClassGroup( gap_class_group )
    
end
export ClassGroup


struct JTorusInvariantDivisorGroup
           GapTorusInvariantDivisorGroup
end
export JTorusInvariantDivisorGroup

function TorusInvariantDivisorGroup( v::JToricVariety )
    
    gap_TorusInvariantDivisorGroup = GAP.Globals.TorusInvariantDivisorGroup( v.GapToricVariety )
    return JTorusInvariantDivisorGroup( gap_TorusInvariantDivisorGroup )
    
end
export TorusInvariantDivisorGroup


struct JMapFromCharacterToPrincipalDivisor
           GapMapFromCharacterToPrincipalDivisor
end
export JMapFromCharacterToPrincipalDivisor

function MapFromCharacterToPrincipalDivisor( v::JToricVariety )
    
    gap_MapFromCharacterToPrincipalDivisor = GAP.Globals.MapFromCharacterToPrincipalDivisor( v.GapToricVariety )
    return JMapFromCharacterToPrincipalDivisor( gap_MapFromCharacterToPrincipalDivisor )
    
end
export MapFromCharacterToPrincipalDivisor


struct JMapFromWeilDivisorsToClassGroup
           GapMapFromWeilDivisorsToClassGroup
end
export JMapFromWeilDivisorsToClassGroup

function MapFromWeilDivisorsToClassGroup( v::JToricVariety )
    
    gap_MapFromWeilDivisorsToClassGroup = GAP.Globals.MapFromWeilDivisorsToClassGroup( v.GapToricVariety )
    return JMapFromWeilDivisorsToClassGroup( gap_MapFromWeilDivisorsToClassGroup )
    
end
export MapFromWeilDivisorsToClassGroup


function Dimension( v::JToricVariety )
    
    return UInt128( GAP.Globals.Dimension( v.GapToricVariety ) )
    
end
export Dimension


function DimensionOfTorusfactor( v::JToricVariety )
    
    return UInt128( GAP.Globals.DimensionOfTorusfactor( v.GapToricVariety ) )
    
end
export DimensionOfTorusfactor


struct JCoordinateRingOfTorus
           GapCoordinateRingOfTorus
end
export JCoordinateRingOfTorus

function CoordinateRingOfTorus( v::JToricVariety )
    
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapToricVariety )
    return JCoordinateRingOfTorus( gap_CoordinateRingOfTorus )
    
end
export CoordinateRingOfTorus


function ListOfVariablesOfCoordinateRingOfTorus(v::JToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoordinateRingOfTorus(v.GapToricVariety)
    return Vector{String}(vars)
end
export ListOfVariablesOfCoordinateRingOfTorus


function IsProductOf(v::JToricVariety)
    factors = CapAndHomalg.GAP.Globals.IsProductOf(v.GapToricVariety)
    return [ JToricVariety( f ) for f in factors ]
end
export IsProductOf


struct JCharacterLattice
           GapCharacterLattice
end
export JCharacterLattice

function CharacterLattice( v::JToricVariety )
    
    gap_CharacterLattice = GAP.Globals.CharacterLattice( v.GapToricVariety )
    return JCharacterLattice( gap_CharacterLattice )
    
end
export CharacterLattice


function TorusInvariantPrimeDivisors( v::JToricVariety )
    divisors = GAP.Globals.TorusInvariantPrimeDivisors( v.GapToricVariety )
    return [ JToricDivisor( d ) for d in divisors ]    
end
export TorusInvariantPrimeDivisors


struct JIrrelevantIdeal
           GapIrrelevantIdeal
end
export JIrrelevantIdeal

function IrrelevantIdeal( v::JToricVariety )
    
    gap_IrrelevantIdeal = GAP.Globals.IrrelevantIdeal( v.GapToricVariety )
    return JIrrelevantIdeal( gap_IrrelevantIdeal )
    
end
export IrrelevantIdeal


struct JSRIdeal
           GapSRIdeal
end
export JSRIdeal

function SRIdeal( v::JToricVariety )
    
    gap_SRIdeal = GAP.Globals.SRIdeal( v.GapToricVariety )
    return JSRIdeal( gap_SRIdeal )
    
end
export SRIdeal


struct JMorphismFromCoxVariety
           GapMorphismFromCoxVariety
end
export JMorphismFromCoxVariety

function MorphismFromCoxVariety( v::JToricVariety )
    
    gap_MorphismFromCoxVariety = GAP.Globals.MorphismFromCoxVariety( v.GapToricVariety )
    return JMorphismFromCoxVariety( gap_MorphismFromCoxVariety )
    
end
export MorphismFromCoxVariety


function CoxVariety( v::JToricVariety )
    
    gap_CoxVariety = GAP.Globals.CoxVariety( v.GapToricVariety )
    return JToricVariety( gap_CoxVariety )
    
end
export CoxVariety


struct JuliaFan
           gap_fan
           rays::Vector{Vector{Int64}}
           cones::Vector{Vector{Int64}}
end
export JuliaFan

function FanOfVariety( v::JToricVariety )

    # collect data
    gap_fan = GAP.Globals.FanOfVariety( v.GapToricVariety )
    rays = Vector{Vector{Int64}}( GAP.Globals.RayGenerators( gap_fan ) )
    cones = Vector{Vector{Int64}}( GAP.Globals.RaysInMaximalCones( gap_fan ) )
    cones = [ findall( x -> x == 1, c ) for c in cones ]
    
    # return the fan
    return JuliaFan( gap_fan, rays, cones )
    
end
export FanOfVariety


function Fan( v::JToricVariety )
    
    return FanOfVariety( v )
    
end
export Fan


struct JCartierTorusInvariantDivisorGroup
           GapCartierTorusInvariantDivisorGroup
end
export JCartierTorusInvariantDivisorGroup

function CartierTorusInvariantDivisorGroup( v::JToricVariety )
    
    gap_CartierTorusInvariantDivisorGroup = GAP.Globals.CartierTorusInvariantDivisorGroup( v.GapToricVariety )
    return JCartierTorusInvariantDivisorGroup( gap_CartierTorusInvariantDivisorGroup )
    
end
export CartierTorusInvariantDivisorGroup


struct JPicardGroup
           GapPicardGroup
end
export JPicardGroup

function PicardGroup( v::JToricVariety )
    
    gap_PicardGroup = GAP.Globals.PicardGroup( v.GapToricVariety )
    return JPicardGroup( gap_PicardGroup )
    
end
export PicardGroup


function NameOfVariety( v::JToricVariety )
    
    if ! ( Bool( GAP.Globals.HasNameOfVariety( v.GapToricVariety ) ) )
            return "No name set for this variety"
    end
    
    return String( GAP.Globals.NameOfVariety( v.GapToricVariety ) )
    
end
export NameOfVariety


function SetNameOfVariety( v::JToricVariety, s::String )
    
    GAP.Globals.SetNameOfVariety( v.GapToricVariety, GapObj( s ) )
    return true
    
end
export SetNameOfVariety


struct JZariskiCotangentSheaf
           GapZariskiCotangentSheaf
end
export JZariskiCotangentSheaf

function ZariskiCotangentSheaf( v::JToricVariety )
    
    gap_ZariskiCotangentSheaf = GAP.Globals.ZariskiCotangentSheaf( v.GapToricVariety )
    return JZariskiCotangentSheaf( gap_ZariskiCotangentSheaf )
    
end
export ZariskiCotangentSheaf


struct JCotangentSheaf
           GapCotangentSheaf
end
export JCotangentSheaf

function CotangentSheaf( v::JToricVariety )
    
    gap_CotangentSheaf = GAP.Globals.CotangentSheaf( v.GapToricVariety )
    return JCotangentSheaf( gap_CotangentSheaf )
    
end
export CotangentSheaf


function EulerCharacteristic( v::JToricVariety )
    
    return Int128( GAP.Globals.EulerCharacteristic( v.GapToricVariety ) )
    
end
export EulerCharacteristic


#struct JUnderlyingSheaf
#          GapUnderlyingSheaf
#end
#export JUnderlyingSheaf

#function UnderlyingSheaf( v::JToricVariety )
    
#   gap_Underlying = GAP.Globals.UnderlyingSheaf( v.GapToricVariety )
#    return JUnderlyingSheaf( gap_UnderlyingSheaf )
    
#end
#export UnderlyingSheaf


######################
# 2: Methods of ToricVarieties
######################

function CoordinateRingOfTorus( v::JToricVariety, names::Vector{String} )
    
    gap_names = [ GapObj( names[ i ] ) for i in 1 : size(names)[1] ]
    gap_names = GapObj( gap_names )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapToricVariety, gap_names )
    return JCoordinateRingOfTorus( gap_CoordinateRingOfTorus )
    
end
export CoordinateRingOfTorus


function CoxRing( v::JToricVariety, names::String )
    
    gap_names = GapObj( names )
    gap_CoxRing = GAP.Globals.CoxRing( v.GapToricVariety, gap_names )
    return JCoxRing( gap_CoxRing )
    
end
export CoordinateRingOfTorus


function Base.:*( v::JToricVariety, w::JToricVariety )
    
    gap_ToricVariety = v.GapToricVariety * w.GapToricVariety
    return JToricVariety( gap_ToricVariety )
    
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
