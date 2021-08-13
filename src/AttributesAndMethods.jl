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


function ListOfVariablesOfCoxRing( v::JToricVariety )
    
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
