include("ToricDivisors.jl")

######################
# 1: Attributes of ToricVarieties
######################

function affine_open_covering( v::NormalToricVariety )
    gap_cover = GAP.Globals.AffineOpenCovering( v.GapNTV )
    return [ NormalToricVariety( v ) for v in gap_cover ]
end
export affine_open_covering


struct cox_ring
           GapCoxRing::GapObj
end
export cox_ring

function cox_ring( v::NormalToricVariety )
    gap_ring = GAP.Globals.CoxRing( v.GapNTV )
    return cox_ring( gap_ring )
end
export cox_ring


function list_of_variables_of_cox_ring(v::NormalToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoxRing(v.GapNTV)
    return Vector{String}(vars)
end
export list_of_variables_of_cox_ring


struct class_group
           GapClassGroup::GapObj
end
export class_group

function class_group( v::NormalToricVariety )
    gap_class_group = GAP.Globals.ClassGroup( v.GapNTV )
    return class_group( gap_class_group )
end
export class_group


struct torus_invariant_divisor_group
           GapTorusInvariantDivisorGroup::GapObj
end
export torus_invariant_divisor_group

function torus_invariant_divisor_group( v::NormalToricVariety )
    gap_TorusInvariantDivisorGroup = GAP.Globals.TorusInvariantDivisorGroup( v.GapNTV )
    return torus_invariant_divisor_group( gap_TorusInvariantDivisorGroup )
end
export torus_invariant_divisor_group


struct map_from_character_to_principal_divisor
           GapMapFromCharacterToPrincipalDivisor::GapObj
end
export map_from_character_to_principal_divisor

function map_from_character_to_principal_divisor( v::NormalToricVariety )
    gap_MapFromCharacterToPrincipalDivisor = GAP.Globals.MapFromCharacterToPrincipalDivisor( v.GapNTV )
    return map_from_character_to_principal_divisor( gap_MapFromCharacterToPrincipalDivisor )
end
export map_from_character_to_principal_divisor


struct map_from_weil_divisors_to_class_group
           GapMapFromWeilDivisorsToClassGroup::GapObj
end
export map_from_weil_divisors_to_class_group

function map_from_weil_divisors_to_class_group( v::NormalToricVariety )
    gap_MapFromWeilDivisorsToClassGroup = GAP.Globals.MapFromWeilDivisorsToClassGroup( v.GapNTV )
    return map_from_weil_divisors_to_class_group( gap_MapFromWeilDivisorsToClassGroup )
end
export map_from_weil_divisors_to_class_group


function dimension( v::NormalToricVariety )
    return GAP.Globals.Dimension(v.GapNTV)::Int
end
export dimension


function dimension_of_torusfactor( v::NormalToricVariety )
    return GAP.Globals.DimensionOfTorusfactor( v.GapNTV )::Int
end
export dimension_of_torusfactor


struct coordinate_ring_of_torus
           GapCoordinateRingOfTorus::GapObj
end
export coordinate_ring_of_torus

function coordinate_ring_of_torus( v::NormalToricVariety )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapNTV )
    return coordinate_ring_of_torus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


function list_of_variables_of_coordinate_ring_of_torus(v::NormalToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoordinateRingOfTorus(v.GapNTV)
    return Vector{String}(vars)
end
export list_of_variables_of_coordinate_ring_of_torus


function is_product_of(v::NormalToricVariety)
    factors = GAP.Globals.IsProductOf(v.GapNTV)
    return [ NormalToricVariety( f ) for f in factors ]
end
export is_product_of


struct character_lattice
           GapCharacterLattice::GapObj
end
export character_lattice

function character_lattice( v::NormalToricVariety )
    gap_CharacterLattice = GAP.Globals.CharacterLattice( v.GapNTV )
    return character_lattice( gap_CharacterLattice )
end
export character_lattice


function torus_invariant_prime_divisors( v::NormalToricVariety )
    divisors = GAP.Globals.TorusInvariantPrimeDivisors( v.GapNTV )
    return [ toric_divisor( d ) for d in divisors ]    
end
export torus_invariant_prime_divisors


struct irrelevant_ideal
           GapIrrelevantIdeal::GapObj
end
export irrelevant_ideal

function irrelevant_ideal( v::NormalToricVariety )
    gap_IrrelevantIdeal = GAP.Globals.IrrelevantIdeal( v.GapNTV )
    return irrelevant_ideal( gap_IrrelevantIdeal )
end
export irrelevant_ideal


struct stanley_reisner_ideal
           GapSRIdeal::GapObj
end
export stanley_reisner_ideal

function stanley_reisner_ideal( v::NormalToricVariety )
    gap_SRIdeal = GAP.Globals.SRIdeal( v.GapNTV )
    return stanley_reisner_ideal( gap_SRIdeal )
end
export stanley_reisner_ideal


struct morphism_from_cox_variety
           GapMorphismFromCoxVariety::GapObj
end
export morphism_from_cox_variety

function morphism_from_cox_variety( v::NormalToricVariety )
    gap_MorphismFromCoxVariety = GAP.Globals.MorphismFromCoxVariety( v.GapNTV )
    return morphism_from_cox_variety( gap_MorphismFromCoxVariety )
end
export morphism_from_cox_variety


function cox_variety( v::NormalToricVariety )
    gap_CoxVariety = GAP.Globals.CoxVariety( v.GapNTV )
    return NormalToricVariety( gap_CoxVariety )
end
export cox_variety


struct fan
           gap_fan::GapObj
           rays::Vector{Vector{Int}}
           cones::Vector{Vector{Int}}
end
export fan

function fan_of_variety( v::NormalToricVariety )
    # collect data
    gap_fan = GAP.Globals.FanOfVariety( v.GapNTV )
    rays = Vector{Vector{Int}}( GAP.Globals.RayGenerators( gap_fan ) )
    cones = Vector{Vector{Int}}( GAP.Globals.RaysInMaximalCones( gap_fan ) )
    cones = [ findall( x -> x == 1, c ) for c in cones ]
    
    # return the fan
    return fan( gap_fan, rays, cones )
end
export fan_of_variety


function fan( v::NormalToricVariety )
    return fan_of_variety( v )
end
export fan


struct cartier_torus_invariant_divisor_group
           GapCartierTorusInvariantDivisorGroup::GapObj
end
export cartier_torus_invariant_divisor_group

function cartier_torus_invariant_divisor_group( v::NormalToricVariety )
    gap_CartierTorusInvariantDivisorGroup = GAP.Globals.CartierTorusInvariantDivisorGroup( v.GapNTV )
    return cartier_torus_invariant_divisor_group( gap_CartierTorusInvariantDivisorGroup )
end
export cartier_torus_invariant_divisor_group


struct picard_group
           GapPicardGroup::GapObj
end
export picard_group

function picard_group( v::NormalToricVariety )
    gap_PicardGroup = GAP.Globals.PicardGroup( v.GapNTV )
    return picard_group( gap_PicardGroup )
end
export picard_group


function name_of_variety( v::NormalToricVariety )
    if ! ( Bool( GAP.Globals.HasNameOfVariety( v.GapNTV ) ) )
            return "No name set for this variety"
    end
    
    return String( GAP.Globals.NameOfVariety( v.GapNTV ) )
end
export name_of_variety


function set_name_of_variety( v::NormalToricVariety, s::String )
    GAP.Globals.SetNameOfVariety( v.GapNTV, GapObj( s ) )
    return true
end
export set_name_of_variety


struct zariski_cotangent_sheaf
           GapZariskiCotangentSheaf::GapObj
end
export zariski_cotangent_sheaf

function zariski_cotangent_sheaf( v::NormalToricVariety )
    gap_ZariskiCotangentSheaf = GAP.Globals.ZariskiCotangentSheaf( v.GapNTV )
    return zariski_cotangent_sheaf( gap_ZariskiCotangentSheaf )
end
export zariski_cotangent_sheaf


struct cotangent_sheaf
           GapCotangentSheaf::GapObj
end
export cotangent_sheaf

function cotangent_sheaf( v::NormalToricVariety )
    gap_CotangentSheaf = GAP.Globals.CotangentSheaf( v.GapNTV )
    return cotangent_sheaf( gap_CotangentSheaf )
end
export cotangent_sheaf


function euler_characteristic( v::NormalToricVariety )
    return GAP.Globals.EulerCharacteristic( v.GapNTV )::Int
end
export euler_characteristic


#struct underlying_sheaf
#          GapUnderlyingSheaf::GapObj
#end
#export underlying_sheaf

#function underlying_sheaf( v::NormalToricVariety )
#   gap_Underlying = GAP.Globals.UnderlyingSheaf( v.GapNTV )
#    return underlying_sheaf( gap_UnderlyingSheaf )
#end
#export underlying_sheaf


######################
# 2: Methods of ToricVarieties
######################

function coordinate_ring_of_torus( v::NormalToricVariety, names::Vector{String} )
    gap_names = [ GapObj( names[ i ] ) for i in 1 : size(names)[1] ]
    gap_names = GapObj( gap_names )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapNTV, gap_names )
    return coordinate_ring_of_torus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


function cox_ring( v::NormalToricVariety, names::String )
    gap_names = GapObj( names )
    gap_CoxRing = GAP.Globals.CoxRing( v.GapNTV, gap_names )
    return cox_ring( gap_CoxRing )
end
export cox_ring


function Base.:*( v::NormalToricVariety, w::NormalToricVariety )
    gap_NormalToricVariety = v.GapNTV * w.GapNTV
    return NormalToricVariety( gap_NormalToricVariety )
end
export *


struct character_to_rational_function
           gap_CharacterToRationalFunction::GapObj
end
export character_to_rational_function

function character_to_rational_function( l::Vector{Int}, v::NormalToricVariety )
    gap_CharacterToRationalFunction = GAP.Globals.CharacterToRationalFunction( GapObj( l ), v.GapNTV )
    return character_to_rational_function( gap_CharacterToRationalFunction )
end
export character_to_rational_function


function weil_divisors_of_variety( v::NormalToricVariety )
    gap_divisors = GAP.Globals.WeilDivisorsOfVariety( v.GapNTV )
    return [ toric_divisor( d ) for d in gap_divisors ]
end
export weil_divisors_of_variety


function factors( v::NormalToricVariety )
    gap_factors = GAP.Globals.Factors( v.GapNTV )
    return [ NormalToricVariety( f ) for f in gap_factors ]
end
export factors


function blowup_on_ith_minimal_torus_orbit( v::NormalToricVariety, i::Int )
    gap_blowup_variety = GAP.Globals.BlowUpOnIthMinimalTorusOrbit( v.GapNTV, GapObj( i ) )
    return NormalToricVariety( gap_blowup_variety )
end
export blowup_on_ith_minimal_torus_orbit

#= 

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
                  [ IsNormalToricVariety, IsInt ] );

#! @Description
#!  
#! @Arguments vari, p
DeclareOperation( "NrOfqRationalPoints",
                  [ IsNormalToricVariety, IsInt ] );

=#
