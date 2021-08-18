include("ToricDivisors.jl")

######################
# 1: Attributes of ToricVarieties
######################

function affine_open_covering( v::toric_variety )
    gap_cover = GAP.Globals.AffineOpenCovering( v.GapToricVariety )
    return [ toric_variety( v ) for v in gap_cover ]
end
export affine_open_covering


struct cox_ring
           GapCoxRing::GapObj
end
export cox_ring

function cox_ring( v::toric_variety )
    gap_ring = GAP.Globals.CoxRing( v.GapToricVariety )
    return cox_ring( gap_ring )
end
export cox_ring


function list_of_variables_of_cox_ring(v::toric_variety)
    vars = GAP.Globals.ListOfVariablesOfCoxRing(v.GapToricVariety)
    return Vector{String}(vars)
end
export list_of_variables_of_cox_ring


struct class_group
           GapClassGroup::GapObj
end
export class_group

function class_group( v::toric_variety )
    gap_class_group = GAP.Globals.ClassGroup( v.GapToricVariety )
    return class_group( gap_class_group )
end
export class_group


struct torus_invariant_divisor_group
           GapTorusInvariantDivisorGroup::GapObj
end
export torus_invariant_divisor_group

function torus_invariant_divisor_group( v::toric_variety )
    gap_TorusInvariantDivisorGroup = GAP.Globals.TorusInvariantDivisorGroup( v.GapToricVariety )
    return torus_invariant_divisor_group( gap_TorusInvariantDivisorGroup )
end
export torus_invariant_divisor_group


struct map_from_character_to_principal_divisor
           GapMapFromCharacterToPrincipalDivisor::GapObj
end
export map_from_character_to_principal_divisor

function map_from_character_to_principal_divisor( v::toric_variety )
    gap_MapFromCharacterToPrincipalDivisor = GAP.Globals.MapFromCharacterToPrincipalDivisor( v.GapToricVariety )
    return map_from_character_to_principal_divisor( gap_MapFromCharacterToPrincipalDivisor )
end
export map_from_character_to_principal_divisor


struct map_from_weil_divisors_to_class_group
           GapMapFromWeilDivisorsToClassGroup::GapObj
end
export map_from_weil_divisors_to_class_group

function map_from_weil_divisors_to_class_group( v::toric_variety )
    gap_MapFromWeilDivisorsToClassGroup = GAP.Globals.MapFromWeilDivisorsToClassGroup( v.GapToricVariety )
    return map_from_weil_divisors_to_class_group( gap_MapFromWeilDivisorsToClassGroup )
end
export map_from_weil_divisors_to_class_group


function dimension( v::toric_variety )
    return GAP.Globals.Dimension(v.GapToricVariety)::Int
end
export dimension


function dimension_of_torusfactor( v::toric_variety )
    return GAP.Globals.DimensionOfTorusfactor( v.GapToricVariety )::Int
end
export dimension_of_torusfactor


struct coordinate_ring_of_torus
           GapCoordinateRingOfTorus::GapObj
end
export coordinate_ring_of_torus

function coordinate_ring_of_torus( v::toric_variety )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapToricVariety )
    return coordinate_ring_of_torus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


function list_of_variables_of_coordinate_ring_of_torus(v::toric_variety)
    vars = GAP.Globals.ListOfVariablesOfCoordinateRingOfTorus(v.GapToricVariety)
    return Vector{String}(vars)
end
export list_of_variables_of_coordinate_ring_of_torus


function is_product_of(v::toric_variety)
    factors = GAP.Globals.IsProductOf(v.GapToricVariety)
    return [ toric_variety( f ) for f in factors ]
end
export is_product_of


struct character_lattice
           GapCharacterLattice::GapObj
end
export character_lattice

function character_lattice( v::toric_variety )
    gap_CharacterLattice = GAP.Globals.CharacterLattice( v.GapToricVariety )
    return character_lattice( gap_CharacterLattice )
end
export character_lattice


function torus_invariant_prime_divisors( v::toric_variety )
    divisors = GAP.Globals.TorusInvariantPrimeDivisors( v.GapToricVariety )
    return [ toric_divisor( d ) for d in divisors ]    
end
export torus_invariant_prime_divisors


struct irrelevant_ideal
           GapIrrelevantIdeal::GapObj
end
export irrelevant_ideal

function irrelevant_ideal( v::toric_variety )
    gap_IrrelevantIdeal = GAP.Globals.IrrelevantIdeal( v.GapToricVariety )
    return irrelevant_ideal( gap_IrrelevantIdeal )
end
export irrelevant_ideal


struct stanley_reisner_ideal
           GapSRIdeal::GapObj
end
export stanley_reisner_ideal

function stanley_reisner_ideal( v::toric_variety )
    gap_SRIdeal = GAP.Globals.SRIdeal( v.GapToricVariety )
    return stanley_reisner_ideal( gap_SRIdeal )
end
export stanley_reisner_ideal


struct morphism_from_cox_variety
           GapMorphismFromCoxVariety::GapObj
end
export morphism_from_cox_variety

function morphism_from_cox_variety( v::toric_variety )
    gap_MorphismFromCoxVariety = GAP.Globals.MorphismFromCoxVariety( v.GapToricVariety )
    return morphism_from_cox_variety( gap_MorphismFromCoxVariety )
end
export morphism_from_cox_variety


function cox_variety( v::toric_variety )
    gap_CoxVariety = GAP.Globals.CoxVariety( v.GapToricVariety )
    return toric_variety( gap_CoxVariety )
end
export cox_variety


struct fan
           gap_fan::GapObj
           rays::Vector{Vector{Int}}
           cones::Vector{Vector{Int}}
end
export fan

function fan_of_variety( v::toric_variety )
    # collect data
    gap_fan = GAP.Globals.FanOfVariety( v.GapToricVariety )
    rays = Vector{Vector{Int}}( GAP.Globals.RayGenerators( gap_fan ) )
    cones = Vector{Vector{Int}}( GAP.Globals.RaysInMaximalCones( gap_fan ) )
    cones = [ findall( x -> x == 1, c ) for c in cones ]
    
    # return the fan
    return fan( gap_fan, rays, cones )
end
export fan_of_variety


function fan( v::toric_variety )
    return fan_of_variety( v )
end
export fan


struct cartier_torus_invariant_divisor_group
           GapCartierTorusInvariantDivisorGroup::GapObj
end
export cartier_torus_invariant_divisor_group

function cartier_torus_invariant_divisor_group( v::toric_variety )
    gap_CartierTorusInvariantDivisorGroup = GAP.Globals.CartierTorusInvariantDivisorGroup( v.GapToricVariety )
    return cartier_torus_invariant_divisor_group( gap_CartierTorusInvariantDivisorGroup )
end
export cartier_torus_invariant_divisor_group


struct picard_group
           GapPicardGroup::GapObj
end
export picard_group

function picard_group( v::toric_variety )
    gap_PicardGroup = GAP.Globals.PicardGroup( v.GapToricVariety )
    return picard_group( gap_PicardGroup )
end
export picard_group


function name_of_variety( v::toric_variety )
    if ! ( Bool( GAP.Globals.HasNameOfVariety( v.GapToricVariety ) ) )
            return "No name set for this variety"
    end
    
    return String( GAP.Globals.NameOfVariety( v.GapToricVariety ) )
end
export name_of_variety


function set_name_of_variety( v::toric_variety, s::String )
    GAP.Globals.SetNameOfVariety( v.GapToricVariety, GapObj( s ) )
    return true
end
export set_name_of_variety


struct zariski_cotangent_sheaf
           GapZariskiCotangentSheaf::GapObj
end
export zariski_cotangent_sheaf

function zariski_cotangent_sheaf( v::toric_variety )
    gap_ZariskiCotangentSheaf = GAP.Globals.ZariskiCotangentSheaf( v.GapToricVariety )
    return zariski_cotangent_sheaf( gap_ZariskiCotangentSheaf )
end
export zariski_cotangent_sheaf


struct cotangent_sheaf
           GapCotangentSheaf::GapObj
end
export cotangent_sheaf

function cotangent_sheaf( v::toric_variety )
    gap_CotangentSheaf = GAP.Globals.CotangentSheaf( v.GapToricVariety )
    return cotangent_sheaf( gap_CotangentSheaf )
end
export cotangent_sheaf


function euler_characteristic( v::toric_variety )
    return GAP.Globals.EulerCharacteristic( v.GapToricVariety )::Int
end
export euler_characteristic


#struct underlying_sheaf
#          GapUnderlyingSheaf::GapObj
#end
#export underlying_sheaf

#function underlying_sheaf( v::toric_variety )
#   gap_Underlying = GAP.Globals.UnderlyingSheaf( v.GapToricVariety )
#    return underlying_sheaf( gap_UnderlyingSheaf )
#end
#export underlying_sheaf


######################
# 2: Methods of ToricVarieties
######################

function coordinate_ring_of_torus( v::toric_variety, names::Vector{String} )
    gap_names = [ GapObj( names[ i ] ) for i in 1 : size(names)[1] ]
    gap_names = GapObj( gap_names )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapToricVariety, gap_names )
    return coordinate_ring_of_torus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


function cox_ring( v::toric_variety, names::String )
    gap_names = GapObj( names )
    gap_CoxRing = GAP.Globals.CoxRing( v.GapToricVariety, gap_names )
    return cox_ring( gap_CoxRing )
end
export cox_ring


function Base.:*( v::toric_variety, w::toric_variety )
    gap_ToricVariety = v.GapToricVariety * w.GapToricVariety
    return toric_variety( gap_ToricVariety )
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
