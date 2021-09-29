######################
# 1: Attributes of ToricVarieties
######################

"""
    affine_open_covering( v::AbstractNormalToricVariety )

Computes an affine open cover of the normal toric variety `v`, i.e. returns a list of affine toric varieties.
"""
function affine_open_covering( v::AbstractNormalToricVariety )
    gap_cover = GAP.Globals.AffineOpenCovering( v.GapNTV )
    return [ AbstractNormalToricVariety( v ) for v in gap_cover ]
end
export affine_open_covering


struct CoxRing
           GapCoxRing::GapObj
end
export CoxRing


"""
    cox_ring( v::AbstractNormalToricVariety )

Computes the Cox ring of the normal toric variety `v`.
"""
function cox_ring( v::AbstractNormalToricVariety )
    gap_ring = GAP.Globals.CoxRing( v.GapNTV )
    return CoxRing( gap_ring )
end
export cox_ring


"""
    list_of_variables_of_cox_ring( v::AbstractNormalToricVariety )

Computes the list of homogeneous variables of the Cox ring of the normal toric variety `v`.
"""
function list_of_variables_of_cox_ring(v::AbstractNormalToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoxRing(v.GapNTV)
    return Vector{String}(vars)
end
export list_of_variables_of_cox_ring


struct ClassGroup
           GapClassGroup::GapObj
end
export ClassGroup


"""
    class_group( v::AbstractNormalToricVariety )

Computes the class group of the normal toric variety `v`.
"""
function class_group( v::AbstractNormalToricVariety )
    gap_class_group = GAP.Globals.ClassGroup( v.GapNTV )
    return ClassGroup( gap_class_group )
end
export class_group


struct TorusInvariantDivisorGroup
           GapTorusInvariantDivisorGroup::GapObj
end
export TorusInvariantDivisorGroup


"""
    torus_invariant_divisor_group( v::AbstractNormalToricVariety )

Computes the torus invariant divisor class group of the normal toric variety `v`.
"""
function torus_invariant_divisor_group( v::AbstractNormalToricVariety )
    gap_TorusInvariantDivisorGroup = GAP.Globals.TorusInvariantDivisorGroup( v.GapNTV )
    return TorusInvariantDivisorGroup( gap_TorusInvariantDivisorGroup )
end
export torus_invariant_divisor_group


struct MapFromCharacterToPrincipalDivisor
           GapMapFromCharacterToPrincipalDivisor::GapObj
end
export MapFromCharacterToPrincipalDivisor


"""
    map_from_character_to_principal_divisor( v::AbstractNormalToricVariety )

Computes the map from the character lattice to the principal divisors of the normal toric variety `v`.
"""
function map_from_character_to_principal_divisor( v::AbstractNormalToricVariety )
    gap_MapFromCharacterToPrincipalDivisor = GAP.Globals.MapFromCharacterToPrincipalDivisor( v.GapNTV )
    return MapFromCharacterToPrincipalDivisor( gap_MapFromCharacterToPrincipalDivisor )
end
export map_from_character_to_principal_divisor


struct MapFromWeilDivisorsToClassGroup
           GapMapFromWeilDivisorsToClassGroup::GapObj
end
export MapFromWeilDivisorsToClassGroup


"""
    map_from_weil_divisors_to_class_group( v::AbstractNormalToricVariety )

Computes the map from the Weil divisors to the Class group of the normal toric variety `v`.
"""
function map_from_weil_divisors_to_class_group( v::AbstractNormalToricVariety )
    gap_MapFromWeilDivisorsToClassGroup = GAP.Globals.MapFromWeilDivisorsToClassGroup( v.GapNTV )
    return MapFromWeilDivisorsToClassGroup( gap_MapFromWeilDivisorsToClassGroup )
end
export map_from_weil_divisors_to_class_group


"""
    dim( v::AbstractNormalToricVariety )

Computes the dimension of the normal toric variety `v`.
"""
function dim( v::AbstractNormalToricVariety )
    return GAP.Globals.Dimension(v.GapNTV)::Int
end
export dim


"""
    dim_of_torusfactor( v::AbstractNormalToricVariety )

Computes the dimension of the torus factor of the normal toric variety `v`.
"""
function dim_of_torusfactor( v::AbstractNormalToricVariety )
    return GAP.Globals.DimensionOfTorusfactor( v.GapNTV )::Int
end
export dim_of_torusfactor


struct CoordinateRingOfTorus
           GapCoordinateRingOfTorus::GapObj
end
export CoordinateRingOfTorus


"""
    coordinate_ring_of_torus( v::AbstractNormalToricVariety )

Computes the coordinate ring of the torus of the normal toric variety `v`.
"""
function coordinate_ring_of_torus( v::AbstractNormalToricVariety )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapNTV )
    return CoordinateRingOfTorus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


"""
    list_of_variables_of_coordinate_ring_of_torus( v::AbstractNormalToricVariety )

Computes the list of homogeneous coordinates of the coordinate ring of the torus of the normal toric variety `v`.
"""
function list_of_variables_of_coordinate_ring_of_torus(v::AbstractNormalToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoordinateRingOfTorus(v.GapNTV)
    return Vector{String}(vars)
end
export list_of_variables_of_coordinate_ring_of_torus


"""
    is_product_of( v::AbstractNormalToricVariety )

Identifies the factors from which the normal toric variety `v` has been constructed in GAP.
"""
function is_product_of(v::AbstractNormalToricVariety)
    factors = GAP.Globals.IsProductOf(v.GapNTV)
    return [ NormalToricVariety( f ) for f in factors ]
end
export is_product_of


"""
    factors( v::AbstractNormalToricVariety )

Identifies the factors from which the normal toric variety `v` has been constructed in GAP.
"""
function factors( v::AbstractNormalToricVariety )
    gap_factors = GAP.Globals.Factors( v.GapNTV )
    return [ NormalToricVariety( f ) for f in gap_factors ]
end
export factors


struct CharacterLattice
           GapCharacterLattice::GapObj
end
export CharacterLattice


"""
    character_lattice( v::AbstractNormalToricVariety )

Computes the character lattice of the normal toric variety `v`.
"""
function character_lattice( v::AbstractNormalToricVariety )
    gap_CharacterLattice = GAP.Globals.CharacterLattice( v.GapNTV )
    return CharacterLattice( gap_CharacterLattice )
end
export character_lattice


"""
    torus_invariant_prime_divisors( v::AbstractNormalToricVariety )

Computes the torus invariant prime divisors of the normal toric variety `v`.
"""
function torus_invariant_prime_divisors( v::AbstractNormalToricVariety )
    divisors = GAP.Globals.TorusInvariantPrimeDivisors( v.GapNTV )
    return [ ToricDivisor( d ) for d in divisors ]    
end
export torus_invariant_prime_divisors


struct IrrelevantIdeal
           GapIrrelevantIdeal::GapObj
end
export IrrelevantIdeal


"""
    irrelevant_ideal( v::AbstractNormalToricVariety )

Computes the irrelevant ideal of the normal toric variety `v`.
"""
function irrelevant_ideal( v::AbstractNormalToricVariety )
    gap_IrrelevantIdeal = GAP.Globals.IrrelevantIdeal( v.GapNTV )
    return IrrelevantIdeal( gap_IrrelevantIdeal )
end
export irrelevant_ideal


struct StanleyReisnerIdeal
           GapSRIdeal::GapObj
end
export StanleyReisnerIdeal


"""
    stanley_reisner_ideal( v::AbstractNormalToricVariety )

Computes the Stanley-Reisner ideal of the normal toric variety `v`.
"""
function stanley_reisner_ideal( v::AbstractNormalToricVariety )
    gap_SRIdeal = GAP.Globals.SRIdeal( v.GapNTV )
    return StanleyReisnerIdeal( gap_SRIdeal )
end
export stanley_reisner_ideal


struct MorphismFromCoxVariety
           GapMorphismFromCoxVariety::GapObj
end
export MorphismFromCoxVariety


"""
    morphism_from_cox_variety( v::AbstractNormalToricVariety )

Computes the morphism from the Cox variety of the normal toric variety `v`.
"""
function morphism_from_cox_variety( v::AbstractNormalToricVariety )
    gap_MorphismFromCoxVariety = GAP.Globals.MorphismFromCoxVariety( v.GapNTV )
    return MorphismFromCoxVariety( gap_MorphismFromCoxVariety )
end
export morphism_from_cox_variety


"""
    cox_variety( v::AbstractNormalToricVariety )

Computes the Cox variety of the normal toric variety `v`.
"""
function cox_variety( v::AbstractNormalToricVariety )
    gap_CoxVariety = GAP.Globals.CoxVariety( v.GapNTV )
    return NormalToricVariety( gap_CoxVariety )
end
export cox_variety


struct Fan
           gap_fan::GapObj
           rays::Vector{Vector{Int}}
           cones::Vector{Vector{Int}}
end
export Fan


"""
    fan_of_variety( v::AbstractNormalToricVariety )

Computes the fan of the normal toric variety `v`.
"""
function fan_of_variety( v::AbstractNormalToricVariety )
    # collect data
    gap_fan = GAP.Globals.FanOfVariety( v.GapNTV )
    rays = Vector{Vector{Int}}( GAP.Globals.RayGenerators( gap_fan ) )
    cones = Vector{Vector{Int}}( GAP.Globals.RaysInMaximalCones( gap_fan ) )
    cones = [ findall( x -> x == 1, c ) for c in cones ]
    
    # return the fan
    return Fan( gap_fan, rays, cones )
end
export fan_of_variety


"""
    fan( v::AbstractNormalToricVariety )

A convenience method for `fan_of_variety( v::AbstractNormalToricVariety )`.
"""
function fan( v::AbstractNormalToricVariety )
    return fan_of_variety( v::AbstractNormalToricVariety )
end
export fan


struct CartierTorusInvariantDivisorGroup
           GapCartierTorusInvariantDivisorGroup::GapObj
end
export CartierTorusInvariantDivisorGroup


"""
    cartier_torus_invariant_divisor_group( v::AbstractNormalToricVariety )

Computes the group of Cartier and torus invariant divisors of the normal toric variety `v`.
"""
function cartier_torus_invariant_divisor_group( v::AbstractNormalToricVariety )
    gap_CartierTorusInvariantDivisorGroup = GAP.Globals.CartierTorusInvariantDivisorGroup( v.GapNTV )
    return CartierTorusInvariantDivisorGroup( gap_CartierTorusInvariantDivisorGroup )
end
export cartier_torus_invariant_divisor_group


struct PicardGroup
           GapPicardGroup::GapObj
end
export PicardGroup


"""
    picard_group( v::AbstractNormalToricVariety )

Computes the Picard group of the normal toric variety `v`.
"""
function picard_group( v::AbstractNormalToricVariety )
    gap_PicardGroup = GAP.Globals.PicardGroup( v.GapNTV )
    return PicardGroup( gap_PicardGroup )
end
export picard_group


"""
    name_of_variety( v::AbstractNormalToricVariety )

Returns the name of the normal toric variety `v`, if set. Otherwise returns "No name set for this variety".
"""
function name_of_variety( v::AbstractNormalToricVariety )
    if ! GAP.Globals.HasNameOfVariety( v.GapNTV )
            return "No name set for this variety"
    end
    
    return String( GAP.Globals.NameOfVariety( v.GapNTV ) )
end
export name_of_variety


struct ZariskiCotangentSheaf
           GapZariskiCotangentSheaf::GapObj
end
export ZariskiCotangentSheaf


"""
    zariski_cotangent_sheaf( v::AbstractNormalToricVariety )

Returns the Zariski cotangent sheaf of the normal toric variety `v`.
"""
function zariski_cotangent_sheaf( v::AbstractNormalToricVariety )
    gap_ZariskiCotangentSheaf = GAP.Globals.ZariskiCotangentSheaf( v.GapNTV )
    return ZariskiCotangentSheaf( gap_ZariskiCotangentSheaf )
end
export zariski_cotangent_sheaf


struct CotangentSheaf
           GapCotangentSheaf::GapObj
end
export CotangentSheaf


"""
    cotangent_sheaf( v::AbstractNormalToricVariety )

Returns the cotangent sheaf of the normal toric variety `v`.
"""
function cotangent_sheaf( v::AbstractNormalToricVariety )
    gap_CotangentSheaf = GAP.Globals.CotangentSheaf( v.GapNTV )
    return CotangentSheaf( gap_CotangentSheaf )
end
export cotangent_sheaf


"""
    euler_characteristic( v::AbstractNormalToricVariety )

Computes the Euler characteristic of the normal toric variety `v`.
"""
function euler_characteristic( v::AbstractNormalToricVariety )
    return GAP.Globals.EulerCharacteristic( v.GapNTV )::Int
end
export euler_characteristic


"""
    weil_divisors_of_variety( v::AbstractNormalToricVariety )

Compute the Weil divisors of the normal toric variety `v`.
"""
function weil_divisors_of_variety( v::AbstractNormalToricVariety )
    gap_divisors = GAP.Globals.WeilDivisorsOfVariety( v.GapNTV )
    return [ ToricDivisor( d ) for d in gap_divisors ]
end
export weil_divisors_of_variety


struct ZariskiCotangentSheafViaEulerSequence
           GapZariskiCotangentSheafViaEulerSequence::GapObj
end
export ZariskiCotangentSheafViaEulerSequence


"""
    zariski_cotangent_sheaf_via_euler_sequence( v::AbstractNormalToricVariety )

Computes the Zariski cotangent sheaf of the normal toric variety `v` via the Euler sequence.
"""
function zariski_cotangent_sheaf_via_euler_sequence( v::AbstractNormalToricVariety )
    gap_ZariskiCotangentSheafViaEulerSequence = GAP.Globals.ZariskiCotangentSheafViaEulerSequence( v.GapNTV )
    return ZariskiCotangentSheafViaEulerSequence( gap_ZariskiCotangentSheafViaEulerSequence )
end
export zariski_cotangent_sheaf_via_euler_sequence


struct ZariskiCotangentSheafViaPoincareResidueMap
           GapZariskiCotangentSheafViaPoincareResidueMap::GapObj
end
export ZariskiCotangentSheafViaPoincareResidueMap


"""
    zariski_cotangent_sheaf_via_poincare_residue_map( v::AbstractNormalToricVariety )

Computes the Zariski cotangent sheaf of the normal toric variety `v` via the Poincare residue map.
"""
function zariski_cotangent_sheaf_via_poincare_residue_map( v::AbstractNormalToricVariety )
    gap_ZariskiCotangentSheafViaPoincareResidueMap = GAP.Globals.ZariskiCotangentSheafViaPoincareResidueMap( v.GapNTV )
    return ZariskiCotangentSheafViaPoincareResidueMap( gap_ZariskiCotangentSheafViaPoincareResidueMap )
end
export zariski_cotangent_sheaf_via_poincare_residue_map


#struct UnderlyingSheaf
#          GapUnderlyingSheaf::GapObj
#end
#export UnderlyingSheaf

#function underlying_sheaf( v::AbstractNormalToricVariety )
#   gap_Underlying = GAP.Globals.UnderlyingSheaf( v.GapNTV )
#    return UnderlyingSheaf( gap_UnderlyingSheaf )
#end
#export underlying_sheaf


struct NefCone
           polymakeNefCone::Polymake.BigObject
end
export NefCone

"""
    nef_cone( v::NormalToricVariety )

Computes the nef cone of the normal toric variety `v`.
"""
function nef_cone( v::NormalToricVariety )
    return NefCone( v.polymakeNTV.NEF_CONE )
end
export nef_cone


struct MoriCone
           polymakeNefCone::Polymake.BigObject
end
export MoriCone

"""
    mori_cone( v::NormalToricVariety )

Computes the mori cone of the normal toric variety `v`.
"""
function mori_cone( v::NormalToricVariety )
    return MoriCone( v.polymakeNTV.MORI_CONE )
end
export mori_cone


######################
# 2: Methods of ToricVarieties
######################

"""
    set_name_of_variety( v::AbstractNormalToricVariety, name::String )

Sets the name of the normal toric variety `v` to `name`.
"""
function set_name_of_variety( v::AbstractNormalToricVariety, s::String )
    GAP.Globals.SetNameOfVariety( v.GapNTV, GapObj( s ) )
    return true
end
export set_name_of_variety


"""
    coordinate_ring_of_torus( v::AbstractNormalToricVariety, names::Vector{String} )

Compute the coordinate ring of the torus factor of the normal toric variety `v`, using `names` as label for the homogeneous coordinates.
"""
function coordinate_ring_of_torus( v::AbstractNormalToricVariety, names::Vector{String} )
    gap_names = [ GapObj( names[ i ] ) for i in 1 : size(names)[1] ]
    gap_names = GapObj( gap_names )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapNTV, gap_names )
    return CoordinateRingOfTorus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


"""
    cox_ring( v::AbstractNormalToricVariety, name::String )

Compute the Cox ring of the normal toric variety `v`, using `name` as label for the homogeneous coordinates.
"""
function cox_ring( v::AbstractNormalToricVariety, names::String )
    gap_names = GapObj( names )
    gap_CoxRing = GAP.Globals.CoxRing( v.GapNTV, gap_names )
    return CoxRing( gap_CoxRing )
end
export cox_ring


"""
    Base.:*( v::AbstractNormalToricVariety, w::AbstractNormalToricVariety )

Compute the direct product of the normal toric varieties `v` and `w`.

# Examples
```julia-repl
julia> projective_space( 2 ) * projective_space( 2 )
NormalToricVariety(GAP: <A projective toric variety of dimension 4 which is a product of 2 toric varieties>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x00005645a1a00930))
```
"""
function Base.:*( v::AbstractNormalToricVariety, w::AbstractNormalToricVariety )
    gap_NormalToricVariety = v.GapNTV * w.GapNTV
    return NormalToricVariety( gap_NormalToricVariety )
end
export *


struct CharacterToRationalFunction
           gap_CharacterToRationalFunction::GapObj
end
export CharacterToRationalFunction


"""
    character_to_rational_function( l::Vector{Int}, v::AbstractNormalToricVariety )

Turn the character `l` of the normal toric variety `v` into a rational function.
"""
function character_to_rational_function( l::Vector{Int}, v::AbstractNormalToricVariety )
    gap_CharacterToRationalFunction = GAP.Globals.CharacterToRationalFunction( GapObj( l ), v.GapNTV )
    return CharacterToRationalFunction( gap_CharacterToRationalFunction )
end
export character_to_rational_function


"""
    blowup_on_ith_minimal_torus_orbit( v::AbstractNormalToricVariety, i::Int )

Compute the blow up of the normal toric variety `v` on the i-th minimal torus orbit.
"""
function blowup_on_ith_minimal_torus_orbit( v::AbstractNormalToricVariety, i::Int )
    gap_blowup_variety = GAP.Globals.BlowUpOnIthMinimalTorusOrbit( v.GapNTV, GapObj( i ) )
    return NormalToricVariety( gap_blowup_variety )
end
export blowup_on_ith_minimal_torus_orbit


"""
    ith_betti_number( v::AbstractNormalToricVariety, i::Int )

Compute the i-th Betti number of the normal toric variety `v`.
"""
function ith_betti_number( v::AbstractNormalToricVariety, i::Int )
    return GAP.Globals.ithBettiNumber( v.GapNTV, GapObj( i ) )::Int
end
export ith_betti_number


"""
    nr_of_q_rational_points( v::AbstractNormalToricVariety, i::Int )

Compute the number of q-rational points of the normal toric variety `v`.
"""
function nr_of_q_rational_points( v::AbstractNormalToricVariety, i::Int )
    return GAP.Globals.NrOfqRationalPoints( v.GapNTV, GapObj( i ) )::Int
end
export nr_of_q_rational_points
