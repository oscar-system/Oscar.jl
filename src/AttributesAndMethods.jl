######################
# 1: Attributes of ToricVarieties
######################

"""
    affine_open_covering( v::normalToricVariety )

Computes an affine open cover of the normal toric variety `v`, i.e. returns a list of affine toric varieties.
"""
function affine_open_covering( v::normalToricVariety )
    gap_cover = GAP.Globals.AffineOpenCovering( v.GapNTV )
    return [ NormalToricVariety( v ) for v in gap_cover ]
end
export affine_open_covering


struct coxRing
           GapCoxRing::GapObj
end
export coxRing


"""
    cox_ring( v::normalToricVariety )

Computes the Cox ring of the normal toric variety `v`.
"""
function cox_ring( v::normalToricVariety )
    gap_ring = GAP.Globals.CoxRing( v.GapNTV )
    return coxRing( gap_ring )
end
export cox_ring


"""
    list_of_variables_of_cox_ring( v::normalToricVariety )

Computes the list of homogeneous variables of the Cox ring of the normal toric variety `v`.
"""
function list_of_variables_of_cox_ring(v::normalToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoxRing(v.GapNTV)
    return Vector{String}(vars)
end
export list_of_variables_of_cox_ring


struct classGroup
           GapClassGroup::GapObj
end
export classGroup


"""
    class_group( v::normalToricVariety )

Computes the class group of the normal toric variety `v`.
"""
function class_group( v::normalToricVariety )
    gap_class_group = GAP.Globals.ClassGroup( v.GapNTV )
    return classGroup( gap_class_group )
end
export class_group


struct torusInvariantDivisorGroup
           GapTorusInvariantDivisorGroup::GapObj
end
export torusInvariantDivisorGroup


"""
    torus_invariant_divisor_group( v::normalToricVariety )

Computes the torus invariant divisor class group of the normal toric variety `v`.
"""
function torus_invariant_divisor_group( v::normalToricVariety )
    gap_TorusInvariantDivisorGroup = GAP.Globals.TorusInvariantDivisorGroup( v.GapNTV )
    return torusInvariantDivisorGroup( gap_TorusInvariantDivisorGroup )
end
export torus_invariant_divisor_group


struct mapFromCharacterToPrincipalDivisor
           GapMapFromCharacterToPrincipalDivisor::GapObj
end
export mapFromCharacterToPrincipalDivisor


"""
    map_from_character_to_principal_divisor( v::normalToricVariety )

Computes the map from the character lattice to the principal divisors of the normal toric variety `v`.
"""
function map_from_character_to_principal_divisor( v::normalToricVariety )
    gap_MapFromCharacterToPrincipalDivisor = GAP.Globals.MapFromCharacterToPrincipalDivisor( v.GapNTV )
    return mapFromCharacterToPrincipalDivisor( gap_MapFromCharacterToPrincipalDivisor )
end
export map_from_character_to_principal_divisor


struct mapFromWeilDivisorsToClassGroup
           GapMapFromWeilDivisorsToClassGroup::GapObj
end
export mapFromWeilDivisorsToClassGroup


"""
    map_from_weil_divisors_to_class_group( v::normalToricVariety )

Computes the map from the Weil divisors to the Class group of the normal toric variety `v`.
"""
function map_from_weil_divisors_to_class_group( v::normalToricVariety )
    gap_MapFromWeilDivisorsToClassGroup = GAP.Globals.MapFromWeilDivisorsToClassGroup( v.GapNTV )
    return mapFromWeilDivisorsToClassGroup( gap_MapFromWeilDivisorsToClassGroup )
end
export map_from_weil_divisors_to_class_group


"""
    dim( v::normalToricVariety )

Computes the dimension of the normal toric variety `v`.
"""
function dim( v::normalToricVariety )
    return GAP.Globals.Dimension(v.GapNTV)::Int
end
export dim


"""
    dim_of_torusfactor( v::normalToricVariety )

Computes the dimension of the torus factor of the normal toric variety `v`.
"""
function dim_of_torusfactor( v::normalToricVariety )
    return GAP.Globals.DimensionOfTorusfactor( v.GapNTV )::Int
end
export dim_of_torusfactor


struct coordinateRingOfTorus
           GapCoordinateRingOfTorus::GapObj
end
export coordinateRingOfTorus


"""
    coordinate_ring_of_torus( v::normalToricVariety )

Computes the coordinate ring of the torus of the normal toric variety `v`.
"""
function coordinate_ring_of_torus( v::normalToricVariety )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapNTV )
    return coordinateRingOfTorus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


"""
    list_of_variables_of_coordinate_ring_of_torus( v::normalToricVariety )

Computes the list of homogeneous coordinates of the coordinate ring of the torus of the normal toric variety `v`.
"""
function list_of_variables_of_coordinate_ring_of_torus(v::normalToricVariety)
    vars = GAP.Globals.ListOfVariablesOfCoordinateRingOfTorus(v.GapNTV)
    return Vector{String}(vars)
end
export list_of_variables_of_coordinate_ring_of_torus


"""
    is_product_of( v::normalToricVariety )

Identifies the factors from which the normal toric variety `v` has been constructed in GAP.
"""
function is_product_of(v::normalToricVariety)
    factors = GAP.Globals.IsProductOf(v.GapNTV)
    return [ NormalToricVariety( f ) for f in factors ]
end
export is_product_of


"""
    factors( v::normalToricVariety )

Identifies the factors from which the normal toric variety `v` has been constructed in GAP.
"""
function factors( v::normalToricVariety )
    gap_factors = GAP.Globals.Factors( v.GapNTV )
    return [ NormalToricVariety( f ) for f in gap_factors ]
end
export factors


struct character_lattice
           GapCharacterLattice::GapObj
end
export character_lattice


"""
    character_lattice( v::normalToricVariety )

Computes the character lattice of the normal toric variety `v`.
"""
function character_lattice( v::normalToricVariety )
    gap_CharacterLattice = GAP.Globals.CharacterLattice( v.GapNTV )
    return character_lattice( gap_CharacterLattice )
end
export character_lattice


"""
    torus_invariant_prime_divisors( v::normalToricVariety )

Computes the torus invariant prime divisors of the normal toric variety `v`.
"""
function torus_invariant_prime_divisors( v::normalToricVariety )
    divisors = GAP.Globals.TorusInvariantPrimeDivisors( v.GapNTV )
    return [ toricDivisor( d ) for d in divisors ]    
end
export torus_invariant_prime_divisors


struct irrelevant_ideal
           GapIrrelevantIdeal::GapObj
end
export irrelevant_ideal


"""
    irrelevant_ideal( v::normalToricVariety )

Computes the irrelevant ideal of the normal toric variety `v`.
"""
function irrelevant_ideal( v::normalToricVariety )
    gap_IrrelevantIdeal = GAP.Globals.IrrelevantIdeal( v.GapNTV )
    return irrelevant_ideal( gap_IrrelevantIdeal )
end
export irrelevant_ideal


struct stanley_reisner_ideal
           GapSRIdeal::GapObj
end
export stanley_reisner_ideal


"""
    stanley_reisner_ideal( v::normalToricVariety )

Computes the Stanley-Reisner ideal of the normal toric variety `v`.
"""
function stanley_reisner_ideal( v::normalToricVariety )
    gap_SRIdeal = GAP.Globals.SRIdeal( v.GapNTV )
    return stanley_reisner_ideal( gap_SRIdeal )
end
export stanley_reisner_ideal


struct morphism_from_cox_variety
           GapMorphismFromCoxVariety::GapObj
end
export morphism_from_cox_variety


"""
    morphism_from_cox_variety( v::normalToricVariety )

Computes the morphism from the Cox variety of the normal toric variety `v`.
"""
function morphism_from_cox_variety( v::normalToricVariety )
    gap_MorphismFromCoxVariety = GAP.Globals.MorphismFromCoxVariety( v.GapNTV )
    return morphism_from_cox_variety( gap_MorphismFromCoxVariety )
end
export morphism_from_cox_variety


"""
    cox_variety( v::normalToricVariety )

Computes the Cox variety of the normal toric variety `v`.
"""
function cox_variety( v::normalToricVariety )
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


"""
    fan_of_variety( v::normalToricVariety )

Computes the fan of the normal toric variety `v`.
"""
function fan_of_variety( v::normalToricVariety )
    # collect data
    gap_fan = GAP.Globals.FanOfVariety( v.GapNTV )
    rays = Vector{Vector{Int}}( GAP.Globals.RayGenerators( gap_fan ) )
    cones = Vector{Vector{Int}}( GAP.Globals.RaysInMaximalCones( gap_fan ) )
    cones = [ findall( x -> x == 1, c ) for c in cones ]
    
    # return the fan
    return fan( gap_fan, rays, cones )
end
export fan_of_variety


"""
    fan( v::normalToricVariety )

A convenience method for `fan_of_variety( v::normalToricVariety )`.
"""
function fan( v::normalToricVariety )
    return fan_of_variety( v::normalToricVariety )
end
export fan


struct cartier_torus_invariant_divisor_group
           GapCartierTorusInvariantDivisorGroup::GapObj
end
export cartier_torus_invariant_divisor_group


"""
    cartier_torus_invariant_divisor_group( v::normalToricVariety )

Computes the group of Cartier and torus invariant divisors of the normal toric variety `v`.
"""
function cartier_torus_invariant_divisor_group( v::normalToricVariety )
    gap_CartierTorusInvariantDivisorGroup = GAP.Globals.CartierTorusInvariantDivisorGroup( v.GapNTV )
    return cartier_torus_invariant_divisor_group( gap_CartierTorusInvariantDivisorGroup )
end
export cartier_torus_invariant_divisor_group


struct picard_group
           GapPicardGroup::GapObj
end
export picard_group


"""
    picard_group( v::normalToricVariety )

Computes the Picard group of the normal toric variety `v`.
"""
function picard_group( v::normalToricVariety )
    gap_PicardGroup = GAP.Globals.PicardGroup( v.GapNTV )
    return picard_group( gap_PicardGroup )
end
export picard_group


"""
    name_of_variety( v::normalToricVariety )

Returns the name of the normal toric variety `v`, if set. Otherwise returns "No name set for this variety".
"""
function name_of_variety( v::normalToricVariety )
    if ! GAP.Globals.HasNameOfVariety( v.GapNTV )
            return "No name set for this variety"
    end
    
    return String( GAP.Globals.NameOfVariety( v.GapNTV ) )
end
export name_of_variety


struct zariski_cotangent_sheaf
           GapZariskiCotangentSheaf::GapObj
end
export zariski_cotangent_sheaf


"""
    zariski_cotangent_sheaf( v::normalToricVariety )

Returns the Zariski cotangent sheaf of the normal toric variety `v`.
"""
function zariski_cotangent_sheaf( v::normalToricVariety )
    gap_ZariskiCotangentSheaf = GAP.Globals.ZariskiCotangentSheaf( v.GapNTV )
    return zariski_cotangent_sheaf( gap_ZariskiCotangentSheaf )
end
export zariski_cotangent_sheaf


struct cotangent_sheaf
           GapCotangentSheaf::GapObj
end
export cotangent_sheaf


"""
    cotangent_sheaf( v::normalToricVariety )

Returns the cotangent sheaf of the normal toric variety `v`.
"""
function cotangent_sheaf( v::normalToricVariety )
    gap_CotangentSheaf = GAP.Globals.CotangentSheaf( v.GapNTV )
    return cotangent_sheaf( gap_CotangentSheaf )
end
export cotangent_sheaf


"""
    euler_characteristic( v::normalToricVariety )

Computes the Euler characteristic of the normal toric variety `v`.
"""
function euler_characteristic( v::normalToricVariety )
    return GAP.Globals.EulerCharacteristic( v.GapNTV )::Int
end
export euler_characteristic


"""
    weil_divisors_of_variety( v::normalToricVariety )

Compute the Weil divisors of the normal toric variety `v`.
"""
function weil_divisors_of_variety( v::normalToricVariety )
    gap_divisors = GAP.Globals.WeilDivisorsOfVariety( v.GapNTV )
    return [ toricDivisor( d ) for d in gap_divisors ]
end
export weil_divisors_of_variety


struct zariski_cotangent_sheaf_via_euler_sequence
           GapZariskiCotangentSheafViaEulerSequence::GapObj
end
export zariski_cotangent_sheaf_via_euler_sequence


"""
    zariski_cotangent_sheaf_via_euler_sequence( v::normalToricVariety )

Computes the Zariski cotangent sheaf of the normal toric variety `v` via the Euler sequence.
"""
function zariski_cotangent_sheaf_via_euler_sequence( v::normalToricVariety )
    gap_ZariskiCotangentSheafViaEulerSequence = GAP.Globals.ZariskiCotangentSheafViaEulerSequence( v.GapNTV )
    return zariski_cotangent_sheaf_via_euler_sequence( gap_ZariskiCotangentSheafViaEulerSequence )
end
export zariski_cotangent_sheaf_via_euler_sequence


struct zariski_cotangent_sheaf_via_poincare_residue_map
           GapZariskiCotangentSheafViaPoincareResidueMap::GapObj
end
export zariski_cotangent_sheaf_via_poincare_residue_map


"""
    zariski_cotangent_sheaf_via_poincare_residue_map( v::normalToricVariety )

Computes the Zariski cotangent sheaf of the normal toric variety `v` via the Poincare residue map.
"""
function zariski_cotangent_sheaf_via_poincare_residue_map( v::normalToricVariety )
    gap_ZariskiCotangentSheafViaPoincareResidueMap = GAP.Globals.ZariskiCotangentSheafViaPoincareResidueMap( v.GapNTV )
    return zariski_cotangent_sheaf_via_poincare_residue_map( gap_ZariskiCotangentSheafViaPoincareResidueMap )
end
export zariski_cotangent_sheaf_via_poincare_residue_map


#struct underlying_sheaf
#          GapUnderlyingSheaf::GapObj
#end
#export underlying_sheaf

#function underlying_sheaf( v::normalToricVariety )
#   gap_Underlying = GAP.Globals.UnderlyingSheaf( v.GapNTV )
#    return underlying_sheaf( gap_UnderlyingSheaf )
#end
#export underlying_sheaf


######################
# 2: Methods of ToricVarieties
######################

"""
    set_name_of_variety( v::normalToricVariety, name::String )

Sets the name of the normal toric variety `v` to `name`.
"""
function set_name_of_variety( v::normalToricVariety, s::String )
    GAP.Globals.SetNameOfVariety( v.GapNTV, GapObj( s ) )
    return true
end
export set_name_of_variety


"""
    coordinate_ring_of_torus( v::normalToricVariety, names::Vector{String} )

Compute the coordinate ring of the torus factor of the normal toric variety `v`, using `names` as label for the homogeneous coordinates.
"""
function coordinate_ring_of_torus( v::normalToricVariety, names::Vector{String} )
    gap_names = [ GapObj( names[ i ] ) for i in 1 : size(names)[1] ]
    gap_names = GapObj( gap_names )
    gap_CoordinateRingOfTorus = GAP.Globals.CoordinateRingOfTorus( v.GapNTV, gap_names )
    return coordinateRingOfTorus( gap_CoordinateRingOfTorus )
end
export coordinate_ring_of_torus


"""
    cox_ring( v::normalToricVariety, name::String )

Compute the Cox ring of the normal toric variety `v`, using `name` as label for the homogeneous coordinates.
"""
function cox_ring( v::normalToricVariety, names::String )
    gap_names = GapObj( names )
    gap_CoxRing = GAP.Globals.CoxRing( v.GapNTV, gap_names )
    return coxRing( gap_CoxRing )
end
export cox_ring


"""
    Base.:*( v::normalToricVariety, w::normalToricVariety )

Compute the direct product of the normal toric varieties `v` and `w`.

# Examples
```julia-repl
julia> projective_space( 2 ) * projective_space( 2 )
NormalToricVariety(GAP: <A projective toric variety of dimension 4 which is a product of 2 toric varieties>, Polymake.BigObjectAllocated(Ptr{Nothing} @0x00005645a1a00930))
```
"""
function Base.:*( v::normalToricVariety, w::normalToricVariety )
    gap_NormalToricVariety = v.GapNTV * w.GapNTV
    return NormalToricVariety( gap_NormalToricVariety )
end
export *


struct character_to_rational_function
           gap_CharacterToRationalFunction::GapObj
end
export character_to_rational_function


"""
    character_to_rational_function( l::Vector{Int}, v::normalToricVariety )

Turn the character `l` of the normal toric variety `v` into a rational function.
"""
function character_to_rational_function( l::Vector{Int}, v::normalToricVariety )
    gap_CharacterToRationalFunction = GAP.Globals.CharacterToRationalFunction( GapObj( l ), v.GapNTV )
    return character_to_rational_function( gap_CharacterToRationalFunction )
end
export character_to_rational_function


"""
    blowup_on_ith_minimal_torus_orbit( v::normalToricVariety, i::Int )

Compute the blow up of the normal toric variety `v` on the i-th minimal torus orbit.
"""
function blowup_on_ith_minimal_torus_orbit( v::normalToricVariety, i::Int )
    gap_blowup_variety = GAP.Globals.BlowUpOnIthMinimalTorusOrbit( v.GapNTV, GapObj( i ) )
    return NormalToricVariety( gap_blowup_variety )
end
export blowup_on_ith_minimal_torus_orbit


"""
    ith_betti_number( v::normalToricVariety, i::Int )

Compute the i-th Betti number of the normal toric variety `v`.
"""
function ith_betti_number( v::normalToricVariety, i::Int )
    return GAP.Globals.ithBettiNumber( v.GapNTV, GapObj( i ) )::Int
end
export ith_betti_number


"""
    nr_of_q_rational_points( v::normalToricVariety, i::Int )

Compute the number of q-rational points of the normal toric variety `v`.
"""
function nr_of_q_rational_points( v::normalToricVariety, i::Int )
    return GAP.Globals.NrOfqRationalPoints( v.GapNTV, GapObj( i ) )::Int
end
export nr_of_q_rational_points
