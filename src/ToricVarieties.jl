######################
# 1: The Julia type for ToricVarieties
######################

struct NormalToricVariety
           GapNTV::GapObj
           polymakeNTV::Polymake.BigObject
end
export NormalToricVariety


######################
# 2: Generic constructors
######################

function NormalToricVariety( rays::Matrix{Int}, cones::Vector{Vector{Int}} )
    # construct the toric variety in GAP
    gap_rays = GapObj( rays, recursive = true )
    gap_cones = GapObj( cones, recursive = true )
    fan = GAP.Globals.Fan( gap_rays, gap_cones )
    variety = GAP.Globals.ToricVariety( fan )

    Incidence = Oscar.IncidenceMatrix(cones)
    arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)

    pmntv = Polymake.fulton.NormalToricVariety(
        RAYS = Oscar.matrix_for_polymake(rays),
        MAXIMAL_CONES = arr,
    )

    # wrap it into a struct and return
    return NormalToricVariety( variety, pmntv )
end

function NormalToricVariety(GapNTV::GapObj)
   pmNTV = ntv_gap2polymake(GapNTV)
   return NormalToricVariety(GapNTV, pmNTV)
end

function ntv_gap2polymake(GapNTV::GapObj)
    ff = GAP.Globals.Fan(GapNTV)
    R = GAP.Globals.RayGenerators(ff)
    MC = GAP.Globals.RaysInMaximalCones(ff)
    rays = Matrix{Int}(R)
    cones = [findall(x->x!=0, vec) for vec in [[e for e in mc] for mc in MC]]
    Incidence = Oscar.IncidenceMatrix(cones)
    arr = @Polymake.convert_to Array{Set{Int}} Polymake.common.rows(Incidence.pm_incidencematrix)
    pmntv = Polymake.fulton.NormalToricVariety(
        RAYS = Oscar.matrix_for_polymake(rays),
        MAXIMAL_CONES = arr,
    )
    return pmntv
end
export ntv_gap2polymake

function ntv_polymake2gap(polymakeNTV::Polymake.BigObject)
    rays = Matrix{Int}(polymakeNTV.RAYS)
    gap_rays = GapObj( rays, recursive = true )
    cones = [findall(x->x!=0, v) for v in eachrow(polymakeNTV.MAXIMAL_CONES)]
    gap_cones = GapObj( cones, recursive = true )
    fan = GAP.Globals.Fan( gap_rays, gap_cones )
    variety = GAP.Globals.ToricVariety( fan )
    return variety
end
export ntv_polymake2gap


######################
# 3: Special constructors
######################

function projective_space( x )
    # construct the projective space in gap
    variety = GAP.Globals.ProjectiveSpace( x )
    f = Polymake.fan.normal_fan(Polymake.polytope.simplex(x))
    pmntv = Polymake.fulton.NormalToricVariety(f)
    
    # wrap it and return
    return NormalToricVariety( variety, pmntv )
end
export projective_space


######################
# 4: Properties
######################

function is_normal_variety( v::NormalToricVariety )
    return true
end
export is_normal_variety


function is_affine( v::NormalToricVariety )
    return v.polymakeNTV.AFFINE
end
export is_affine


function is_projective( v::NormalToricVariety )
    return v.polymakeNTV.PROJECTIVE
end
export is_projective


function is_smooth( v::NormalToricVariety )
    return v.polymakeNTV.SMOOTH
end
export is_smooth


function is_complete( v::NormalToricVariety )
    return v.polymakeNTV.COMPLETE
end
export is_complete


function has_torusfactor( v::NormalToricVariety )
    return GAP.Globals.HasTorusfactor( v.GapNTV )::Bool
end
export has_torusfactor


function is_orbifold( v::NormalToricVariety )
    return GAP.Globals.IsOrbifold( v.GapNTV )::Bool
end
export is_orbifold


function is_simplicial( v::NormalToricVariety )
    return v.polymakeNTV.SIMPLICIAL
end
export is_simplicial


function is_isomorphic_to_projective_space( v::NormalToricVariety )
    return GAP.Globals.IsIsomorphicToProjectiveSpace( v.GapNTV )::Bool
end
export is_isomorphic_to_projective_space


function is_direct_product_of_projective_spaces( v::NormalToricVariety )
    return GAP.Globals.IsDirectProductOfPNs( v.GapNTV )::Bool
end
export is_direct_product_of_projective_spaces
