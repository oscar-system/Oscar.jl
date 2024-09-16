include("Types.jl")
include("CoveredScheme.jl")
include("ProjectiveModules.jl")
include("SpaceGerms.jl")
include("Tjurina.jl")
include("CoveredProjectiveSchemes.jl")

include("Auxiliary.jl")
include("BlowupMorphism.jl")
include("duValSing.jl")
include("elliptic_surface.jl")
include("critical_locus.jl")

include("ToricIdealSheaves/auxiliary.jl")
include("ToricIdealSheaves/constructors.jl")
include("ToricIdealSheaves/attributes.jl")
include("ToricIdealSheaves/methods.jl")

include("ToricDivisors/constructors.jl")
include("ToricDivisors/attributes.jl")

include("NormalToricVarieties/attributes.jl")

include("ToricBlowups/types.jl")
include("ToricBlowups/constructors.jl")
include("ToricBlowups/attributes.jl")
include("ToricBlowups/methods.jl")

include("DerivedPushforward.jl")
include("Resolution_structure.jl")

# Exports
export CompleteIntersectionGerm
export HypersurfaceGerm
export MorphismFromRationalFunctions
export SpaceGerm

export ambient_germ
export basis_representation
export complete_intersection_germ
export defining_ring_element
export defining_ring_elements
export elliptic_parameter
export germ_at_point
export has_du_val_singularities
export hypersurface_germ
export ideal_sheaf
export index_of_new_ray
export is_isolated_singularity
export milnor_algebra
export milnor_number
export point
export rational_point_coordinates
export standard_covering
export total_transform
export two_neighbor_step

export tjurina_algebra
export tjurina_number
export order_as_series
export is_finitely_determined
export determinacy_bound
export sharper_determinacy_bound
export is_contact_equivalent


# Deprecated after 0.15
Base.@deprecate_binding _compute_inherited_glueing _compute_inherited_gluing
Base.@deprecate_binding base_glueing base_gluing
Base.@deprecate_binding inherit_glueings! inherit_gluings!
Base.@deprecate_binding AbsProjectiveGlueing AbsProjectiveGluing
Base.@deprecate_binding CoveredProjectiveGlueingData CoveredProjectiveGluingData
Base.@deprecate_binding InheritGlueingData InheritGluingData
Base.@deprecate_binding LazyProjectiveGlueing LazyProjectiveGluing
Base.@deprecate_binding ProjectiveGlueing ProjectiveGluing
Base.@deprecate_binding ProjectiveGlueingData ProjectiveGluingData
