include("Types.jl")
include("CoveredScheme.jl")
include("ProjectiveModules.jl")
include("SpaceGerms.jl")
include("Tjurina.jl")

include("Auxiliary.jl")
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
include("Resolution_tools.jl")


# Exports
export CompleteIntersectionGerm
export HypersurfaceGerm
export MorphismFromRationalFunctions
export SpaceGerm

export ambient_germ
export basis_representation
export blow_up_along_minimal_supercone_coordinates
export complete_intersection_germ
export cox_ring_module_homomorphism
export defining_ring_element
export defining_ring_elements
export elliptic_parameter
export germ_at_point
export has_du_val_singularities
export hypersurface_germ
export ideal_sheaf
export index_of_exceptional_ray
export is_isolated_singularity
export intersection_matrix
export milnor_algebra
export milnor_number
export minimal_supercone_coordinates_of_exceptional_ray
export point
export rational_point_coordinates
export standard_covering
export standard_coordinates
export strict_transform
export strict_transform_with_index
export total_transform
export two_neighbor_step

export tjurina_algebra
export tjurina_number
export order_as_series
export is_finitely_determined
export determinacy_bound
export sharper_determinacy_bound
export is_contact_equivalent
export tjurina_module


# Deprecated after 0.15
Base.@deprecate_binding base_glueing base_gluing
Base.@deprecate_binding inherit_glueings! inherit_gluings!
Base.@deprecate_binding AbsProjectiveGlueing AbsProjectiveGluing
Base.@deprecate_binding CoveredProjectiveGlueingData CoveredProjectiveGluingData
Base.@deprecate_binding InheritGlueingData InheritGluingData
Base.@deprecate_binding LazyProjectiveGlueing LazyProjectiveGluing
Base.@deprecate_binding ProjectiveGlueing ProjectiveGluing
Base.@deprecate_binding ProjectiveGlueingData ProjectiveGluingData
