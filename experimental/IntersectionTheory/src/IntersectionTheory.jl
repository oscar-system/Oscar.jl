module IntersectionTheory
using ..Oscar

import Base: +, -, *, ^, ==, div, zero, one, parent
import ..Oscar: AffAlgHom, Ring, MPolyDecRingElem, symmetric_power, exterior_power, pullback, canonical_bundle, graph, euler_characteristic, pullback
import ..Oscar: basis, betti, chow_ring, codomain, degree, det, dim, domain, dual, gens, hilbert_polynomial, hom, integral, rank, signature, partitions
import ..Oscar.AbstractAlgebra: combinations
import ..Oscar.AbstractAlgebra.Generic: FunctionalMap
import..Oscar: pullback, pushforward, base, OO, product, compose
import ..Oscar: trivial_line_bundle
import ..Oscar: chern_class

export a_hat_genus
export abstract_bundle
export abstract_flag_bundle
export abstract_flag_variety
export abstract_grassmannian
export abstract_hirzebruch_surface
export abstract_point
export abstract_projective_bundle
export abstract_projective_space
export abstract_variety
export base
export betti
export blowup
export blowup_points
export bundles
export canonical_bundle
export canonical_class
export chern_character
export chern_class
export chern_number
export chern_numbers
export chow_ring
export complete_intersection
export compose
export cotangent_bundle
export degeneracy_locus
export dual_basis
export euler
export euler_pairing
export graph
export hyperplane_class
export intersection_matrix
export l_genus
export linear_subspaces_on_hypersurface
export line_bundle
export OO
export point_class
export pontryagin_class
export present_finite_extension_ring
export product
export pullback
export pushforward
export schubert_class
export schubert_classes
export schur_functor
export segre_class
export structure_map
export tangent_bundle
export tautological_bundles
export tn_flag_variety
export tn_grassmannian
export todd_class
export top_chern_class
export total_chern_class
export total_pontryagin_class
export trivial_line_bundle
export zero_locus_section

export MPolyDecRingOrQuo
export AbstractVariety
export AbstractVarietyMap
export AbstractBundle
export TnBundle
export TnVariety

include("Types.jl")
include("Misc.jl")

include("Bott.jl")   # integration using Bott's formula
include("Main.jl")   # basic constructors and functionality
include("blowup.jl") # blowup
include("schubert.jl") # Schubert calculus
# include("Moduli.jl") # moduli of matrices, twisted cubics
# include("Weyl.jl")   # weyl groups

end # module

using .IntersectionTheory

export a_hat_genus
export abstract_bundle
export abstract_flag_bundle
export abstract_flag_variety
export abstract_grassmannian
export abstract_hirzebruch_surface
export abstract_point
export abstract_projective_bundle
export abstract_projective_space
export abstract_variety
export base
export betti
export blowup
export blowup_points
export bundles
export canonical_bundle
export canonical_class
export chern_character
export chern_class
export chern_number
export chern_numbers
export chow_ring
export complete_intersection
export compose
export cotangent_bundle
export degeneracy_locus
export dual_basis
export euler
export euler_pairing
export graph
export hyperplane_class
export intersection_matrix
export l_genus
export linear_subspaces_on_hypersurface
export line_bundle
export OO
export point_class
export pontryagin_class
export present_finite_extension_ring
export product
export pullback
export pushforward
export schubert_class
export schubert_classes
export schur_functor
export segre_class
export structure_map
export tangent_bundle
export tautological_bundles
export tn_flag_variety
export tn_grassmannian
export todd_class
export top_chern_class
export total_chern_class
export total_pontryagin_class
export trivial_line_bundle
export zero_locus_section

export MPolyDecRingOrQuo
export AbstractVariety
export AbstractVarietyMap
export AbstractBundle
export TnBundle
export TnVariety
