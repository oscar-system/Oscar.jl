module IntersectionTheory
using ..Oscar

import Base: +, -, *, ^, ==, div, zero, one, parent
import ..Oscar: AffAlgHom, Ring, MPolyDecRingElem, symmetric_power, exterior_power, pullback, canonical_bundle, graph, euler_characteristic, pullback
import ..Oscar: basis, betti_numbers, chow_ring, codomain, degree, det, dim, domain, dual, gens, hilbert_polynomial, hom, integral, rank, signature, partitions, blow_up
import ..Oscar: pullback, pushforward, base, OO, product, compose, identity_map, map, fixed_points
import ..Oscar: trivial_line_bundle
import ..Oscar: intersection_matrix
import ..Oscar: chern_class
import ..Oscar: IntegerUnion
import ..Oscar: localization
import ..AbstractAlgebra: polynomial


export a_hat_genus
export abstract_bundle
export flag_bundle
export abstract_flag_variety
export abstract_grassmannian
export abstract_hirzebruch_surface
export abstract_point
export abstract_projective_space
export abstract_variety
export base
export betti_numbers
export blow_up
export blow_up_points
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
export det
export dual
export dual_basis
export euler_number
export euler_pairing
export extend_inclusion
export fixed_points
export gromov_witten_invariant
export graph
export hom
export identity_map
export instanton_number
export intersection_matrix
export kontsevich_moduli_space
export l_genus
export lines_on_hypersurface
export linear_subspaces_on_hypersurface
export line_bundle
export localization
export map
export OO
export point_class
export polarization
### export r_polynomial
export polynomial
export pontryagin_class
export product
export projective_bundle
export pullback
export pushforward
export schubert_class
export schubert_classes
export schur_functor
export total_segre_class
export segre_class
export structure_map
export tangent_bundle
export tautological_bundles
export tn_bundle
export tn_flag_variety
export tn_grassmannian
export tn_representation
export tn_variety
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
export TnBundleChern
export TnRep
export TnVariety


include("Types.jl")
include("Misc.jl")

include("Bott.jl")   # integration using Bott's formula
include("Main.jl")   # basic constructors and functionality
include("blowup.jl") # blow_up
include("schubert.jl") # Schubert calculus
include("kontsevich.jl") # Kontsevich moduli spaces
# include("Moduli.jl") # moduli of matrices, twisted cubics
# include("Weyl.jl")   # weyl groups

end # module

using .IntersectionTheory

export a_hat_genus
export abstract_bundle
export flag_bundle
export abstract_flag_variety
export abstract_grassmannian
export abstract_hirzebruch_surface
export abstract_point
export abstract_projective_space
export abstract_variety
export base
export betti_numbers
export blow_up
export blow_up_points
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
export det
export dual
export dual_basis
export euler_number
export euler_pairing
export extend_inclusion
export fixed_points
export graph
export gromov_witten_invariant
export hom
export identity_map
export instanton_number
export intersection_matrix
export kontsevich_moduli_space
export l_genus
export lines_on_hypersurface
export linear_subspaces_on_hypersurface
export line_bundle
export localization
export map
export OO
export point_class
export polarization
### export r_polynomial
export polynomial
export pontryagin_class
export product
export projective_bundle
export pullback
export pushforward
export schubert_class
export schubert_classes
export schur_functor
export total_segre_class
export segre_class
export structure_map
export tangent_bundle
export tautological_bundles
export tn_bundle
export tn_flag_variety
export tn_grassmannian
export tn_representation
export tn_variety
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
export TnBundleChern
export TnRep
export TnVariety
