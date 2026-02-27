module IntersectionTheory
using ..Oscar

import Base: +, -, *, ^, ==, div, zero, one, parent
import ..Oscar: AffAlgHom, Ring, MPolyDecRingElem, symmetric_power, exterior_power,
  pullback, canonical_bundle, graph, euler_characteristic, pullback
import ..Oscar: basis, betti_numbers, chow_ring, codomain, degree, det, dim, domain, dual,
  gens, hilbert_polynomial, hom, integral, rank, signature, partitions
import ..Oscar: pullback, pushforward, base, OO, product, compose, identity_map, map,
  fixed_points, number_of_fixed_points
import ..Oscar: trivial_line_bundle
import ..Oscar: intersection_matrix
import ..Oscar: chern_class
import ..Oscar: IntegerUnion
import ..Oscar: localization
import ..Oscar: _homogeneous_components
import ..AbstractAlgebra: polynomial


export a_hat_genus
export abstract_bundle
export abstract_cayley_grassmannian
export abstract_cayley_plane
export abstract_curve
export abstract_K3_surface
export abstract_quadric
export adams
export cannibalistic
export flag_bundle
export abstract_flag_variety
export abstract_grassmannian
export abstract_hirzebruch_surface
export abstract_point
export abstract_projective_space
export abstract_variety
export base
export betti_numbers
export blowup
export blowup_points
export canonical_bundle
export canonical_class
export chern_character
export chern_class
export chern_number
export chern_numbers
export chi
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
export kernel_bundle
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
export pullback_morphism
export pushforward
export pushforward_morphism
export schubert_class
export schubert_classes
export schur_functor
export section_class
export total_segre_class
export segre_class
export set_point_class
export set_tangent_bundle
export set_polarization
export set_tautological_bundles
export set_structure_map
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
export trim!
export trivial_line_bundle
export zero_locus_section

export MPolyDecRingOrQuo
export MPolyDecRingOrQuoElem
export AbstractVariety
export AbstractVarietyMap
export AbstractBundle
export TnBundle
export TnBundleChern
export TnRep
export TnVariety

include("types.jl")
include("misc.jl")

include("bott.jl")   # integration using Bott's formula
include("main.jl")   # basic constructors and functionality
include("blowup.jl") # blowup
include("schubert.jl") # Schubert calculus
include("kontsevich.jl") # Kontsevich moduli spaces
# include("Moduli.jl") # moduli of matrices, twisted cubics
# include("Weyl.jl")   # weyl groups

end # module

using .IntersectionTheory

export a_hat_genus
export abstract_bundle
export abstract_cayley_grassmannian
export abstract_cayley_plane
export abstract_curve
export abstract_K3_surface
export abstract_quadric
export adams
export cannibalistic
export flag_bundle
export abstract_flag_variety
export abstract_grassmannian
export abstract_hirzebruch_surface
export abstract_point
export abstract_projective_space
export abstract_variety
export base
export betti_numbers
export blowup
export blowup_points
export canonical_bundle
export canonical_class
export chern_character
export chern_class
export chern_number
export chern_numbers
export chi
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
export kernel_bundle
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
export pullback_morphism
export pushforward
export pushforward_morphism
export schubert_class
export schubert_classes
export schur_functor
export section_class
export total_segre_class
export segre_class
export set_point_class
export set_tangent_bundle
export set_polarization
export set_tautological_bundles
export set_structure_map
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
export trim!
export trivial_line_bundle
export zero_locus_section

export MPolyDecRingOrQuo
export MPolyDecRingOrQuoElem
export AbstractVariety
export AbstractVarietyMap
export AbstractBundle
export TnBundle
export TnBundleChern
export TnRep
export TnVariety
