module IntersectionTheory
using ..Oscar

import Base: +, -, *, ^, ==, div, zero, one, parent
import ..Oscar: AffAlgHom, Ring, MPolyDecRingElem, symmetric_power, exterior_power, OO, pullback, canonical_bundle, graph, euler_characteristic, pullback
import ..Oscar: basis, chi, codomain, degree, det, dim, domain, dual, gens, hilbert_polynomial, hom, integral, rank, signature
import ..Oscar.AbstractAlgebra: combinations
import ..Oscar.AbstractAlgebra.Generic: FunctionalMap, Partition, partitions

export betti, euler, graph
export a_hat_genus, l_genus, pontryagin, chern_number, chern_numbers
export intersection_matrix, dual_basis
export bundles, tangent_bundle, cotangent_bundle, canonical_bundle, canonical_class
export schur_functor
export complete_intersection, zero_locus_section, degeneracy_locus
export abstract_variety, bundle
export schubert_class, schubert_classes
export abstract_variety, abstract_bundle
export abstract_projective_variety
export abstract_projective_bundle
export abstract_grassmannian
export abstract_flag_variety
export chern_character
export total_chern_class
export chern_class
export top_chern_class
export segre_class
export pontryagin_class
export total_pontryagin_class
export todd_class
export euler_pairing
export line_bundle
export trivial_line_bundle

include("Types.jl")
include("Misc.jl")

include("Bott.jl")   # integration using Bott's formula
include("Main.jl")   # basic constructions for Schubert calculus
# include("Blowup.jl") # blowup
# include("Moduli.jl") # moduli of matrices, twisted cubics
# include("Weyl.jl")   # weyl groups

end # module
using .IntersectionTheory
