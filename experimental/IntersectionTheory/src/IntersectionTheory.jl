module IntersectionTheory
using ..Oscar

import Base: +, -, *, ^, ==, div, zero, one, parent
import ..Oscar: AffAlgHom, Ring, MPolyDecRingElem, symmetric_power, exterior_power, pullback, canonical_bundle, graph, euler_characteristic, pullback
import ..Oscar: basis, betti, codomain, degree, det, dim, domain, dual, gens, hilbert_polynomial, hom, integral, rank, signature
import ..Oscar.AbstractAlgebra: combinations
import ..Oscar.AbstractAlgebra.Generic: FunctionalMap

export a_hat_genus
export abstract_bundle
export abstract_flag_variety
export abstract_grassmannian
export abstract_projective_bundle
export abstract_projective_space
export abstract_variety
export betti
export bundles
export canonical_bundle
export canonical_class
export chern_character
export chern_class
export chern_number
export chern_numbers
export complete_intersection
export cotangent_bundle
export degeneracy_locus
export dual_basis
export euler
export euler_pairing
export graph
export intersection_matrix
export l_genus
export line_bundle
export pontryagin_class
export schubert_class
export schubert_classes
export schur_functor
export segre_class
export tangent_bundle
export todd_class
export top_chern_class
export total_chern_class
export total_pontryagin_class
export trivial_line_bundle
export zero_locus_section

include("Types.jl")
include("Misc.jl")

include("Bott.jl")   # integration using Bott's formula
include("Main.jl")   # basic constructions for Schubert calculus
# include("Blowup.jl") # blowup
# include("Moduli.jl") # moduli of matrices, twisted cubics
# include("Weyl.jl")   # weyl groups

end # module

using .IntersectionTheory
