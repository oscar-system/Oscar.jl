module LinearAlgebraicGroups

using ..Oscar

import Random

using AbstractAlgebra.PrettyPrinting

import ..Oscar: gen
import ..Oscar: gens
import ..Oscar: has_gens
import ..Oscar: number_of_generators
import ..Oscar: isfinite
import ..Oscar: is_subgroup
import ..Oscar: order
import ..Oscar: root_system
import ..Oscar: degree
import ..Oscar: elem_type
import ..Oscar: parent
import ..Oscar: parent_type
import ..Oscar: one

export LinearAlgebraicGroup, LinearAlgebraicGroupElem
export linear_algebraic_group, linear_algebraic_group_elem
export root_subgroup
export maximal_torus
export torus_element
export apply_root_to_torus_element
export representative_of_root_in_group
export borel_subgroup
export bruhat_cell
export bruhat_decomposition
export bruhat_cell_representative

include("Types.jl")
include("LinearAlgebraicGroup.jl")

end

using .LinearAlgebraicGroups

export LinearAlgebraicGroup, LinearAlgebraicGroupElem
export linear_algebraic_group, linear_algebraic_group_elem
export root_subgroup
export maximal_torus
export torus_element
export apply_root_to_torus_element
export representative_of_root_in_group
export borel_subgroup
export bruhat_cell
export bruhat_decomposition
export bruhat_cell_representative
