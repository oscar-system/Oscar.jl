function weight end
function word end

module LieAlgebras

using ..Oscar

import Oscar: GAPWrap, IntegerUnion, MapHeader

# not importet in Oscar
using AbstractAlgebra: CacheDictType, ProductIterator, get_cached!, ordinal_number_string

using AbstractAlgebra.PrettyPrinting

# functions with new methods
import ..Oscar:
  _iso_oscar_gap,
  action,
  basis_matrix,
  basis,
  canonical_injection,
  canonical_injections,
  canonical_projection,
  canonical_projections,
  center,
  centralizer,
  character,
  characteristic,
  coeff,
  coefficient_ring,
  coefficients,
  compose,
  derived_series,
  dim,
  direct_sum,
  dual,
  elem_type,
  expressify,
  exterior_power,
  gen,
  gens,
  hom,
  hom_tensor,
  ideal,
  identity_map,
  image,
  inv,
  is_abelian,
  is_exterior_power,
  is_isomorphism,
  is_nilpotent,
  is_perfect,
  is_simple,
  is_solvable,
  is_tensor_product,
  is_welldefined,
  kernel,
  lower_central_series,
  matrix,
  ngens,
  normalizer,
  parent_type,
  rank,
  root,
  roots,
  sub,
  symbols,
  symmetric_power,
  tensor_product,
  weyl_vector,
  word,
  zero_map,
  ⊕,
  ⊗

import Base: getindex, deepcopy_internal, hash, issubset, iszero, parent, zero

export AbstractLieAlgebra, AbstractLieAlgebraElem
export LieAlgebra, LieAlgebraElem
export LieAlgebraHom
export LieAlgebraIdeal
export LieAlgebraModule, LieAlgebraModuleElem
export LieAlgebraModuleHom
export LieSubalgebra
export LinearLieAlgebra, LinearLieAlgebraElem
export RootSpaceElem
export RootSystem
export WeightLatticeElem
export WeylGroup, WeylGroupElem

export abelian_lie_algebra
export abstract_module
export base_lie_algebra
export base_module
export base_modules
export bracket
export cartan_bilinear_form
export cartan_matrix
export cartan_symmetrizer
export cartan_type
export cartan_type_with_ordering
export chevalley_basis
export coefficient_vector
export coerce_to_lie_algebra_elem
export combinations
export conjugate_dominant_weight
export coxeter_matrix
export derived_algebra
export dim_of_simple_module
export dominant_character
export exterior_power
export fundamental_weight
export fundamental_weights
export general_linear_lie_algebra
export hom_direct_sum
export hom_power
export is_cartan_matrix
export is_cartan_type
export is_direct_sum
export is_dual
export is_negative_root_with_index
export is_positive_root_with_index
export is_root_with_index
export is_self_normalizing
export is_simple_root_with_index
export is_standard_module
export is_symmetric_power
export is_tensor_power
export lie_algebra
export lmul!
export longest_element
export lower_central_series
export matrix_repr_basis
export multicombinations
export negative_root
export negative_roots
export num_positive_roots
export num_roots, nroots
export num_simple_roots
export permutations
export permutations_with_sign
export positive_root
export positive_roots
export reduced_expressions
export reflect, reflect!
export root_system_type, has_root_system_type
export root_system, has_root_system
export show_dynkin_diagram
export simple_module
export simple_root
export simple_roots
export special_linear_lie_algebra
export special_orthogonal_lie_algebra
export standard_module
export symmetric_power
export tensor_power
export tensor_product_decomposition
export trivial_module
export universal_enveloping_algebra
export weyl_group
export word

include("Combinatorics.jl")
include("CartanMatrix.jl")
include("CoxeterGroup.jl")
include("RootSystem.jl")
include("DynkinDiagram.jl")
include("WeylGroup.jl")

include("Util.jl")
include("LieAlgebra.jl")
include("AbstractLieAlgebra.jl")
include("LinearLieAlgebra.jl")
include("LieSubalgebra.jl")
include("LieAlgebraIdeal.jl")
include("LieAlgebraHom.jl")
include("LieAlgebraModule.jl")
include("LieAlgebraModuleHom.jl")
include("iso_oscar_gap.jl")
include("iso_gap_oscar.jl")
include("GapWrapper.jl")

end

using .LieAlgebras

export AbstractLieAlgebra, AbstractLieAlgebraElem
export LieAlgebra, LieAlgebraElem
export LieAlgebraHom
export LieAlgebraIdeal
export LieAlgebraModule, LieAlgebraModuleElem
export LieAlgebraModuleHom
export LieSubalgebra
export LinearLieAlgebra, LinearLieAlgebraElem
export RootSpaceElem
export RootSystem
export WeightLatticeElem
export WeylGroup, WeylGroupElem

export abelian_lie_algebra
export abstract_module
export base_lie_algebra
export base_module
export base_modules
export bracket
export cartan_bilinear_form
export cartan_matrix
export cartan_symmetrizer
export cartan_type
export cartan_type_with_ordering
export chevalley_basis
export coerce_to_lie_algebra_elem
export conjugate_dominant_weight
export coxeter_matrix
export derived_algebra
export dim_of_simple_module
export dominant_character
export exterior_power
export fundamental_weight
export fundamental_weights
export general_linear_lie_algebra
export hom_direct_sum
export hom_power
export is_cartan_matrix
export is_cartan_type
export is_direct_sum
export is_dual
export is_negative_root_with_index
export is_positive_root_with_index
export is_root_with_index
export is_self_normalizing
export is_simple_root_with_index
export is_standard_module
export is_symmetric_power
export is_tensor_power
export is_tensor_product
export lie_algebra
export lmul!
export longest_element
export lower_central_series
export matrix_repr_basis
export matrix_repr_basis
export negative_root
export negative_roots
export num_positive_roots
export num_roots, nroots
export num_simple_roots
export positive_root
export positive_roots
export reduced_expressions
export reflect, reflect!
export root
export root_system_type, has_root_system_type
export root_system, has_root_system
export roots
export show_dynkin_diagram
export simple_module
export simple_root
export simple_roots
export special_linear_lie_algebra
export special_orthogonal_lie_algebra
export standard_module
export symmetric_power
export tensor_power
export tensor_product_decomposition
export trivial_module
export universal_enveloping_algebra
export weyl_group
export word
